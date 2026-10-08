import pandas as pd,numpy as np,json,re
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold,LeaveOneOut
from sklearn.metrics import roc_auc_score,average_precision_score,f1_score
d=pd.read_csv('/tmp/emb/Data S2.tsv',sep='\t'); meta=d.iloc[:,:3]; X=d.iloc[:,3:].astype(np.float32); cols=list(X.columns); A=X.values.T
if A.max()>100: A=np.log2(A+1)
pre=np.array([re.split(r'[ _#]',c)[0] for c in cols]); arr=np.array([c.startswith('Arrested') and 'Resveratrol' not in c for c in cols]); res=np.array(['Resveratrol' in c for c in cols])
typ=np.array([re.search(r'Type (\w+)',c).group(1) if re.search(r'Type (\w+)',c) else '' for c in cols])
R={'d1':dict(n_samples=len(cols),n_features=A.shape[1],arrested=int(arr.sum()),resveratrol=int(res.sum()),arrest_types={k:int((typ[arr]==k).sum()) for k in set(typ[arr])},groups={k:int((pre==k).sum()) for k in sorted(set(pre))},feature_classes=meta.iloc[:,2].value_counts().head(6).to_dict())}
def topvar(Atr,k): return np.argsort(-Atr.var(0))[:k]
def fitpred(Atr,ytr,Ate,k=2000,C=0.05):
    ix=topvar(Atr,k); mu=Atr[:,ix].mean(0); sd=Atr[:,ix].std(0)+1e-6; m=LogisticRegression(C=C,max_iter=400,class_weight='balanced').fit((Atr[:,ix]-mu)/sd,ytr); return m.decision_function((Ate[:,ix]-mu)/sd)
def cv(sel,y,k=2000,rep=3,shuffle=False,feats=None):
    M=A[sel] if feats is None else A[sel][:,feats]; yy=y[sel]; au=[];ap=[]
    for r in range(rep):
        o=np.zeros(len(yy)); rng=np.random.default_rng(r)
        for tr,te in StratifiedKFold(5,shuffle=True,random_state=r).split(M,yy):
            ytr=rng.permutation(yy[tr]) if shuffle else yy[tr]; o[te]=fitpred(M[tr],ytr,M[te],k)
        au.append(roc_auc_score(yy,o)); ap.append(average_precision_score(yy,o))
    return dict(auroc=round(float(np.mean(au)),4),auprc=round(float(np.mean(ap)),4),n=int(len(yy)),pos=int(yy.sum()))
def save(): json.dump(R,open('results/results.json','w'),indent=1,default=float)
clv=np.isin(pre,['2C','4C','8C','E3']); sel=arr|clv; y=arr.astype(int)
R['d2']=cv(sel,y); R['d2_gate_ge_0.90']=bool(R['d2']['auroc']>=.9); print('d2',R['d2'],flush=True); save()
c3=np.isin(pre,['E3','8C','4C'])&~arr; y3=(pre=='E3').astype(int)
R['d3_E3_vs_8C4C_same_stage_cross_source']=cv(c3,y3); R['d3_batch_dominated']=bool(R['d3_E3_vs_8C4C_same_stage_cross_source']['auroc']>=R['d2']['auroc']-.05); print('d3',R['d3_E3_vs_8C4C_same_stage_cross_source'],flush=True); save()
# d4
s4=np.where(arr&np.isin(typ,['I','II','III']))[0]; yt=typ[s4]; pr=[]
for i in range(len(s4)):
    tr=np.delete(np.arange(len(s4)),i); ix=topvar(A[s4][tr],500); mu=A[s4][tr][:,ix].mean(0); sd=A[s4][tr][:,ix].std(0)+1e-6
    pr.append(LogisticRegression(C=.05,max_iter=300).fit((A[s4][tr][:,ix]-mu)/sd,yt[tr]).predict(((A[s4][i][ix]-mu)/sd).reshape(1,-1))[0])
R['d4']=dict(n=len(s4),counts={k:int((yt==k).sum()) for k in set(yt)},macroF1=round(float(f1_score(yt,pr,average='macro')),3),majority_acc=round(float(pd.Series(yt).value_counts(normalize=True).iloc[0]),3)); print('d4',R['d4'],flush=True)
# d5
genes=['DUX4','ZSCAN4','BTG4','ZNF280A','LEUTX','TPRX1','KLF17','DPPA3']; nm=meta.iloc[:,1].values; gi=[np.where(nm==g)[0][0] for g in genes if (nm==g).any()]
sc=A[:,gi].mean(1); R['d5']=dict(genes_found=[genes[i] for i,g in enumerate(genes) if (nm==g).any()],auroc_mean_score_arrested_vs_cleavage=round(float(roc_auc_score(arr[sel].astype(int),sc[sel])),4),auroc_sign_free=round(float(max(roc_auc_score(arr[sel].astype(int),sc[sel]),1-roc_auc_score(arr[sel].astype(int),sc[sel]))),4)); save()
R['d6']={str(k):cv(sel,y,k=k,rep=2) for k in (50,500,2000,5000)}; print('d6',R['d6'],flush=True); save()
R['d7_shuffle']=cv(sel,y,rep=3,shuffle=True); print('d7',R['d7_shuffle'],flush=True); save()
tt=meta.iloc[:,2].values; R['d8']={}
for nm_,mk in (('protein_coding',tt=='protein_coding'),('non_protein_coding',tt!='protein_coding')):
    R['d8'][nm_]=cv(sel,y,rep=2,feats=np.where(mk)[0]) if mk.sum()>100 else None
save()
s9=arr|res; y9=res.astype(int); R['d9_resveratrol_vs_untreated_arrested']=cv(arr|res,y9,k=500,rep=3); save()
nn=~(pre=='Late'); R['d10_arrested_vs_all_nonarrested_stagemismatched']=cv(np.ones(len(cols),bool)&~res,y,rep=2); save(); print(json.dumps(R,indent=1,default=float)[:3000])
