import json
src=open('src/run.py').read().split("clv=np.isin(pre")[0]
exec(src)
from scipy.stats import mannwhitneyu
clv=np.isin(pre,['2C','4C','8C','E3']); sel=arr|clv; y=arr.astype(int)
R={}
# Q1
s9=arr|res; y9=res.astype(int); obs=cv(s9,y9,k=500,rep=1)['auroc']
def perm_cv(sel_,yy,k,seed):
    M=A[sel_]; y_=yy[sel_]; rng=np.random.default_rng(1000+seed); yp=rng.permutation(y_); o=np.zeros(len(y_))
    for tr,te in StratifiedKFold(5,shuffle=True,random_state=0).split(M,yp): o[te]=fitpred(M[tr],yp[tr],M[te],k)
    return roc_auc_score(yp,o)
# observed with identical single-repeat CV
M=A[s9]; y_=y9[s9]; o=np.zeros(len(y_))
for tr,te in StratifiedKFold(5,shuffle=True,random_state=0).split(M,y_): o[te]=fitpred(M[tr],y_[tr],M[te],500)
obs=roc_auc_score(y_,o); null=[perm_cv(s9,y9,500,i) for i in range(200)]
R['Q1']=dict(observed=float(obs),null_mean=float(np.mean(null)),null_max=float(np.max(null)),p=float((1+sum(n>=obs for n in null))/201)); print('Q1',R['Q1'],flush=True)
# Q2
null2=[perm_cv(sel,y,2000,i) for i in range(60)]
R['Q2']=dict(null_mean=float(np.mean(null2)),null_sd=float(np.std(null2)),null_min=float(np.min(null2)),null_max=float(np.max(null2))); print('Q2',R['Q2'],flush=True)
# Q3
R['Q3']={}
for nm_,grp in (('E3_vs_4C',['E3','4C']),('E3_vs_8C',['E3','8C'])):
    m=np.isin(pre,grp)&~arr; R['Q3'][nm_]=cv(m,(pre=='E3').astype(int),rep=2)
for g_ in ['2C','4C','8C','E3']:
    m=arr|((pre==g_)&~arr); R['Q3']['arrested_vs_'+g_]=cv(m,y,rep=2)
print('Q3',R['Q3'],flush=True)
# Q4
genes=['DUX4','ZSCAN4','BTG4','ZNF280A','LEUTX','TPRX1','KLF17','DPPA3']; nm=meta.iloc[:,1].values; gi=[np.where(nm==g)[0][0] for g in genes if (nm==g).any()]
sc=A[:,gi].mean(1); a=sc[arr]; c=sc[clv&~arr]; R['Q4']=dict(n_genes=len(gi),mannwhitney_p=float(mannwhitneyu(a,c).pvalue),mean_arrested=float(a.mean()),mean_control=float(c.mean())); print('Q4',R['Q4'])
json.dump(R,open('results/audit.json','w'),indent=1,default=float); print('done')
