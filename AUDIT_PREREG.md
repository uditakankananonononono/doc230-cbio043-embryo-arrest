# AUDIT PREREG - 043 (written before the audit ran)
Questions the original run left open. Same data (SHA256 re-verified), same code (src/run.py functions re-used via exec of its first part).
Q1 D9 (resveratrol-treated vs untreated arrested, n 30, 7 vs 23, AUROC 0.959): is it beyond chance? Permutation test: 200 label permutations, same CV, p = (1 + #null >= observed)/(201). Pre-specified: p < 0.05 -> signal beyond chance (but treated and untreated arrested embryos may still differ by batch/date, which the data cannot test).
Q2 D7 anomaly (shuffle AUROC 0.309): estimate the null distribution of the D2 CV AUROC under label permutation (60 permutations, 1 repeat each). Pre-specified: if the null mean is < 0.45 the D7 value is explained as a known negative bias of permuted-label CV; otherwise unexplained.
Q3 D3 control validity: report D3 separately for E3 vs 4C and E3 vs 8C (two stage-matched contrasts), and arrested vs each control group alone (2C, 4C, 8C, E3). Pre-specified: if every contrast >= 0.95, "arrest is not separable from source" stands as written.
Q4 D5 (ZGA score) significance: Mann-Whitney U p-value, arrested vs cleavage controls, sign-free.
No change to original gates or verdicts unless an audit result contradicts them; any change will be written into RESULTS and README.
