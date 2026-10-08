# CBIO043 results - partial: transcriptome only; batch-dominated. (results/results.json, src/run.py, results/run.log)
Data: PLOS Biol 3001682 normalized matrix (figshare 19775992 Data S2): 1,027 samples x 17,648 features (15,189 protein-coding, 1,426 lncRNA, 1,033 TE); 23 arrested (Type I 10, II 5, III 8) + 7 resveratrol-treated arrested; controls from other public datasets.
- D2 arrested vs cleavage-stage controls (2C/4C/8C/E3; n 95, 23 arrested): AUROC 1.000 (gate >= 0.90 met).
- D3 batch control, E3 vs 4C/8C (same biological stage, different source): AUROC 1.000, so source alone is perfectly separable. By the preregistered rule D2 is BATCH-DOMINATED: a perfect arrested-vs-control score cannot be read as biology.
- D4 arrest type (I/II/III) leave-one-out macro-F1 1.00 (n 23; majority 0.43). The types were defined in the source paper by clustering the same expression, so this is circular, not independent validation.
- D5 preregistered ZGA/maternal-clearance score (8 genes): AUROC 0.666 sign-free (0.334 in the stated direction, i.e. the score is lower in arrested, opposite to a simple "ZGA failure up" reading); weak.
- D6 feature count 50/500/2000/5000: AUROC 1.000 at all sizes. D8 protein-coding only and non-coding/TE only: both 1.000.
- D7 label shuffle: AUROC 0.309 (below 0.5, a known artifact of permuted labels in tiny, class-balanced CV; not above chance).
- D9 resveratrol-treated vs untreated arrested: AUROC 0.959 (n 30, 7 treated).
- D10 arrested vs all non-arrested incl. blastocysts (stage-mismatched): 1.000, expected inflated.
Honest reading: arrest cannot be separated from dataset of origin here, so these results show source separability, not an embryo-arrest predictor. A valid test needs arrested and normal embryos from the same lab and protocol, which this public matrix does not provide. "Multi-omic" not tested (transcriptome only). 0 novel findings.

## Audit (second pass) - verdict unchanged: batch-dominated, do not count D2 as an arrest predictor
Prereg AUDIT_PREREG.md (committed before running); src/audit.py, results/audit.json, results/audit.log. Data SHA256 re-verified.
- Q1 D9 permutation test (200 permutations, single-repeat CV): observed AUROC 0.969, null mean 0.472 (max 0.925), p = 0.005. Beyond chance by the pre-specified rule. Caveat: treated and untreated arrested embryos may still differ by batch or date; the matrix cannot test that. n 7 vs 23.
- Q2 D7 null (60 permutations of the D2 CV): null mean 0.454, sd 0.093, range 0.21-0.66. The observed D7 0.309 lies inside the null range (about 1.6 sd below the null mean). By my pre-specified cutoff (null mean < 0.45) the negative bias is NOT fully established (0.454); the honest statement is "inside the permutation null, mild negative bias of permuted-label CV; not above chance".
- Q3 D3 control: E3 vs 4C 1.000, E3 vs 8C 1.000; arrested vs each of 2C, 4C, 8C, E3 alone all 1.000. All >= 0.95, so "arrest is not separable from source of origin here" stands as written.
- Q4 D5 ZGA score: Mann-Whitney p = 0.017 (not multiplicity-corrected), score lower in arrested (mean 8.38 vs 8.98), opposite to a "ZGA failure up" reading.
No gate or verdict changed.
