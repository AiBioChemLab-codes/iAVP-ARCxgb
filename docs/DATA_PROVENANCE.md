# Data provenance

The iAVP-ARCXGB benchmark dataset was assembled from public peptide databases and published datasets described in the manuscript and Supporting Information.

## Source composition recorded in the Supporting Information

| Source | AVPs | non-AVPs | Total |
|---|---:|---:|---:|
| AVPdb | 2683 | 0 | 2683 |
| HIPdb | 981 | 0 | 981 |
| APD3 | 252 | 0 | 252 |
| CAMPR3 | 117 | 0 | 117 |
| DBAASP | 1521 | 0 | 1521 |
| DRAMP | 2004 | 0 | 2004 |
| LAMP2 | 4361 | 0 | 4361 |
| YADAMP | 3459 | 0 | 3459 |
| AVPpred | 604 | 1056 | 1660 |
| AVPIden | 2662 | 10095 | 12757 |
| AI4AVP | 2641 | 293 | 2934 |
| Guan et al. | 2662 | 2662 | 5324 |
| Beltran et al. | 1337 | 1337 | 2674 |
| Vishnepolsky et al. | 246 | 246 | 492 |

After sequence curation, length filtering, duplicate removal, and CD-HIT redundancy reduction at 30% sequence similarity, the final processed dataset contained 292 AVPs and 1,042 non-AVPs (1,334 sequences in total).

The released `data/train.fasta` and `data/test.fasta` files contain the final processed sequences used by the project. The manuscript and Supporting Information provide the bibliographic references for each source database and published dataset.
