# Spearman-s-rho
Repository contains Supplementary Material to "The Exact Null Distribution of Spearman’s ρ for n = 26(1)30", written by J.W.D.Eaton and W.F.Scott. 
## Repositary Contents
### Scripts
FORTRAN scripts used to compute frequency tables of S statistic. Note, scripts must be altered to personal file path when writing CSV's
- perm_mod.f90, computes frequency tables of S modulo m for given range of n.
- perm_crt.f90, performs CRT on modulo frequencies to restore full tables

### Frequency Tables
Frequency tables of S and CDF of rho for n = 18 to 30.
- cdf_ni.csv, frequency table for n = i
