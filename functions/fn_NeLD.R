# ----------- NeLD formulas ------------------------
# Effective population size (Ne) from Linkge Disequilibrium (LD)

# formula for calculating how much r2 is expected to come from sampling error due to sample size  (R2S)
Calc.r2S <- function(n) if (n < 30) {    # with sample size either <30 or > 30
  0.0018 + (0.907/n) + (4.44/n^2)
} else {
  (1/n) + (3.19/n^2)
}

# (0.9/11) + (4/11^2)
#(1/11) + (0.907/(11^2))

# formula for calculating Ne with sample size <30 or > 30
Calc.NeLD <- function(r2_adj, n) {
  if (n < 30) {
    (0.308 + sqrt(0.308^2 - 2.08*r2_adj))/(2*r2_adj)
  } else {
    (1/3 + sqrt(1/9 - 2.76*r2_adj))/(2*r2_adj)  ## check brackets of formulas
  }
}

# formula for correcting for finite number of chromasomes -  Bettongs have 11n chromasomes (2n = 22) (from Waples et al 2016)
ChromAdj.Ne <- function(Ne_raw) {
  Ne <-  Ne_raw/ (0.098 + 0.219*log(11))
}
