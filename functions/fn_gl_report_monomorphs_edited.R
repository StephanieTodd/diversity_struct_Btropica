gl.report.monomorphs.edited <- function (gl) 
{
  x <- gl
  cat("Identifying monomorphic loci\n")
  a <- vector(mode = "logical", length = nLoc(x))
  for (i in 1:nLoc(x)) {
    a[i] <- NA
  }
  b <- vector(mode = "logical", length = nLoc(x))
  for (i in 1:nLoc(x)) {
    b[i] <- NA
  }
  c <- vector(mode = "logical", length = nLoc(x))
  for (i in 1:nLoc(x)) {
    c[i] <- NA
  }
  d <- vector(mode = "logical", length = nLoc(x))
  for (i in 1:nLoc(x)) {
    d[i] <- NA
  }
  
  xmat <- as.matrix(x)
  for (i in (1:nLoc(x))) {
    if (all(is.na(xmat[, i]))) {
      d[i] <- TRUE
      a[i] <- FALSE
      b[i] <- FALSE
      c[i] <- FALSE
    }
    else {
      a[i] <- all(xmat[, i] == 0, na.rm = TRUE)
      b[i] <- all(xmat[, i] == 2, na.rm = TRUE)
      c[i] <- all(xmat[, i] == 1, na.rm = TRUE)
      d[i] <- FALSE
    }
    
  }
  s1 <- sum(a, na.rm = TRUE) + sum(b, na.rm = TRUE) + sum(c, 
                                                          na.rm = TRUE)
  s2 <- s1 - sum(d, na.rm = TRUE)
  polym <- nLoc(x) - s2
  cat("\nBreakdown of", nLoc(x), "loci\n")
  cat("  Polymorphic loci:", polym, "\n  Monomorphic loci:", 
      s1, "\n  Loci with no scores (all NA):", sum(d), 
      "\n")
}
