gl.report.hwe.edited <- function (gl, p = 0.05, subset = "each") 
{
  x <- gl
  flag <- 0
  if (class(x) != "genlight") {
    cat("Fatal Error: genlight object required!\n")
    stop()
  }
  if (subset[1] == "all") {
    gl2 <- x
    flag <- 1
  }
  else if (subset[1] == "each") {
    poplist <- seppop(x)
  }
  else if (nPop(x[pop(x) %in% subset])) {
    flag <- 1
    gl2 <- x[pop(x) %in% subset]
  }
  else {
    cat("Fatal Error: subset parameter must be \"each\", \"all\", or a list of populations\n")
    stop()
  }
  if (flag == 1) {
    result <- utils.hwe(gl2, prob = p)
  }
  else {
    count <- 0
    for (i in poplist) {
      count <- count + 1
      if (count == 1) {
        result <- utils.hwe(i, prob = p)
        Population <- rep(names(poplist)[count], nrow(result))
        result <- cbind(Population, result)
      }
      else {
        r <- utils.hwe(i, prob = p)
        Population <- rep(names(poplist)[count], nrow(r))
        r <- cbind(Population, r)
        result <- rbind(result, r)
      }
    }
  }
  rprob <- as.numeric(as.character(result$Prob))
  result <- result[(rprob > 0 & rprob <= p), ]
  result <- result[order(result$Locus), ]
  cat("Reporting significant departures from Hardy-Weinberg Equilibrium\n")
  if (nrow(result) == 0) {
    cat("No significant departures\n")
  }
  else {
    cat("NB: Departures significant at the alpha level of", 
        p, "are listed\n")
    if (p > 0.05) {
      cat("ns --", p, "< p < 0.05; * -- 0.05 < p < 0.01; ** -- 0.01 < p < 0.001; etc\n")
    }
    else {
      cat("ns -- p > 0.05; * -- 0.05 < p < 0.01; ** -- 0.01 < p < 0.001; etc\n")
    }
    cat("Critical values for significance of Bonferroni Corrected significance vary with sample size\n\n")
    print(result, row.names = FALSE)
  }
  silent(result)
}