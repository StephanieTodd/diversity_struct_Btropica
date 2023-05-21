# Edited gl.report.heterozygosity:
gl <- mappinggl
gl.report.heterozygosity.edited <- function (gl) {
  x <- gl
  sgl <- seppop(x)
  sgl <- sgl[c(3,1,2,4)]  # reorder pops to correct order!
  hs <- data.frame(lapply(sgl, function(x) mean(colMeans(as.matrix(x) == 1, na.rm = TRUE), na.rm = TRUE))) # MOVED na.rm = TRUE to arg for colMeans() instead of as.matrix
  sums <- data.frame(lapply(sgl, function(x) mean(colSums(as.matrix(x) == 1, na.rm = TRUE), na.rm = TRUE))) # as above
  
  n <- t(sums/hs)
  n2 <- as.integer(n) # convert to whole indiviudals 
  n <- cbind(row.names(n), n)
  d1 <- data.frame(names = names(hs), hs = as.numeric(hs[1,]))
  d2 <- data.frame(n)
  d2$X2 <- n2 # overwrites floating point 
  names(d2) <- c("names", "freq")
  dd <- join(d1, d2, by = "names")
  # dd <- dd[order(dd$hs), ]
  par(mai = c(2.5, 1, 0.5, 0.2))
  barplot(dd$hs, names.arg = paste(dd$names, dd$freq, sep = " | "), 
          las = 2, cex.names = 1, space = 0, border = F, col = rainbow(nrow(dd)), 
          ylab = "Overall mean Observed Heterozygosity",
          ylim = c(0,0.2))
  return(dd)
}

gl.report.heterozygosity.edited(mappinggl)
