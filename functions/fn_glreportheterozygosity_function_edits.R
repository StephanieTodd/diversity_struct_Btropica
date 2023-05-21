# Original :
 function (gl)  {
  x <- gl
  sgl <- seppop(x)
  hs <- data.frame(lapply(sgl, function(x) mean(colMeans(as.matrix(x, 
                                                                   na.rm = TRUE) == 1), na.rm = TRUE)))  # output function function to 
  sums <- data.frame(lapply(sgl, function(x) mean(colSums(as.matrix(x, 
                                                                    na.rm = TRUE) == 1), na.rm = TRUE)))
  n <- t(sums/hs)
  n <- cbind(row.names(n), n)
  d1 <- data.frame(names = names(hs), hs = as.numeric(hs[1, 
                                                         ]))
  d2 <- data.frame(n)
  names(d2) <- c("names", "freq")
  dd <- join(d1, d2, by = "names")
  dd <- dd[order(dd$hs), ]
  par(mai = c(2.5, 1, 0.5, 0.2))
  barplot(dd$hs, names.arg = paste(dd$names, dd$freq, sep = " | "), 
          las = 2, cex.names = 1, space = 0, border = F, col = rainbow(nrow(dd)), 
          main = "Observed Heterozygosity")
  return(dd)
  
}

 # Edited:
 
gl.report.heterozygosity.edited <- function (gl) 
 {
   x <- gl
   sgl <- seppop(x)
   hs <- data.frame(lapply(sgl, function(x) mean(colMeans(as.matrix(x, 
                                                                    na.rm = TRUE) == 1, na.rm = TRUE), na.rm = TRUE))) # ADDED na.rm = TRUE as argfor colMeans()
   sums <- data.frame(lapply(sgl, function(x) mean(colSums(as.matrix(x, 
                                                                     na.rm = TRUE) == 1, na.rm = TRUE), na.rm = TRUE)))
   n <- t(sums/hs)
   n2 <- as.integer(n)
   n <- cbind(row.names(n), n)
   d1 <- data.frame(names = names(hs), hs = as.numeric(hs[1, 
                                                          ]))
   d2 <- data.frame(n)
   d2$X2 <- n2
   names(d2) <- c("names", "freq")
   dd <- join(d1, d2, by = "names")
   dd <- dd[order(dd$hs), ]
   par(mai = c(2.5, 1, 0.5, 0.2))
   barplot(dd$hs, names.arg = paste(dd$names, dd$freq, sep = " | "), 
           las = 2, cex.names = 1, space = 0, border = F, col = rainbow(nrow(dd)), 
           main = "Observed Heterozygosity")
   return(dd)
}



gl.report.heterozygosity.edited(mapping.gl)

# potential other way to get n
ntest <- t(data.frame(samplesize))
colnames(ntest) <- Population



# break it down

test <- as.matrix(mapping.gl, 
          na.rm = TRUE)  # NA.Rm doest really work in a matrix - you would have to have all alleles in and individual NA or all indivudals at a locus NA

test2 <- test ==1  # creates analouge matrix of TRUE where genotype == 1 (i.e. T/F for heterozygotes)

test3 <- colSums(test2, na.rm = TRUE) # calculates the sum of TRUE (i.e. counts total number of heterozygotes by locus) 
                        # loci that are either a) monomorphic in all pops or b) fixed in all pops have ZERO heterozygotes, and are listed as NA!!
                        # therefore by removing NAs when you take the mean you are removing loci!

test4 <- mean(test3, na.rm = TRUE) # average heterozygosity across loci

# mean < colsums < as.matrix < ==1
# 1. convert gl to matrix
# 2. matrix == 1
# 3. colSums(T/F matrix) or colMeans()
# 4. average heterozygosity (across loci)

# separate functions out from rest

    # testfn_colsums(gl)   ## the mean number of bettongs (sequenced or heterozygous) per locus
    # 2 changes posotion of , na.rm = TRUE
    testfn_colsums2 <- function(gl) {
      mean(colSums(as.matrix(gl) == 1, na.rm = TRUE), na.rm = TRUE)
    } 
    
    # testfn_colmeans(gl)  ## the overall mean heterozygosity (mean of  heterozygotes/total for each locus)
    testfn_colmeans2 <- function(gl) {
      mean(colMeans(as.matrix(gl) == 1, na.rm = TRUE), na.rm = TRUE)
    } 

testsgl <- seppop(mapping.gl)

hs2 <- data.frame(lapply(testsgl, testfn_colmeans2))
sums2 <- data.frame(lapply(testsgl, testfn_colsums2))

# plot:
  
  n <- t(sums2/hs2)  # t= transpose
  n2 <- as.integer(n)
  n <- cbind(row.names(n), n)
  
  
  d1 <- data.frame(names = names(hs2), hs2 = as.numeric(hs2[1, 
                                                         ]))
  d2 <- data.frame(n)
  d2$X2 <- n2
  names(d2) <- c("names", "freq")
  dd <- join(d1, d2, by = "names")
  dd <- dd[order(dd$hs2), ]
  par(mai = c(2.5, 1, 0.5, 0.2))
  barplot(dd$hs2, names.arg = paste(dd$names, dd$freq, sep = " | "), 
          las = 2, cex.names = 1, space = 0, border = F, col = rainbow(nrow(dd)), 
          main = "Observed Heterozygosity")
  
  return(dd)  # only needed for function??


# 3 without mean colmeans (returns dataframe)

testfn_colmeans3 <- function(gl) {
  colMeans(as.matrix(gl) == 1, na.rm = TRUE)
} 

hs3 <- data.frame(lapply(testsgl, testfn_colmeans3))

boxplot(hs3)

hs3_long <-reshape2:: melt(hs3)

colnames(hs3_long) <- c("Population", "Hs")


hs3_aov <- aov(Hs ~ Population, hs3_long)
hs3_Tukey <- TukeyHSD(hs3_aov)

# 4 with mean and se colmeans

testfn_colmeans4 <- function(gl) {
  meanandse(colMeans(as.matrix(gl) == 1, na.rm = TRUE))
} 

hs4 <- data.frame(lapply(testsgl, testfn_colmeans4))





