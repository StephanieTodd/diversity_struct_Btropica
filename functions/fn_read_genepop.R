
# source("C:/Users/steph/OneDrive - James Cook University/Desktop (manual backup)/Rfunctions.R")
# genepop = "data/mapping_genpop_6213_110_HWE2.gen"
# popname = "data/popnames.txt"
read.GENEPOP_edited <- function (genepop, popname = NULL) {
  all_lines <- scan(genepop, what = character(), quiet = TRUE, 
                    sep = "\n",blank.lines.skip = FALSE)  # 
  gp_title <- all_lines[1]
  cline <- gsub(" ", "", all_lines)
  cline <- gsub("\t", "", cline)
  poploc <- which(toupper(cline) == "POP")
  MarkerList <- cline[2:(poploc[1] - 1)] # equivalent to  
  # MarkerList <- gsub(",", "\b", MarkerList)
  # MarkerList <- unlist(strsplit(MarkerList, "\n"))
  #MarkerList <- unlist(strsplit(MarkerList, ",")[[1]])
  MarkerList <- unlist(strsplit(MarkerList, ","))
  MarkerList <- unique(MarkerList) # no idea why but was 2x markers
  
  MarkerList <- MarkerList[! (MarkerList %in% "52395188-7-A/G16016638-52-C/A")]
  
  #write_csv(MarkerList)
  numMarker <- length(MarkerList)
  numPop   <- length(poploc)
  popstart <- poploc + 1
  popend <- c(poploc[-1] - 1, length(all_lines))
  numInd <- popend - popstart + 1
  numIndAll <- sum(numInd)
  PopID <- rep(1:numPop, numInd)
  rm(cline)
  gc()
  if (is.null(popname)) {
    pop.names <- paste0("pop", 1:numPop)
  } else {
    pop.names <- scan(popname, what = character(), quiet = TRUE, 
                      blank.lines.skip = TRUE)
  }
  gtdata <- all_lines[-poploc]  # remove "POP" lines
  gtdata <- gtdata[-(1:(poploc[1] - 1))]  # remove locus names line and header
  gtdata <- unlist(strsplit(gtdata, ",")) # commas separate individuals from their genotypes
  # each individual's genotype is one long string of alleles. 
  # gtdata <- unlist(strsplit(gtdata, " "))

  IndID <- gtdata[c(TRUE, FALSE)]  # ind id is every second item
  
  # IndID <- gsub(" ", "", IndID)  # dont need as my ind ids dont have spaces or tabs
  # IndID <- gsub("\t", "", IndID)

  # gtdata <- gsub(" ", "\t", gtdata[c(FALSE, TRUE)]) # pastes "\t" instead of tab - dont need anyway can split list on spaces
  
  gtdata <- gtdata[c(FALSE, TRUE)] #  # remove Ind IDs
  
  gtdata <- matrix(unlist(strsplit(gtdata, " ")), nrow = numIndAll, byrow = TRUE)
 
  gtdata <- gtdata[, -which(colSums(gtdata == "") != 0)]  ## # remove columns with zero genotypes (first col is blank)
 
  
  gp_digit <- as.integer(nchar(as.character(gtdata[1, 1]))/2)
  gp_na <- paste(rep("0", gp_digit), collapse = "")
  rm(all_lines)
  gc()
  htdata1 <- substr(gtdata, 1, gp_digit)
  htdata2 <- substr(gtdata, gp_digit + 1, gp_digit * 2)
  htdata1 <- gsub(gp_na, NA, htdata1)
  htdata2 <- gsub(gp_na, NA, htdata2)
  Haplotype <- list()
  Ho <- rep(0, numPop)
  for (cpop in 1:numPop) {   # cpop = 'current' pop
    cpopind <- PopID == cpop  
    cnpop <- numInd[cpop]  # number of individuals in pop
    
    chaplo <- array(data = "",  dim = c(cnpop, numMarker, 2))  # haplo = haploid??  # numMarker =1 could be issue
    h1 <- htdata1[cpopind, , drop = FALSE]
    h2 <- htdata2[cpopind, , drop = FALSE]
    # cnpop*numMarker =454200
    # 454200 / 127155
    chaplo[, , 1] <- h1 # this is the memory error!!
    chaplo[, , 2] <- h2  
    Ho[cpop] <- mean(colMeans(h1 != h2, na.rm = TRUE), na.rm = TRUE)
    Haplotype[[cpop]] <- chaplo
  }
  call <- !is.na(htdata1)
  CallRate.loci <- colMeans(call)
  CallRate.ind <- list()
  call.ind <- rowMeans(call)
  for (cp in 1:numPop) {
    CallRate.ind[[cp]] <- call.ind[PopID == cp]
  }
  htdata <- rbind(htdata1, htdata2)
  rm(htdata1, htdata2, gtdata)
  gc()
  AlleleCount <- list()
  AlleleFreq <- list()
  IndObs <- list()
  numAlleles <- rep(0, numMarker)
  AlleleList <- list()
  for (cm in 1:numMarker) {
    cgt <- table(htdata[, cm], c(PopID, PopID), useNA = "no")
    colnames(cgt) <- NULL
    numAlleles[cm] <- nrow(cgt)
    AlleleList[[cm]] <- rownames(cgt)
    cgtnum <- colSums(cgt)
    AlleleCount[[cm]] <- cgt
    numcall <- as.integer(cgtnum/2)
    IndObs[[cm]] <- numcall
    AlleleFreq[[cm]] <- t(t(cgt)/cgtnum)
  }
  He <- sapply(AlleleFreq, function(x) {
    1 - colSums(x^2)
  })
  He <- rowMeans(He, na.rm = TRUE)
  IndNames <- list()
  for (cpop in 1:numPop) {
    IndNames[[cpop]] <- IndID[PopID == cpop]
  }
  
  return(list(genotype = Haplotype, allele_count = AlleleCount,
              allele_freq = AlleleFreq, ind_count = IndObs, num_pop = numPop,
              pop_sizes = numInd, pop_names = pop.names, ind_names = IndNames,
              num_loci = numMarker, loci_names = MarkerList, num_allele = numAlleles,
              allele_list = AlleleList, call_rate_loci = CallRate.loci,
              call_rate_ind = CallRate.ind, He = He, Ho = Ho))
}


read.GENEPOP_edited2 <- function (genepop, popname = NULL) {
  all_lines <- scan(genepop, what = character(), quiet = TRUE, 
                    sep = "\n",blank.lines.skip = FALSE)  # 
  gp_title <- all_lines[1]
  cline <- gsub(" ", "", all_lines)
  cline <- gsub("\t", "", cline)
  poploc <- which(toupper(cline) == "POP")
  #MarkerList <- cline[2:(poploc[1] - 1)] # equivalent to  
  # MarkerList <- gsub(",", "\b", MarkerList)
  # MarkerList <- unlist(strsplit(MarkerList, "\n"))
  #MarkerList <- unlist(strsplit(MarkerList, ",")[[1]])
  MarkerList <- unlist(strsplit(MarkerList, ","))
  MarkerList <- unique(MarkerList) # no idea why but was 2x markers
  
  MarkerList <- MarkerList[! (MarkerList %in% "52395188-7-A/G16016638-52-C/A")]
  
  #write_csv(MarkerList)
  numMarker <- length(MarkerList)
  numPop   <- length(poploc)
  popstart <- poploc + 1
  popend <- c(poploc[-1] - 1, length(all_lines))
  numInd <- popend - popstart + 1
  numIndAll <- sum(numInd)
  PopID <- rep(1:numPop, numInd)
  rm(cline)
  gc()
  if (is.null(popname)) {
    pop.names <- paste0("pop", 1:numPop)
  } else {
    pop.names <- scan(popname, what = character(), quiet = TRUE, 
                      blank.lines.skip = TRUE)
  }
  gtdata <- all_lines[-poploc]  # remove "POP" lines
  gtdata <- gtdata[-(1:(poploc[1] - 1))]  # remove locus names line and header
  gtdata <- unlist(strsplit(gtdata, ",")) # commas separate individuals from their genotypes
  # each individual's genotype is one long string of alleles. 
  # gtdata <- unlist(strsplit(gtdata, " "))
  
  IndID <- gtdata[c(TRUE, FALSE)]  # ind id is every second item
  
  # IndID <- gsub(" ", "", IndID)  # dont need as my ind ids dont have spaces or tabs
  # IndID <- gsub("\t", "", IndID)
  
  # gtdata <- gsub(" ", "\t", gtdata[c(FALSE, TRUE)]) # pastes "\t" instead of tab - dont need anyway can split list on spaces
  
  gtdata <- gtdata[c(FALSE, TRUE)] #  # remove Ind IDs
  
  gtdata <- matrix(unlist(strsplit(gtdata, " ")), nrow = numIndAll, byrow = TRUE)
  
  gtdata <- gtdata[, -which(colSums(gtdata == "") != 0)]  ## # remove columns with zero genotypes (first col is blank)
  
  
  gp_digit <- as.integer(nchar(as.character(gtdata[1, 1]))/2)
  gp_na <- paste(rep("0", gp_digit), collapse = "")
  rm(all_lines)
  gc()
  htdata1 <- substr(gtdata, 1, gp_digit)
  htdata2 <- substr(gtdata, gp_digit + 1, gp_digit * 2)
  htdata1 <- gsub(gp_na, NA, htdata1)
  htdata2 <- gsub(gp_na, NA, htdata2)
  Haplotype <- list()
  Ho <- rep(0, numPop)
  for (cpop in 1:numPop) {   # cpop = 'current' pop
    cpopind <- PopID == cpop  
    cnpop <- numInd[cpop]  # number of individuals in pop
    
    chaplo <- array(data = "",  dim = c(cnpop, numMarker, 2))  # haplo = haploid??  # numMarker =1 could be issue
    h1 <- htdata1[cpopind, , drop = FALSE]
    h2 <- htdata2[cpopind, , drop = FALSE]
    # cnpop*numMarker =454200
    # 454200 / 127155
    chaplo[, , 1] <- h1 # this is the memory error!!
    chaplo[, , 2] <- h2  
    Ho[cpop] <- mean(colMeans(h1 != h2, na.rm = TRUE), na.rm = TRUE)
    Haplotype[[cpop]] <- chaplo
  }
  call <- !is.na(htdata1)
  CallRate.loci <- colMeans(call)
  CallRate.ind <- list()
  call.ind <- rowMeans(call)
  for (cp in 1:numPop) {
    CallRate.ind[[cp]] <- call.ind[PopID == cp]
  }
  htdata <- rbind(htdata1, htdata2)
  rm(htdata1, htdata2, gtdata)
  gc()
  AlleleCount <- list()
  AlleleFreq <- list()
  IndObs <- list()
  numAlleles <- rep(0, numMarker)
  AlleleList <- list()
  for (cm in 1:numMarker) {
    cgt <- table(htdata[, cm], c(PopID, PopID), useNA = "no")
    colnames(cgt) <- NULL
    numAlleles[cm] <- nrow(cgt)
    AlleleList[[cm]] <- rownames(cgt)
    cgtnum <- colSums(cgt)
    AlleleCount[[cm]] <- cgt
    numcall <- as.integer(cgtnum/2)
    IndObs[[cm]] <- numcall
    AlleleFreq[[cm]] <- t(t(cgt)/cgtnum)
  }
  He <- sapply(AlleleFreq, function(x) {
    1 - colSums(x^2)
  })
  He <- rowMeans(He, na.rm = TRUE)
  IndNames <- list()
  for (cpop in 1:numPop) {
    IndNames[[cpop]] <- IndID[PopID == cpop]
  }
  
  return(list(genotype = Haplotype, allele_count = AlleleCount,
              allele_freq = AlleleFreq, ind_count = IndObs, numPop = numPop,
              pop_sizes = numInd, pop_names = pop.names, ind_names = IndNames,
              numMarker = numMarker, loci_names = MarkerList, numAllele = numAlleles,
              allele_list = AlleleList, call_rate_loci = CallRate.loci,
              call_rate_ind = CallRate.ind, He = He, Ho = Ho))
}
#numPop, numMarker, numAllele
  