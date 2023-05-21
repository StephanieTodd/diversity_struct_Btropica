# Run this after Ne_Estimator
# replaces missing bits - opy down to fill columns one & two e.g. "1:3_18/11/16    75" 
# and reads in data tables
# datcode = "6213_110_HWE"
source("functions/fn_NeLD.R")
source("scripts/02_global_variables_dartR.R")
if(!exists("gls")){
  gls <- list("allinds" = "allinds","nomigrants" =  "nomigrants")
}

#popabv <- c("Sp", "Dav", "Emu", "Tin")


for (i in 1:length(gls)) {
  for (j in 1:length(popabv)) {
    
    l <- ifelse(j == 1, 1, i)  # dont need need to run with no migrants for Spurgeon 
    # if (j  %in% 1:2) {  # I dont know why this doesnt work??
    #   i=1
    # }
    # 
    
    ## do with this file
    LDxLD_file <- paste0("Ne_Estimator/output/", popabv[j], "_", datcode, "_",  names(gls)[l], "_genpopLDxLD.txt")
    
    # read in the results table
    LDxLD <- readLines( LDxLD_file) # , n = 19 
    # print(LDxLD)
    
    # fill in the blank rows
    for (k in c(18, 19)) {
      substr(LDxLD[k], start = 1, stop = 18) <- substr(LDxLD[17], start = 1, stop = 18)
    }
    
    # write the entire table back to the file
    # print(LDxLD)
    writeLines( LDxLD , LDxLD_file )
    
    # done!
    table <- read.table(LDxLD_file,
                        header = FALSE, sep = "", skip = 16, nrows = 3)
    
    names(table) <- c("Population",  "SampSize",  "CritValue", "WeightedHMean", "IndepAlleles", "r2",
                      "Exp(r2)",   "Ne", "ParametricCI_min", "ParametricCI_max", "JackknifeCI_min",
                      "JackknifeCI_max", "Eff.df")
    table$Population <- Population[j]
    Ne_adj <- ChromAdj.Ne(as.matrix(table[,8:12]))
    Ne_adj <- as.data.frame(Ne_adj)
    names(Ne_adj) <- paste0(names(table[,8:12]), ".adj")
    
    table <- cbind.data.frame(table, Ne_adj)
    
    assign(paste0("Ne_estimator_", popabv[j], "_", names(gls)[l]), table)
    
  }
  cat("\n---------------------------\n")
}

# recombine results for populations into single df ------
Neobjs <- ls(pattern = "Ne_estimator_")
Netbls <- list()
for (i in 1:length(Neobjs)) {
  Netbls[[i]] <- get(Neobjs[i]) %>% 
    mutate(version = ifelse(grepl("allinds", Neobjs[i]), "all inds", "no migrants"))

}



Ne_estimator <- reduce(Netbls, rbind.data.frame)


names(Ne_estimator)[14] <- "Ne_adj"
write.csv(Ne_estimator, paste0("output/Ne_estimator", Vsfx, ".csv"))
rm(list = Neobjs)

cat("Ne estimated using the LD method in NeEstimator (Do et al. 2014), from loci with maf >0.02, and with the chromasome adjustment of Waples et al. 2016 applied: \n")
