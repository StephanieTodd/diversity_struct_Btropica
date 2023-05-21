# 19/4/2023
# copied from QC2.Rmd because issues with working directory of .Rmd

### load packages
{
  library(adegenet)
  # library(PopGenReport)  
  # library(boot)
  library(dartR)
  library(StAMPP)
  library(ggplot2)
  # library(plyr) causes problems!
  library(dplyr)
  library(igraph)
  library(tidyr)
  library(reshape2)
  library(SOAR)
}

source("functions/Rfunctions.R")

# QC SETTINGS ---------
descr    = ""
filt.2nds = T    # QC2 filters secondaries, QC3 includes secondaries, 
IND_CR = 0.8     # 80% highest threshold without losing any of the Spurgeon individuals
LOC_CR = 1       # CR = "call rate" threshold
READS_MIN = 10   # minimum read depth 
READS_MAX = 1000 # maximum read depth (find out why too many reads is a problem?)
REP_AVG   = 0.9  # reproducibility
MAF_GLOB = 0     # global minor allele freq filter
filt.HWE = T     # remove loci out of HWE (across all pops), with a Bonferroni correction 
recalc = T       # recalculate 

# 6213 is the new order number
# 110 = 100% loc filt and min 10 reads
# HWE = filtered Bonsig loci
#datcode = "6213_110"
datcode = "6213_110_HWE"  # only need to change string here
# "6213_allinds"  # only need to change string here

Vsfx        <- paste0("_", datcode)    # to ID any csvs exported to global output folder
finalname   <- paste0("mappinggl", datcode) # "mappinggl4QC2_30809095"
Rflder      <- paste0("Rdata", datcode)

dir.create(Rflder)

# save settings
QC_settings <- list("datcode" =datcode,
                    "filt.2nds" = filt.2nds, 
                    "IND_CR" = IND_CR,
                    "LOC_CR" = LOC_CR,          # CR = "call rate" threshold
                    "READS_MIN" = READS_MIN,
                    "READS_MAX" = READS_MAX, # find out why too many reads is a problem?
                    "REP_AVG"   = REP_AVG,
                    "MAF_GLOB" = REP_AVG,     # global minor allele freq filter
                    "filt.HWE" = filt.HWE,
                    "recalc" = recalc)


# Read in data -----------

# Following manual changes to original Dart dataset were made (required for it to import without errors)
# - find & replace ind.names " " with "_"
# - highlight duplicate names (two) and added suffixes to 2nd
# T142_453371_ETC2_02/11/15 >> "_2" 
# 425155_SPC06_22/11/18  >> "_rep"

data_file <- "data/Report_DBet21-6213_2_moreOrders_SNP_mapping_2.csv" # SNP mapping (one line per genotype)
ind_file <-  "data/individual_metadata_6213_v3.csv" # indiviudal metadata v2 has duplicates sorted, v3 resolves individual that was swapped from Tin to Emu 
# relative path is broken for some reason

QC_settings["data_file"] <- data_file
QC_settings["ind_file"]  <- ind_file

mapping.gl <- gl.read.dart(filename = data_file, ind.metafile = ind_file)

mapping.gl.import <- mapping.gl  # mapping.gl is overwritten each time...

#mapping.gl <- mapping.gl.import


# REMOVE DUPLICATE INDIVIDUALS ----------

# # rename id to id2
# mapping.gl <- gl.recode.ind(mapping.gl, ind.recode = "data/recode_ind2.csv")
ind_meta <- read.csv("data/individual_metadata_6213_v3.csv")
inds2keep <- ind_meta %>% 
  filter(KeepInd)
inds2keep <- inds2keep$id
length(inds2keep)

mapping.gl <- gl.keep.ind(mapping.gl, inds2keep,
                          recalc = TRUE, mono.rm = TRUE)


# NB mappinggl623_110 has "390507_DCB_H4_08/08/04" and "2A4D27_DCB_B6_06/08/04" but mappinggl623_110_HWE doesnt 
# (im not sure why they werent filtered in the first place)


## add pre-filtering ind call rate as a column to ind_metadata
# 
# ind.call.rate_prefilt <- 1 - rowSums(is.na(as.matrix(mapping.gl)))/nLoc(mapping.gl)
# 
# ind.call.rate_prefilt <- as.data.frame(ind.call.rate_prefilt)
# ind.call.rate_prefilt$id <- row.names(ind.call.rate_prefilt)
# 
# ind_meta <- merge(ind_meta, ind.call.rate_prefilt, all = T)



# Filter Reproducibility ---------------

# gl.filter.reproducibility & gl.filter.RepAvg() are the same
gl.report.reproducibility(mapping.gl)
mapping.gl <- gl.filter.reproducibility(mapping.gl, 
                                        threshold = REP_AVG)
mapping.gl.RepAvg <- mapping.gl


# Filter Secondaries -----------------
gl.report.secondaries(mapping.gl)
if(filt.2nds){                    
  mapping.gl <- gl.filter.secondaries(mapping.gl, 
                                      method = "best") # "best" filters out loci after ordering the genlight object on based on repeatability, avgPIC in that order
}
mapping.gl         <- gl.recalc.metrics(mapping.gl, verbose =  0)  # need to reclac metric & create mapping.gl2 regardless of filtering 2ndaries
mapping.gl.filt2nd <- mapping.gl

# Read Depth
gl.report.rdepth(mapping.gl)
mapping.gl <- gl.filter.rdepth(mapping.gl, 
                               lower = READS_MIN, 
                               upper = READS_MAX)  # find out why max read depth should be cut off?
mapping.gl.reads <- mapping.gl

# Filter Call rate per indiviudal  ----------------
gl.report.callrate(mapping.gl, method="ind")


# # add post filtering (reproducibility, read depth & secondaries) ind call rate as a column to ind_metadata
# ind.call.rate_postfilt <- 1 - rowSums(is.na(as.matrix(mapping.gl)))/nLoc(mapping.gl)
# 
# ind.call.rate_postfilt <- as.data.frame(ind.call.rate_postfilt)
# ind.call.rate_postfilt$id <- row.names(ind.call.rate_postfilt)
# 
# ind_meta <- merge(ind_meta, ind.call.rate_postfilt, all = T)


## only T065_424474_ETD2_21/05/15 below threshold post filt!!
# write.csv(ind_meta, "data/individual_metadata_6213.csv")

# filter ind call rate
mapping.gl <- gl.filter.callrate(mapping.gl, 
                                 method="ind", 
                                 threshold = IND_CR)
if(recalc){
  mapping.gl <- gl.recalc.metrics(mapping.gl, verbose = 0)
}

mapping.gl.indfilt <- mapping.gl

# No inds removed!
  
# Filter Call rate per locus -------------
gl.report.callrate(mapping.gl, method="loc")
mapping.gl <- gl.filter.callrate(mapping.gl, 
                                 method = "loc",
                                 threshold = LOC_CR)
if(recalc){
  mapping.gl <- gl.recalc.metrics(mapping.gl, verbose = 0)
}

mapping.gl.locfilt <- mapping.gl
# mapping.gl <- mapping.gl.locfilt


# Filter MAF (minor allele freq) ----------
# 18,478 binary SNPs with MAF filt vs 19,162 without! (684 with MAF < 0.01)
# Rick et al. 2019  'To further improve genotyping quality for this study, loci
# with a minor allele frequency (MAF) less than 5%, minimum heterozygosity less
# than 5%, and a maximum heterozygosity greater than 95% were excluded'

# White et al. 2018 'filtered out loci with >25% missing data and minor allele
# frequencies of <0.05. Although removing loci with low minor allele frequencies
# prohibits tracing the loss of rare alleles, we believe this conservative step
# is necessary to avoid incorporating erroneously called SNPs.'

# mapping.gl <- gl.filter.maf(mapping.gl, 
#                             threshold = MAF_GLOB)
# if(recalc){
#      mapping.gl <- gl.recalc.metrics(mapping.gl, verbose = 0)
# }
# mapping.gl.maffilt <- mapping.gl

# causes number of SNPs to drop to 3829 (from 5,176) even though no MAF threshold set



# filter HWE -----------------------

# Rick et al (2019) 'filtered both datasets removing loci that were in linkage
# disequilibrium (LD) and did not conform to Hardy-Weinberg expectations
# following a Bonferroni correction'
if (filt.HWE) {
  # hweloc <- mapping_hwe_Fsh  # cheat if already saved from dartR_notebook
  
  # 21/4/23 this function incorrectly summarises npop 'the number of population where the same locus is significantly out of hwe.'
  # So now keep all loci, and summarise myself
  hweloc_df <- gl.report.hwe(mapping.gl, alpha_val = 0.05, multi_comp = T, multi_comp_method = 'bonferroni', 
                             sig_only = FALSE, subset = "each", method_sig =  'Exact', v = 0) # 21/4/23 sig_only = FALSE

  # add columns for npops_Sig and npops_BonSig
  npops_Sig <- hweloc_df %>% 
    filter(Sig == "sig") %>% 
    group_by(Locus) %>% 
    dplyr::summarise(npops_Sig = n())
      
  npops_BonSig <- hweloc_df %>% 
    filter(Sig.adj == "sig") %>% 
    group_by(Locus) %>% 
    dplyr::summarise(npops_BonSig = n())
  
  hweloc_df %<>% 
    left_join(npops_Sig) %>% 
    left_join(npops_BonSig) %>% 
    data.frame() 
    
  hweloc_df <- replace_nulls_multi(hweloc_df, vars = c("npops_Sig", "npops_BonSig"), replacement = 0)
  write_csv(hweloc_df, suffix = paste0("_prefilt_", datcode))

  
  # summarise all combinations of Sig and BonSig
  summary_HWE_BonSig_loci <- hweloc_df %>% 
    filter(npops_BonSig > 0) %>% 
    distinct(Locus, .keep_all = T) %>% 
    group_by(npops_BonSig, npops_Sig) %>% 
    dplyr::summarise(n_loci = n()) 
  write_csv(summary_HWE_BonSig_loci)
  
  # summarise npops BonSig
  hweloc_df %>% 
    filter(npops_BonSig > 0) %>% 
    distinct(Locus, .keep_all = T) %>% 
    group_by(npops_BonSig) %>% 
    dplyr::summarise(n_loci = n()) 
  
  
  # which pops are responsible for the loci that are BonSig in just one pop? 
  # were checking that they're not all from one pop, which is more likely to have a biological reason 
  single_pop_BonSig <- hweloc_df %>% 
    filter(npops_BonSig == 1, Sig.adj == "sig") %>% 
    group_by(Population) %>% 
    dplyr::summarise(n_loci = n()) 
  write_csv(single_pop_BonSig)
  

  # which pops pairs are responsible for the loci that are BonSig in two pops? 
  two_pop_BonSig <- hweloc_df %>% 
    filter(npops_BonSig == 2, Sig.adj == "sig") %>% 
    group_by(Locus) %>% 
    mutate(pop = 1:n()) %>% 
    pivot_wider(id_cols = Locus, values_from = Population, names_from = pop, names_prefix = "pop") %>% 
    group_by(pop1, pop2) %>% 
    summarise(n_loci = n())
  write_csv(two_pop_BonSig)
    
  
  
  hweloc <- c(unique(hweloc$Locus), "16021466-40-C/T" )  # for some reason this locus now shows up as BonSig in Davies
  QC_settings[["QC_BonSighweloc"]] <- hweloc
  write_csv_if(hweloc, outname = "QC_BonSigHWEloc_", flag = T, folder = Rflder, suffix = paste(Vsfx,Sys.Date(), sep = "_"))
  
  loci2keep <- locNames(mapping.gl)[!locNames(mapping.gl) %in% hweloc]
  
  mapping.gl <- gl.keep.loc(mapping.gl, loci2keep)
  
  # # # NB this is the same as:
  # mapping.gl <- gl.filter.hwe(mapping.gl,  n.pop.threshold = 1,  
  #                          method_sig = "Exact",
  #                          multi_comp = TRUE,
  #                          multi_comp_method = "bonferroni",
  #                          alpha_val = 0.05,
  #                          pvalue_type = "midp",
  #                          cc_val = 0.5,
  #                          min_sample_size = 5,
  #                          verbose = 0)
  # # except for locus  "16021466-40-C/T"
  

}



# Final checks and export ----------

# check for monomorphic loci
gl.report.monomorphs(mapping.gl) 
mappinggl <- mapping.gl # mappinggl is the generic name used in dartR_notebook 
print(mapping.gl)

## EXPORT TO Rdata folder
save(mappinggl, file = paste(Rflder,  "/", finalname, ".Rdata", sep = ""))

# NB: file name and object name can be different. File & path name referr to specific version, object name is generic
# save a copy of .nb.html

# rstudioapi::documentSave(rstudioapi::getActiveDocumentContext()$id)
save_Rdata(QC_settings, withdate = T)
# file.copy(from = "QC.nb.html", to = paste(Rflder, paste0("QCrmd_", datcode, ".nb.html" ), sep ="/"), copy.date = TRUE, overwrite = T)

