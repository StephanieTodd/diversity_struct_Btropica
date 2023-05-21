# ========================================================================================
# Ne_estimator software used separetly to estimate r2 and do most of Ne calculations --------------
# This file convertsand export genlight data to genepop format required by Ne Estimator software
# currently there is an issue with population specificationn in gi2genpop which means I need to run each pop separately
# ============================================================================================

library(dartR)
library(tidyverse)
library(magrittr)

source(file="functions/fn_Export_a_biallelic_SNP_Genind_object_to_a_Genepop_.R")
source("functions/fn_NeLD.R")

if(!exists("datcode")){
  datcode = "6213_110_HWE"
  Rflder    <-  paste0("Rdata", datcode)
  datname   <-  paste0("mappinggl", datcode, ".Rdata")
  Vsfx      <-  paste0("_", datcode)  # suffix for output csvs & mapping_gi
}

# Load data --------------------
load(file = paste(Rflder, datname, sep ="/"))
# load(file = "Rdata4QC2_8050_1090_reassnR/popsgls_reassn_4QC2_8050_1090_allloc.Rdata") # reassn broken after R update



# Optional - if removing migrants for 'local' vs 'm'etapop' Ne estimates --------------------

ancestrydat <- read.csv(paste0("../SambaR/v", datcode, "/SambaR_output/Structure/LEA_ancestry_K4.csv"))

# migrant classes:
# f0 <40%
# f1 40-60%
# f2 60-70%
# These arent really based on anything scientific

THRESH <- 0.6  # define migrant admixture threshold!

ancestrydat %<>% 
  mutate(capturepop_prob = case_when(capturepop == "Spurgeon" ~ Spurgeon_cluster,   # could be done more efficiently using pattern matching but im lazy...
                                     capturepop == "Davies" ~ Davies_cluster,
                                     capturepop == "Emu" ~ Emu_cluster,
                                     capturepop == "Tinaroo" ~ Tinaroo_cluster))

migrantsdat <- ancestrydat %>% 
  filter(capturepop_prob < THRESH) %>% 
  mutate(migrant_cls = case_when(capturepop_prob < 0.4 ~"f0",
                                 between(capturepop_prob, 0.4, 0.6) ~ "f1",
                                 capturepop_prob > 0.6 ~ "f2"))


# # from 'data/migrants.csv' (also see 'SambaR/gmatrix_3.xls')
# migrants <- c("T088_424412_TTA3_26/07/14",
#               "T183_423980_TTB4_23/11/16",
#               "424685_TTT6_03/04/16",
#               "T025_424564_TTB2_01/07/15",
#               "T022_453372_TTB7_31/08/15",
#               "T177_424640_TTH6_10/05/16",
#               "438103_TTB4_31/08/16",
#               "707051_TTH3_30/08/16",
#               "112761_TTB7_23/11/16",
#               "424753_ETD1_07/11/16",
#               "T072_423800_ETF3_20/02/15",
#               "424749_TTE7_01/08/16",
#               "T165_425294_TTE5_30/08/16",
#               "425146_TTH4_28/08/16",
#               "T162_424503_TTF2_28/08/16",
#               "425147_TTB5_28/08/16",
#               "T201_103932_DTG5_03/08/16",
#               "424385_TTF1_28/08/16")
# 
# 

# NEED TO SAVE/CARRY THROUGH MIGRANT THRESHOLD!!!!

migrants2rm <- migrantsdat %>% 
 .[["id"]]
 # c(.[["id"]], "Roadkill")  # doesnt work
  
Ne_migrants_summary <-  migrantsdat %>% 
    group_by(capturepop) %>% 
    summarise(count  =n(),
              threshold = THRESH) 
write_csv(Ne_migrants_summary)


mappinggl_migrantsrm <- gl.drop.ind(mappinggl, c(migrants2rm, "Roadkill"), recalc = T)


# # removemigrants # superceeded by above
# removemigrants <- function(x, ancestrydat){
#   migrants <- ancestrydat %>%
#     filter(Migrant == "f0"|
#              Migrant == "f1")
#   gl.drop.ind(x, ind.list = migrants$id, recalc = TRUE, mono.rm = FALSE, verbose = NULL)
# }
# 
# mappinggl_nomigrants <- removemigrants(mappinggl, ancestrydat)
# temp <- mappinggl
# mappinggl <- mappinggl_nomigrants

# convert data -----------------------------------------------------------
# ## converting mappinggl as a whole ## currently doesnt work in Ne Estimator
# mapping_gi <- gl2gi(mappinggl)
# genind2genepop(mapping_gi, file = "mapping_genpop_migrantsrm.gen")


# convert and export data -------------------
gls <- list("allinds" = mappinggl,"nomigrants" =  mappinggl_migrantsrm)
popabv <- c("Sp", "Dav", "Emu", "Tin")
for (i in 1:length(gls)) {
  cat("\n--------------------------\nrunning with ", names(gls)[i])
  # separate gl into pops
  popsgls <- seppop(gls[[i]])  
  popsgls <- popsgls[c("Spurgeon", "Davies Creek", "Emu Creek", "Tinaroo")] # reorder so in correct pop order
  print(popsgls)

  popgi_ls <-list()
  for (j in 1:length(popsgls)) {
    ## convert to genid
    popgi_ls[[j]] <- gl2gi(popsgls[[j]])  ###
  
    ## convert to genpop and save to file
    genind2genepop(popgi_ls[[j]], file = paste0("Ne_Estimator/data/", popabv[j], "_", datcode, "_",  names(gls)[i], "_genpop.gen"))
    
  }
  cat("\n---------------------------\n")
}


# -------------------------------------------------
# now run NeLR with the following paremeters:
# check 'output files in tabular format for top 3 crit values' box  to create 'LDxLD.txt' output.
# uncheck methods other than NeLD
# use critical values = 0.01,0.02,0.05
# uncheck 'also run without frequency restriction'
# >>> run
# then: open up LDxLD text files and copy down to fill columns one & two e.g. "1:3_18/11/16    75" 

# we set Pcrit in the program
# to screen out alleles at frequency ,0.02; Waples and Do
# (2010) found that this criterion provides a generally good
# balance between maximizing precision and minimizing bias.


# replace missing bits - opy down to fill columns one & two e.g. "1:3_18/11/16    75" 

source("scripts/04_NeLD_output_wrangling.RR")

