library(dartR)

if (!exists("datcode")) {            
  datcode <- "6213_110_HWE"
  dat_descr <- "new data same QC as 100pc filt, duplicates removed"
}

Rflder    <-  paste0("Rdata", datcode)
datname   <-  paste0("mappinggl", datcode, ".Rdata")
Vsfx      <-  paste0("_", datcode)  # suffix for output csvs & mapping_gi


# Load data (from .Rdata files created in 'QC.R'): --------------------
# load(file = paste("../DArT-R_v2", Rflder, datname, sep ="/"))
load(file = paste(Rflder, datname, sep ="/"))


# datcodestbl <- read.table(file = "../DArT-R_v2/datcodes.txt", header = T)
# 
# cat("\n data QC version imported: \n", datcodestbl[as.character(datcodestbl$datcode) == datcode, "descr"], 
#     " (", datcode, ") \n\n")



# === DEFINE GLOBAL VARIABLES - part A =========================

# Mappinggl with Lower Davies pop
mappinggl_ldav     <- mappinggl
mappinggl_ldav$pop <- mappinggl$other$ind.metrics$pop2

# # COMBINED INTO Lamb Range
# mappingglLR <-  gl.recode.pop(mappinggl, pop.recode = "../DArT-R_v2/data/pop_recode_LR.csv")
# popsgls_lr  <- seppop(mappingglLR)

# error: head(mappinggl@other$loc.metrics.flags$monomorphs)


## define list of population names for indexing      
popabv      <- c("sp", "dav", "emu", "tin") # population abbreviations for naming objects
Population  <- c("Spurgeon", "Davies Creek", "Emu Creek", "Tinaroo")   # list of population names
Population2 <- c("Spurgeon", "Davies", "Emu", "Tinaroo")   # list of single word population names

# with lower davies
popabv_ldav  <- c("Sp", "LwDav", "UpDav", "Emu", "Tin")
Population_ldav  <- c("Spurgeon", "Lower Davies", "Upper Davies", "Emu Creek", "Tinaroo")   # list of population names
Population2_ldav <- c("Spurgeon", "LwrDavies", "UprDavies", "Emu", "Tinaroo")   # list of single word population names

cat("\n use '_ldav' suffixes for Davies Creek separated into Lower and Upper Davies populations \n")

# darkcyan
# deepskyblue4
#mediumpurple3
# indianred3
# salmon & salmon2
#mediumpurple4
# indianred1

# population colours for plotting 
popcol       <- c("sandybrown", "plum4", "darkseagreen2", "dodgerblue4")
popcol2      <- c("plum4",  "darkseagreen2", "sandybrown", "dodgerblue4") # for population order of mappinggl
popcol_ldav  <- c("sandybrown", "orchid", "orchid4", "darkseagreen2", "dodgerblue4")
popcol2_ldav <- c("orchid4",  "darkseagreen2", "orchid", "sandybrown", "dodgerblue4") # for population order of mappinggl "Davies Creek" "Emu Creek"    "Lower Davies" "Spurgeon"     "Tinaroo" 
popcol3_ldav <- c("darkseagreen2", "orchid",  "sandybrown", "dodgerblue4", "orchid4")  # for some reason ggplot puts them in a weird order 

names(popcol)  <- Population
names(popcol2) <- Population2[order(Population2)]
popcol3 <- popcol2
names(popcol3) <- Population[order(Population)]

poporder <- c(2,3,1,4)

# levels(mappinggl_ldav$pop)


# ===== SUBSET GENELIGHT INTO POPS =========================================

# using seppop
popsgls <- seppop(mappinggl)  # separate gl into list of pop gls
popsgls <- popsgls[c("Spurgeon", "Davies Creek", "Emu Creek", "Tinaroo")] # reorder so in correct pop order


popsgls_ldav <- seppop(mappinggl_ldav) 
popsgls_ldav <- popsgls_ldav[c("Spurgeon", "Lower Davies", "Upper Davies", "Emu Creek",  "Tinaroo")] # reorder so in correct pop order

names(popsgls_ldav)


# get lookup table of ind numbers and barcodes
indrename2nr_lookup2 <- get_indrename2nr_df()

write.csv(indrename2nr_lookup2, "output/indrename2nr_lookup2.csv")
