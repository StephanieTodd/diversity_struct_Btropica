
# the controller runs the code with different datasets (datcodes)

datcodestbl <- read.table(file = "datcodes.txt", header = T)

## Define version: 4QC2_[IND_CR][LOC_CR]_[READS_MIN][REP_AVG]_[MAF_GLOB]_reassnR -----------------------------
# datcode = "4QC2_8050_1090_reassnR" # <<<<<< only need to change this
# datcode = "4QC2_8050_1090_01" # basic + global MAF filt
#datcode = "4QC2_8010_1090" # 100% filt
# datcode = "4QC2_8050_1090" # basic
# datcode = "4QC3_8050_1090" # basic, 'QC3' = including secondaries

# datcodestbl[5,2]


#datcode = "6213_110"
#datcode = "6213_allinds"

datcode = "6213_110_HWE"
dat_descr = "new data same QC as 100pc filt, duplicates removed, HWE BonSig loci removed"
rmarkdown::render("dartR_notebook2.Rmd", "html_document")

file.copy(from = "dartR_notebook2.html", to = paste("results_writing", paste0("2_dartR_", datcode, ".html" ), sep ="/"), copy.date = TRUE, overwrite = TRUE)




# descriptions <- c("old samples 100pc_filt", "new data same QC as 100pc filt")
# versions_to_run <- c("4QC2_8010_1090", "6213_110")

descriptions <- c("new data same QC as 100pc filt, duplicates removed")
versions_to_run <- c("6213_110")


# 
# datcode = versions_to_run[1]
# dat_descr = descriptions[1]


for (i in 1:length(versions_to_run)) {
  datcode = versions_to_run[i]
  dat_descr = descriptions[i]
  cat("n\ running with ", versions_to_run[i], ".......")
  
  rmarkdown::render("dartR_notebook2.Rmd", "html_document")
  
  file.copy(from = "dartR_notebook2.html", to = paste("results_writing", paste0("dartR_", datcode, ".html" ), sep ="/"), copy.date = TRUE, overwrite = TRUE)
  
}

# save notebook html to Rdata folder

# rstudioapi::documentSave(rstudioapi::getActiveDocumentContext()$id)


file.copy(from = "QC.nb.html", to = paste(Rflder, paste0("QCrmd_", datcode, ".nb.html" ), sep ="/"), copy.date = FALSE, overwrite = T)

popsgls_ldav <- seppop(mappinggl_ldav)
popsgls_ldav
