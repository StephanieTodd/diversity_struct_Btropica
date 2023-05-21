
### Wahuld effect - correlate Fst vs Fis 
# Checking if Heterozygosity defecit could be caused y Wahuld effect
# 
## notes added 19/4/23 in response to reviewers comments
# Fst requires populations to be defined, Fis does not
# previously I was testing for wahuld effects in a pooled population sample -
# which is effectively asking if there is structure between these pops
# a different approach is to calculate Fst between populations, and Fis within
# just one population. This is asking 'how much of population X's heterozygote
# excess is due to immigration from genetically differentiated pop Y?
library(magrittr)
library(tidyverse)
library(dartR)
library(SOAR)

# ------ Plot Fst vs Fis (all pops) -----
Fst_loci <- as.data.frame(Fst_loci)
plot(Fst_loci$Fst, Fst_loci$Fis)



# write.csv(Fst_loci, "output/Fst_loci.csv")

#--- calculate correlation (r2) between Fst and Fis ---
FstFis_lm <- lm(Fis ~ Fst, data = Fst_loci)
summary(FstFis_lm)$r.squared
ggplotRegression(FstFis_lm)

# # ----- just sig loci
# Fst_sigloci <- hwe_fst3[hwe_fst3$sig_bin == "sig",]
# plot(Fst_sigloci$Fst, Fst_sigloci$Fis.bypop)



# =====================================================================================
#### Subset tests for Wahuld Effect - Fst vs Fis correlation investigations
# ===========================================================================

#---- zoom into negative Fis area of plot with funny V pattern---

# negFis_loci <- Fst_loci[Fst_loci$Fis < 0,]
# plot(negFis_loci$Fst, negFis_loci$Fis)


#---- subset test for Wahuld effects between Emu and Tin-----
TinEmu_loci <- mapping_loci[mapping_loci$population == "Tinaroo" | mapping_loci$population == "Emu Creek",] ## rerun, used Daveis my mistake

# calculate Fst and save
if(!file.exists(file = paste(Rflder, "TinEmu_Fst.Rdata", sep = "/"))) {
  TinEmu_Fst <- Fst(TinEmu_loci) # takes a while..
  save(TinEmu_Fst, file = paste(Rflder, "TinEmu_Fst.Rdata", sep = "/"))   # so save to file...
}
load(file = paste(Rflder, "TinEmu_Fst.Rdata", sep = "/")) # and read back in

Tin_loci <- TinEmu_loci[TinEmu_loci$population == "Tinaroo", ]
Tin_Fis  <- Fst(Tin_loci)

Emu_loci <- TinEmu_loci[TinEmu_loci$population == "Emu Creek", ]
Emu_Fis  <- Fst(Emu_loci)

Tin_Fis %<>%  
  data.frame() %>% 
  mutate(locus = row.names(.)) %>% 
  dplyr::rename(Fis_Tin = Fis) %>% 
  dplyr::select(locus, Fis_Tin)
 
 
Emu_Fis %<>%  
  data.frame() %>% 
  mutate(locus = row.names(.)) %>% 
  dplyr::rename(Fis_Emu = Fis) %>% 
  dplyr::select(locus, Fis_Emu)

test <- TinEmu_Fst %>% 
  data.frame(.) %>% 
  dplyr::mutate(locus = row.names(.)) %>% 
  left_join(Tin_Fis) %>% 
  left_join(Emu_Fis)
    
test %<>% 
  filter(!is.na(Fst))

Emu_mod <- lm(Fis_Emu ~ Fst, data = test)


ggplotRegression(Emu_mod)


# convert to df and plot
TinEmu_Fst <- as.data.frame(TinEmu_Fst)

# calculate r2
FstFis_TinEmu_lm <- lm(Fst ~ Fis, data = TinEmu_Fst)
summary(FstFis_TinEmu_lm)$r.squared

ggplot(TinEmu_Fst, aes( x= Fst, y= Fis)) +
  

plot(Fis ~ Fst, TinEmu_Fst)
abline(FstFis_TinEmu_lm)

# --- subset test for Emu and Davies ---
DavEmu_loci <- mapping_loci[mapping_loci$population == "Davies Creek" | mapping_loci$population == "Emu Creek",]

# calculate Fst and save
if(!file.exists(file = paste(Rflder, "DavEmu_Fst.Rdata", sep = "/"))) {
  DavEmu_Fst <- Fst(DavEmu_loci) # takes a while..
  save(DavEmu_Fst, file = paste(Rflder, "DavEmu_Fst.Rdata", sep = "/"))   # so save to file...
}
load(file = paste(Rflder, "DavEmu_Fst.Rdata", sep = "/")) # and read back in

# convert to df and plot
DavEmu_Fst <- as.data.frame(DavEmu_Fst)

# calculate r2
FstFis_DavEmu_lm <- lm(Fst ~ Fis, data = DavEmu_Fst)
summary(FstFis_DavEmu_lm)$r.squared

plot(DavEmu_Fst$Fst, DavEmu_Fst$Fis)
abline(FstFis_DavEmu_lm)



# #---- subset test removing Spurgeon (LR only)-----
LR_loci <- mapping_loci[mapping_loci$population == "Davies Creek" | mapping_loci$population == "Emu Creek"| mapping_loci$population == "Tinaroo",]

# # calculate Fst and save
if(!file.exists(file = paste(Rflder, "LR_Fst.Rdata", sep = "/"))) {
  LR_Fst <- Fst(LR_loci) # takes a while..
  save(LR_Fst, file = paste(Rflder, "LR_Fst.Rdata", sep = "/"))   # so save to file...
}
load(file = paste(Rflder, "LR_Fst.Rdata", sep = "/")) # and read back in

# convert to df and plot
LR_Fst <- as.data.frame(LR_Fst)

# calculate r2
FstFis_LR_lm <- lm(Fst ~ Fis, data = LR_Fst)
summary(FstFis_LR_lm)$r.squared

plot(LR_Fst$Fst, LR_Fst$Fis, xlab = "Lamb Range Fst", ylab = "Lamb Range Fis")
abline(FstFis_LR_lm)

#---- subset test loci with higher Fst (> 0.2) -------
LR_Fst_pt2 <- LR_Fst[LR_Fst$Fst > 0.2,]
plot(LR_Fst_pt2$Fst, LR_Fst_pt2$Fis)
abline(FstFis_LR_lm)

# calculate r2
FstFis_pt2_LR_lm <- lm(Fst ~ Fis, data = LR_Fst_pt2)
summary(FstFis_pt2_LR_lm)$r.squared

# Section 2 ================
# De Meeûs 2018  method - using variance
Attach()
Fst_loci <- data.frame(Fst_loci)

# we want Fst between all 3 LR pops, and for each neighbour pair
LR_loci <- mapping_loci %>% 
  filter(population != "Spurgeon")

LR_loci_pool <- LR_loci %>% 
  dplyr::select( -population)

LR_loci_pool <- LR_loci %>% 
  mutate(population = factor("Lamb Range"))

EmuTin_loci <- mapping_loci %>% 
  filter(population == "Emu Creek" |population == "Tinaroo")

DavEmu_loci <- mapping_loci %>% 
  filter(population == "Emu Creek" |population == "Davies Creek")

# but we want Fis within each pop
Dav_loci <- mapping_loci %>% 
  filter(population == "Davies Creek")

Emu_loci <- mapping_loci %>% 
  filter(population == "Emu Creek")

Tin_loci <- mapping_loci %>% 
  filter(population == "Tinaroo")

loci_popls <- list("LR" = LR_loci, "LR-pool" = LR_loci_pool, "Emu-Tin" = EmuTin_loci, "Davies-Emu" = DavEmu_loci, "Dav" = Dav_loci, "Emu" = Emu_loci, "Tin" = Tin_loci)

# calculate F statistics
Fst_popls <- lapply(loci_popls, Fst)
# Fis_LR_pool <- Fst(LR_loci_pool)
# Fst_popls[[7]] <- Fis_LR_pool
# names(Fst_popls)[7] <- "LR-pool"


# convert to df and rename cols 
Fst_popls <- lapply(Fst_popls, data.frame)
for (i in 1:length(Fst_popls)) {
  Fst_popls[[i]] %<>% 
    rename_with(~ paste0(.x, "_", names(Fst_popls)[i])) %>% 
    mutate(locus = row.names(.))
  
}


Fst_Fis <- reduce(Fst_popls, left_join)

Fst_Fis  %<>%  
  remove_empty_cols(.) # %>% 
  filter(!is.na(Fst_LR))



# correlate
# calculate correlation (r2) between Fst and Fis ---
Emu_lm <- lm(Fis_Emu ~ Fst_LR, data = Fst_Fis)
summary(Emu_lm)$r.squared
ggplotRegression(Emu_lm)

Dav_lm <- lm(Fis_Davies ~ Fst_LR, data = Fst_Fis)
summary(Dav_lm)$r.squared
ggplotRegression(Dav_lm)

# calculate variance 
SE_Fst_LR <- se(Fst_Fis$Fst_LR)
SE_Fis_Emu <- se(Fst_Fis$Fis_Emu)

SE_Fis_Emu > SE_Fst_LR

SD_Fst_LR <- sd(Fst_Fis$Fst_LR)
SD_Fis_Emu <- sd(Fst_Fis$Fis_Emu, na.rm = T)

SD_Fis_Emu > SD_Fst_LR


SE_Fst_LR <- se(Fst_Fis$Fst_LR)
SE_Fis_Emu <- se(Fst_Fis$Fis_Emu)

# homozgote excess
mean(Fst_Fis$Fis_Davies, na.rm = T)

Fst_Fis %>%
  ggplot( aes(x=Fis_Davies)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) +
  ggtitle("Davies Creek distribution of locus homozygote excess (Fis)")

Fst_Fis %>%
  ggplot( aes(x=Fst_LR)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) +
  ggtitle("Lamb Range distribution of locus Fst")
