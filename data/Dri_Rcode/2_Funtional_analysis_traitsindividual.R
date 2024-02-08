#### Code that follows the paper Mixed-species bird flocks enhance the benefits of group aggregation by minimizing variation in some functional traits while maximizing variation in others by GF Dri, NC Caceres, J Della-Flora, CS Dambros
#### Written by GF Dri and CS Dambros 
#### Contact gabriela.franzoi at maine.edu 
#### November 2021

## This code performs functional richness for each trait individually (OBS FRic) at both local(flock) and regional (site) scales 

library(FD)

# Import functional traits table
# Row: species; Column: traits
trait_sp<-read.csv("Data/SpeciesTraitsOK.csv", row.names = 1)
head(trait_sp)

# Import species occurrence table
# Row: communities; Column: species
occ_sp<-read.csv("Data/Occurrence.csv", row.names = 1)
#occ_sp[rowSums(occ_sp)>2,] # just selecting communities with more than 2 species

##### LOCAL ####
##### wings ####
indicesfd_wing<-dbFD(trait_sp[,"wings",drop=FALSE], occ_sp[rowSums(occ_sp)>2,], w.abun = FALSE, stand.x = FALSE, corr = "none", calc.FRic = TRUE, m = "min", stand.FRic = TRUE, scale.RaoQ = FALSE, calc.CWM = FALSE, CWM.type = "all", calc.FDiv = TRUE, print.pco = FALSE, messages = TRUE)

#### beak ####
indicesfd_beak<-dbFD(trait_sp[,"beak_volume",drop=FALSE], occ_sp[rowSums(occ_sp)>2,], w.abun = FALSE, stand.x = FALSE, corr = "none", calc.FRic = TRUE, m = "min", stand.FRic = TRUE, scale.RaoQ = FALSE, calc.CWM = FALSE, CWM.type = "all", calc.FDiv = TRUE, print.pco = FALSE, messages = TRUE)

#### tail ####
indicesfd_retrizes<-dbFD(trait_sp[,"rectrizes",drop=FALSE], occ_sp[rowSums(occ_sp)>2,], w.abun = FALSE, stand.x = FALSE, corr = "none", calc.FRic = TRUE, m = "min", stand.FRic = TRUE, scale.RaoQ = FALSE, calc.CWM = FALSE, CWM.type = "all", calc.FDiv = TRUE, print.pco = FALSE, messages = TRUE)

#### tarsus ####
indicesfd_tarsus<-dbFD(trait_sp[,"rectrizes",drop=FALSE], occ_sp[rowSums(occ_sp)>2,], w.abun = FALSE, stand.x = FALSE, corr = "none", calc.FRic = TRUE, m = "min", stand.FRic = TRUE, scale.RaoQ = FALSE, calc.CWM = FALSE, CWM.type = "all", calc.FDiv = TRUE, print.pco = FALSE, messages = TRUE)

#### diet ####
indicesfd_diet<-dbFD(trait_sp[,"Inset.Plant",drop=FALSE], occ_sp[rowSums(occ_sp)>2,], w.abun = FALSE, stand.x = FALSE, corr = "none", calc.FRic = TRUE, m = "min", stand.FRic = TRUE, scale.RaoQ = FALSE, calc.CWM = FALSE, CWM.type = "all", calc.FDiv = TRUE, print.pco = FALSE, messages = TRUE)

#### foragign strata ####
indicesfd_est<-dbFD(trait_sp[,"Dossel",drop=FALSE], occ_sp[rowSums(occ_sp)>2,], w.abun = FALSE, stand.x = FALSE, corr = "none", calc.FRic = TRUE, m = "min", stand.FRic = TRUE, scale.RaoQ = FALSE, calc.CWM = FALSE, CWM.type = "all", calc.FDiv = TRUE, print.pco = FALSE, messages = TRUE)

#### body mass #####
indicesfd_bmass<-dbFD(trait_sp[,"BodyMass.Value",drop=FALSE], occ_sp[rowSums(occ_sp)>2,], w.abun = FALSE, stand.x = FALSE, corr = "none", calc.FRic = TRUE, m = "min", stand.FRic = TRUE, scale.RaoQ = FALSE, calc.CWM = FALSE, CWM.type = "all", calc.FDiv = TRUE, print.pco = FALSE, messages = TRUE)

#### REGIONAL ####
#### wing ####
indicesfd_site_wing<-dbFD(trait_sp[,"wings",drop=FALSE], occ_sp_site, w.abun = FALSE, stand.x = FALSE, corr = "none", calc.FRic = TRUE, m = "min", stand.FRic = TRUE, scale.RaoQ = FALSE, calc.CWM = FALSE, CWM.type = "all", calc.FDiv = TRUE, print.pco = FALSE, messages = TRUE)

#### beak ####
indicesfd_site_beak<-dbFD(trait_sp[,"beak_volume",drop=FALSE], occ_sp_site, w.abun = FALSE, stand.x = FALSE, corr = "none", calc.FRic = TRUE, m = "min", stand.FRic = TRUE, scale.RaoQ = FALSE, calc.CWM = FALSE, CWM.type = "all", calc.FDiv = TRUE, print.pco = FALSE, messages = TRUE)

#### tail ####
indicesfd_site_retrizes<-dbFD(trait_sp[,"rectrizes",drop=FALSE], occ_sp_site, w.abun = FALSE, stand.x = FALSE, corr = "none", calc.FRic = TRUE, m = "min", stand.FRic = TRUE, scale.RaoQ = FALSE, calc.CWM = FALSE, CWM.type = "all", calc.FDiv = TRUE, print.pco = FALSE, messages = TRUE)

#### tarsus ####
indicesfd_site_tarsus<-dbFD(trait_sp[,"rectrizes",drop=FALSE], occ_sp_site, w.abun = FALSE, stand.x = FALSE, corr = "none", calc.FRic = TRUE, m = "min", stand.FRic = TRUE, scale.RaoQ = FALSE, calc.CWM = FALSE, CWM.type = "all", calc.FDiv = TRUE, print.pco = FALSE, messages = TRUE)

#### diet ####
indicesfd_site_diet<-dbFD(trait_sp[,"Inset.Plant",drop=FALSE], occ_sp_site, w.abun = FALSE, stand.x = FALSE, corr = "none", calc.FRic = TRUE, m = "min", stand.FRic = TRUE, scale.RaoQ = FALSE, calc.CWM = FALSE, CWM.type = "all", calc.FDiv = TRUE, print.pco = FALSE, messages = TRUE)

#### foraging strata ####
indicesfd_site_est<-dbFD(trait_sp[,"Dossel",drop=FALSE], occ_sp_site, w.abun = FALSE, stand.x = FALSE, corr = "none", calc.FRic = TRUE, m = "min", stand.FRic = TRUE, scale.RaoQ = FALSE, calc.CWM = FALSE, CWM.type = "all", calc.FDiv = TRUE, print.pco = FALSE, messages = TRUE)

#### body mass ####
indicesfd_site_bmass<-dbFD(trait_sp[,"BodyMass.Value",drop=FALSE], occ_sp_site, w.abun = FALSE, stand.x = FALSE, corr = "none", calc.FRic = TRUE, m = "min", stand.FRic = TRUE, scale.RaoQ = FALSE, calc.CWM = FALSE, CWM.type = "all", calc.FDiv = TRUE, print.pco = FALSE, messages = TRUE)

#
## End