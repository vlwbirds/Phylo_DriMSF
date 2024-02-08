#### Code that follows the paper Mixed-species bird flocks enhance the benefits of group aggregation by minimizing variation in some functional traits while maximizing variation in others by GF Dri, NC Caceres, J Della-Flora, CS Dambros
#### Written by GF Dri and CS Dambros 
#### Contact gabriela.franzoi at maine.edu 
#### November 2021

## This code performs null functional richness for each trait individually (SES FRic)
## Local scale (flock) analysis

## You will need objects created in the previous codes (1_Functional_analysis_traitscombined and 2_Functional_analysis_traitsindividual)

##### FUNCTIONAL ANALYSIS #####
#### FRic Wing ####
indiceswingfun<-function(i){
  require(FD)

  trait_sp_rand<-trait_sp[sample(nrow(trait_sp)),]
  rownames(trait_sp_rand)<-rownames(trait_sp)
  
  wing_rand<-data.frame(trait_sp_rand[,1], row.names = row.names(trait_sp_rand))
  
  # FRic with randomized traits (SES FRic)
  null_indiceswing_pre<- dbFD(wing_rand, occ_sp[rowSums(occ_sp)>2,], w.abun=FALSE, stand.x = FALSE, corr="none", calc.FRic=TRUE, m="min",
                                   stand.FRic=TRUE, scale.RaoQ=FALSE, 
                                   calc.CWM=FALSE, calc.FDiv= TRUE, print.pco=FALSE, messages = FALSE
  )
  
  class(null_indiceswing_pre)
  
  null_indiceswing_vec<-unlist(null_indiceswing_pre)
  
  # Save results
  dir.create("NullResultIndiceWing",showWarnings = FALSE)
  
  writeLines(as.character(null_indiceswing_vec),paste0("NullResultIndiceWing/resuIndiceWingTEST",i,".txt"))
  
}

# paralleling the analysis to run the function 99 times 
library(snow)
{
  cl<-makeCluster(3)
  clusterExport(cl,list("trait_sp","occ_sp"))
  clusterApply(cl,1:99,indiceswingfun)
  stopCluster(cl)
} 

# importing back to R the randomized functional indices
null.files.wing<-list.files("NullResultIndiceWing/",full.names = TRUE)

check2.wing<-readLines(null.files.wing[1])

# we need to remove one row for each file, so lets find the row and remove it
toRemove.wing<-sum(c(length(indicesfd_wing$nbsp),
                         length(indicesfd_wing$sing.sp),
                         length(indicesfd_wing$FRic),
                         length(indicesfd_wing$qual.FRic)))

checknum2.wing<-as.numeric(check2.wing[-toRemove.wing]) # removing the row 

# organizing the files 
checknum3.wing<-lapply(null.files.wing,function(x){
  check2.wing<-readLines(x)
  checknum2.wing<-as.numeric(check2.wing[-toRemove.wing])
})

# tranforming it as an array (3 dimension matrix)
ver.wing<-do.call(c,checknum3.wing)
null.resu.wing<-array(ver.wing,dim=c(length(checknum3.wing[[1]])/6,6,length(checknum3.wing))) 

# observed trait FRic - average of all randomized (null) SES FRic trait
difobsnull.wing<-(indicesfd_wing$FRic-rowMeans(null.resu.wing[,3,],na.rm = TRUE))

# standard deviation of null FRic trait (SES FRic)
sdnull.wing<-apply(null.resu.wing[,3,],na.rm = TRUE,1,sd)

# Association between observed and expected FRic
summary(lm(indicesfd_wing$FRic~I(difobsnull.wing/sdnull.wing)))
plot(indicesfd_wing$FRic~I(difobsnull.wing/sdnull.wing)) # strong association

# association between SES FRic trait with group size (species richness of a flock)
summary(lm(I(difobsnull.wing/sdnull.wing)~null.resu.wing[,1,1])) # no evidence for the association

# linear mixed effect model to investigate the association between observed and expected FRic of the trait with flock size and enviromental covariates
# SES
summary(lmer(I(difobsnull.wing/sdnull.wing)~null.resu.wing[,1,1]+indicesfd_env$Temp+indicesfd_env$Prec+indicesfd_env$npp+indicesfd_env$Edge+indicesfd_env$Fitofisionomia+indicesfd_env$Ecoregion+(1|indicesfd_env$Localidade)))

# OBS
summary(lmer(indicesfd_wing$FRic~null.resu.wing[,1,1]+indicesfd_env$Temp+indicesfd_env$Prec+indicesfd_env$npp+indicesfd_env$Edge+indicesfd_env$Fitofisionomia+indicesfd_env$Ecoregion+(1|indicesfd_env$Localidade)))

#### FRic Beak ####
indicesbeakfun<-function(i){
  require(FD)

  trait_sp_rand<-trait_sp[sample(nrow(trait_sp)),]
  rownames(trait_sp_rand)<-rownames(trait_sp)
  
  beak_rand<-data.frame(trait_sp_rand[,4], row.names = row.names(trait_sp_rand))
  
  # FRic with randomized traits (SES FRic)
  null_indicesbeak_pre<- dbFD(
    beak_rand, occ_sp[rowSums(occ_sp)>2,], w.abun=FALSE, stand.x = FALSE, corr="none", calc.FRic=TRUE, m="min",
    stand.FRic=TRUE, scale.RaoQ=FALSE, 
    calc.CWM=FALSE, calc.FDiv= TRUE, print.pco=FALSE, messages = FALSE
  )
  
  class(null_indicesbeak_pre)
  
  null_indicesbeak_vec<-unlist(null_indicesbeak_pre)
  
  # Save results
  dir.create("NullResultIndiceBeak",showWarnings = FALSE)
  
  writeLines(as.character(null_indicesbeak_vec),paste0("NullResultIndiceBeak/resuIndiceBeak",i,".txt"))
  
}

# paralleling the analysis to run the function 99 times 
library(snow)
{
  cl<-makeCluster(3)
  clusterExport(cl,list("trait_sp","occ_sp"))
  clusterApply(cl,1:99,indicesbeakfun)
  stopCluster(cl)
} 

# importing back to R the randomized functional indices
null.files.beak<-list.files("NullResultIndiceBeak/",full.names = TRUE)

check2.beak<-readLines(null.files.beak[1])

# we need to remove one row for each file, so lets find the row and remove it
toRemove.beak<-sum(c(length(indicesfd_beak$nbsp),
                     length(indicesfd_beak$sing.sp),
                     length(indicesfd_beak$FRic),
                     length(indicesfd_beak$qual.FRic)))

checknum2.beak<-as.numeric(check2.beak[-toRemove.beak]) # removing the row 

# organizing the files 
checknum3.beak<-lapply(null.files.beak,function(x){
  check2.beak<-readLines(x)
  checknum2.beak<-as.numeric(check2.beak[-toRemove.beak])
})

# tranforming it as an array (3 dimension matrix)
ver.beak<-do.call(c,checknum3.beak)
null.resu.beak<-array(ver.beak,dim=c(length(checknum3.beak[[1]])/6,6,length(checknum3.beak))) 

# observed trait FRic - average of all randomized (null) SES FRic trait
difobsnull.beak<-(indicesfd_beak$FRic-rowMeans(null.resu.beak[,3,],na.rm = TRUE))

# standard deviation of null FRic trait (SES FRic)
sdnull.beak<-apply(null.resu.beak[,3,],na.rm = TRUE,1,sd)

# association between SES FRic trait with group size (species richness of a flock)
summary(lm(indicesfd_beak$FRic~I(difobsnull.beak/sdnull.beak)))
plot(indicesfd_beak$FRic~I(difobsnull.beak/sdnull.beak)) # strong association positiva

# association between SES FRic trait with group size (species richness of a flock)
summary(lm(I(difobsnull.beak/sdnull.beak)~null.resu.beak[,1,1])) # no evidence for the association

# SES
summary(lmer(I(difobsnull.beak/sdnull.beak)~null.resu.beak[,1,1]+indicesfd_env$Temp+indicesfd_env$Prec+indicesfd_env$npp+indicesfd_env$Edge+indicesfd_env$Fitofisionomia+indicesfd_env$Ecoregion+(1|indicesfd_env$Localidade)))

summary(lmer(indicesfd_beak$FRic~null.resu.beak[,1,1]+indicesfd_env$Temp+indicesfd_env$Prec+indicesfd_env$npp+indicesfd_env$Edge+indicesfd_env$Fitofisionomia+indicesfd_env$Ecoregion+(1|indicesfd_env$Localidade)))


#### FRic Tail ####
indicesretrizesfun<-function(i){
  require(FD)

  trait_sp_rand<-trait_sp[sample(nrow(trait_sp)),]
  rownames(trait_sp_rand)<-rownames(trait_sp)
  
  retrizes_rand<-data.frame(trait_sp_rand[,2], row.names = row.names(trait_sp_rand))
  
  # FRic with randomized traits (SES FRic)
  null_indicesretrizes_pre<- dbFD(
    retrizes_rand, occ_sp[rowSums(occ_sp)>2,], w.abun=FALSE, stand.x = FALSE, corr="none", calc.FRic=TRUE, m="min",
    stand.FRic=TRUE, scale.RaoQ=FALSE, 
    calc.CWM=FALSE, calc.FDiv= TRUE, print.pco=FALSE, messages = FALSE
  )
  
  class(null_indicesretrizes_pre)
  
  null_indicesretrizes_vec<-unlist(null_indicesretrizes_pre)
  
  # Save results
  dir.create("NullResultIndiceRetrizes",showWarnings = FALSE)
  
  writeLines(as.character(null_indicesretrizes_vec),paste0("NullResultIndiceRetrizes/resuIndiceRetrizes",i,".txt"))
  
}

# paralleling the analysis to run the function 99 times 
library(snow)
{
  cl<-makeCluster(3)
  clusterExport(cl,list("trait_sp","occ_sp"))
  clusterApply(cl,1:99,indicesretrizesfun)
  stopCluster(cl)
} 
# importing back to R the randomized functional indices
null.files.retrizes<-list.files("NullResultIndiceRetrizes/",full.names = TRUE)

check2.retrizes<-readLines(null.files.retrizes[1])

# we need to remove one row for each file, so lets find the row and remove it
toRemove.retrizes<-sum(c(length(indicesfd_retrizes$nbsp),
                         length(indicesfd_retrizes$sing.sp),
                         length(indicesfd_retrizes$FRic),
                         length(indicesfd_retrizes$qual.FRic)))

checknum2.retrizes<-as.numeric(check2.retrizes[-toRemove.retrizes]) # removing the row 

# organizing the files 
checknum3.retrizes<-lapply(null.files.retrizes,function(x){
  check2.retrizes<-readLines(x)
  checknum2.retrizes<-as.numeric(check2.retrizes[-toRemove.retrizes])
})

# tranforming it as an array (3 dimension matrix)
ver.retrizes<-do.call(c,checknum3.retrizes)
null.resu.retrizes<-array(ver.retrizes,dim=c(length(checknum3.retrizes[[1]])/6,6,length(checknum3.retrizes))) 

# observed trait FRic - average of all randomized (null) SES FRic trait
difobsnull.retrizes<-(indicesfd_retrizes$FRic-rowMeans(null.resu.retrizes[,3,],na.rm = TRUE))

# standard deviation of null FRic trait (SES FRic)
sdnull.retrizes<-apply(null.resu.retrizes[,3,],na.rm = TRUE,1,sd)

# association between FRic trait with group size (species richness of a flock)
summary(lm(indicesfd_retrizes$FRic~I(difobsnull.retrizes/sdnull.retrizes)))
plot(indicesfd_retrizes$FRic~I(difobsnull.retrizes/sdnull.retrizes)) # strong association

# association between SES FRic trait with group size (species richness of a flock)
summary(lm(I(difobsnull.retrizes/sdnull.retrizes)~null.resu.retrizes[,1,1])) # positive association

# SES
summary(lmer(I(difobsnull.retrizes/sdnull.retrizes)~null.resu.retrizes[,1,1]+indicesfd_env$Temp+indicesfd_env$Prec+indicesfd_env$npp+indicesfd_env$Edge+indicesfd_env$Fitofisionomia+indicesfd_env$Ecoregion+(1|indicesfd_env$Localidade)))

summary(lmer(indicesfd_retrizes$FRic~null.resu.retrizes[,1,1]+indicesfd_env$Temp+indicesfd_env$Prec+indicesfd_env$npp+indicesfd_env$Edge+indicesfd_env$Fitofisionomia+indicesfd_env$Ecoregion+(1|indicesfd_env$Localidade)))

#### FRic Tarsus ####
indicestarsusfun<-function(i){
  require(FD)

  trait_sp_rand<-trait_sp[sample(nrow(trait_sp)),]
  rownames(trait_sp_rand)<-rownames(trait_sp)
  
  tarsus_rand<-data.frame(trait_sp_rand[,3], row.names = row.names(trait_sp_rand))
  
  # FRic with randomized traits (SES FRic)
  null_indicestarsus_pre<- dbFD(
    tarsus_rand, occ_sp[rowSums(occ_sp)>2,], w.abun=FALSE, stand.x = FALSE, corr="none", calc.FRic=TRUE, m="min",
    stand.FRic=TRUE, scale.RaoQ=FALSE, 
    calc.CWM=FALSE, calc.FDiv= TRUE, print.pco=FALSE, messages = FALSE
  )
  
  class(null_indicestarsus_pre)
  
  null_indicestarsus_vec<-unlist(null_indicestarsus_pre)
  
  # Save results
  dir.create("NullResultIndiceTarsus",showWarnings = FALSE)
  
  writeLines(as.character(null_indicestarsus_vec),paste0("NullResultIndiceTarsus/resuIndiceTarsus",i,".txt"))
  
}

# paralleling the analysis to run the function 99 times 
library(snow)
{
  cl<-makeCluster(3)
  clusterExport(cl,list("trait_sp","occ_sp"))
  clusterApply(cl,1:99,indicestarsusfun)
  stopCluster(cl)
} 

# importing back to R the randomized functional indices
null.files.tarsus<-list.files("NullResultIndiceTarsus/",full.names = TRUE)

check2.tarsus<-readLines(null.files.tarsus[1])

# we need to remove one row for each file, so lets find the row and remove it
toRemove.tarsus<-sum(c(length(indicesfd_tarsus$nbsp),
                       length(indicesfd_tarsus$sing.sp),
                       length(indicesfd_tarsus$FRic),
                       length(indicesfd_tarsus$qual.FRic)))

checknum2.tarsus<-as.numeric(check2.tarsus[-toRemove.tarsus]) # removing the row 

# organizing the files 
checknum3.tarsus<-lapply(null.files.tarsus,function(x){
  check2.tarsus<-readLines(x)
  checknum2.tarsus<-as.numeric(check2.tarsus[-toRemove.tarsus])
})

# tranforming it as an array (3 dimension matrix)
ver.tarsus<-do.call(c,checknum3.tarsus)
null.resu.tarsus<-array(ver.tarsus,dim=c(length(checknum3.tarsus[[1]])/6,6,length(checknum3.tarsus))) 

# observed trait FRic - average of all randomized (null) SES FRic trait
difobsnull.tarsus<-(indicesfd_tarsus$FRic-rowMeans(null.resu.tarsus[,3,],na.rm = TRUE))

# standard deviation of null FRic trait (SES FRic)
sdnull.tarsus<-apply(null.resu.tarsus[,3,],na.rm = TRUE,1,sd)

# association between FRic trait with group size (species richness of a flock)
summary(lm(indicesfd_tarsus$FRic~I(difobsnull.tarsus/sdnull.tarsus)))
plot(indicesfd_tarsus$FRic~I(difobsnull.tarsus/sdnull.tarsus)) # strong association positiva

# association between SES FRic trait with group size (species richness of a flock)
summary(lm(I(difobsnull.tarsus/sdnull.tarsus)~null.resu.tarsus[,1,1])) # no evidence for the association

# SES
summary(lmer(I(difobsnull.tarsus/sdnull.tarsus)~null.resu.tarsus[,1,1]+indicesfd_env$Temp+indicesfd_env$Prec+indicesfd_env$npp+indicesfd_env$Edge+indicesfd_env$Fitofisionomia+indicesfd_env$Ecoregion+(1|indicesfd_env$Localidade)))

summary(lmer(indicesfd_tarsus$FRic~null.resu.tarsus[,1,1]+indicesfd_env$Temp+indicesfd_env$Prec+indicesfd_env$npp+indicesfd_env$Edge+indicesfd_env$Fitofisionomia+indicesfd_env$Ecoregion+(1|indicesfd_env$Localidade)))

#### FRic Dieta ####
indicesdietfun<-function(i){
  require(FD)

  trait_sp_rand<-trait_sp[sample(nrow(trait_sp)),]
  rownames(trait_sp_rand)<-rownames(trait_sp)
  

  # FRic with randomized traits (SES FRic)
  null_indicesdiet_pre<- dbFD(
    diet_rand, occ_sp[rowSums(occ_sp)>2,], w.abun=FALSE, stand.x = FALSE, corr="none", calc.FRic=TRUE, m="min",
    stand.FRic=TRUE, scale.RaoQ=FALSE, 
    calc.CWM=FALSE, calc.FDiv= TRUE, print.pco=FALSE, messages = FALSE
  )
  
  class(null_indicesdiet_pre)
  
  null_indicesdiet_vec<-unlist(null_indicesdiet_pre)
  
  # Save results
  dir.create("NullResultIndiceDiet",showWarnings = FALSE)
  
  writeLines(as.character(null_indicesdiet_vec),paste0("NullResultIndiceDiet/resuIndiceDiet",i,".txt"))
  
}

# paralleling the analysis to run the function 99 times 
library(snow)
{
  cl<-makeCluster(3)
  clusterExport(cl,list("trait_sp","occ_sp"))
  clusterApply(cl,1:99,indicesdietfun)
  stopCluster(cl)
} 

# importing back to R the randomized functional indices
null.files.diet<-list.files("NullResultIndiceDiet/",full.names = TRUE)

check2.diet<-readLines(null.files.diet[1])

# we need to remove one row for each file, so lets find the row and remove it
toRemove.diet<-sum(c(length(indicesfd_diet$nbsp),
                     length(indicesfd_diet$sing.sp),
                     length(indicesfd_diet$FRic),
                     length(indicesfd_diet$qual.FRic)))

checknum2.diet<-as.numeric(check2.diet[-toRemove.diet]) # removing the row 

# organizing the files 
checknum3.diet<-lapply(null.files.diet,function(x){
  check2.diet<-readLines(x)
  checknum2.diet<-as.numeric(check2.diet[-toRemove.diet])
})

# tranforming it as an array (3 dimension matrix)
ver.diet<-do.call(c,checknum3.diet)
null.resu.diet<-array(ver.diet,dim=c(length(checknum3.diet[[1]])/6,6,length(checknum3.diet))) 

# observed trait FRic - average of all randomized (null) SES FRic trait
difobsnull.diet<-(indicesfd_diet$FRic-rowMeans(null.resu.diet[,3,],na.rm = TRUE))

# standard deviation of null FRic trait (SES FRic)
sdnull.diet<-apply(null.resu.diet[,3,],na.rm = TRUE,1,sd)

# association between FRic trait with group size (species richness of a flock)
summary(lm(indicesfd_diet$FRic~I(difobsnull.diet/sdnull.diet)))
plot(indicesfd_diet$FRic~I(difobsnull.diet/sdnull.diet)) # strong association positiva

# association between SES FRic trait with group size (species richness of a flock)
summary(lm(I(difobsnull.diet/sdnull.diet)~null.resu.diet[,1,1])) # positive association

# SES
summary(lmer(sesfric_diet~null.resu.diet[,1,1]+indicesfd_env$Temp+indicesfd_env$Prec+indicesfd_env$npp+indicesfd_env$Edge+indicesfd_env$Fitofisionomia+indicesfd_env$Ecoregion+(1|indicesfd_env$Localidade)))

summary(lmer(indicesfd_diet$FRic~null.resu.diet[,1,1]+indicesfd_env$Temp+indicesfd_env$Prec+indicesfd_env$npp+indicesfd_env$Edge+indicesfd_env$Fitofisionomia+indicesfd_env$Ecoregion+(1|indicesfd_env$Localidade)))


#### FRic Foraging strata ####
indicesestfun<-function(i){
  require(FD)

  trait_sp_rand<-trait_sp[sample(nrow(trait_sp)),]
  rownames(trait_sp_rand)<-rownames(trait_sp)
  
  # FRic with randomized traits (SES FRic)
  null_indicesest_pre<- dbFD(
    est_rand, occ_sp[rowSums(occ_sp)>2,], w.abun=FALSE, stand.x = FALSE, corr="none", calc.FRic=TRUE, m="min",
    stand.FRic=TRUE, scale.RaoQ=FALSE, 
    calc.CWM=FALSE, calc.FDiv= TRUE, print.pco=FALSE, messages = FALSE
  )
  
  class(null_indicesest_pre)
  
  null_indicesest_vec<-unlist(null_indicesest_pre)
  
  # Save results
  dir.create("NullResultIndiceEst",showWarnings = FALSE)
  
  writeLines(as.character(null_indicesest_vec),paste0("NullResultIndiceEst/resuIndiceEst",i,".txt"))
  
}

# paralleling the analysis to run the function 99 times 
library(snow)
{
  cl<-makeCluster(3)
  clusterExport(cl,list("trait_sp","occ_sp"))
  clusterApply(cl,1:99,indicesestfun)
  stopCluster(cl)
} 

# importing back to R the randomized functional indices
null.files.est<-list.files("NullResultIndiceEst/",full.names = TRUE)

check2.est<-readLines(null.files.est[1])

# we need to remove one row for each file, so lets find the row and remove it
toRemove.est<-sum(c(length(indicesfd_est$nbsp),
                    length(indicesfd_est$sing.sp),
                    length(indicesfd_est$FRic),
                    length(indicesfd_est$qual.FRic)))

checknum2.est<-as.numeric(check2.est[-toRemove.est]) # removing the row 

# organizing the files 
checknum3.est<-lapply(null.files.est,function(x){
  check2.est<-readLines(x)
  checknum2.est<-as.numeric(check2.est[-toRemove.est])
})

# tranforming it as an array (3 dimension matrix)
ver.est<-do.call(c,checknum3.est)
null.resu.est<-array(ver.est,dim=c(length(checknum3.est[[1]])/6,6,length(checknum3.est))) 

# observed trait FRic - average of all randomized (null) SES FRic trait
difobsnull.est<-(indicesfd_est$FRic-rowMeans(null.resu.est[,3,],na.rm = TRUE))

# standard deviation of null FRic trait (SES FRic)
sdnull.est<-apply(null.resu.est[,3,],na.rm = TRUE,1,sd)

# association between FRic trait with group size (species richness of a flock)
summary(lm(indicesfd_est$FRic~I(difobsnull.est/sdnull.est)))
plot(indicesfd_est$FRic~I(difobsnull.est/sdnull.est)) # strong association positiva

# association between SES FRic trait with group size (species richness of a flock)
summary(lm(I(difobsnull.est/sdnull.est)~null.resu.est[,1,1])) # no evidence for the association

# SES
summary(lmer(I(difobsnull.est/sdnull.est)~null.resu.est[,1,1]+indicesfd_env$Temp+indicesfd_env$Prec+indicesfd_env$npp+indicesfd_env$Edge+indicesfd_env$Fitofisionomia+indicesfd_env$Ecoregion+(1|indicesfd_env$Localidade)))

summary(lmer(indicesfd_est$FRic~null.resu.est[,1,1]+indicesfd_env$Temp+indicesfd_env$Prec+indicesfd_env$npp+indicesfd_env$Edge+indicesfd_env$Fitofisionomia+indicesfd_env$Ecoregion+(1|indicesfd_env$Localidade)))

#### FRic Body Mass ####
indicesbmassfun<-function(i){
  require(FD)

  trait_sp_rand<-trait_sp[sample(nrow(trait_sp)),]
  rownames(trait_sp_rand)<-rownames(trait_sp)
  
  body_mass_rand<-data.frame(trait_sp_rand[,30], row.names = row.names(trait_sp_rand))
  
  # FRic with randomized traits (SES FRic)
  null_indicesbmass_pre<- dbFD(
    body_mass_rand, occ_sp[rowSums(occ_sp)>2,], w.abun=FALSE, stand.x = FALSE, corr="none", calc.FRic=TRUE, m="min",
    stand.FRic=TRUE, scale.RaoQ=FALSE, 
    calc.CWM=FALSE, calc.FDiv= TRUE, print.pco=FALSE, messages = FALSE
  )
  
  class(null_indicesbmass_pre)
  
  null_indicesbmass_vec<-unlist(null_indicesbmass_pre)
  
  # Save results
  dir.create("NullResultIndiceBMass",showWarnings = FALSE)
  
  writeLines(as.character(null_indicesbmass_vec),paste0("NullResultIndiceBMass/resuIndiceBMass",i,".txt"))
  
}

# paralleling the analysis to run the function 99 times 
library(snow)
{
  cl<-makeCluster(3)
  clusterExport(cl,list("trait_sp","occ_sp"))
  clusterApply(cl,1:99,indicesbmassfun)
  stopCluster(cl)
} 

# importing back to R the randomized functional indices
null.files.bmass<-list.files("NullResultIndiceBMass/",full.names = TRUE)

check2.bmass<-readLines(null.files.bmass[1])

# we need to remove one row for each file, so lets find the row and remove it
toRemove.bmass<-sum(c(length(indicesfd_bmass$nbsp),
                      length(indicesfd_bmass$sing.sp),
                      length(indicesfd_bmass$FRic),
                      length(indicesfd_bmass$qual.FRic)))

checknum2.bmass<-as.numeric(check2.bmass[-toRemove.bmass]) # removing the row 

# organizing the files 
checknum3.bmass<-lapply(null.files.bmass,function(x){
  check2.bmass<-readLines(x)
  checknum2.bmass<-as.numeric(check2.bmass[-toRemove.bmass])
})

# tranforming it as an array (3 dimension matrix)
ver.bmass<-do.call(c,checknum3.bmass)
null.resu.bmass<-array(ver.bmass,dim=c(length(checknum3.bmass[[1]])/6,6,length(checknum3.bmass))) 

# observed trait FRic - average of all randomized (null) SES FRic trait
difobsnull.bmass<-(indicesfd_bmass$FRic-rowMeans(null.resu.bmass[,3,],na.rm = TRUE))

# standard deviation of null FRic trait (SES FRic)
sdnull.bmass<-apply(null.resu.bmass[,3,],na.rm = TRUE,1,sd)

# association between SES FRic trait with group size (species richness of a flock)
summary(lm(I(difobsnull.bmass/sdnull.bmass)~null.resu.bmass[,1,1])) # negative association

# association between FRic trait with group size (species richness of a flock)
summary(lm(indicesfd_bmass$FRic~null.resu.bmass[,1,1])) # positive association

# SES
summary(lmer(I(difobsnull.bmass/sdnull.bmass)~null.resu.bmass[,1,1]+indicesfd_env$Temp+indicesfd_env$Prec+indicesfd_env$npp+indicesfd_env$Edge+indicesfd_env$Fitofisionomia+indicesfd_env$Ecoregion+(1|indicesfd_env$Localidade)))

summary(lmer(indicesfd_bmass$FRic~null.resu.bmass[,1,1]+indicesfd_env$Temp+indicesfd_env$Prec+indicesfd_env$npp+indicesfd_env$Edge+indicesfd_env$Fitofisionomia+indicesfd_env$Ecoregion+(1|indicesfd_env$Localidade)))

#
## End