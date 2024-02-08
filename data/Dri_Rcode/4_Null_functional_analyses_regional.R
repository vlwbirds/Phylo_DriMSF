#### Code that follows the paper Mixed-species bird flocks enhance the benefits of group aggregation by minimizing variation in some functional traits while maximizing variation in others by GF Dri, NC Caceres, J Della-Flora, CS Dambros
#### Written by GF Dri and CS Dambros 
#### Contact gabriela.franzoi at maine.edu 
#### November 2021

## This code performs null functional richness for each trait individually (SES FRic)
## Regional scale (flock) analysis

## You will need objects created in the previous codes (1_Functional_analysis_traitscombined and 2_Functional_analysis_traitsindividual)

##### FUNCTIONAL ANALYSIS #####
#### FRic Wing ####
indiceswingfunsite<-function(i){
  require(FD)

  trait_sp_rand<-trait_sp[sample(nrow(trait_sp)),]
  rownames(trait_sp_rand)<-rownames(trait_sp)
  
  wing_rand<-data.frame(trait_sp_rand[,1], row.names = row.names(trait_sp_rand))
  
   # FRic with randomized traits (SES FRic)
  null_indiceswing_site_pre<- dbFD(wing_rand, occ_sp_site, w.abun=FALSE, stand.x = FALSE, corr="none", calc.FRic=TRUE, m="min",
                                   stand.FRic=TRUE, scale.RaoQ=FALSE, 
                                   calc.CWM=FALSE, calc.FDiv= TRUE, print.pco=FALSE, messages = FALSE
  )
  
  class(null_indiceswing_site_pre)
  
  null_indiceswing_site_vec<-unlist(null_indiceswing_site_pre)
  
  # Save results
  dir.create("NullResultIndiceWingSite",showWarnings = FALSE)
  
  writeLines(as.character(null_indiceswing_site_vec),paste0("NullResultIndiceWingSite/resuIndiceWingSite",i,".txt"))
  
}

# paralleling the analysis to run the function 99 times 
library(snow)
{
  cl<-makeCluster(3)
  clusterExport(cl,list("trait_sp","occ_sp_site"))
  clusterApply(cl,1:99,indiceswingfunsite)
  stopCluster(cl)
} 

# importing back to R the randomized functional indices
null.files.wing.site<-list.files("NullResultIndiceWingSite/",full.names = TRUE)

check2.wingsite<-readLines(null.files.wing.site[1])

# we need to remove one row for each file, so lets find the row and remove it
toRemove.wingSite<-sum(c(length(indicesfd_site_wing$nbsp),
                         length(indicesfd_site_wing$sing.sp),
                         length(indicesfd_site_wing$FRic),
                         length(indicesfd_site_wing$qual.FRic)))

checknum2.wingsite<-as.numeric(check2.wingsite[-toRemove.wingSite]) # removing the row 

# organizing the files 
checknum3.wingsite<-lapply(null.files.wing.site,function(x){
  check2.wingsite<-readLines(x)
  checknum2.wingsite<-as.numeric(check2.wingsite[-toRemove.wingSite])
})

# tranforming it as an array (3 dimension matrix)
ver.wingsite<-do.call(c,checknum3.wingsite)
null.resu.wing.site<-array(ver.wingsite,dim=c(length(checknum3.wingsite[[1]])/6,6,length(checknum3.wingsite))) 

# observed trait FRic - average of all randomized (null) SES FRic trait
difobsnull.wingsite<-(indicesfd_site_wing$FRic-rowMeans(null.resu.wing.site[,3,],na.rm = TRUE))

# standard deviation of null FRic trait (SES FRic)
sdnull.wingsite<-apply(null.resu.wing.site[,3,],na.rm = TRUE,1,sd)

# Association between observed and expected FRic
summary(lm(indicesfd_site_wing$FRic~I(difobsnull.wingsite/sdnull.wingsite)))
plot(indicesfd_site_wing$FRic~I(difobsnull.wingsite/sdnull.wingsite)) # strong association

# association between SES FRic trait with group size (species richness of a flock)
summary(lm(I(difobsnull.wingsite/sdnull.wingsite)~null.resu.wing.site[,1,1])) # no evidence for the association

# linear regression  model to investigate the association between observed and expected FRic of the trait with group  size and enviromental covariates
# SES
summary(lm(I(difobsnull.wingsite/sdnull.wingsite)~null.resu.wing.site[,1,1]+indicesfd_site_env$Temp+indicesfd_site_env$Prec+indicesfd_site_env$npp+indicesfd_site_env$Edge+indicesfd_site_env$Fitofisionomia+indicesfd_site_env$Ecoregion))

# OBS
summary(lm(indicesfd_site_wing$FRic~null.resu.wing.site[,1,1]+indicesfd_site_env$Temp+indicesfd_site_env$Prec+indicesfd_site_env$npp+indicesfd_site_env$Edge+indicesfd_site_env$Fitofisionomia+indicesfd_site_env$Ecoregion))

#### FRic Beak ####
indicesbeakfunsite<-function(i){
  require(FD)

    trait_sp_rand<-trait_sp[sample(nrow(trait_sp)),]
  rownames(trait_sp_rand)<-rownames(trait_sp)
  
  beak_rand<-data.frame(trait_sp_rand[,4], row.names = row.names(trait_sp_rand))
  
   # FRic with randomized traits (SES FRic)
  null_indicesbeak_site_pre<- dbFD(
    beak_rand, occ_sp_site, w.abun=FALSE, stand.x = FALSE, corr="none", calc.FRic=TRUE, m="min",
    stand.FRic=TRUE, scale.RaoQ=FALSE, 
    calc.CWM=FALSE, calc.FDiv= TRUE, print.pco=FALSE, messages = FALSE
  )
  
  class(null_indicesbeak_site_pre)
  
  null_indicesbeak_site_vec<-unlist(null_indicesbeak_site_pre)
  
  # Save results
  dir.create("NullResultIndiceBeakSite",showWarnings = FALSE)
  
  writeLines(as.character(null_indicesbeak_site_vec),paste0("NullResultIndiceBeakSite/resuIndiceBeakSite",i,".txt"))
  
}

# paralleling the analysis to run the function 99 times 
library(snow)
{
  cl<-makeCluster(3)
  clusterExport(cl,list("trait_sp","occ_sp_site"))
  clusterApply(cl,1:99,indicesbeakfunsite)
  stopCluster(cl)
} 

# importing back to R the randomized functional indices
null.files.beak.site<-list.files("NullResultIndiceBeakSite/",full.names = TRUE)

check2.beaksite<-readLines(null.files.beak.site[1])

# we need to remove one row for each file, so lets find the row and remove it
toRemove.beakSite<-sum(c(length(indicesfd_site_beak$nbsp),
                         length(indicesfd_site_beak$sing.sp),
                         length(indicesfd_site_beak$FRic),
                         length(indicesfd_site_beak$qual.FRic)))

checknum2.beaksite<-as.numeric(check2.beaksite[-toRemove.beakSite]) # removing the row 

# organizing the files 
checknum3.beaksite<-lapply(null.files.beak.site,function(x){
  check2.beaksite<-readLines(x)
  checknum2.beaksite<-as.numeric(check2.beaksite[-toRemove.beakSite])
})

# tranforming it as an array (3 dimension matrix)
ver.beaksite<-do.call(c,checknum3.beaksite)
null.resu.beak.site<-array(ver.beaksite,dim=c(length(checknum3.beaksite[[1]])/6,6,length(checknum3.beaksite))) 

dim(null.resu.beak.site)
class(null.resu.beak.site)
str(null.resu.beak.site)

# observed trait FRic - average of all randomized (null) SES FRic trait
difobsnull.beaksite<-(indicesfd_site_beak$FRic-rowMeans(null.resu.beak.site[,3,],na.rm = TRUE))

# standard deviation of null FRic trait (SES FRic)
sdnull.beaksite<-apply(null.resu.beak.site[,3,],na.rm = TRUE,1,sd)

# Association between observed and expected FRic
summary(lm(indicesfd_site_beak$FRic~I(difobsnull.beaksite/sdnull.beaksite)))
plot(indicesfd_site_beak$FRic~I(difobsnull.beaksite/sdnull.beaksite)) # strong association

# association between SES FRic trait with group size (species richness of a flock)
summary(lm(I(difobsnull.beaksite/sdnull.beaksite)~null.resu.beak.site[,1,1])) # no evidence for the association

# linear regression  model to investigate the association between observed and expected FRic of the trait with group  size and enviromental covariates
# SES
summary(lm(I(difobsnull.beaksite/sdnull.beaksite)~null.resu.beak.site[,1,1]+indicesfd_site_env$Temp+indicesfd_site_env$Prec+indicesfd_site_env$npp+indicesfd_site_env$Edge+indicesfd_site_env$Fitofisionomia+indicesfd_site_env$Ecoregion))

# OBS
summary(lm(indicesfd_site_beak$FRic~null.resu.beak.site[,1,1]+indicesfd_site_env$Temp+indicesfd_site_env$Prec+indicesfd_site_env$npp+indicesfd_site_env$Edge+indicesfd_site_env$Fitofisionomia+indicesfd_site_env$Ecoregion))

#### FRic Tail ####
indicesretrizesfunsite<-function(i){
  require(FD)
 
  trait_sp_rand<-trait_sp[sample(nrow(trait_sp)),]
  rownames(trait_sp_rand)<-rownames(trait_sp)
  
  retrizes_rand<-data.frame(trait_sp_rand[,2], row.names = row.names(trait_sp_rand))
  
   # FRic with randomized traits (SES FRic)
  null_indicesretrizes_site_pre<- dbFD(
    retrizes_rand, occ_sp_site, w.abun=FALSE, stand.x = FALSE, corr="none", calc.FRic=TRUE, m="min",
    stand.FRic=TRUE, scale.RaoQ=FALSE, 
    calc.CWM=FALSE, calc.FDiv= TRUE, print.pco=FALSE, messages = FALSE
  )
  
  class(null_indicesretrizes_site_pre)
  
  null_indicesretrizes_site_vec<-unlist(null_indicesretrizes_site_pre)
  
  # Save results
  dir.create("NullResultIndiceRetrizesSite",showWarnings = FALSE)
  
  writeLines(as.character(null_indicesretrizes_site_vec),paste0("NullResultIndiceRetrizesSite/resuIndiceRetrizesSite",i,".txt"))
  
}

# paralleling the analysis to run the function 99 times 
library(snow)
{
  cl<-makeCluster(3)
  clusterExport(cl,list("trait_sp","occ_sp_site"))
  clusterApply(cl,1:99,indicesretrizesfunsite)
  stopCluster(cl)
} 

# importing back to R the randomized functional indices
null.files.retrizes.site<-list.files("NullResultIndiceRetrizesSite/",full.names = TRUE)

check2.retrizessite<-readLines(null.files.retrizes.site[1])

# we need to remove one row for each file, so lets find the row and remove it
toRemove.retrizesSite<-sum(c(length(indicesfd_site_retrizes$nbsp),
                             length(indicesfd_site_retrizes$sing.sp),
                             length(indicesfd_site_retrizes$FRic),
                             length(indicesfd_site_retrizes$qual.FRic)))

checknum2.retrizessite<-as.numeric(check2.retrizessite[-toRemove.retrizesSite]) # removing the row 

# organizing the files 
checknum3.retrizessite<-lapply(null.files.retrizes.site,function(x){
  check2.retrizessite<-readLines(x)
  checknum2.retrizessite<-as.numeric(check2.retrizessite[-toRemove.retrizesSite])
})

# tranforming it as an array (3 dimension matrix)
ver.retrizessite<-do.call(c,checknum3.retrizessite)
null.resu.retrizes.site<-array(ver.retrizessite,dim=c(length(checknum3.retrizessite[[1]])/6,6,length(checknum3.retrizessite))) 

# observed trait FRic - average of all randomized (null) SES FRic trait
difobsnull.retrizessite<-(indicesfd_site_retrizes$FRic-rowMeans(null.resu.retrizes.site[,3,],na.rm = TRUE))

# standard deviation of null FRic trait (SES FRic)
sdnull.retrizessite<-apply(null.resu.retrizes.site[,3,],na.rm = TRUE,1,sd)

# Association between observed and expected FRic
summary(lm(indicesfd_site_retrizes$FRic~I(difobsnull.retrizessite/sdnull.retrizessite)))
plot(indicesfd_site_retrizes$FRic~I(difobsnull.retrizessite/sdnull.retrizessite)) # strong association

# association between SES FRic trait with group size (species richness of a flock)
summary(lm(I(difobsnull.retrizessite/sdnull.retrizessite)~null.resu.retrizes.site[,1,1])) # no evidence for the association

# linear regression  model to investigate the association between observed and expected FRic of the trait with group  size and enviromental covariates
# SES
summary(lm(I(difobsnull.retrizessite/sdnull.retrizessite)~null.resu.retrizes.site[,1,1]+indicesfd_site_env$Temp+indicesfd_site_env$Prec+indicesfd_site_env$npp+indicesfd_site_env$Edge+indicesfd_site_env$Fitofisionomia+indicesfd_site_env$Ecoregion))

# OBS
summary(lm(indicesfd_site_retrizes$FRic~null.resu.retrizes.site[,1,1]+indicesfd_site_env$Temp+indicesfd_site_env$Prec+indicesfd_site_env$npp+indicesfd_site_env$Edge+indicesfd_site_env$Fitofisionomia+indicesfd_site_env$Ecoregion))


#### FRic Tarsus ####
indicestarsusfunsite<-function(i){
  require(FD)

  trait_sp_rand<-trait_sp[sample(nrow(trait_sp)),]
  rownames(trait_sp_rand)<-rownames(trait_sp)
  
  tarsus_rand<-data.frame(trait_sp_rand[,9], row.names = row.names(trait_sp_rand))
  
   # FRic with randomized traits (SES FRic)
  null_indicestarsus_site_pre<- dbFD(
    tarsus_rand, occ_sp_site, w.abun=FALSE, stand.x = FALSE, corr="none", calc.FRic=TRUE, m="min",
    stand.FRic=TRUE, scale.RaoQ=FALSE, 
    calc.CWM=FALSE, calc.FDiv= TRUE, print.pco=FALSE, messages = FALSE
  )
  
  class(null_indicestarsus_site_pre)
  
  null_indicestarsus_site_vec<-unlist(null_indicestarsus_site_pre)
  
  # Save results
  dir.create("NullResultIndiceTarsusSite",showWarnings = FALSE)
  
  writeLines(as.character(null_indicestarsus_site_vec),paste0("NullResultIndiceTarsusSite/resuIndiceTarsusSite",i,".txt"))
  
}

# paralleling the analysis to run the function 99 times 
library(snow)
{
  cl<-makeCluster(3)
  clusterExport(cl,list("trait_sp","occ_sp_site"))
  clusterApply(cl,1:99,indicestarsusfunsite)
  stopCluster(cl)
} 

# importing back to R the randomized functional indices
null.files.tarsus.site<-list.files("NullResultIndiceTarsusSite/",full.names = TRUE)

check2.tarsussite<-readLines(null.files.tarsus.site[1])

# we need to remove one row for each file, so lets find the row and remove it
toRemove.tarsusSite<-sum(c(length(indicesfd_site_tarsus$nbsp),
                           length(indicesfd_site_tarsus$sing.sp),
                           length(indicesfd_site_tarsus$FRic),
                           length(indicesfd_site_tarsus$qual.FRic)))

checknum2.tarsussite<-as.numeric(check2.tarsussite[-toRemove.tarsusSite]) # removing the row 

# organizing the files 
checknum3.tarsussite<-lapply(null.files.tarsus.site,function(x){
  check2.tarsussite<-readLines(x)
  checknum2.tarsussite<-as.numeric(check2.tarsussite[-toRemove.tarsusSite])
})

# tranforming it as an array (3 dimension matrix)
ver.tarsussite<-do.call(c,checknum3.tarsussite)
null.resu.tarsus.site<-array(ver.tarsussite,dim=c(length(checknum3.tarsussite[[1]])/6,6,length(checknum3.tarsussite))) 

# observed trait FRic - average of all randomized (null) SES FRic trait
difobsnull.tarsussite<-(indicesfd_site_tarsus$FRic-rowMeans(null.resu.tarsus.site[,3,],na.rm = TRUE))

# standard deviation of null FRic trait (SES FRic)
sdnull.tarsussite<-apply(null.resu.tarsus.site[,3,],na.rm = TRUE,1,sd)

# Association between observed and expected FRic
summary(lm(indicesfd_site_tarsus$FRic~I(difobsnull.tarsussite/sdnull.tarsussite)))
plot(indicesfd_site_tarsus$FRic~I(difobsnull.tarsussite/sdnull.tarsussite)) # strong association

# association between SES FRic trait with group size (species richness of a flock)
summary(lm(I(difobsnull.tarsussite/sdnull.tarsussite)~null.resu.tarsus.site[,1,1])) # no evidence for the association

# linear regression  model to investigate the association between observed and expected FRic of the trait with group  size and enviromental covariates
# SES
summary(lm(I(difobsnull.tarsussite/sdnull.tarsussite)~null.resu.tarsus.site[,1,1]+indicesfd_site_env$Temp+indicesfd_site_env$Prec+indicesfd_site_env$npp+indicesfd_site_env$Edge+indicesfd_site_env$Fitofisionomia+indicesfd_site_env$Ecoregion))

# OBS
summary(lm(indicesfd_site_tarsus$FRic~null.resu.tarsus.site[,1,1]+indicesfd_site_env$Temp+indicesfd_site_env$Prec+indicesfd_site_env$npp+indicesfd_site_env$Edge+indicesfd_site_env$Fitofisionomia+indicesfd_site_env$Ecoregion))

#### FRic Diet ####
indicesdietfunsite<-function(i){
  require(FD)
  
  trait_sp_rand<-trait_sp[sample(nrow(trait_sp)),]
  rownames(trait_sp_rand)<-rownames(trait_sp)
  
  diet_rand<-data.frame(trait_sp_rand[,7], row.names = row.names(trait_sp_rand))
  
   # FRic with randomized traits (SES FRic)
  null_indicesdiet_site_pre<- dbFD(diet_rand, occ_sp_site, w.abun=FALSE, stand.x = FALSE, corr="none", calc.FRic=TRUE, m="min",
    stand.FRic=TRUE, scale.RaoQ=FALSE, 
    calc.CWM=FALSE, calc.FDiv= TRUE, print.pco=FALSE, messages = FALSE
  )
  
  class(null_indicesdiet_site_pre)
  
  null_indicesdiet_site_vec<-unlist(null_indicesdiet_site_pre)
  
  # Save results
  dir.create("NullResultIndiceDietSite",showWarnings = FALSE)
  
  writeLines(as.character(null_indicesdiet_site_vec),paste0("NullResultIndiceDietSite/resuIndiceDietSite",i,".txt"))
  
}

# paralleling the analysis to run the function 99 times 
library(snow)
{
  cl<-makeCluster(3)
  clusterExport(cl,list("trait_sp","occ_sp_site"))
  clusterApply(cl,1:99,indicesdietfunsite)
  stopCluster(cl)
} 

# importing back to R the randomized functional indices
null.files.diet.site<-list.files("NullResultIndiceDietSite/",full.names = TRUE)

check2.dietsite<-readLines(null.files.diet.site[1])

# we need to remove one row for each file, so lets find the row and remove it
toRemove.dietSite<-sum(c(length(indicesfd_site_diet$nbsp),
                         length(indicesfd_site_diet$sing.sp),
                         length(indicesfd_site_diet$FRic),
                         length(indicesfd_site_diet$qual.FRic)))

checknum2.dietsite<-as.numeric(check2.dietsite[-toRemove.dietSite]) # removing the row 

# organizing the files 
checknum3.dietsite<-lapply(null.files.diet.site,function(x){
  check2.dietsite<-readLines(x)
  checknum2.dietsite<-as.numeric(check2.dietsite[-toRemove.dietSite])
})

# tranforming it as an array (3 dimension matrix)
ver.dietsite<-do.call(c,checknum3.dietsite)
null.resu.diet.site<-array(ver.dietsite,dim=c(length(checknum3.dietsite[[1]])/6,6,length(checknum3.dietsite))) 

# observed trait FRic - average of all randomized (null) SES FRic trait
difobsnull.dietsite<-(indicesfd_site_diet$FRic-rowMeans(null.resu.diet.site[,3,],na.rm = TRUE))

# standard deviation of null FRic trait (SES FRic)
sdnull.dietsite<-apply(null.resu.diet.site[,3,],na.rm = TRUE,1,sd)

# Association between observed and expected FRic
summary(lm(indicesfd_site_diet$FRic~I(difobsnull.dietsite/sdnull.dietsite)))
plot(indicesfd_site_diet$FRic~I(difobsnull.dietsite/sdnull.dietsite)) # strong association

# association between SES FRic trait with group size (species richness of a flock)
summary(lm(I(difobsnull.dietsite/sdnull.dietsite)~null.resu.diet.site[,1,1])) # no evidence for the association

# linear regression  model to investigate the association between observed and expected FRic of the trait with group  size and enviromental covariates
# SES
summary(lm(I(difobsnull.dietsite/sdnull.dietsite)~null.resu.diet.site[,1,1]+indicesfd_site_env$Temp+indicesfd_site_env$Prec+indicesfd_site_env$npp+indicesfd_site_env$Edge+indicesfd_site_env$Fitofisionomia+indicesfd_site_env$Ecoregion))

# OBS
summary(lm(indicesfd_site_diet$FRic~null.resu.diet.site[,1,1]+indicesfd_site_env$Temp+indicesfd_site_env$Prec+indicesfd_site_env$npp+indicesfd_site_env$Edge+indicesfd_site_env$Fitofisionomia+indicesfd_site_env$Ecoregion))

#### FRic Foraging strata ####
indicesestfunsite<-function(i){
  require(FD)

  trait_sp_rand<-trait_sp[sample(nrow(trait_sp)),]
  rownames(trait_sp_rand)<-rownames(trait_sp)
  
  est_fuz_rand<-prep.fuzzy(trait_sp_rand[,25:29],col.blocks = 5)
  strata_ktab_rand<-ktab.list.df(list(est_fuz_rand))
  strata_dist_rand<-dist.ktab(x=strata_ktab_rand, type = "F", option = "scaledBYrange")
  
  pcoa_est_rand<-cmdscale(strata_dist_rand, add=T)
  pcoa1_est_rand<-scores(pcoa_est_rand)[,1]
  
   # FRic with randomized traits (SES FRic)
  null_indicesest_site_pre<- dbFD(
    pcoa1_est_rand, occ_sp_site, w.abun=FALSE, stand.x = FALSE, corr="none", calc.FRic=TRUE, m="min",
    stand.FRic=TRUE, scale.RaoQ=FALSE, 
    calc.CWM=FALSE, calc.FDiv= TRUE, print.pco=FALSE, messages = FALSE
  )
  
  class(null_indicesest_site_pre)
  
  null_indicesest_site_vec<-unlist(null_indicesest_site_pre)
  
  # Save results
  dir.create("NullResultIndiceForagSite",showWarnings = FALSE)
  
  writeLines(as.character(null_indicesest_site_vec),paste0("NullResultIndiceForagSite/resuIndiceForagSite",i,".txt"))
  
}

# paralleling the analysis to run the function 99 times 
library(snow)
{
  cl<-makeCluster(3)
  clusterExport(cl,list("trait_sp","occ_sp_site"))
  clusterApply(cl,1:99,indicesestfunsite)
  stopCluster(cl)
} 

# importing back to R the randomized functional indices
null.files.forag.site<-list.files("NullResultIndiceForagSite/",full.names = TRUE)

check2.foragsite<-readLines(null.files.forag.site[1])

# we need to remove one row for each file, so lets find the row and remove it
toRemove.estSite<-sum(c(length(indicesfd_site_est$nbsp),
                        length(indicesfd_site_est$sing.sp),
                        length(indicesfd_site_est$FRic),
                        length(indicesfd_site_est$qual.FRic)))

checknum2.estsite<-as.numeric(check2.foragsite[-toRemove.estSite]) # removing the row 

# organizing the files 
checknum3.estsite<-lapply(null.files.forag.site,function(x){
  check2.estsite<-readLines(x)
  checknum2.estsite<-as.numeric(check2.estsite[-toRemove.estSite])
})

# tranforming it as an array (3 dimension matrix)
ver.estsite<-do.call(c,checknum3.estsite)
null.resu.est.site<-array(ver.estsite,dim=c(length(checknum3.estsite[[1]])/6,6,length(checknum3.estsite))) 

# observed trait FRic - average of all randomized (null) SES FRic trait
difobsnull.estsite<-(indicesfd_site_est$FRic-rowMeans(null.resu.est.site[,3,],na.rm = TRUE))

# standard deviation of null FRic trait (SES FRic)
sdnull.estsite<-apply(null.resu.est.site[,3,],na.rm = TRUE,1,sd)

# Association between observed and expected FRic
summary(lm(indicesfd_site_est$FRic~I(difobsnull.estsite/sdnull.estsite)))
plot(indicesfd_site_bmass$FRic~I(difobsnull.estsite/sdnull.estsite)) # strong association

# association between SES FRic trait with group size (species richness of a flock)
summary(lm(I(difobsnull.estsite/sdnull.estsite)~null.resu.est.site[,1,1])) # no evidence for the association

# linear regression  model to investigate the association between observed and expected FRic of the trait with group  size and enviromental covariates
# SES
summary(lm(I(difobsnull.estsite/sdnull.estsite)~null.resu.est.site[,1,1]+indicesfd_site_env$Temp+indicesfd_site_env$Prec+indicesfd_site_env$npp+indicesfd_site_env$Edge+indicesfd_site_env$Fitofisionomia+indicesfd_site_env$Ecoregion))

# OBS
summary(lm(indicesfd_site_est$FRic~null.resu.est.site[,1,1]+indicesfd_site_env$Temp+indicesfd_site_env$Prec+indicesfd_site_env$npp+indicesfd_site_env$Edge+indicesfd_site_env$Fitofisionomia+indicesfd_site_env$Ecoregion))

#### FRic Body Mass ####
indicesbmassfunsite<-function(i){
  require(FD)
  
  trait_sp_rand<-trait_sp[sample(nrow(trait_sp)),]
  rownames(trait_sp_rand)<-rownames(trait_sp)
  
  body_mass_rand<-data.frame(trait_sp_rand[,30], row.names = row.names(trait_sp_rand))
  
   # FRic with randomized traits (SES FRic)
  null_indicesbmass_site_pre<- dbFD(
    body_mass_rand, occ_sp_site, w.abun=FALSE, stand.x = FALSE, corr="none", calc.FRic=TRUE, m="min",
    stand.FRic=TRUE, scale.RaoQ=FALSE, 
    calc.CWM=FALSE, calc.FDiv= TRUE, print.pco=FALSE, messages = FALSE
  )
  
  class(null_indicesbmass_site_pre)
  
  null_indicesbmass_site_vec<-unlist(null_indicesbmass_site_pre)
  
  # Save results
  dir.create("NullResultIndiceBMassSite",showWarnings = FALSE)
  
  writeLines(as.character(null_indicesbmass_site_vec),paste0("NullResultIndiceBMassSite/resuIndiceBMassSite",i,".txt"))
  
}

# paralleling the analysis to run the function 99 times 
library(snow)
{
  cl<-makeCluster(3)
  clusterExport(cl,list("trait_sp","occ_sp_site"))
  clusterApply(cl,1:99,indicesbmassfunsite)
  stopCluster(cl)
} 

# importing back to R the randomized functional indices
null.files.bmass.site<-list.files("NullResultIndiceBMassSite/",full.names = TRUE)

check2.bmasssite<-readLines(null.files.bmass.site[1])

# we need to remove one row for each file, so lets find the row and remove it
toRemove.bmassSite<-sum(c(length(indicesfd_site_bmass$nbsp),
                          length(indicesfd_site_bmass$sing.sp),
                          length(indicesfd_site_bmass$FRic),
                          length(indicesfd_site_bmass$qual.FRic)))

checknum2.bmasssite<-as.numeric(check2.bmasssite[-toRemove.bmassSite]) # removing the row 

# organizing the files 
checknum3.bmasssite<-lapply(null.files.bmass.site,function(x){
  check2.bmasssite<-readLines(x)
  checknum2.bmasssite<-as.numeric(check2.bmasssite[-toRemove.bmassSite])
})

# tranforming it as an array (3 dimension matrix)
ver.bmasssite<-do.call(c,checknum3.bmasssite)
null.resu.bmass.site<-array(ver.bmasssite,dim=c(length(checknum3.bmasssite[[1]])/6,6,length(checknum3.bmasssite))) 

# observed trait FRic - average of all randomized (null) SES FRic trait
difobsnull.bmasssite<-(indicesfd_site_bmass$FRic-rowMeans(null.resu.bmass.site[,3,],na.rm = TRUE))

# standard deviation of null FRic trait (SES FRic)
sdnull.bmasssite<-apply(null.resu.bmass.site[,3,],na.rm = TRUE,1,sd)

# Association between observed and expected FRic
summary(lm(indicesfd_site_bmass$FRic~I(difobsnull.bmasssite/sdnull.bmasssite)))
plot(indicesfd_site_bmass$FRic~I(difobsnull.bmasssite/sdnull.bmasssite)) # strong association

# association between SES FRic trait with group size (species richness of a flock)
summary(lm(I(difobsnull.bmasssite/sdnull.bmasssite)~null.resu.bmass.site[,1,1])) # no evidence for the association

# linear regression  model to investigate the association between observed and expected FRic of the trait with group  size and enviromental covariates
# SES
summary(lm(I(difobsnull.bmasssite/sdnull.bmasssite)~null.resu.bmass.site[,1,1]+indicesfd_site_env$Temp+indicesfd_site_env$Prec+indicesfd_site_env$npp+indicesfd_site_env$Edge+indicesfd_site_env$Fitofisionomia+indicesfd_site_env$Ecoregion))

# OBS
summary(lm(indicesfd_site_bmass$FRic~null.resu.bmass.site[,1,1]+indicesfd_site_env$Temp+indicesfd_site_env$Prec+indicesfd_site_env$npp+indicesfd_site_env$Edge+indicesfd_site_env$Fitofisionomia+indicesfd_site_env$Ecoregion))

#
## End