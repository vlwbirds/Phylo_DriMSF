#### Code that follows the paper Mixed-species bird flocks enhance the benefits of group aggregation by minimizing variation in some functional traits while maximizing variation in others by GF Dri, NC Caceres, J Della-Flora, CS Dambros
#### Written by GF Dri and CS Dambros 
#### Contact gabriela.franzoi at maine.edu 
#### November 2021

## This code performs linear mixed-effect models at local scale (flock) and linear regression models at regional scale (site)

## You will need objects created in the previous codes (1_Functional_analysis_traitscombined; 2_Functional_analysis_traitsindividual; 3_Null_functional_analyses_local and 4_Null_functional_analyses_regional)

library(MuMIn)
library(lme4)
library(lmerTest)

### LOCAL ####
#### EXPECTED ####
#### all traits combined ####
data<-data.frame(scale(ses_mf),Temp=scale(indicesfd_env$Temp),npp=scale(indicesfd_env$npp),Fitofisionomia=indicesfd_env$Fitofisionomia,Ecoregion=indicesfd_env$Ecoregion,Edge=scale(indicesfd_env$Edge),log.nbsp=scale(log(indicesfd_env$nbsp)),Localidade=indicesfd_env$Localidade)

summary(M1<-lmer(ses_mf~Temp+npp+Fitofisionomia+Ecoregion+log.nbsp+(1|Localidade),data=data))

summary(M2<-lmer(ses_mf~Temp+npp+Ecoregion+log.nbsp+(1|Localidade),data=data))

summary(M3<-lmer(ses_mf~Temp+npp+Fitofisionomia+log.nbsp+(1|Localidade),data=data))

anova(M1,M2)
anova(M1,M3)

r.squaredGLMM(M1)

##### diet ####
summary(M1<-lmer(scale(sesfric_diet)~scale(log(null.resu.diet[,1,1]))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Ecoregion+indicesfd_env$Fitofisionomia+(1|indicesfd_env$Localidade)))

summary(M2<-lmer(scale(sesfric_diet)~scale(log(null.resu.diet[,1,1]))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Ecoregion+(1|indicesfd_env$Localidade)))

summary(M3<-lmer(scale(sesfric_diet)~scale(log(null.resu.diet[,1,1]))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Fitofisionomia+(1|indicesfd_env$Localidade)))

anova(M1,M2)
anova(M1,M3)

r.squaredGLMM(M1)

#### wings ####
sesfric_wing<-difobsnull.wing/sdnull.wing
summary(M1<-lmer(scale(sesfric_wing)~scale(log(null.resu.wing[,1,1]))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Ecoregion+indicesfd_env$Fitofisionomia+(1|indicesfd_env$Localidade)))

summary(M2<-lmer(scale(sesfric_wing)~scale(log(null.resu.wing[,1,1]))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Ecoregion+(1|indicesfd_env$Localidade)))

summary(M3<-lmer(scale(sesfric_wing)~scale(log(null.resu.wing[,1,1]))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Fitofisionomia+(1|indicesfd_env$Localidade)))

anova(M1,M2)
anova(M1,M3)
r.squaredGLMM(M2)

#### tail ####
sesfric_rectrizes<-difobsnull.retrizes/sdnull.retrizes

summary(M1<-lmer(scale(sesfric_rectrizes)~scale(log(null.resu.retrizes[,1,1]))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Ecoregion+indicesfd_env$Fitofisionomia+(1|indicesfd_env$Localidade)))

summary(M2<-lmer(scale(sesfric_rectrizes)~scale(log(null.resu.retrizes[,1,1]))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Ecoregion+(1|indicesfd_env$Localidade)))

summary(M3<-lmer(scale(sesfric_rectrizes)~scale(log(null.resu.retrizes[,1,1]))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Fitofisionomia+(1|indicesfd_env$Localidade)))

anova(M1,M2)
anova(M1,M3)
r.squaredGLMM(M2)


#### body mass ####
sesfric_bmass<-(difobsnull.bmass/sdnull.bmass)

summary(M1<-lmer(scale(sesfric_bmass)~scale(log(null.resu.bmass[,1,1]))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Ecoregion+indicesfd_env$Fitofisionomia+(1|indicesfd_env$Localidade)))

summary(M2<-lmer(scale(sesfric_bmass)~scale(log(null.resu.bmass[,1,1]))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Ecoregion+(1|indicesfd_env$Localidade)))

summary(M3<-lmer(scale(sesfric_bmass)~scale(log(null.resu.bmass[,1,1]))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Fitofisionomia+(1|indicesfd_env$Localidade)))

anova(M1,M2)
anova(M1,M3)

r.squaredGLMM(M1)

#### beak ####
sesfric_beak<-difobsnull.beak/sdnull.beak

summary(M1<-lmer(scale(sesfric_beak)~scale(log(null.resu.beak[,1,1]))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Ecoregion+indicesfd_env$Fitofisionomia+(1|indicesfd_env$Localidade)))

summary(M2<-lmer(scale(sesfric_beak)~scale(log(null.resu.beak[,1,1]))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Ecoregion+(1|indicesfd_env$Localidade)))

summary(M3<-lmer(scale(sesfric_beak)~scale(log(null.resu.beak[,1,1]))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Fitofisionomia+(1|indicesfd_env$Localidade)))

anova(M1,M2)
anova(M1,M3)

r.squaredGLMM(M2)

#### foraging strata ####
sesfric_est<-difobsnull.est/sdnull.est

summary(M1<-lmer(scale(sesfric_est)~scale(log(null.resu.est[,1,1]))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Ecoregion+indicesfd_env$Fitofisionomia+(1|indicesfd_env$Localidade)))

summary(M2<-lmer(scale(sesfric_est)~scale(log(null.resu.est[,1,1]))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Ecoregion+(1|indicesfd_env$Localidade)))

summary(M3<-lmer(scale(sesfric_est)~scale(log(null.resu.est[,1,1]))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Fitofisionomia+(1|indicesfd_env$Localidade)))

anova(M1,M2)
anova(M1,M3)

r.squaredGLMM(M2)

#### tarsus ####
sesfric_tarsus<-difobsnull.tarsus/sdnull.tarsus

summary(M1<-lmer(scale(sesfric_tarsus)~scale(log(null.resu.tarsus[,1,1]))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Ecoregion+indicesfd_env$Fitofisionomia+(1|indicesfd_env$Localidade)))

summary(M2<-lmer(scale(sesfric_tarsus)~scale(log(null.resu.tarsus[,1,1]))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Ecoregion+(1|indicesfd_env$Localidade)))

summary(M3<-lmer(scale(sesfric_tarsus)~scale(log(null.resu.tarsus[,1,1]))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Fitofisionomia+(1|indicesfd_env$Localidade)))

anova(M1,M2)
anova(M1,M3)
r.squaredGLMM(M2)

#### OBSERVED ####
#### all traits combined ####
summary(M1<-lmer(scale(indicesfd_env$FRic)~scale(log(indicesfd_env$nbsp))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Ecoregion+indicesfd_env$Fitofisionomia+(1|indicesfd_env$Localidade)))

summary(M2<-lmer(scale(indicesfd_env$FRic)~scale(log(indicesfd_env$nbsp))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Ecoregion+(1|indicesfd_env$Localidade)))

summary(M3<-lmer(scale(indicesfd_env$FRic)~scale(log(indicesfd_env$nbsp))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Fitofisionomia+(1|indicesfd_env$Localidade)))

anova(M1,M2)
anova(M1,M3)

r.squaredGLMM(M1)

#### diet ####
summary(M2<-lmer(scale(indicesfd_diet$FRic)~scale(log(indicesfd_diet$nbsp))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Ecoregion+indicesfd_env$Fitofisionomia+(1|indicesfd_env$Localidade)))

summary(M2.1<-lmer(scale(indicesfd_diet$FRic)~scale(log(indicesfd_diet$nbsp))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Ecoregion+(1|indicesfd_env$Localidade)))

summary(M2.2<-lmer(scale(indicesfd_diet$FRic)~scale(log(indicesfd_diet$nbsp))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Fitofisionomia+(1|indicesfd_env$Localidade)))

anova(M2,M2.1)
anova(M2,M2.2)

r.squaredGLMM(M2)

#### body mass ####
summary(M1<-lmer(scale(indicesfd_bmass$FRic)~scale(log(indicesfd_bmass$nbsp))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Ecoregion+indicesfd_env$Fitofisionomia+(1|indicesfd_env$Localidade)))

summary(M2<-lmer(scale(indicesfd_bmass$FRic)~scale(log(indicesfd_bmass$nbsp))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Ecoregion+(1|indicesfd_env$Localidade)))

summary(M3<-lmer(scale(indicesfd_bmass$FRic)~scale(log(indicesfd_bmass$nbsp))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Fitofisionomia+(1|indicesfd_env$Localidade)))

anova(M1,M2)
anova(M1,M3)
r.squaredGLMM(M1)

#### beak #####
summary(M1<-lmer(scale(indicesfd_beak$FRic)~scale(log(indicesfd_beak$nbsp))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Ecoregion+indicesfd_env$Fitofisionomia+(1|indicesfd_env$Localidade)))

summary(M2<-lmer(scale(indicesfd_beak$FRic)~scale(log(indicesfd_beak$nbsp))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Ecoregion+(1|indicesfd_env$Localidade)))

summary(M3<-lmer(scale(indicesfd_beak$FRic)~scale(log(indicesfd_beak$nbsp))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Fitofisionomia+(1|indicesfd_env$Localidade)))

anova(M1,M2)
anova(M1,M3)
r.squaredGLMM(M2)

##### tail ####
summary(M1<-lmer(scale(indicesfd_retrizes$FRic)~scale(log(indicesfd_retrizes$nbsp))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Ecoregion+indicesfd_env$Fitofisionomia+(1|indicesfd_env$Localidade)))

summary(M2<-lmer(scale(indicesfd_retrizes$FRic)~scale(log(indicesfd_retrizes$nbsp))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Ecoregion+(1|indicesfd_env$Localidade)))

summary(M3<-lmer(scale(indicesfd_retrizes$FRic)~scale(log(indicesfd_retrizes$nbsp))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Fitofisionomia+(1|indicesfd_env$Localidade)))

anova(M1,M2)
anova(M1,M3)
r.squaredGLMM(M2)

#### wing ####
summary(M1<-lmer(scale(indicesfd_wing$FRic)~scale(log(indicesfd_wing$nbsp))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Ecoregion+indicesfd_env$Fitofisionomia+(1|indicesfd_env$Localidade)))

summary(M2<-lmer(scale(indicesfd_wing$FRic)~scale(log(indicesfd_wing$nbsp))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Ecoregion+(1|indicesfd_env$Localidade)))

summary(M3<-lmer(scale(indicesfd_wing$FRic)~scale(log(indicesfd_wing$nbsp))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Fitofisionomia+(1|indicesfd_env$Localidade)))

anova(M1,M2)
anova(M1,M3)
r.squaredGLMM(M2)

#### foraging strata #### 
summary(M1<-lmer(scale(indicesfd_est$FRic)~scale(log(indicesfd_est$nbsp))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Ecoregion+indicesfd_env$Fitofisionomia+(1|indicesfd_env$Localidade)))

summary(M2<-lmer(scale(indicesfd_est$FRic)~scale(log(indicesfd_est$nbsp))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Ecoregion+(1|indicesfd_env$Localidade)))

summary(M3<-lmer(scale(indicesfd_est$FRic)~scale(log(indicesfd_est$nbsp))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Fitofisionomia+(1|indicesfd_env$Localidade)))

anova(M1,M2)
anova(M1,M3)
r.squaredGLMM(M2)

#### tarsus ####
summary(M1<-lmer(scale(indicesfd_tarsus$FRic)~scale(log(indicesfd_tarsus$nbsp))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Ecoregion+indicesfd_env$Fitofisionomia+(1|indicesfd_env$Localidade)))

summary(M2<-lmer(scale(indicesfd_tarsus$FRic)~scale(log(indicesfd_tarsus$nbsp))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Ecoregion+(1|indicesfd_env$Localidade)))

summary(M3<-lmer(scale(indicesfd_tarsus$FRic)~scale(log(indicesfd_tarsus$nbsp))+scale(indicesfd_env$Temp)+scale(indicesfd_env$npp)+indicesfd_env$Fitofisionomia+(1|indicesfd_env$Localidade)))

anova(M1,M2)
anova(M1,M3)

r.squaredGLMM(M2)

#### REGIONAL ####
# Creating SES objects for each FRic
sesfric_wingsite<-difobsnull.wingsite/sdnull.wingsite
sesfric_beaksite<-difobsnull.beaksite/sdnull.beaksite
sesfric_tarsussite<-difobsnull.tarsussite/sdnull.tarsussite
sesfric_bmasssite<-difobsnull.bmasssite/sdnull.bmasssite
sesfric_dietsite<-difobsnull.dietsite/sdnull.dietsite
sesfric_estsite<-difobsnull.estsite/sdnull.estsite
sesfric_rectrizessite<-difobsnull.retrizessite/sdnull.retrizessite

#### EXPECTED ####
#### all traits combined ####
summary(M1<-lm(scale(ses_site)~scale(indicesfd_site_env$Temp)+scale(indicesfd_site_env$npp)+indicesfd_site_env$Fitofisionomia+indicesfd_site_env$Ecoregion+scale(log(null.resu.site[,1,1]))))

summary(M2<-lm(scale(ses_site)~scale(indicesfd_site_env$Temp)+scale(indicesfd_site_env$npp)+indicesfd_site_env$Ecoregion+scale(log(null.resu.site[,1,1]))))

summary(M3<-lm(scale(ses_site)~scale(indicesfd_site_env$Temp)+scale(indicesfd_site_env$npp)+indicesfd_site_env$Fitofisionomia+scale(log(null.resu.site[,1,1]))))

anova(M1,M2)
anova(M1,M3)

#### beak ####
summary(M1<-lm(scale(sesfric_beaksite)~scale(indicesfd_site_env$Temp)+scale(indicesfd_site_env$npp)+indicesfd_site_env$Fitofisionomia+indicesfd_site_env$Ecoregion+scale(log(null.resu.beak.site[,1,1]))))

summary(M2<-lm(scale(sesfric_beaksite)~scale(indicesfd_site_env$Temp)+scale(indicesfd_site_env$npp)+indicesfd_site_env$Ecoregion+scale(log(null.resu.beak.site[,1,1]))))

summary(M3<-lm(scale(sesfric_beaksite)~scale(indicesfd_site_env$Temp)+scale(indicesfd_site_env$npp)+indicesfd_site_env$Fitofisionomia+scale(log(null.resu.beak.site[,1,1]))))

anova(M1,M2)
anova(M1,M3)

#### wings ####
summary(M1<-lm(scale(sesfric_wingsite)~scale(indicesfd_site_env$Temp)+scale(indicesfd_site_env$npp)+indicesfd_site_env$Fitofisionomia+indicesfd_site_env$Ecoregion+scale(log(null.resu.wing.site[,1,1]))))

summary(M2<-lm(scale(sesfric_wingsite)~scale(indicesfd_site_env$Temp)+scale(indicesfd_site_env$npp)+indicesfd_site_env$Ecoregion+scale(log(null.resu.wing.site[,1,1]))))

summary(M3<-lm(scale(sesfric_wingsite)~scale(indicesfd_site_env$Temp)+scale(indicesfd_site_env$npp)+indicesfd_site_env$Fitofisionomia+scale(log(null.resu.wing.site[,1,1]))))

#### diet ####
summary(M1<-lm(scale(sesfric_dietsite)~scale(log(null.resu.diet.site[,1,1]))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Fitofisionomia+indicesfd_site_env$Ecoregion))

summary(M2<-lm(scale(sesfric_dietsite)~scale(log(null.resu.diet.site[,1,1]))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Ecoregion))

summary(M3<-lm(scale(sesfric_dietsite)~scale(log(null.resu.diet.site[,1,1]))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Fitofisionomia))

anova(M1,M2)
anova(M1,M3)

#### body mass ####
summary(M1<-lm(scale(sesfric_bmasssite)~scale(log(null.resu.bmass.site[,1,1]))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Fitofisionomia+indicesfd_site_env$Ecoregion))

summary(M2<-lm(scale(sesfric_bmasssite)~scale(log(null.resu.bmass.site[,1,1]))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Ecoregion))

summary(M3<-lm(scale(sesfric_bmasssite)~scale(log(null.resu.bmass.site[,1,1]))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Fitofisionomia))

anova(M1,M2)
anova(M1,M3)

#### foraging strata ####
summary(M1<-lm(scale(sesfric_estsite)~scale(log(null.resu.est.site[,1,1]))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Fitofisionomia+indicesfd_site_env$Ecoregion))

summary(M2<-lm(scale(sesfric_estsite)~scale(log(null.resu.est.site[,1,1]))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Ecoregion))

summary(M3<-lm(scale(sesfric_estsite)~scale(log(null.resu.est.site[,1,1]))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Fitofisionomia))

anova(M1,M2)
anova(M1,M3)

#### tail ####
summary(M1<-lm(scale(sesfric_rectrizessite)~scale(log(null.resu.retrizes.site[,1,1]))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Fitofisionomia+indicesfd_site_env$Ecoregion))

summary(M2<-lm(scale(sesfric_rectrizessite)~scale(log(null.resu.retrizes.site[,1,1]))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Ecoregion))

summary(M3<-lm(scale(sesfric_rectrizessite)~scale(log(null.resu.retrizes.site[,1,1]))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Fitofisionomia))

anova(M1,M2)
anova(M1,M3)

#### tarsus ####
summary(M1<-lm(scale(sesfric_tarsussite)~scale(log(null.resu.tarsus.site[,1,1]))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Fitofisionomia+indicesfd_site_env$Ecoregion))

summary(M2<-lm(scale(sesfric_tarsussite)~scale(log(null.resu.tarsus.site[,1,1]))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Ecoregion))

summary(M3<-lm(scale(sesfric_tarsussite)~scale(log(null.resu.tarsus.site[,1,1]))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Fitofisionomia))

anova(M1,M2)
anova(M1,M3)

#### OBSERVED ####
#### all traits combined ####
summary(M1<-lm(scale(indicesfd_site_env$FRic)~scale(log(indicesfd_site_env$nbsp))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Fitofisionomia+indicesfd_site_env$Ecoregion))

summary(M2<-lm(scale(indicesfd_site_env$FRic)~scale(log(indicesfd_site_env$nbsp))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Ecoregion))

summary(M3<-lm(scale(indicesfd_site_env$FRic)~scale(log(indicesfd_site_env$nbsp))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Fitofisionomia))

anova(M1,M2)
anova(M1,M3)

#### beak ####
summary(M1<-lm(scale(indicesfd_site_beak$FRic)~scale(log(indicesfd_site_beak$nbsp))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Fitofisionomia+indicesfd_site_env$Ecoregion))

summary(M2<-lm(scale(indicesfd_site_beak$FRic)~scale(log(indicesfd_site_beak$nbsp))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Ecoregion))

summary(M3<-lm(scale(indicesfd_site_beak$FRic)~scale(log(indicesfd_site_beak$nbsp))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Fitofisionomia))

anova(M1,M2)
anova(M1,M3)

#### wings ####
summary(M1<-lm(scale(indicesfd_site_wing$FRic)~scale(log(indicesfd_site_wing$nbsp))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Fitofisionomia+indicesfd_site_env$Ecoregion))

summary(M2<-lm(scale(indicesfd_site_wing$FRic)~scale(log(indicesfd_site_wing$nbsp))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Ecoregion))

summary(M3<-lm(scale(indicesfd_site_wing$FRic)~scale(log(indicesfd_site_wing$nbsp))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Fitofisionomia))

anova(M1,M2)
anova(M1,M3)

#### diet ####
summary(M1<-lm(scale(indicesfd_site_diet$FRic)~scale(log(indicesfd_site_diet$nbsp))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Fitofisionomia+indicesfd_site_env$Ecoregion))

summary(M2<-lm(scale(indicesfd_site_diet$FRic)~scale(log(indicesfd_site_diet$nbsp))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Ecoregion))

summary(M3<-lm(scale(indicesfd_site_diet$FRic)~scale(log(indicesfd_site_diet$nbsp))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Fitofisionomia))

anova(M1,M2)
anova(M1,M3)

#### body mass ####
summary(M1<-lm(scale(indicesfd_site_bmass$FRic)~scale(log(indicesfd_site_bmass$nbsp))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Fitofisionomia+indicesfd_site_env$Ecoregion))

summary(M2<-lm(scale(indicesfd_site_bmass$FRic)~scale(log(indicesfd_site_bmass$nbsp))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Ecoregion))

summary(M3<-lm(scale(indicesfd_site_bmass$FRic)~scale(log(indicesfd_site_bmass$nbsp))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Fitofisionomia))

anova(M1,M2)
anova(M1,M3)

#### foraging strata ####
summary(M1<-lm(scale(indicesfd_site_est$FRic)~scale(log(indicesfd_site_est$nbsp))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Fitofisionomia+indicesfd_site_env$Ecoregion))

summary(M2<-lm(scale(indicesfd_site_est$FRic)~scale(log(indicesfd_site_est$nbsp))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Ecoregion))

summary(M3<-lm(scale(indicesfd_site_est$FRic)~scale(log(indicesfd_site_est$nbsp))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Fitofisionomia))

anova(M1,M2)
anova(M1,M3)

#### tarsus ####
summary(M1<-lm(scale(indicesfd_site_tarsus$FRic)~scale(log(indicesfd_site_tarsus$nbsp))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Fitofisionomia+indicesfd_site_env$Ecoregion))

summary(M2<-lm(scale(indicesfd_site_tarsus$FRic)~scale(log(indicesfd_site_tarsus$nbsp))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Ecoregion))

summary(M3<-lm(scale(indicesfd_site_tarsus$FRic)~scale(log(indicesfd_site_tarsus$nbsp))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Fitofisionomia))

anova(M1,M2)
anova(M1,M3)

#### tail ####
summary(M1<-lm(scale(indicesfd_site_retrizes$FRic)~scale(log(indicesfd_site_retrizes$nbsp))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Fitofisionomia+indicesfd_site_env$Ecoregion))

summary(M2<-lm(scale(indicesfd_site_retrizes$FRic)~scale(log(indicesfd_site_retrizes$nbsp))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Ecoregion))

summary(M3<-lm(scale(indicesfd_site_retrizes$FRic)~scale(log(indicesfd_site_retrizes$nbsp))+indicesfd_site_env$Temp+indicesfd_site_env$npp+indicesfd_site_env$Fitofisionomia))

anova(M1,M2)
anova(M1,M3)

#
## End