#### Code that follows the paper Mixed-species bird flocks enhance the benefits of group aggregation by minimizing variation in some functional traits while maximizing variation in others by GF Dri, NC Caceres, J Della-Flora, CS Dambros
#### Written by GF Dri and CS Dambros 
#### Contact gabriela.franzoi at maine.edu 
#### Novemeber 2021

# This code calculates observed and expected functional richness at both local and regional scales for all traits combined 

## Load required packages
library(snow)
library(parallel)
library(FD)
library(ecodist)
library(cluster)
library(MuMIn)
library(lme4)
library(ade4)

# Import functional traits table
# Row: species; Column: traits
trait_sp<-read.csv("Data/SpeciesTraitsOK.csv", row.names = 1)
head(trait_sp)

# Import species occurrence table
# Row: communities; Column: species
occ_sp<-read.csv("Data/Occurrence.csv", row.names = 1)
#occ_sp[rowSums(occ_sp)>2,] # just selecting communities with more than 2 species

# Environmental variables table
# Site
env_site<-read.csv("Data/Env_site.csv", row.names = 1)
head(env_site)
nrow(env_site)

# Flock
env_mf<-read.csv("Data/Env_flock.csv", row.names = 1)
head(env_mf)
nrow(env_mf)

####  LOCAL ANALYSIS ####
#### Observed Functional Richness ####
# Distance matrix for functional traits
traits_dist<-daisy(trait_sp,metric = "euclidean")

# Check for euclidean properties
is.euclid(traits_dist)

# Calculate functional richness from Villeger et al. 2008
indicesfd_mf<-dbFD(traits_dist, occ_sp[rowSums(occ_sp)>2,], w.abun = FALSE, stand.x = FALSE, corr = "none", calc.FRic = TRUE, m = "max", stand.FRic = TRUE, scale.RaoQ = FALSE, calc.CWM = FALSE, CWM.type = "all", calc.FDiv = TRUE, print.pco = FALSE, messages = TRUE)

# Aggregate functional diversity indices with environmental variables to exclude NA from flocks that have less than 2 species
indicesmf_env<-cbind(indicesfd_mf, env_mf)
indicesfd_env<-na.omit(indicesmf_env)

#### Null Functional Richness ####
# We create a function "NullFunIndices" to calculate the null functional indices from Villeger et al. 2008 
trait_list<-trait_sp

# Calculate null functional richness
# There is two options: using lapply function or paralleling to improve computational effort 

# Option 1
# Not paralleling the analysis
#NullFunIndices(trait_list,occ_sp,trait_type = c("Q", "Q", "Q", "Q", "Q", "F", "F"),i = 1)
lapply(1:3,NullFunIndices,trait_list=trait_list,occ_sp=occ_sp)

# Option 2
# Paralleling the analysis
# Detect computer cores for parallel runs (use all cores minus 1)
ncores<-detectCores()-1

# Create functions to run null models in parallel
{
  cl<-makeCluster(ncores)
  clusterExport(cl,list("trait_sp","occ_sp"))
  clusterApply(cl,1:99,NullFunIndices,trait_list=trait_list,occ_sp=occ_sp)
  stopCluster(cl)
} 
# Null Functional Richness Results #

# Import the 99 runs
null.files<-list.files("NullResultIndice/",full.names = TRUE)

# Reorganize the imported files
# There is an index that has length 1 (not 355 like all other indices).Therefore, we need to exclude such index
check2<-readLines(null.files[1])

# Check the number of the index that will be removed
toRemove<-sum(c(length(indicesfd_mf$nbsp),
                length(indicesfd_mf$sing.sp),
                length(indicesfd_mf$FRic),
                length(indicesfd_mf$qual.FRic)))

# Removing the index with different length
checknum2<-as.numeric(check2[-toRemove]) 

# Organize the files
checknum3<-lapply(null.files,function(x){
  check2<-readLines(x)
  checknum2<-as.numeric(check2[-toRemove])
})

# Transform into an array
ver<-do.call(c,checknum3)
null.resu<-array(ver,dim=c(length(checknum3[[1]])/7,7,length(checknum3))) # Rows: communities; Columns: functional indices; 3D: replicates 

## Calculating SES for Functional Richness ##
# SES = (Observed value - Mean null value)/Standard deviation of null values

# Difference between observed and null values 
difobsnull<-(indicesfd_mf$FRic-rowMeans(null.resu[,3,],na.rm = TRUE))

# Standard deviation of null values
sdnull<-apply(null.resu[,3,],na.rm = TRUE,1,sd)

# plotting
plot(rnorm(null.resu[,1,1],null.resu[,1,1],0.1),rnorm(indicesfd_mf$FRic,indicesfd_mf$FRic,0.005),pch=21,bg="grey70",col="grey70",las=1,xlab="Species richness",ylab="Corrected functional richness",ylim=c(0,1))
legend(1,1,"Expected",lwd=2,col=1,lty = 1,cex=0.8,bty="n")
#points(null.resu.bmass[,1,1],indicesfd_bmass$FRic,pch=21,bg=2)
lines(ms<-spline(null.resu[,1,1],rowMeans(null.resu[,3,],na.rm = TRUE),n = 12),lwd=2)
#lines(spline(null.resu.diet[,1,1],indicesfd_diet$FRic,n = 7),lwd=3,lty=3,col=2)
m1<-lm(I(indicesfd_mf$FRic-rowMeans(null.resu[,3,],na.rm = TRUE))~null.resu[,1,1])
lines(ms$x,ms$y+(coef(m1)[1]+coef(m1)[2]*ms$x),lwd=2,lty=2,col=2)
legend(1,0.9,"Observed",lwd=2,col=2,lty = 2,cex=0.8,bty="n")

#### REGIONAL ANALYSIS ####

# Aggregate species by site
occ_sp_site_raw<-aggregate(occ_sp[rowSums(occ_sp)>2,],list(env_mf$Localidade),sum)
occ_sp_site<-(occ_sp_site_raw[,-1]>0)*1
rownames(occ_sp_site)<-occ_sp_site_raw[,1]
dim(occ_sp_site)

rownames(trait_sp)==colnames(occ_sp_site)

#### Observed Functional Richness ####

# Calculate functional richness from Villeger et al. 2008
indicesfd_site<-dbFD(traits_dist, occ_sp_site, w.abun = FALSE, stand.x = FALSE, corr = "none", calc.FRic = TRUE, m = "max", stand.FRic = TRUE, scale.RaoQ = FALSE, calc.CWM = FALSE, CWM.type = "all", calc.FDiv = TRUE, print.pco = FALSE, messages = TRUE)


# Aggregate functional diversity indices with environmental variables to exclude NA from flocks that have less than 2 species
indicesfd_site_env<-cbind(indicesfd_site, env_site)
indicesfd_site_env<-na.omit(indicesfd_site_env)
head(indicesfd_site_env)

#### Null Functional Richness ####
# We create a function "NullFunIndices" to calculate the null functional indices from Villeger et al. 2008 
trait_list<-trait_sp

# Calculate null functional richness
# There is two options: using lapply function or paralleling to improve computational effort 

# Option 1
# Not paralleling the analysis
#NullFunIndices(trait_list,occ_sp,trait_type = c("Q", "Q", "Q", "Q", "Q", "F", "F"),i = 1)
lapply(1:3,NullFunIndices,trait_list=trait_list,occ_sp=occ_sp_site)

# Option 2
# Paralleling the analysis
# Detect computer cores for parallel runs (use all cores minus 1)
ncores<-detectCores()-1

# Create functions to run null models in parallel
{
  cl<-makeCluster(ncores)
  clusterExport(cl,list("trait_sp","occ_sp"))
  clusterApply(cl,1:99,NullFunIndices,trait_list=trait_list,occ_sp=occ_sp_site)
  stopCluster(cl)
} 

# Null Functional Richness Results #

# Import the 99 runs
null.files.site<-list.files("NullResultIndiceSite/",full.names = TRUE)

# Reorganize the imported files
# There is an index that has length 1 (not 355 like all other indices).Therefore, we need to exclude such index
check2<-readLines(null.files.site[1])

# Check the number of the index that will be removed
toRemove<-sum(c(length(indicesfd_site$nbsp),
                length(indicesfd_site$sing.sp),
                length(indicesfd_site$FRic),
                length(indicesfd_site$qual.FRic)))

# Removing the index with different length
checknum2<-as.numeric(check2[-toRemove]) 

# Organize the files
checknum3<-lapply(null.files.site,function(x){
  check2<-readLines(x)
  checknum2<-as.numeric(check2[-toRemove])
})

# Transform into an array
ver<-do.call(c,checknum3)
null.resu.siteOK<-array(ver,dim=c(length(checknum3[[1]])/7,7,length(checknum3))) # Rows: communities; Columns: functional indices; 3D: replicates 

## Calculating SES for Functional Richness ##
# SES = (Observed value - Mean null value)/Standard deviation of null values

# Difference between observed and null values 
difobsnullsite<-(indicesfd_site$FRic-rowMeans(null.resu.site[,3,],na.rm = TRUE))

# Standard deviation of null values
sdnullsite<-apply(null.resu.site[,3,],na.rm = TRUE,1,sd)

# Plotting
plot(rnorm(null.resu.site[,1,1],null.resu.site[,1,1],0.1),rnorm(indicesfd_site$FRic,indicesfd_site$FRic,0.005),pch=21,bg="grey70",col="grey70",las=1,xlab="Species richness",ylab="Corrected functional richness",ylim=c(0,1))
legend(16,1,"Expected",lwd=2,col=1,lty = 1,cex=0.8,bty="n")
#points(null.resu.bmass[,1,1],indicesfd_bmass$FRic,pch=21,bg=2)
lines(ms<-spline(null.resu.site[,1,1],rowMeans(null.resu.site[,3,],na.rm = TRUE),n = 4),lwd=2)
m1<-lm(I(indicesfd_site$FRic-rowMeans(null.resu.site[,3,],na.rm = TRUE))~null.resu.site[,1,1])
lines(ms$x,ms$y+(coef(m1)[1]+coef(m1)[2]*ms$x),lwd=2,lty=2,col=2)
legend(16,0.9,"Observed",lwd=2,col=2,lty = 2,cex=0.8,bty="n")

#
## End