#### Code that follows the paper Mixed-species bird flocks enhance the benefits of group aggregation by minimizing variation in some functional traits while maximizing variation in others by GF Dri, NC Caceres, J Della-Flora, CS Dambros
#### Written by GF Dri and CS Dambros 
#### Contact gabriela.franzoi at maine.edu 
#### November 2021

# This code creates the plots used in the paper

## You will need objects created in the previous codes (1_Functional_analysis_traitscombined; 2_Functional_analysis_traitsindividual; 3_Null_functional_analyses_local and 4_Null_functional_analyses_regional)

### All traits combined ####
## Local
plot(rnorm(null.resu[,1,1],null.resu[,1,1],0.1),rnorm(indicesfd_mf$FRic,indicesfd_mf$FRic,0.005),pch=21,bg="grey70",col="grey70",las=1,xlab="Species richness",ylab="Corrected functional richness",ylim=c(0,1))
legend(1,1,"Expected",lwd=2,col=1,lty = 1,cex=0.8,bty="n")
#points(null.resu.bmass[,1,1],indicesfd_bmass$FRic,pch=21,bg=2)
lines(ms<-spline(null.resu[,1,1],rowMeans(null.resu[,3,],na.rm = TRUE),n = 12),lwd=2)
#lines(spline(null.resu.diet[,1,1],indicesfd_diet$FRic,n = 7),lwd=3,lty=3,col=2)
m1<-lm(I(indicesfd_mf$FRic-rowMeans(null.resu[,3,],na.rm = TRUE))~null.resu[,1,1])
lines(ms$x,ms$y+(coef(m1)[1]+coef(m1)[2]*ms$x),lwd=2,lty=2,col=2)
legend(1,0.9,"Observed",lwd=2,col=2,lty = 2,cex=0.8,bty="n")

## Regional
plot(rnorm(null.resu.site[,1,1],null.resu.site[,1,1],0.1),rnorm(indicesfd_site$FRic,indicesfd_site$FRic,0.005),pch=21,bg="grey70",col="grey70",las=1,xlab="Species richness",ylab="Corrected functional richness",ylim=c(0,1))
legend(16,1,"Expected",lwd=2,col=1,lty = 1,cex=0.8,bty="n")
#points(null.resu.bmass[,1,1],indicesfd_bmass$FRic,pch=21,bg=2)
lines(ms<-spline(null.resu.site[,1,1],rowMeans(null.resu.site[,3,],na.rm = TRUE),n = 4),lwd=2)
m1<-lm(I(indicesfd_site$FRic-rowMeans(null.resu.site[,3,],na.rm = TRUE))~null.resu.site[,1,1])
lines(ms$x,ms$y+(coef(m1)[1]+coef(m1)[2]*ms$x),lwd=2,lty=2,col=2)
legend(16,0.9,"Observed",lwd=2,col=2,lty = 2,cex=0.8,bty="n")

#### DIET ####
## Local
plot(rnorm(null.resu.diet[,1,1],null.resu.diet[,1,1],0.1),rnorm(indicesfd_diet$FRic,indicesfd_diet$FRic,0.005),pch=21,bg="grey70",col="grey70",las=1,xlab="Species richness",ylab="Diet amplitude (FRic for diet)")
legend(20,0.25,c("Expected","Observed"),lwd=c(2,2),col=c(1,2),lty = c(1,2),bty="n",cex=0.8)
#points(null.resu.bmass[,1,1],indicesfd_bmass$FRic,pch=21,bg=2)
lines(ms<-spline(null.resu.diet[,1,1],rowMeans(null.resu.diet[,3,],na.rm = TRUE),n = 12),lwd=2)
#lines(spline(null.resu.diet[,1,1],indicesfd_diet$FRic,n = 7),lwd=3,lty=3,col=2)
m1<-lm(I(indicesfd_diet$FRic-rowMeans(null.resu.diet[,3,],na.rm = TRUE))~null.resu.diet[,1,1])
lines(ms$x,ms$y+(coef(m1)[1]+coef(m1)[2]*ms$x),lwd=2,lty=2,col=2)

## Regional
plot(rnorm(null.resu.diet.site[,1,1],null.resu.diet.site[,1,1],0.1),rnorm(indicesfd_site_diet$FRic,indicesfd_site_diet$FRic,0.005),pch=21,bg="grey70",col="grey70",las=1,xlab="Species richness",ylab="Diet amplitude (FRic for diet)",ylim=c(0,1))
legend(50,0.25,"Expected",lwd=2,col=1,lty =1,bty="n",cex=0.8)
#points(null.resu.bmass[,1,1],indicesfd_bmass$FRic,pch=21,bg=2)
lines(ms<-spline(null.resu.diet.site[,1,1],rowMeans(null.resu.diet.site[,3,],na.rm = TRUE),n = 4),lwd=2)

#### TAIL ####
## Local
plot(rnorm(null.resu.retrizes[,1,1],null.resu.retrizes[,1,1],0.1),rnorm(indicesfd_retrizes$FRic,indicesfd_retrizes$FRic,0.005),pch=21,bg="grey70",col="grey70",las=1,xlab="Species richness",ylab="retrizes amplitude (FRic for retrizes)")
legend(21,0.24,c("Expected","Observed"),lwd=c(2,2),col=c(1,2),lty = c(1,2),cex=0.8,bty="n")
#points(null.resu.bmass[,1,1],indicesfd_bmass$FRic,pch=21,bg=2)
lines(ms<-spline(null.resu.retrizes[,1,1],rowMeans(null.resu.retrizes[,3,],na.rm = TRUE),n = 12),lwd=2)
#lines(spline(null.resu.retrizes[,1,1],indicesfd_retrizes$FRic,n = 7),lwd=3,lty=3,col=2)
m1<-lm(I(indicesfd_retrizes$FRic-rowMeans(null.resu.retrizes[,3,],na.rm = TRUE))~null.resu.retrizes[,1,1])
lines(ms$x,ms$y+(coef(m1)[1]+coef(m1)[2]*ms$x),lwd=2,lty=2,col=2)

## Regional
plot(rnorm(null.resu.retrizes.site[,1,1],null.resu.retrizes.site[,1,1],0.1),rnorm(indicesfd_site_retrizes$FRic,indicesfd_site_retrizes$FRic,0.005),pch=21,bg="grey70",col="grey70",las=1,xlab="Species richness",ylab="retrizes amplitude (FRic for retrizes)",ylim=c(0,1))
legend(53,0.2,"Expected",lwd=2,col=1,lty =1,bty="n",cex=0.8)
#points(null.resu.bmass[,1,1],indicesfd_bmass$FRic,pch=21,bg=2)
lines(ms<-spline(null.resu.retrizes.site[,1,1],rowMeans(null.resu.retrizes.site[,3,],na.rm = TRUE),n = 4),lwd=2)

#### BODY MASS ####
## Local
plot(rnorm(null.resu.bmass[,1,1],null.resu.bmass[,1,1],0.1),rnorm(indicesfd_bmass$FRic,indicesfd_bmass$FRic,0.005),pch=21,bg="grey70",col="grey70",xlab="Species richness",ylab="Body mass amplitude (FRic for body mass)",las=1)
legend(21,0.23,c("Expected","Observed"),lwd=c(2,2),col=c(1,2),lty = c(1,2),bty="n",cex=0.8)
#points(null.resu.bmass[,1,1],indicesfd_bmass$FRic,pch=21,bg=2)
lines(ms<-spline(null.resu.bmass[,1,1],rowMeans(null.resu.bmass[,3,],na.rm = TRUE),n = 12),lwd=2)
#lines(spline(null.resu.diet[,1,1],indicesfd_diet$FRic,n = 7),lwd=3,lty=3,col=2)
m1<-lm(I(indicesfd_bmass$FRic-rowMeans(null.resu.bmass[,3,],na.rm = TRUE))~null.resu.bmass[,1,1])
lines(ms$x,ms$y+(coef(m1)[1]+coef(m1)[2]*ms$x),lwd=2,lty=2,col=2)

## Regional
plot(rnorm(null.resu.bmass.site[,1,1],null.resu.bmass.site[,1,1],0.1),rnorm(indicesfd_site_bmass$FRic,indicesfd_site_bmass$FRic,0.005),pch=21,bg="grey70",col="grey70",las=1,xlab="Species richness",ylab="Body mass amplitude (FRic for body mass)",ylim=c(0,1))
legend(55,0.12,"Expected",lwd=2,col=1,lty =1,bty="n",cex=0.8)
#points(null.resu.bmass[,1,1],indicesfd_bmass$FRic,pch=21,bg=2)
lines(ms<-spline(null.resu.bmass.site[,1,1],rowMeans(null.resu.bmass.site[,3,],na.rm = TRUE),n = 4),lwd=2)
#lines(spline(null.resu.diet[,1,1],indicesfd_diet$FRic,n = 7),lwd=3,lty=3,col=2)

#
## End