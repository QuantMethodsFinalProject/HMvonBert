#rm(list=ls())

#### Load libraries
library(mgcv)
library(MCMCpack) # rwish function
library(R2jags)


#read in data
setwd("S:/FlatheadOtolithProject/Data/Rprojects/MetaAnalysis/HMvonBert")
dat <- read.csv("FHCGrowthMetaData.csv")

#remove prob pop's
dat<-dat[dat$River!="North Raccoon", ]
#length(unique(dat$River)) #check J=24
dat<-dat[dat$River!="Lake Lyndon B. Johnson", ]

#remove old data from analysis
dat <- subset(dat, Year >= 2000) 
sort(unique(dat$Year))


#now lets grab only river pops
riv<- subset(dat, Waterbody=="River")
length(unique(riv$River)) #should be 17
riv <- droplevels(riv) #drop before ordering!

#river name column
riv <- riv[order(riv$River), ]
riv$River_name <- riv$River
riv$River <- as.numeric(riv$River)


#final check:
unique(riv$River) #should be in a consecutive numeric order
unique(riv$River_name) #should not include removed pops above AS WELL AS Lake pops

###################################ADD PREDICTORS
#stream order
hist(unique(riv$StreamOrder))#variety isnt ideal
so <- as.numeric(by(riv$StreamOrder, riv$River, mean)) 
z_so <- as.numeric(scale(so))

# Categorical predictor ('Native' is the reference cell, which means is is set to zero)
riv$Status_binary <- ifelse(riv$Status == 'introduced', 1, 0)
status <- as.numeric(by(riv$Status_binary, riv$River, mean)) 
#######################################
sink("MetaAnalysisHMvonBmodel_River_pred.txt")
cat("
    model{
    for(i in 1:n){
    y[i] ~ dnorm(y.hat[i], tau.y)
    y.hat[i] <- Linf[g[i]] * (1-exp(-k[g[i]] * (age[i] - t0[g[i]] )))
    }
    
    tau.y <- pow(sigma.y,-2)
    sigma.y ~ dunif(0,100)
    
    # Parameters modeled on log-scale
    for(j in 1:J){
    Linf[j] <- exp(BB[j,1])
    k[j] <- exp(BB[j,2])
    t0[j] <- exp(BB[j,3])-10 # A constant of 10 is added (subtracted) to t0 to ensure that negative values are possible, becuase t0 is estimated on log-scale
    BB[j,1:K] ~ dmnorm (BB.hat[j,], Tau.B[,]) # Multivariate normal dist'n;  Tau.B is a precision 
    BB.hat[j,1] <- mu.Linf + Lgamma1 * status[j]+ Lgamma2 * z_so[j] + Lgamma3 * z_so[j] * status[j]
    BB.hat[j,2] <- mu.k + Kgamma1 * status[j] + Kgamma2 * z_so[j] + Kgamma3 * z_so[j] * status[j]
    BB.hat[j,3] <- mu.t0
    }
    
    # Priors for population-average parameters
    mu.Linf ~ dnorm(0,.0001)
    mu.k ~ dnorm(0,.0001)
    mu.t0 ~ dnorm(0,.0001)
    
    Lgamma1 ~ dnorm(0,.0001)
    Kgamma1 ~ dnorm(0,.0001)
    
    Lgamma2 ~ dnorm(0,.0001)
    Kgamma2 ~ dnorm(0,.0001)
    
    Lgamma3 ~ dnorm(0,.0001)
    Kgamma3 ~ dnorm(0,.0001)
    
    # Model variance-covariance
    Tau.B[1:K,1:K] ~ dwish(W[,], df)
    df <- K+1
    Sigma.B[1:K,1:K] <- inverse(Tau.B[,])
    for (k in 1:K){
    for (k.prime in 1:K){
    rho.B[k,k.prime] <- Sigma.B[k,k.prime]/
    sqrt(Sigma.B[k,k]*Sigma.B[k.prime,k.prime])
    }
    sigma.B[k] <- sqrt(Sigma.B[k,k])
    }
    
    } # end model
    ",fill=TRUE)
sink()

################
J <- length(unique(riv$River))


# Create identity matrix for Wishart dist'n
#!!Number of varying parameters to estimate (K)
K <- 3

# Create identity matrix for Wishart dist'n
W <- diag(3)

# load data
data <- list(y = riv$TL, age = riv$Age, g = riv$River, n = dim(riv)[1],
             J = J, W=W, K=K, status=status , z_so=z_so )

# Initial values
inits <- function(){list(mu.Linf = rnorm(1,3,0.001), mu.k = rnorm(1,1,0.001), mu.t0 = rnorm(1,0.7,0.001),Lgamma1=rnorm(1),Kgamma1=rnorm(1), 
                         sigma.y = runif(1,1,10), Lgamma2=rnorm(1),Kgamma2=rnorm(1), Lgamma3=rnorm(1),Kgamma3=rnorm(1), 
                         BB=array(c(rep(log(950) +rnorm(1,0.01,0.01),J),rep(log(0.04)+rnorm(1,0.001,0.1),J),rep(log(-2+10)+rnorm(1,0.01,0.1),J)),
                                  c(J,K)), Tau.B=rwish(K+1,diag(K)) ) }

# Parameters monitored
params1 <- c("mu.Linf", "mu.k", "mu.t0", "sigma.y","BB","Sigma.B","rho.B",
             "Lgamma1","Kgamma1", "Lgamma2","Kgamma2", "Lgamma3","Kgamma3")


# MCMC settings
ni <- 100000
nt <- 3
nb <- 50000
nc <- 3
############################################
start.time = Sys.time()  

riv.pred.out <- jags(data = data, inits = inits, parameters.to.save = params1, 
                 model.file = "MetaAnalysisHMvonBmodel_River_pred.txt", n.chains = nc, 
                 n.thin = nt, n.iter = ni, n.burnin = nb)

#Calculate computation time
end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 
##############################
#Saving the JAGS output
#DO NOT RUN "SAVE" BY ACCIDENT IF NOT IN YOUR ENVIRONMENT, WILL HAVE TO RE-RUN MODEL IF SO 
#saveRDS(riv.pred.out, compress=TRUE, "S:/FlatheadOtolithProject/Data/Rprojects/MetaAnalysis/HMvonBert/JAGS.rds/riv.pred.out.rds")

#remove from envirn when not using
#rm(riv.pred.out)

#read in model output
riv.pred.out <- readRDS("S:/FlatheadOtolithProject/Data/Rprojects/MetaAnalysis/HMvonBert/JAGS.rds/riv.pred.out.rds")
identical(riv.pred.out, riv.pred.out, ignore.environment = TRUE)
##############################
print(riv.pred.out, digits = 2)

# Find which parameters, if any, have Rhat > 1.1
which(riv.pred.out$BUGSoutput$summary[, c("Rhat")] > 1.1) #three increase iterations

traceplot(riv.pred.out)
#####################
############Relationships
# Does Linf and K differ between native and introduced?
#Linf
L_natslope <- riv.pred.out$BUGSoutput$sims.list$Lgamma2 #ref cell
L_introslope <- riv.pred.out$BUGSoutput$sims.list$Lgamma2 + riv.pred.out$BUGSoutput$sims.list$Lgamma3
L_Diffslope <- L_natslope - L_introslope 
quantile(L_Diffslope, c(0.025, 0.975))
quantile(L_Diffslope, c(0.05, 0.95)) #overlaps zero, so no difference
mean(L_introslope)
mean(L_natslope)

#K coef
K_natslope <- riv.pred.out$BUGSoutput$sims.list$Kgamma2 #ref cell
K_introslope <- riv.pred.out$BUGSoutput$sims.list$Lgamma2 +  riv.pred.out$BUGSoutput$sims.list$Lgamma3
K_Diffslope <- K_natslope - K_introslope 
quantile(K_Diffslope, c(0.025, 0.975))
quantile(K_Diffslope, c(0.05, 0.95)) #overlaps zero, so no difference
mean(K_introslope)
mean(K_natslope)
###############
# Does Linf and K differ by stream ordder?
#Linf
#NATIVE POPS
mean(L_natslope)
# What is the probability that the effect of stream order is negative
mean(L_natslope < 0)
quantile(L_natslope,c(0.05, 0.95)) #does it overlap zero? YES= NOT different
quantile(L_natslope,c(0.025, 0.975)) 
#INTRO POPS
mean(L_introslope)
# What is the probability that the effect of stream order is negative
mean(L_introslope < 0)
quantile(L_introslope,c(0.05, 0.95))# does it overlap zero? YES= NOT different
quantile(L_introslope,c(0.025, 0.975))

#K coef
#NATIVE POPS
mean(K_natslope)
mean(K_natslope < 0)
quantile(K_natslope,c(0.025, 0.975))
quantile(K_natslope,c(0.05, 0.95)) # does it overlap zero? YES= NOT different
#INTRO POPS
mean(K_introslope)
mean(K_introslope <0)
quantile(K_introslope,c(0.025, 0.975))
quantile(K_introslope,c(0.05, 0.95)) # does it overlap zero? YES= NOT different

########################################################################################################
##########################################PLOTS#########################################################
######### plot stream order by mean growth parameters
#############calcs for k
#separate z_so by status
natives <- subset(riv, Status=="native")
nat_so <- as.numeric(by(natives$StreamOrder, natives$River, mean)) 
z_natso <- as.numeric(scale(nat_so))

intros <- subset(riv, Status=="introduced")
intro_so <- as.numeric(by(intros$StreamOrder, intros$River, mean)) 
z_introso <- as.numeric(scale(intro_so)) 

# Select random slopes  
meanK.nat <- riv.pred.out$BUGSoutput$mean$BB[c(3:6,8,10,17),2] #select only native pops
meanK.intro <- riv.pred.out$BUGSoutput$mean$BB[c(1:2, 7,9,11:16),2] #select only intro pops
meanK.all <- riv.pred.out$BUGSoutput$mean$BB[,2]

#fake pred values
fake.nat <- seq(min(z_natso), max(z_natso), length=50)
fake.intro <- seq(min(z_introso), max(z_introso), length=50)

# Create container
K_est.lineNat <- matrix(NA, ncol=length(fake.nat), nrow=riv.pred.out$BUGSoutput$n.sims)
K_est.lineIntro <- matrix(NA, ncol=length(fake.intro), nrow=riv.pred.out$BUGSoutput$n.sims)

# calc predicted values for natives
for(i in 1:riv.pred.out$BUGSoutput$n.sims){
  for(j in 1:length(fake.nat) ){
    K_est.lineNat[i, j] <- riv.pred.out$BUGSoutput$sims.list$mu.k[i] +  riv.pred.out$BUGSoutput$sims.list$Kgamma2[i] * fake.nat[j] 
  }
}


# Derive introduced pops ints and slopes
K_IntroInt <- riv.pred.out$BUGSoutput$sims.list$mu.k + riv.pred.out$BUGSoutput$sims.list$Kgamma1
K_IntroSlope <- riv.pred.out$BUGSoutput$sims.list$Kgamma2 + riv.pred.out$BUGSoutput$sims.list$Kgamma3

#calc predicted values for intros
for(i in 1:riv.pred.out$BUGSoutput$n.sims){
  for(j in 1:length(fake.intro) ){
    K_est.lineIntro[i, j] <- K_IntroInt[i] + K_IntroSlope[i] * fake.intro[j]
  }
}

### Obtain posterior mean fitted line and associated 
# 95% CRIs (upper and lower CRIs for predicted values in the matrix est.line
#Natives
K_fitted.mean1 <-  apply(K_est.lineNat, 2, mean )
K_upper.CRI1 <- apply(K_est.lineNat, 2, quantile, probs=c(0.975))
K_lower.CRI1 <- apply(K_est.lineNat, 2, quantile, probs=c(0.025))
#Intro
K_fitted.mean2 <-  apply(K_est.lineIntro, 2, mean )
K_upper.CRI2 <- apply(K_est.lineIntro, 2, quantile, probs=c(0.975))
K_lower.CRI2 <- apply(K_est.lineIntro, 2, quantile, probs=c(0.025))


## Grab 95% CIs for slopes of each status
#natives
natpops <- c(3:6,8,10,17) #remember: these river #'s do not match up to table1 (two pops removed)
u.Knat <- numeric(length(meanK.nat))
l.Knat <- numeric(length(meanK.nat))
for(i in 1:length(meanK.nat)){ 
  u.Knat[i] <- quantile(riv.pred.out$BUGSoutput$sims.list$BB[,natpops[i],2],probs=c(0.975))
  l.Knat[i] <-quantile(riv.pred.out$BUGSoutput$sims.list$BB[,natpops[i],2], probs=c(0.025))
}

#intros
intropops <- c(1:2, 7,9,11:16) #remember: these river #'s do not match up to table1 (two pops removed)
u.Kintro <- numeric(length(meanK.intro))
l.Kintro <- numeric(length(meanK.intro))
for(i in 1:length(meanK.intro)){ 
  u.Kintro[i] <- quantile(riv.pred.out$BUGSoutput$sims.list$BB[,intropops[i],2],probs=c(0.975))
  l.Kintro[i] <-quantile(riv.pred.out$BUGSoutput$sims.list$BB[,intropops[i],2], probs=c(0.025))
}

#############calcs for Linf
#separate z_so by status 
#already done with k
# natives <- subset(riv, Status=="native")
# nat_so <- as.numeric(by(natives$StreamOrder, natives$River, mean)) 
# z_natso <- as.numeric(scale(nat_so))
# 
# intros <- subset(riv, Status=="introduced")
# intro_so <- as.numeric(by(intros$StreamOrder, intros$River, mean)) 
# z_introso <- as.numeric(scale(intro_so)) 

# Select random slopes  
meanLinf.nat <- riv.pred.out$BUGSoutput$mean$BB[c(3:6,8,10,17),1] #select only native pops
meanLinf.intro <- riv.pred.out$BUGSoutput$mean$BB[c(1:2, 7,9,11:16),1] #select only intro pops
meanLinf.all <- riv.pred.out$BUGSoutput$mean$BB[,1]

#fake pred values
#already done with k
# fake.nat <- seq(min(z_natso), max(z_natso), length=50)
# fake.intro <- seq(min(z_introso), max(z_introso), length=50)

# Create container
Linf_est.lineNat <- matrix(NA, ncol=length(fake.nat), nrow=riv.pred.out$BUGSoutput$n.sims)
Linf_est.lineIntro <- matrix(NA, ncol=length(fake.intro), nrow=riv.pred.out$BUGSoutput$n.sims)

# calc predicted values for natives
for(i in 1:riv.pred.out$BUGSoutput$n.sims){
  for(j in 1:length(fake.nat) ){
    Linf_est.lineNat[i, j] <- riv.pred.out$BUGSoutput$sims.list$mu.Linf[i] +  riv.pred.out$BUGSoutput$sims.list$Lgamma2[i] * fake.nat[j] 
  }
}


# Derive introduced pops ints and slopes
Linf_IntroInt <- riv.pred.out$BUGSoutput$sims.list$mu.Linf + riv.pred.out$BUGSoutput$sims.list$Lgamma1
Linf_IntroSlope <- riv.pred.out$BUGSoutput$sims.list$Lgamma2 + riv.pred.out$BUGSoutput$sims.list$Lgamma3

#calc predicted values for intros
for(i in 1:riv.pred.out$BUGSoutput$n.sims){
  for(j in 1:length(fake.intro) ){
    Linf_est.lineIntro[i, j] <- Linf_IntroInt[i] + Linf_IntroSlope[i] * fake.intro[j]
  }
}

### Obtain posterior mean fitted line and associated 
# 95% CRIs (upper and lower CRIs for predicted values in the matrix est.line
#Natives
Linf_fitted.mean1 <-  apply(Linf_est.lineNat, 2, mean )
Linf_upper.CRI1 <- apply(Linf_est.lineNat, 2, quantile, probs=c(0.975))
Linf_lower.CRI1 <- apply(Linf_est.lineNat, 2, quantile, probs=c(0.025))
#Intro
Linf_fitted.mean2 <-  apply(Linf_est.lineIntro, 2, mean )
Linf_upper.CRI2 <- apply(Linf_est.lineIntro, 2, quantile, probs=c(0.975))
Linf_lower.CRI2 <- apply(Linf_est.lineIntro, 2, quantile, probs=c(0.025))


## Grab 95% CIs for slopes of each status
#natives
natpops <- c(3:6,8,10,17) #remember: these river #'s do not match up to table1 (two pops removed)
u.Linfnat <- numeric(length(meanLinf.nat))
l.Linfnat <- numeric(length(meanLinf.nat))
for(i in 1:length(meanLinf.nat)){ 
  u.Linfnat[i] <- quantile(riv.pred.out$BUGSoutput$sims.list$BB[,natpops[i],1],probs=c(0.975))
  l.Linfnat[i] <-quantile(riv.pred.out$BUGSoutput$sims.list$BB[,natpops[i],1], probs=c(0.025))
}

#intros
intropops <- c(1:2, 7,9,11:16) #remember: these river #'s do not match up to table1 (two pops removed)
u.Linfintro <- numeric(length(meanLinf.intro))
l.Linfintro <- numeric(length(meanLinf.intro))
for(i in 1:length(meanLinf.intro)){ 
  u.Linfintro[i] <- quantile(riv.pred.out$BUGSoutput$sims.list$BB[,intropops[i],1],probs=c(0.975))
  l.Linfintro[i] <-quantile(riv.pred.out$BUGSoutput$sims.list$BB[,intropops[i],1], probs=c(0.025))
}
#################################################################################
######################CREATE PANEL PLOT##########################################
#################################################################################
res <- 6
name_figure <- "Riv_StreamOrder_GrowthParams.png"
jpeg(filename = name_figure, height =500*res, width = 500*res, res=72*res)
def.par <- par(no.readonly = TRUE) 		# save default, for resetting...

par(mar=c(0.0,0.1,0.1,0.1),oma=c(3,4,0,1),mai=c(0.15,0.15,0.15,0.15) )

nf <- layout(matrix(c(1:4),nrow=2,ncol=2,byrow=TRUE),  TRUE) 
layout.show(nf)

size.labels = 1
size.text = 1.0

x.label <- 'Standarized Stream Order'
y.label <- expression(paste(log[e],'(Growth Coefficent (K))'))

x.label2 <- 'Standarized Stream Order'
y.label2 <- expression(paste(log[e],'(L'[infinity],")"))


# Add transparency to points
source('AddTransFunction.R')

#K plots
#natives
plot(meanK.nat ~ z_natso, data=dat, type="n", axes=F,xlab='',ylab='', ylim=c(-2.5, -1.0) ,xlim=c(-2.5,2))

points(z_natso, meanK.nat, pch=21, cex=1.5, col='black',bg='blue')

#add segments
segments(x0= z_natso, x1= z_natso, y0=l.Knat, y1=u.Knat, col='black',lwd=1)

axis(side=1,cex.axis=size.text, tck=-0.01, mgp=c(3,0.05,0),labels=F )
axis(side=2,cex.axis=size.text, las=1)

# Add fitted line
lines(fake.nat,K_fitted.mean1, lwd = 3, col="blue", lty = 1)
# Add credible region
#natives
i.for <- order(fake.nat)
i.back <- order(fake.nat, decreasing = TRUE )
x.polygon <- c( fake.nat[i.for] , fake.nat[i.back] )
y.polygon <- c( K_lower.CRI1[i.for] , K_upper.CRI1[i.back] )
polygon( x.polygon , y.polygon , col = addTrans("blue",100)  , border = NA )
#axis labs
mtext(y.label, line = 3, side = 2, cex = size.text)
box()

#intros
plot(meanK.intro ~ z_introso, data=dat, type="n", axes=F,xlab='',ylab='',ylim=c(-2.5, -1.0) ,xlim=c(-2.5,2) )

points(z_introso, meanK.intro, pch=21, cex=1.5, col='black',bg='gray')

#add segments
segments(x0= z_introso, x1= z_introso, y0=l.Kintro, y1=u.Kintro, col='black',lwd=1)

axis(side=1,cex.axis=size.text, tck=-0.01, mgp=c(3,0.05,0),labels=F )
axis(side=2,cex.axis=size.text, las=1, labels=F)

# Add fitted line
lines(fake.intro,K_fitted.mean2, lwd = 3, col="gray", lty = 1)
# Add credible region
i.for <- order(fake.intro)
i.back <- order(fake.intro, decreasing = TRUE )
x.polygon <- c( fake.intro[i.for] , fake.intro[i.back] )
y.polygon <- c( K_lower.CRI2[i.for] , K_upper.CRI2[i.back] )
polygon( x.polygon , y.polygon , col = addTrans("gray",100)  , border = NA )
#axis labs
#none
box()


####Linf plots
#natives
plot(meanLinf.nat ~ z_natso, data=dat, type="n", axes=F,xlab='',ylab='' , xlim=c(-2.5,2), ylim=c(6.2,7.5))

points(z_natso, meanLinf.nat, pch=21, cex=1.5, col='black',bg='blue')

#add segments
segments(x0= z_natso, x1= z_natso, y0=l.Linfnat, y1=u.Linfnat, col='black',lwd=1)


axis(side=1,cex.axis=size.text, tck=-0.01, mgp=c(3,0.05,0))
axis(side=2,cex.axis=size.text, las=1)

# Add fitted line
lines(fake.nat,Linf_fitted.mean1, lwd = 3, col="blue", lty = 1)
# Add credible region
i.for <- order(fake.nat)
i.back <- order(fake.nat, decreasing = TRUE )
x.polygon <- c( fake.nat[i.for] , fake.nat[i.back] )
y.polygon <- c( Linf_lower.CRI1[i.for] , Linf_upper.CRI1[i.back] )
polygon( x.polygon , y.polygon , col = addTrans("blue",100)  , border = NA )
#axis labs
mtext(y.label2, line = 3, side = 2, cex = size.text)
box()

#intros
plot(meanLinf.intro ~ z_introso, data=dat, type="n", axes=F,xlab='',ylab="", xlim=c(-2.5,2), ylim=c(6.2,7.5))

points(z_introso, meanLinf.intro, pch=21, cex=1.5, col='black',bg='gray') 
#add segments
segments(x0= z_introso, x1= z_introso, y0=l.Linfintro, y1=u.Linfintro, col='black',lwd=1)


axis(side=1,cex.axis=size.text, tck=-0.01, mgp=c(3,0.05,0))
axis(side=2,cex.axis=size.text, las=1, labels=F)

# Add fitted line
lines(fake.intro,Linf_fitted.mean2, lwd = 3, col="gray", lty = 1)
# Add credible region
i.for <- order(fake.intro)
i.back <- order(fake.intro, decreasing = TRUE )
x.polygon <- c( fake.intro[i.for] , fake.intro[i.back] )
y.polygon <- c( Linf_lower.CRI2[i.for] , Linf_upper.CRI2[i.back] )
polygon( x.polygon , y.polygon , col = addTrans("gray",100)  , border = NA )
#axis labs
#none
box()

mtext(x.label , line = 1, side = 1, cex = size.text, outer = TRUE)
#####################
par(def.par)
dev.off()
