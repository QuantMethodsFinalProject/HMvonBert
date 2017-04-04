#rm(list=ls())


#### Load libraries
library(mgcv)
library(MCMCpack) # rwish function
library(R2jags)


#read in data
setwd("S:/FlatheadOtolithProject/Data/Rprojects/MetaAnalysis/HMvonBert")
dat <- read.csv("FHCGrowthMetaData.csv")

#remove old data from analysis
dat <- subset(dat, Year >= 2000) 
sort(unique(dat$Year)) #check


#only interested in intro pops
intro <-subset(dat, Status=="introduced")
#removing lakes from intro pops did not change relationship, so keeping them in model
# intro2<-intro[intro$River!="Lake Marion", ]
# intro2<-intro2[intro2$River!="Lake Moultrie", ]


#must drop levels before continuing
intro <- droplevels(intro)
# intro2 <- droplevels(intro2)


#river name column
intro <- intro[order(intro$River), ]
intro$River_name <- intro$River
intro$River <- as.numeric(intro$River)

#final check:
length(unique(intro$River)) #should be 12 for intro pops
unique(intro$River) #should be in a consecutive numeric order
unique(intro$River_name) #does not include removed pops above or NATIVES
sort(unique(intro$TimeSinceEstab))
mean(intro$TimeSinceEstab)
###################################ADD PREDICTORS
# Site-level predictors
estab <- as.numeric(by(intro$TimeSinceEstab, intro$River, mean)) 
z_estab <- as.numeric(scale(estab))

#################################################################
########## BUGS CODE ############################################
#################################################################
sink("MetaAnalysisHMvonBmodel_intro_pred.txt")
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
    BB.hat[j,1] <- mu.Linf + Lgamma1 * time[j] 
    BB.hat[j,2] <- mu.k + Kgamma1 * time[j] 
    BB.hat[j,3] <- mu.t0
    }
    
    # Priors for population-average parameters
    mu.Linf ~ dnorm(0,.0001)
    mu.k ~ dnorm(0,.0001)
    mu.t0 ~ dnorm(0,.0001)
    Lgamma1 ~ dnorm(0,.0001)
    Kgamma1 ~ dnorm(0,.0001)
    
    
    
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
J <- length(unique(intro$River))


# Create identity matrix for Wishart dist'n
#!!Number of varying parameters to estimate (K)
K <- 3

# Create identity matrix for Wishart dist'n
W <- diag(3)

# load data
data <- list(y = intro$TL, age = intro$Age, g = intro$River, n = dim(intro)[1],
             J = J, W=W, K=K, time=z_estab)

# Initial values
inits <- function(){list(mu.Linf = rnorm(1,3,0.001), mu.k = rnorm(1,1,0.001), mu.t0 = rnorm(1,0.7,0.001),Lgamma1=rnorm(1),Kgamma1=rnorm(1), 
                         sigma.y = runif(1,1,10),  
                         BB=array(c(rep(log(950) +rnorm(1,0.01,0.01),J),rep(log(0.04)+rnorm(1,0.001,0.1),J),rep(log(-2+10)+rnorm(1,0.01,0.1),J)),
                                  c(J,K)), Tau.B=rwish(K+1,diag(K)) ) }

# Parameters monitored
params1 <- c("mu.Linf", "mu.k", "mu.t0", "sigma.y","BB","Sigma.B","rho.B",
             "Lgamma1","Kgamma1" )


# MCMC settings
ni <- 100000
nt <- 3
nb <- 50000
nc <- 3


############################################
start.time = Sys.time()  

intro.pred.out <- jags(data = data, inits = inits, parameters.to.save = params1, 
            model.file = "MetaAnalysisHMvonBmodel_intro_pred.txt", n.chains = nc, 
            n.thin = nt, n.iter = ni, n.burnin = nb)

#Calculate computation time
end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 
##############################
#Saving the JAGS output
#DO NOT RUN "SAVE" BY ACCIDENT IF NOT IN YOUR ENVIRONMENT, WILL HAVE TO RE-RUN MODEL IF SO 

#saveRDS(intro.pred.out, compress=TRUE, "S:/FlatheadOtolithProject/Data/Rprojects/MetaAnalysis/HMvonBert/JAGS.rds/intro.pred.out.rds")

#remove from envirn when not using
#rm(intro.pred.out) 

#read in model output
intro.pred.out <- readRDS("S:/FlatheadOtolithProject/Data/Rprojects/MetaAnalysis/HMvonBert/JAGS.rds/intro.pred.out.rds")
identical(intro.pred.out, intro.pred.out, ignore.environment = TRUE)
############################
# Summarize the result
print(intro.pred.out, digits = 2)
# str(intro.pred.out)

# Find which parameters, if any, have Rhat > 1.1
which(intro.pred.out$BUGSoutput$summary[, c("Rhat")] > 1.1) #three increase iterations

#check traceplots
traceplot(intro.pred.out)

########################
# Does time since estab effect Linf and k in intro pops?
#Linf
quantile(intro.pred.out$BUGSoutput$sims.list$Lgamma1, c(0.025, 0.975))
quantile(intro.pred.out$BUGSoutput$sims.list$Lgamma1, c(0.05, 0.95))
#overlaps zero, so no difference


#K coef
quantile(intro.pred.out$BUGSoutput$sims.list$Kgamma1, c(0.025, 0.975))
quantile(intro.pred.out$BUGSoutput$sims.list$Kgamma1, c(0.05, 0.95))
#overlaps zero, so no difference



# Prob that slope is positive
#prob that param is in the direction of posterior mean estimate
#Linf
mean(intro.pred.out$BUGSoutput$sims.list$Lgamma1 > 0)*100

#k coef #prob that in direction of posterior mean estimate
mean(intro.pred.out$BUGSoutput$sims.list$Kgamma1 < 0)*100

########################################################################################################
##########################################PLOTS#########################################################
#calc for Linf
# Select random slopes  
mean.Linf <- intro.pred.out$BUGSoutput$mean$BB[,1] #Linf is in column 1

# Fake data to predict
fake1 <- seq(min(z_estab), max(z_estab), length=100)

# Obtain parameters of interest
Linf <- intro.pred.out$BUGSoutput$summary[c("mu.Linf", "Lgamma1"),1]
                    
# Fitted lines
fit1 <- Linf[1] + Linf[2]*fake1

# Obtain 90% CIs for fitted line
est.lineLinf <- matrix(NA, ncol=length(fake1), nrow=intro.pred.out$BUGSoutput$n.sims) #container for predicted values


for(i in 1:intro.pred.out$BUGSoutput$n.sims){
  for(t in 1:length(fake1)){
    est.lineLinf[i,t] <- intro.pred.out$BUGSoutput$sims.list$mu.Linf[i] + intro.pred.out$BUGSoutput$sims.list$Lgamma1[i] * fake1[t] 
  }
}

# CIs for fitted values
upper.CILinf <- apply(est.lineLinf, 2, quantile, probs=c(0.975))
lower.CILinf <- apply(est.lineLinf, 2, quantile, probs=c(0.025))

## Grab 95% CIs for beta's
u.Linf <- numeric(length(mean.Linf) )
l.Linf <- numeric(length(mean.Linf) )
for(i in 1:length(mean.Linf)) { 
  u.Linf[i] <- quantile(intro.pred.out$BUGSoutput$sims.list$BB[,i,1],probs=c(0.975) )
  l.Linf[i] <- quantile(intro.pred.out$BUGSoutput$sims.list$BB[,i,1],probs=c(0.025) )
}

############
#calc for k

# Select random slopes  
mean.K<- intro.pred.out$BUGSoutput$mean$BB[,2] #k is in column 2

# Fake data to predict
fake1 <- seq(min(z_estab), max(z_estab), length=100)

# Obtain parameters of interest
K <- intro.pred.out$BUGSoutput$summary[c("mu.k", "Kgamma1"),1]

# Fitted lines
fit2 <- K[1] + K[2]*fake1

# Obtain 95% CIs for fitted line
est.lineK <- matrix(NA, ncol=length(fake1), nrow=intro.pred.out$BUGSoutput$n.sims) #container for predicted values


for(i in 1:intro.pred.out$BUGSoutput$n.sims){
  for(t in 1:length(fake1)){
    est.lineK[i,t] <- intro.pred.out$BUGSoutput$sims.list$mu.k[i] + intro.pred.out$BUGSoutput$sims.list$Kgamma1[i] * fake1[t] 
  }
}

# CIs for fitted values
upper.CIK <- apply(est.lineK, 2, quantile, probs=c(0.975) )
lower.CIK<- apply(est.lineK, 2, quantile, probs=c(0.025) )

## Grab 95% CIs for beta's
u.K <- numeric(length(mean.K) )
l.K <- numeric(length(mean.K) )
for(i in 1:length(mean.K)) { 
  u.K[i] <- quantile(intro.pred.out$BUGSoutput$sims.list$BB[,i,2],probs=c(0.975) )
  l.K[i] <- quantile(intro.pred.out$BUGSoutput$sims.list$BB[,i,2],probs=c(0.025) )
}



########create panel plot
res <- 6
name_figure <- "Intro_GrowthParams_estabtime.png"
png(filename = name_figure, height = 500*res, width = 700*res, res=72*res)
def.par <- par(no.readonly = TRUE)

par(mar=c(0.0,0.1,0.1,0.1),oma=c(3,4,0,1),mai=c(0.15,0.15,0.15,0.15) )

size.labels = 1
size.text = 1

# Add transparency to points
source('AddTransFunction.R')

nf <- layout(matrix(c(1:2),nrow=2,ncol=1,byrow=TRUE),  TRUE) 
layout.show(nf)

size.labels = 1
size.text = 1.0

x.label2 <- 'Standarized Time Since Establishment'
y.label2 <- expression(paste(log[e],'(L'[infinity],")"))

x.label <- 'Standarized Time Since Establishment'
y.label <- expression(paste(log[e],'(Growth rate (K))'))

#k plot
plot(mean.K ~ z_estab ,pch=16,axes=F, xlab='',ylab='',cex=0.8,type='n',ylim=c(min(l.K), max(u.K)) )

points(z_estab, mean.K,pch=21,cex=2,col='black', bg='gray'  )

segments(x0=z_estab, x1=z_estab,y0=l.K, y1=u.K, col="black",lwd=2)

lines(fake1,fit2, lwd = 5, col="gray", lty = 1)

axis(side=1,cex.axis=size.text, mgp=c(0,0.5,0),tck= -0.01, labels=F )
axis(side=2,cex.axis=size.text,  mgp=c(0,0.5,0),tck= -0.01, las=1 )

i.for <- order(fake1)
i.back <- order(fake1, decreasing = TRUE )
x.polygon <- c( fake1[i.for] , fake1[i.back] )
y.polygon <- c( lower.CIK[i.for] , upper.CIK[i.back] )
polygon( x.polygon , y.polygon , col = addTrans('gray', 100) , border = NA )

mtext(y.label, line = 1.9, side = 2, cex = size.text)

box()


#Linf plot
plot(mean.Linf ~ z_estab,pch=16,axes=F, xlab='',ylab='',cex=0.8,type='n' , ylim=c(6.5, 7.4))

points(z_estab, mean.Linf,pch=21,cex=2,col='black', bg='gray' )

segments(x0=z_estab, x1=z_estab, y0=l.Linf, y1=u.Linf, col='black',lwd=2)

lines(fake1,fit1, lwd =5, col="gray",lty = 1 )

axis(side=1,cex.axis=size.text, mgp=c(0,0.5,0),tck= -0.01 ) 
axis(side=2,cex.axis=size.text,  mgp=c(0,0.5,0),tck= -0.01, las=1 ) 

i.for <- order(fake1)
i.back <- order(fake1, decreasing = TRUE )
x.polygon <- c( fake1[i.for] , fake1[i.back] )
y.polygon <- c( lower.CILinf[i.for] , upper.CILinf[i.back] )
polygon( x.polygon , y.polygon , col = addTrans("gray",100), border = NA )

mtext(x.label2, line = 1.5, side = 1, cex = size.text)
mtext(y.label2, line = 1.7, side = 2, cex = size.text)

box()
#####################
par(def.par)
dev.off()

