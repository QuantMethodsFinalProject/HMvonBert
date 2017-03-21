#rm(list=ls())


#### Load libraries
library(FSA)
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

intro2<-intro[intro$River!="Lake Marion", ]
intro2<-intro2[intro2$River!="Lake Moultrie", ]

#length(unique(intro$River)) #check J=12 
length(unique(intro2$River)) #check J=10 (without Lake pops)
#unique(intro$TimeSinceEstab)
unique(intro2$TimeSinceEstab)
mean(intro2$TimeSinceEstab)

#must drop levels before continuing
intro2 <- droplevels(intro2)
intro <- droplevels(intro)

#river name column
# intro <- intro[order(intro$River), ]
# intro$River_name <- intro$River
# intro$River <- as.numeric(intro$River)
intro2 <- intro2[order(intro2$River), ]
intro2$River_name <- intro2$River
intro2$River <- as.numeric(intro2$River)

#final check:
# length(unique(intro$River)) #should be 12 for intro pops
# unique(intro$River) #should be in a consecutive numeric order
# unique(intro$River_name) #does not include removed pops above or NATIVES
# sort(unique(intro$TimeSinceEstab))
length(unique(intro2$River)) #should be 12 for intro pops
unique(intro2$River) #should be in a consecutive numeric order
unique(intro2$River_name) #does not include removed pops above or NATIVES
sort(unique(intro2$TimeSinceEstab))
###################################ADD PREDICTORS
# Site-level predictors
estab <- as.numeric(by(intro$TimeSinceEstab, intro$River, mean)) 
z_estab <- as.numeric(scale(estab))

estab2 <- as.numeric(by(intro2$TimeSinceEstab, intro2$River, mean)) 
z_estab2 <- as.numeric(scale(estab2))
 #remember z estab and zestab2 will have different values for same rivers because means are different 
#################################################################
########## BUGS CODE ############################################
#################################################################
sink("MetaAnalysisHMvonBmodel_INTRO2pred.txt")
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
J <- length(unique(intro2$River))


# Create identity matrix for Wishart dist'n
#!!Number of varying parameters to estimate (K)
K <- 3

# Create identity matrix for Wishart dist'n
W <- diag(3)

# load data
data <- list(y = intro2$TL, age = intro2$Age, g = intro2$River, n = dim(intro2)[1],
             J = J, W=W, K=K, time=z_estab2)

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



intro2.out <- jags(data = data, inits = inits, parameters.to.save = params1, 
            model.file = "MetaAnalysisHMvonBmodel_INTRO2pred.txt", n.chains = nc, 
            n.thin = nt, n.iter = ni, n.burnin = nb)


#Calculate computation time
end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 


# Summarize the result
print(intro2.out, digits = 2)
# str(out)

# Find which parameters, if any, have Rhat > 1.1
which(intro2.out$BUGSoutput$summary[, c("Rhat")] > 1.1) #three increase iterations

#check traceplots
#code may be wrong
traceplot(intro.out)

########################
# Does time since estab effect Linf and k in intro pops?
#Linf
quantile(intro2.out$BUGSoutput$sims.list$Lgamma1, c(0.025, 0.975))
quantile(intro2.out$BUGSoutput$sims.list$Lgamma1, c(0.05, 0.95))
#overlaps zero, so no difference


#K coef
quantile(intro2.out$BUGSoutput$sims.list$Kgamma1, c(0.025, 0.975))
quantile(intro2.out$BUGSoutput$sims.list$Kgamma1, c(0.05, 0.95))
#overlaps zero, so no difference



# Prob that slope is positive
#prob that param is in the direction of posterior mean estimate
#Linf
mean(intro2.out$BUGSoutput$sims.list$Lgamma1 > 0)*100

#k coef #prob that in direction of posterior mean estimate
mean(intro2.out$BUGSoutput$sims.list$Kgamma1 > 0)*100
#mean(intro.out$BUGSoutput$sims.list$Kgamma1 < 0)*100
##################################################################
#plot relationship of L and K per pop against posterior mean esimate of that param

###########################LINF HM REG PLOT

# Select random slopes  
mean.beta <- intro2.out$BUGSoutput$mean$BB[,1] #Linf is in column 1

# Fake data to predict
fake1 <- seq(min(z_estab2), max(z_estab2), length=50)

# Obtain parameters of interest
ests2 <- intro2.out$BUGSoutput$summary[c("mu.Linf", "Lgamma1"),1]
                    
# Fitted lines
fit1 <- ests2[1] + ests2[2]*fake1

# Obtain 90% CIs for fitted line
est.lineA <- matrix(NA, ncol=length(fake1), nrow=intro2.out$BUGSoutput$n.sims) #container for predicted values


for(i in 1:intro2.out$BUGSoutput$n.sims){
  for(t in 1:length(fake1)){
    est.lineA[i,t] <- intro2.out$BUGSoutput$sims.list$mu.Linf[i] + intro2.out$BUGSoutput$sims.list$Lgamma1[i] * fake1[t] 
  }
}

# CIs for fitted values
upper.CIA <- apply(est.lineA, 2, quantile, probs=c(0.975) )
lower.CIA <- apply(est.lineA, 2, quantile, probs=c(0.025) )

## Grab 95% CIs for beta's
u.alpha <- numeric(length(mean.beta) )
l.alpha <- numeric(length(mean.beta) )
for(i in 1:length(mean.beta)) { 
  u.alpha[i] <- quantile(intro2.out$BUGSoutput$sims.list$BB[,i,1],probs=c(0.975) )
  l.alpha[i] <- quantile(intro2.out$BUGSoutput$sims.list$BB[,i,1],probs=c(0.025) )
}

###########################################
####### FIGURE WITH CRI's
res <- 6
name_figure <- "Intro2_log(Linf)_estabtime.png"
png(filename = name_figure, height = 500*res, width = 800*res, res=72*res)
def.par <- par(no.readonly = TRUE)
size.labels = 1
size.text = 1


par(mar=c(0.0,0.1,0.1,0.1),oma=c(3,3,0,1),mai=c(0.0,0.05,0.05,0) )

x.label <- 'Standarized Time Since Establishment'
y.label <- expression(paste(log[e],'(L'[infinity],")"))

#xlab1 <- c(-3,-2,-1,0,1)
#xlab2 <- xlab1* sd(lat1) + mean(lat1)

plot(mean.beta ~ z_estab2,pch=16,axes=F, xlab='',ylab='',cex=0.8,type='n' , ylim=c(6.5, 7.4))

axis(side=1,cex.axis=size.text, mgp=c(0,0.5,0),tck= -0.01 ) #, at=xlab1, labels=round(xlab2,2)
axis(side=2,cex.axis=size.text,  mgp=c(0,0.5,0),tck= -0.01, las=1 ) # at=ylab1, labels=round(ylab2,2)

i.for <- order(fake1)
i.back <- order(fake1, decreasing = TRUE )
x.polygon <- c( fake1[i.for] , fake1[i.back] )
y.polygon <- c( lower.CIA[i.for] , upper.CIA[i.back] )
polygon( x.polygon , y.polygon , col = "lightgray" , border = NA )

points(z_estab2, mean.beta,pch=16,cex=0.8)

segments(x0=z_estab2, x1=z_estab2,
         y0=l.alpha, y1=u.alpha, col='black',lwd=1)


lines(fake1,fit1, lwd = 3, col="black", lty = 1)

mtext(x.label, line = 1.5, side = 1, cex = size.text)
mtext(y.label, line = 1.7, side = 2, cex = size.text)

box()
par(def.par)
dev.off()
####### END


###########################K HM REG PLOT

# Select random slopes  
mean.beta <- intro2.out$BUGSoutput$mean$BB[,2] #k is in column 2

# Fake data to predict
fake1 <- seq(min(z_estab2), max(z_estab2), length=50)

# Obtain parameters of interest
ests2 <- intro2.out$BUGSoutput$summary[c("mu.k", "Kgamma1"),1]

# Fitted lines
fit1 <- ests2[1] + ests2[2]*fake1

# Obtain 95% CIs for fitted line
est.lineA <- matrix(NA, ncol=length(fake1), nrow=intro2.out$BUGSoutput$n.sims) #container for predicted values


for(i in 1:intro2.out$BUGSoutput$n.sims){
  for(t 
      in 1:length(fake1)){
    est.lineA[i,t] <- intro2.out$BUGSoutput$sims.list$mu.k[i] + intro2.out$BUGSoutput$sims.list$Kgamma1[i] * fake1[t] 
  }
}

# CIs for fitted values
upper.CIA <- apply(est.lineA2, 2, quantile, probs=c(0.975) )
lower.CIA <- apply(est.lineA2, 2, quantile, probs=c(0.025) )

## Grab 95% CIs for beta's
u.alpha <- numeric(length(mean.beta) )
l.alpha <- numeric(length(mean.beta) )
for(i in 1:length(mean.beta)) { 
  u.alpha[i] <- quantile(intro2.out$BUGSoutput$sims.list$BB[,i,2],probs=c(0.975) )
  l.alpha[i] <- quantile(intro2.out$BUGSoutput$sims.list$BB[,i,2],probs=c(0.025) )
}

###########################################
####### FIGURE WITH CRI's
res <- 6
name_figure <- "Intro2_log(k)_estabtime.png"
png(filename = name_figure, height = 500*res, width = 800*res, res=72*res)
def.par <- par(no.readonly = TRUE)
size.labels = 1
size.text = 1


par(mar=c(0.0,0.1,0.1,0.1),oma=c(3,4,0,1),mai=c(0.0,0.05,0.05,0) )

x.label <- 'Standarized Time Since Establishment'
y.label <- expression(paste(log[e],'(Growth rate (K))'))

#xlab1 <- c(-3,-2,-1,0,1)
#xlab2 <- xlab1* sd(lat1) + mean(lat1)

plot(mean.beta ~ z_estab2 ,pch=16,axes=F, xlab='',ylab='',cex=0.8,type='n',ylim=c(min(l.alpha), max(u.alpha)) )

axis(side=1,cex.axis=size.text, mgp=c(0,0.5,0),tck= -0.01 ) #, at=xlab1, labels=round(xlab2,2)
axis(side=2,cex.axis=size.text,  mgp=c(0,0.5,0),tck= -0.01, las=1 ) # at=ylab1, labels=round(ylab2,2)

i.for <- order(fake1)
i.back <- order(fake1, decreasing = TRUE )
x.polygon <- c( fake1[i.for] , fake1[i.back] )
y.polygon <- c( lower.CIA[i.for] , upper.CIA[i.back] )
polygon( x.polygon , y.polygon , col = "lightgray" , border = NA )

points(z_estab2, mean.beta,pch=16,cex=0.8)

segments(x0=z_estab2, x1=z_estab2,
         y0=l.alpha, y1=u.alpha, col='black',lwd=1)


lines(fake1,fit1, lwd = 3, col="black", lty = 1)

mtext(x.label, line = 1.5, side = 1, cex = size.text)
mtext(y.label, line = 1.9, side = 2, cex = size.text)

box()
par(def.par)
dev.off()
####### END

