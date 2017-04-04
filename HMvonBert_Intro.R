#### Load libraries
library(mgcv)
library(MCMCpack) # rwish function
library(R2jags)


#read in data
setwd("S:/FlatheadOtolithProject/Data/Rprojects/MetaAnalysis/HMvonBert")
dat <- read.csv("FHCGrowthMetaData.csv")

####clean up dat
#remove prob pop's
intro<-dat[dat$River!="North Raccoon", ] #not really nec since its a native pop but oh well
intro<-intro[intro$River!="Lake Lyndon B. Johnson", ] #not really nec since its a native pop but oh well

#remove old data from analysis
intro <- subset(intro, Year >= 2000)
sort(unique(intro$Year))
#how many rows were removed?
#nrow(newdata) #4632 rows removed WOW!

#subset out only introduced pops
intro <- subset(intro, Status =="introduced")
length(unique(intro$River)) #should be 12

#drop levels BEFORE ordering rivers
intro <- droplevels(intro)

#ordering river names, ordering rivers and turning into numeric
intro <- intro[order(intro$River), ]
intro$River_name <- intro$River
intro$River <- as.numeric(intro$River)


#################################################################
########## BUGS CODE ############################################
#################################################################
sink("MetaAnalysisHMvonBmodel_intro.txt")
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
    BB.hat[j,1] <- mu.Linf
    BB.hat[j,2] <- mu.k
    BB.hat[j,3] <- mu.t0
    }
    
    # Priors for population-average parameters
    mu.Linf ~ dnorm(0,.0001)
    mu.k ~ dnorm(0,.0001)
    mu.t0 ~ dnorm(0,.0001)
    
    
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
             J = J, W=W, K=K )

# Initial values
inits <- function(){list(mu.Linf = rnorm(1,3,0.001), mu.k = rnorm(1,1,0.001), mu.t0 = rnorm(1,0.7,0.001),
                         sigma.y = runif(1,1,10), 
                         BB=array(c(rep(log(950) +rnorm(1,0.01,0.01),J),rep(log(0.04)+rnorm(1,0.001,0.1),J),rep(log(-2+10)+rnorm(1,0.01,0.1),J)),
                                  c(J,K)), Tau.B=rwish(K+1,diag(K)) ) }

# Parameters monitored
params1 <- c("mu.Linf", "mu.k", "mu.t0", "sigma.y","BB","Sigma.B","rho.B")


# MCMC settings
ni <- 100000
nt <- 3
nb <- 50000
nc <- 3


############################################
start.time = Sys.time()  



intro.out <- jags(data = data, inits = inits, parameters.to.save = params1, 
               model.file = "MetaAnalysisHMvonBmodel_intro.txt", n.chains = nc, 
               n.thin = nt, n.iter = ni, n.burnin = nb)


#Calculate computation time
end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 
##############################
#Saving the JAGS output
#DO NOT RUN "SAVE" BY ACCIDENT IF NOT IN YOUR ENVIRONMENT, WILL HAVE TO RE-RUN MODEL IF SO 
#saveRDS(intro.out, compress=TRUE, "S:/FlatheadOtolithProject/Data/Rprojects/MetaAnalysis/HMvonBert/JAGS.rds/intro.out.rds")

#remove from envirn when not using
#rm(intro.out) 

#read in model output
intro.out <- readRDS("S:/FlatheadOtolithProject/Data/Rprojects/MetaAnalysis/HMvonBert/JAGS.rds/intro.out.rds")
identical(intro.out, intro.out, ignore.environment = TRUE)
##########################
# Summarize the result
print(intro.out, digits = 2)
# str(newout)

# Find which parameters, if any, have Rhat > 1.1
which(intro.out$BUGSoutput$summary[, c("Rhat")] > 1.1) #three increase iterations

#check traceplots
#code may be wrong
traceplot(intro.out)

########################################################################################################
##########################################PLOTS#########################################################
##### Plot Intro Popl'n average effect

# fake data to predict
predX <- seq(min(intro$Age), max(intro$Age), length=100) 

# Extract pop'n average coefficents

PopAve <- intro.out$BUGSoutput$summary[c("mu.Linf", "mu.k", "mu.t0"),1]

# Re-transform
PopAve[1] <- exp(PopAve[1])
PopAve[2] <- exp(PopAve[2])
PopAve[3] <- exp(PopAve[3])-10

GroupCoef <- matrix(intro.out$BUGSoutput$summary[1:(3*(length(unique(intro$River)) )),1], c(length(unique(intro$River)),3), byrow=F)

# Re-transform
GroupCoef[,1] <- exp(GroupCoef[,1])
GroupCoef[,2] <- exp(GroupCoef[,2])
GroupCoef[,3] <- exp(GroupCoef[,3])-10


# Generate fake y-axis for creating plot
z <- seq(min(intro$Age),max(intro$Age),length=100) 

# Create Popl'n average mean fitted line
y.pred <- PopAve[1] * (1-exp(-PopAve[2]  * (predX -  PopAve[3] )))


##### Create Figure
res <- 6
name_figure <- "Intro_PopulationHM_vonBert.png"
jpeg(filename = name_figure, height = 500*res, width = 500*res, res=72*res)
def.par <- par(no.readonly = TRUE)

par(mar=c(0.0,0.1,0.1,0.1),oma=c(3,3,0,1),mai=c(0.0,0.05,0.05,0) )

size.labels = 1
size.text = 1
x.label = 'Age (years)'
y.label = 'Length (mm)'


plot(predX,z, ylim=c(min(intro$TL),max(intro$TL)),xlim=c(min(intro$Age),max(intro$Age)), axes=F, ylab='', xlab='', type='n')
points(jitter(intro$Age),intro$TL, cex=0.8, pch=16)

# Add group-specific fits
for(i in 1:J){
  y.pred2 <- GroupCoef[i,1] * (1-exp(-GroupCoef[i,2]  * (predX -  GroupCoef[i,3] )))
  lines(predX, y.pred2,lwd=0.5, col='black',lty=1)
}

# Add popl'n average fit
lines(predX, y.pred, lwd = 5, col="blue", lty = 1)

axis(side=1, cex.axis=size.text, tck=-0.01, mgp=c(0,0.5,0) ) 
axis(side=2,cex.axis=size.text,font=1 ,tck=-0.01, mgp=c(0,0.5,0), las=1) 
mtext(x.label, line = 1.3, side = 1, cex = size.text,outer=T)
mtext(y.label, line = 1.8, side = 2, cex = size.text,outer=T)
box()

# text(25, 350, paste("mu.Linf =", round(PopAve[1], 2), "\n mu.k =", 
#                     round(PopAve[2], 2), "\n mu.t0 =", round(PopAve[3],2)))   #adds coef to plot w/ titles        


par(def.par)
dev.off()

#########################################################
##### Plot Intro Group Specific growth

# fake data to predict
predX <- seq(min(intro$Age), max(intro$Age), length=100) 


# Container for predicted values
est.lineB <- array(NA, c(intro.out$BUGSoutput$n.sims,length(predX),J) )

# Put each groups MCMC draws for all parameters in its own list
group.params <- list()
for(m in 1:J){
  group.params[[m]] <- intro.out$BUGSoutput$sims.list$BB[,m,]
}

# Re-transform: look at - str(group.params)
for(j in 1:J){
  group.params[[j]][,1] <- exp(group.params[[j]][,1] )
  group.params[[j]][,2] <- exp(group.params[[j]][,2] )
  group.params[[j]][,3] <- exp(group.params[[j]][,3] )-10
}



for(k in 1:J){ # loop over groups (J)
  for(i in 1:intro.out$BUGSoutput$n.sims ){  
    for(t in 1:length(predX)){
      # est.lineB[i,t,k] <-  group.params[[k]][i,1] + group.params[[k]][i,2] * predX[t]
      est.lineB[i,t,k] <-  group.params[[k]][i,1] * (1 -exp(-group.params[[k]][i,2] * (predX[t] - group.params[[k]][i,3])))
    }	  
  }
}

groupMean <- array(NA, c(1,length(predX),J) )
upper.CIB <- array(NA, c(1,length(predX),J) )
lower.CIB <- array(NA, c(1,length(predX),J) )

for(i in 1:J){
  # Means
  groupMean[,,i] <- apply(est.lineB[,,i], 2, mean )
  # 95% CIs for fitted values
  upper.CIB[,,i] <- apply(est.lineB[,,i], 2, quantile, probs=c(0.975) )
  lower.CIB[,,i] <- apply(est.lineB[,,i], 2, quantile, probs=c(0.025) )
}



############actual plot
names <- levels(intro$River_name)
res <- 6
name_figure <- "Intro_GroupHM_vonBert.png"
png(filename = name_figure, height = 500*res, width = 800*res, res=72*res)
def.par <- par(no.readonly = TRUE)

size.labels = 1
size.text = 1
axissize <- 1
x.label = 'Age (yrs)'
y.label = "Length (mm)"

nf <- layout(matrix(c(1:12),nrow=3,ncol=4,byrow=TRUE),  TRUE) 
layout.show(nf)
par(mar=c(0.0,0.1,0.1,0.1),oma=c(3,3,0,1),mai=c(0.0,0.05,0.05,0) )

# Add transparency to points
source('AddTransFunction.R')


# Group-specific plot
for(i in 1:J){
  
  plot(predX,z, ylim=c(min(intro$TL,na.rm=T),max(intro$TL,na.rm=T)),
       xlim=c(min(intro$Age),max(intro$Age)), axes=F, ylab='', xlab='', type='n')
  
  points(jitter(intro$Age[intro$River==i]), intro$TL[intro$River==i], cex=0.8, pch=16,col='black' )
  
  
  if( i <=8){
    axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, labels=F ) 
  } else {
    axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01)
  }	
  if( i ==1 | i==5 |  i==9 ){
    axis(side=2,cex.axis=axissize , mgp=c(0,0.3,0),tck= -0.01, las=1)
  } else {
    axis(side=2,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, labels=F)
  }	
  
  
  # Add credible region
  i.for <- order(predX)
  i.back <- order(predX, decreasing = TRUE )
  x.polygon <- c( predX[i.for] , predX[i.back] )
  y.polygon <- c( lower.CIB[,,i][i.for] , upper.CIB[,,i][i.back] )
  polygon( x.polygon , y.polygon , col = addTrans('black',100) , border = NA )
  
  # Add posterior means
  lines(predX, groupMean[,,i],lwd=1, col='black',lty=1)
  
  text(15, 200, names[i], cex=1)
  
  box()
  
}

mtext(y.label, line =1.9, side = 2, cex = size.text,outer=T)
mtext(x.label, line = 1.5, side = 1, cex = size.text, outer=T)


par(def.par)
dev.off()

