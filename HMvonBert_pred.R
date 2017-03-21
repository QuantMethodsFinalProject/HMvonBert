#rm(list=ls())

#### Load libraries
library(FSA)
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
#J =23
#must drop levels before ordering river in next step
dat <- droplevels(dat) 


#remove old data from analysis
dat <- subset(dat, Year >= 2000) 
#newdata <- subset(dat, Year >= 2007)
sort(unique(dat$Year))
#how many rows were removed?
#7770 - nrow(dat)
#should be 172


#must remove empty levels before continuing 
dat <- droplevels(dat)


#river name column
dat <- dat[order(dat$River), ]
dat$River_name <- dat$River
dat$River <- as.numeric(dat$River)
dat <- droplevels(dat)


#final check:
length(unique(dat$River)) #should be 23
unique(dat$River) #should be in a consecutive numeric order
unique(dat$River_name) #does not include removed pops above

###################################ADD PREDICTORS
# Site-level predictors
lat <- as.numeric(by(dat$Lat, dat$River, mean)) 
z_lat <- as.numeric(scale(lat))
#not used until we have lats for all sites

# Categorical predictor ('Native' is the reference cell, which means is is set to zero)
dat$Status_binary <- ifelse(dat$Status == 'introduced', 1, 0)
status <- as.numeric(by(dat$Status_binary, dat$River, mean)) # Don't need to standardize this categorical predictor variable

#################################################################
########## BUGS CODE ############################################
#################################################################
sink("MetaAnalysisHMvonBmodel_pred.txt")
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
    BB.hat[j,1] <- mu.Linf + Lgamma1 * status[j] + Lgamma2 * z_lat[j] + Lgamma3 * z_lat[j] * status[j]
    BB.hat[j,2] <- mu.k + Kgamma1 * status[j] + Kgamma2 * z_lat[j] + Kgamma3 * z_lat[j] * status[j]
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
J <- length(unique(dat$River))


# Create identity matrix for Wishart dist'n
#!!Number of varying parameters to estimate (K)
K <- 3

# Create identity matrix for Wishart dist'n
W <- diag(3)

# load data
data <- list(y = dat$TL, age = dat$Age, g = dat$River, n = dim(dat)[1],
             J = J, W=W, K=K, status=status , z_lat=z_lat)

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



pred.out <- jags(data = data, inits = inits, parameters.to.save = params1, 
            model.file = "MetaAnalysisHMvonBmodel_pred.txt", n.chains = nc, 
            n.thin = nt, n.iter = ni, n.burnin = nb)


#Calculate computation time
end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 


# Summarize the result
print(pred.out, digits = 2)
# str(out)

# Find which parameters, if any, have Rhat > 1.1
which(pred.out$BUGSoutput$summary[, c("Rhat")] > 1.1) #three increase iterations

#check traceplots
#code may be wrong
traceplot(pred.out)
############Relationships
# Does Linf and K differ between native and introduced?
#Linf
L_natslope <- pred.out$BUGSoutput$sims.list$Lgamma2 #ref cell
L_introslope <- pred.out$BUGSoutput$sims.list$Lgamma2 + pred.out$BUGSoutput$sims.list$Lgamma3
L_Diffslope <- L_natslope - L_introslope 
quantile(L_Diffslope, c(0.025, 0.975))
quantile(L_Diffslope, c(0.05, 0.95)) #overlaps zero, so no difference

#K coef
K_natslope <- pred.out$BUGSoutput$sims.list$Kgamma2 #ref cell
K_introslope <- pred.out$BUGSoutput$sims.list$Lgamma2 +  pred.out$BUGSoutput$sims.list$Lgamma3
K_Diffslope <- K_natslope - K_introslope 
quantile(K_Diffslope, c(0.025, 0.975))
quantile(K_Diffslope, c(0.05, 0.95)) #overlaps zero, so no difference


###############
# Does Linf and K differ by lattitude?
#Linf
#NATIVE POPS
mean(L_natslope)
# What is the probability that the effect of latitude is negative
mean(L_natslope < 0)
quantile(L_natslope,c(0.05, 0.95))
quantile(L_natslope,c(0.025, 0.975))# does it overlap zero? NO= they are different
#INTRO POPS
mean(L_introslope)
# What is the probability that the effect of latitude is negative
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


# Example plot
X1 <- seq(min(z_lat),max(z_lat), length=20)
y.hat <- pred.out$BUGSoutput$mean$mu.Linf + pred.out$BUGSoutput$mean$Lgamma2 * X1
plot(X1,y.hat, type='l')
plot(X1, exp(y.hat), type='l')

###########################################
##### Plot Popl'n average effect

# fake data to predict
predX <- seq(min(dat$Age), max(dat$Age), length=100) 

# Extract pop'n average coefficents

PopAve <- pred.out$BUGSoutput$summary[c("mu.Linf", "mu.k", "mu.t0"),1]

# Re-transform
PopAve[1] <- exp(PopAve[1])
PopAve[2] <- exp(PopAve[2])
PopAve[3] <- exp(PopAve[3])-10

GroupCoef <- matrix(pred.out$BUGSoutput$summary[1:(3*(length(unique(dat$River)) )),1], c(length(unique(dat$River)),3), byrow=F)

# Re-transform
GroupCoef[,1] <- exp(GroupCoef[,1])
GroupCoef[,2] <- exp(GroupCoef[,2])
GroupCoef[,3] <- exp(GroupCoef[,3])-10


# Generate fake y-axis for creating plot
z <- seq(min(dat$Age),max(dat$Age),length=100) 

# Create Popl'n average mean fitted line
y.pred <- PopAve[1] * (1-exp(-PopAve[2]  * (predX -  PopAve[3] )))
 

##### Create Figure
res <- 6
name_figure <- "HMMetavonBert(J=23)Predictors2.png"
jpeg(filename = name_figure, height = 500*res, width = 500*res, res=72*res)
def.par <- par(no.readonly = TRUE)

par(mar=c(0.0,0.1,0.1,0.1),oma=c(3,3,0,1),mai=c(0.0,0.05,0.05,0) )

size.labels = 1
size.text = 1
x.label = 'Age (years)'
y.label = 'Length (mm)'


plot(predX,z, ylim=c(min(dat$TL),max(dat$TL)),xlim=c(min(dat$Age),max(dat$Age)), axes=F, ylab='', xlab='', type='n')
points(jitter(dat$Age),dat$TL, cex=0.8, pch=16)

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

text(25, 350, paste("mu.Linf =", round(PopAve[1], 2), "\n mu.k =", 
                    round(PopAve[2], 2), "\n mu.t0 =", round(PopAve[3],2)))   #adds coef to plot w/ titles        


par(def.par)
dev.off()

#########################################################
#Plot of group mod with CI region

# fake data to predict
predX <- seq(min(dat$Age), max(dat$Age), length=100) 


# Container for predicted values
est.lineB <- array(NA, c(out$BUGSoutput$n.sims,length(predX),J) )

# Put each groups MCMC draws for all parameters in its own list
group.params <- list()
for(m in 1:J){
  group.params[[m]] <- out$BUGSoutput$sims.list$BB[,m,]
}

# Re-transform: look at - str(group.params)
for(j in 1:J){
  group.params[[j]][,1] <- exp(group.params[[j]][,1] )
  group.params[[j]][,2] <- exp(group.params[[j]][,2] )
  group.params[[j]][,3] <- exp(group.params[[j]][,3] )-10
}



for(k in 1:J){ # loop over groups (J)
  for(i in 1:out$BUGSoutput$n.sims ){  
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



############PLOT group params w/ CI region
names <- levels(dat$River_name)
res <- 6
name_figure <- "Group_Specific_HMMetavonBert(J=23)Predictors.png"
png(filename = name_figure, height = 500*res, width = 800*res, res=72*res)
def.par <- par(no.readonly = TRUE)

size.labels = 1
size.text = 1
axissize <- 1
x.label = 'Age (yrs)'
y.label = "Length (mm)"

nf <- layout(matrix(c(1:24),nrow=6,ncol=4,byrow=TRUE),  TRUE) 
layout.show(nf)
par(mar=c(0.0,0.1,0.1,0.1),oma=c(3,3,0,1),mai=c(0.0,0.05,0.05,0) )

# Add transparency to points
source('AddTransFunction.R')


# Group-specific plot
for(i in 1:J){
  
  plot(predX,z, ylim=c(min(dat$TL,na.rm=T),max(dat$TL,na.rm=T)),
       xlim=c(min(dat$Age),max(dat$Age)), axes=F, ylab='', xlab='', type='n')
  
  points(jitter(dat$Age[dat$River==i]), dat$TL[dat$River==i], cex=0.8, pch=16,col='black' )
  
  
  if( i <=23){
    axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, labels=F ) 
  } else {
    axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01)
  }	
  
  if( i ==1 | i==5 |  i==9 |  i==13 |  i==17|  i==21 ){
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


###############################################

#Growth parameter comparison per site
#Linf
L_Diff <- matrix(NA,nrow=J,ncol=J)
for (j in 1:J){
  for (jj in 1:J){
    
    diff <- out$BUGSoutput$sims.list$BB[,j,1] - out$BUGSoutput$sims.list$BB[,jj,1]
    diffCI <- quantile(diff,c(0.025,0.975))
    L_Diff[j,jj] <- (diffCI[1] * diffCI[2]) > 0 # Does the 95% CI overlap zero? TRUE = signifcant difference
  }
}
rownames(L_Diff) <- names
colnames(L_Diff) <- names
View(L_Diff)




#k
k_Diff <- matrix(NA,nrow=J,ncol=J)
for (j in 1:J){
  for (jj in 1:J){
    
    diff <- out$BUGSoutput$sims.list$BB[,j,2] - out$BUGSoutput$sims.list$BB[,jj,2]
    diffCI <- quantile(diff,c(0.025,0.975))
    k_Diff[j,jj] <- (diffCI[1] * diffCI[2]) > 0 # Does the 95% CI overlap zero? TRUE = signifcant difference
  }
}
rownames(k_Diff) <- names
colnames(k_Diff) <- names
View(k_Diff)



##########################
#plot relationship of native on introduced
#(show difference in intercepts)



