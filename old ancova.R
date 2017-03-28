######### Calc values for K 
#separate z_lat by status
natives <- subset(dat, Status=="native")
nat_lat <- as.numeric(by(natives$Lat, natives$River, mean)) 
z_natlat <- as.numeric(scale(nat_lat)) #no this method gives you a diff mean and thereofore diff z lat values

intros <- subset(dat, Status=="introduced")
intro_lat <- as.numeric(by(intros$Lat, intros$River, mean)) 
z_introlat <- as.numeric(scale(intro_lat)) #no this method gives you a diff mean and thereofore diff z lat values

# Select random slopes  
meanK.nat <- pred.out$BUGSoutput$mean$BB[c(3:5,8:9,11,13, 14:16,23),2] #select only native pops
meanK.intro <- pred.out$BUGSoutput$mean$BB[c(1:2,6,7,10,12,15,18:22),2] #select only intro pops
meanK.all <- pred.out$BUGSoutput$mean$BB[,2]

#fake pred values
fake.nat <- seq(min(z_natlat), max(z_natlat), length=50)
fake.intro <- seq(min(z_introlat), max(z_introlat), length=50)

# Create container
#Kmeanest <- pred.out$BUGSoutput$summary[c("mu.k", "Kgamma1", "Kgamma2", "Kgamma3"),1]
K_est.lineNat <- matrix(NA, ncol=length(fake.nat), nrow=pred.out$BUGSoutput$n.sims)
K_est.lineIntro <- matrix(NA, ncol=length(fake.intro), nrow=pred.out$BUGSoutput$n.sims)

# calc predicted values for natives
for(i in 1:pred.out$BUGSoutput$n.sims){
  for(j in 1:length(fake.nat) ){
    K_est.lineNat[i, j] <- pred.out$BUGSoutput$sims.list$mu.k[i] +  pred.out$BUGSoutput$sims.list$Kgamma2[i] * fake.nat[j] 
  }
}


# Derive introduced pops ints and slopes
K_IntroInt <- pred.out$BUGSoutput$sims.list$mu.k + pred.out$BUGSoutput$sims.list$Kgamma1
K_IntroSlope <- pred.out$BUGSoutput$sims.list$Kgamma2 + pred.out$BUGSoutput$sims.list$Kgamma3

#calc predicted values for natives
for(i in 1:pred.out$BUGSoutput$n.sims){
  for(j in 1:length(fake.x) ){
    K_est.lineIntro[i, j] <- K_IntroInt[i] + K_IntroSlope[i] * fake.intro[j]
  }
}

### Obtain posterior mean fitted line and associated 
# 95% CRIs (upper and lower CRIs for predicted values in the matrix est.line
#Natives
K_fitted.mean1 <-  apply(K_est.lineNat, 2, mean )
K_upper.CRI1 <- apply(K_est.lineNat, 2, quantile, probs=c(0.95))
K_lower.CRI1 <- apply(K_est.lineNat, 2, quantile, probs=c(0.05))
#Intro
K_fitted.mean2 <-  apply(K_est.lineIntro, 2, mean )
K_upper.CRI2 <- apply(K_est.lineIntro, 2, quantile, probs=c(0.95))
K_lower.CRI2 <- apply(K_est.lineIntro, 2, quantile, probs=c(0.05))


## Grab 95% CIs for slopes of each status??? 
#natives
u.Knat <- numeric(length(meanK.nat))
l.Knat <- numeric(length(meanK.nat))
#natives.slopes <- pred.out$BUGSoutput$mean$BB[c(3:5,8:9,11,13, 14:16,23),2]
for(i in 1:length(pred.out$BUGSoutput$mean$BB[c(3:5,8:9,11,13, 14:16,23),2])){ 
  u.Knat[i] <- quantile(pred.out$BUGSoutput$mean$BB[c(3:5,8:9,11,13, 14:16,23),2],probs=c(0.95))
  l.Knat[i] <-quantile(pred.out$BUGSoutput$mean$BB[c(3:5,8:9,11,13, 14:16,23),2], probs=c(0.5))
}

#intros
u.Kintro <- numeric(length(meanK.intro))
l.Kintro <- numeric(length(meanK.intro))
for(i in 1:length(meanK.intro)) { 
  u.Kintro [i] <- quantile(meanK.intro[[i]],probs=c(0.95))
  l.Kintro [i] <-quantile(meanK.intro[[i]],probs=c(0.5))
}



######### Calc values for Linf 

# Select random slopes  
meanLinf.nat <- pred.out$BUGSoutput$mean$BB[c(3:5,8:9,11,13, 14:16,23),1] #select only native pops
meanLinf.intro <- pred.out$BUGSoutput$mean$BB[c(1:2,6,7,10,12,15,18:22),1] #select only intro pops
meanLinf.all <- pred.out$BUGSoutput$mean$BB[,1]


# Create container
#Linfmeanest <- pred.out$BUGSoutput$summary[c("mu.Linf", "Lgamma1", "Lgamma2", "Lgamma3"),1]
Linf_est.lineNat <- matrix(NA, ncol=length(fake.nat), nrow=pred.out$BUGSoutput$n.sims)
Linf_est.lineIntro <- matrix(NA, ncol=length(fake.intro), nrow=pred.out$BUGSoutput$n.sims)

# calc predicted values for natives
for(i in 1:pred.out$BUGSoutput$n.sims){
  for(j in 1:length(fake.nat) ){
    Linf_est.lineNat[i, j] <- pred.out$BUGSoutput$sims.list$mu.Linf[i] +  pred.out$BUGSoutput$sims.list$Lgamma2[i] * fake.nat[j] 
  }
}

# Derive introduced pops ints and slopes
Linf_IntroInt <- pred.out$BUGSoutput$sims.list$mu.Linf + pred.out$BUGSoutput$sims.list$Lgamma1
Linf_IntroSlope <- pred.out$BUGSoutput$sims.list$Lgamma2 + pred.out$BUGSoutput$sims.list$Lgamma3

#calc predicted values for natives
for(i in 1:pred.out$BUGSoutput$n.sims){
  for(j in 1:length(fake.x) ){
    Linf_est.lineIntro[i, j] <- Linf_IntroInt[i] + Linf_IntroSlope[i] * fake.intro[j]
  }
}

### Obtain posterior mean fitted line and associated 
# 95% CRIs (upper and lower CRIs for predicted values in the matrix est.line
#Natives
Linf_fitted.mean1 <-  apply(Linf_est.lineNat, 2, mean )
Linf_upper.CRI1 <- apply(Linf_est.lineNat, 2, quantile, probs=c(0.95) )
Linf_lower.CRI1 <- apply(Linf_est.lineNat, 2, quantile, probs=c(0.05) )
#Intros
Linf_fitted.mean2 <-  apply(Linf_est.lineIntro, 2, mean )
Linf_upper.CRI2 <- apply(Linf_est.lineIntro, 2, quantile, probs=c(0.95) )
Linf_lower.CRI2 <- apply(Linf_est.lineIntro, 2, quantile, probs=c(0.05) )



#################################################################################
######################CREATE PANEL PLOT##########################################
#################################################################################
res <- 6
name_figure <- "anocoav_panel_plot .jpg"
jpeg(filename = name_figure, height = 500*res, width = 500*res, res=72*res)
def.par <- par(no.readonly = TRUE) 		# save default, for resetting...

nf <- layout(matrix(c(1:4),nrow=2,ncol=2,byrow=TRUE),  TRUE) 
layout.show(nf)

size.labels = 1
size.text = 1.0

x.label <- 'Standarized Latitude'
y.label <- expression(paste(log[e],'(Growth Coefficent (K))'))

x.label2 <- 'Standarized Latitude'
y.label2 <- expression(paste(log[e],'(L'[infinity],")"))


# Add transparency to points
source('AddTransFunction.R')

#K plots
#natives
plot(meanK.nat ~ z_natlat, data=dat, type="n", axes=F,xlab='',ylab='', ylim=c(-2.5, -1.0), xlim=c(-1.5,2))

points(z_natlat, meanK.nat, pch=21, cex=1.5, col='black',bg='blue')

#add segments
segments(x0= z_natlat, x1= z_natlat, y0=l.Knat, y1=u.Knat, col='black',lwd=1)


axis(side=1,cex.axis=size.text, tck=-0.01, mgp=c(0,0,0))
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
plot(meanK.intro ~ z_introlat, data=dat, type="n", axes=F,xlab='',ylab='',ylim=c(-2.5, -1.0), xlim=c(-1.5,2) )

points(z_introlat, meanK.intro, pch=21, cex=1.5, col='black',bg='gray')

axis(side=1,cex.axis=size.text, tck=-0.01, mgp=c(0,0,0))
axis(side=2,cex.axis=size.text, las=1)

# Add fitted line
lines(fake.intro,K_fitted.mean2, lwd = 3, col="gray", lty = 1)
# Add credible region
#natives
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
plot(meanLinf.nat ~ z_natlat, data=dat, type="n", axes=F,xlab='',ylab='',ylim=c(6.4, 7.5), xlim=c(-1.5,2) )

points(z_natlat, meanLinf.nat, pch=21, cex=1.5, col='black',bg='blue')

axis(side=1,cex.axis=size.text, tck=-0.01, mgp=c(0,0,0))
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
plot(meanLinf.intro ~ z_introlat, data=dat, type="n", axes=F,xlab='',ylab='', ylim=c(6.4, 7.5), xlim=c(-1.5,2) )

points(z_introlat, meanLinf.intro, pch=21, cex=1.5, col='black',bg='gray')

axis(side=1,cex.axis=size.text, tck=-0.01, mgp=c(0,0,0))
axis(side=2,cex.axis=size.text, las=1)

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
#####################
par(def.par)
dev.off()