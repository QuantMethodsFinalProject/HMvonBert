
# Select random slopes  
mean.beta <- pred.out$BUGSoutput$mean$BB[,2] #k is in column 2

#fake pred values
fake.x <- seq(min(z_lat), max(z_lat), length=50)

# Obtain parameters of interest
#meanest <- pred.out$BUGSoutput$summary[c("mu.k", "Kgamma1", "Kgamma2", "Kgamma3"),1]
est.lineNat <- matrix(NA, ncol=length(fake.x), nrow=pred.out$BUGSoutput$n.sims)
est.lineIntro <- matrix(NA, ncol=length(fake.x), nrow=pred.out$BUGSoutput$n.sims)

# Predict response variable for all MCMC samples and for every value of fake.x
#Natives
for(i in 1:pred.out$BUGSoutput$n.sims){
  for(j in 1:length(fake.x) ){
    est.lineNat[i, j] <- pred.out$BUGSoutput$sims.list$mu.k[i] + pred.out$BUGSoutput$sims.list$Lgamma1[i] * 0 + pred.out$BUGSoutput$sims.list$Lgamma2[i] * fake.x[j] +
                         pred.out$BUGSoutput$sims.list$Lgamma3[i] * fake.x[j] * 0
  }
}

#Intros
for(i in 1:pred.out$BUGSoutput$n.sims){
  for(j in 1:length(fake.x) ){
    est.lineIntro[i, j] <- pred.out$BUGSoutput$sims.list$mu.k[i] + pred.out$BUGSoutput$sims.list$Lgamma1[i] * 1 + pred.out$BUGSoutput$sims.list$Lgamma2[i] * fake.x[j] +
                           pred.out$BUGSoutput$sims.list$Lgamma3[i] * fake.x[j] * 1
  }
}

### Obtain posterior mean fitted line and associated 
# 95% CRIs (upper and lower CRIs for predicted values in the matrix est.line
#Natives
fitted.mean1 <-  apply(est.lineNat, 2, mean )
upper.CRI1 <- apply(est.lineNat, 2, quantile, probs=c(0.95) )
lower.CRI1 <- apply(est.lineNat, 2, quantile, probs=c(0.05) )
#Intro
fitted.mean2 <-  apply(est.lineIntro, 2, mean )
upper.CRI2 <- apply(est.lineIntro, 2, quantile, probs=c(0.95) )
lower.CRI2 <- apply(est.lineIntro, 2, quantile, probs=c(0.05) )

###Actual plot
res <- 6
name_figure <- "HM ANCOVA .jpg"
jpeg(filename = name_figure, height = 500*res, width = 500*res, res=72*res)
def.par <- par(no.readonly = TRUE) 		# save default, for resetting...

nf <- layout(matrix(c(1:1),nrow=1,ncol=1,byrow=TRUE),  TRUE) 
layout.show(nf)
par(mar=c(1,1,1,1), oma=c(3,3,0,0) )

size.labels = 1
size.text = 1.0

x.label <- 'Standarized Latitude'
y.label <- expression(paste(log[e],'(Growth rate (K))'))

# Add transparency to points
source('AddTransFunction.R')


plot(mean.beta ~ z_lat,pch=16,axes=F, xlab='',ylab='',cex=0.8,type='n',ylim=c(min(lower.CRI2), max(upper.CRI2)) )

points(z_lat[dat$Status=="native"], mean.beta[dat$Status=="native"],pch=16,cex=0.8, col="black")

points(z_lat[dat$Status=="introduced"], mean.beta[dat$Status=="introduced"],pch=16,cex=0.8, col="red")


axis(side=1,cex.axis=size.text, tck=-0.01, mgp=c(0,0.5,0))
axis(side=2,cex.axis=size.text, las=1)

# Add fitted line
lines(fake.x,fitted.mean1, lwd = 3, col="black", lty = 1)
lines(fake.x,fitted.mean2, lwd = 3, col="red", lty = 1)

# Add credible region
#natives
i.for <- order(fake.x)
i.back <- order(fake.x, decreasing = TRUE )
x.polygon <- c( fake.x[i.for] , fake.x[i.back] )
y.polygon <- c( lower.CRI1[i.for] , upper.CRI1[i.back] )
polygon( x.polygon , y.polygon , col = addTrans("black",100)  , border = NA )

#intro
i.for <- order(fake.x)
i.back <- order(fake.x, decreasing = TRUE )
x.polygon <- c( fake.x[i.for] , fake.x[i.back] )
y.polygon <- c( lower.CRI2[i.for] , upper.CRI2[i.back] )
polygon( x.polygon , y.polygon , col = addTrans("red",100)  , border = NA )

#axis labs
mtext(x.label, line = 3, side = 1, cex = size.text)
mtext(y.label, line = 3, side = 2, cex = size.text)


box()
par(def.par)
dev.off()



