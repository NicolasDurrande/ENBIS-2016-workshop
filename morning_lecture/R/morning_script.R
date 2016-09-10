library(shiny)
library(rgl)
source('../../catapult/catapultFunctions.R')
source('../../catapult/catapultSettings.R')

######################
## run catapult simulator
runApp(appDir='../../catapult')

################################
## linear regression in 2D

x = c(0,0.33,0.66,1)
X = as.matrix(expand.grid(x,x))

## DoE
par(mar=c(4.5,5.1,1.5,1.5),cex.axis=1.5,cex.lab=2) 
plot(X,xlab="x1",ylab="x2",pch=4, cex=1,lwd=3)

## Output
Xm = cbind(X,0*X+0.5)
Y = apply(Xm, 1, runExperiment)
Y = matrix(Y[1,])

open3d()
spheres3d(X[,1],X[,2],Y,radius=1)
aspect3d(1,1,1)
axes3d(box=TRUE)
title3d(xlab="x1",ylab="x2",zlab="y",cex=1.5)

## model
XY <- data.frame(x1=X[,1],x2=X[,2],y=Y)
mr <- lm(y ~ 1 + x1 + x2, data = XY)

xg = seq(0,1,length=20)
Xnew = expand.grid(x1=xg,x2=xg)
Yp = predict(mr,Xnew)

persp3d(xg,xg,matrix(Yp,20),col='red',xlab='', ylab='', zlab='',cex=1.5,alpha=0.8)
spheres3d(X[,1],X[,2],Y,radius =1)
title3d(xlab="x1",ylab="x2",zlab="m(x)",cex=1.5)







