library(tikzDevice)
#library(MASS)
library(rgl)
source('../../../catapult/catapultFunctions.R')
source('../../../catapult/catapultSettings.R')

x = c(0,1)
X = as.matrix(expand.grid(x,x,x,x))

Y = apply(X, 1, runExperiment)
Y = matrix(Y[1,])

round(Y,1)


################################
## linear regression
# 
# XY <- data.frame(x1=X[,1],x2=X[,2],x3=X[,3],x4 = X[,4],y=Y)
# mrt <- lm(y ~ ., data = XY)   
# summary(mrt)
# 
# ################################
# ## prediction 
# Xnew = data.frame(x1=.5, x2=.5, x3=.5, x4=.5)
# 
# predict(mr,Xnew)
# runExperiment(as.matrix(Xnew,nrow=1))

##################################################################
## better model in 2D

x = c(0,0.33,0.66,1)
X = as.matrix(expand.grid(x,x))

## DoE
fileName <- 'linReg_DoE.tex'
tikz(fileName, standAlone = TRUE, width=5, height=5)
par(mar=c(4.5,5.1,1.5,1.5),cex.axis=1.5,cex.lab=2) 
plot(X,xlab="$x_1$",ylab="$x_2$",pch=4, cex=1,lwd=3)
dev.off()
tools::texi2dvi(fileName,pdf=T,clean=TRUE)
file.remove(fileName)

## Output
Xm = cbind(X,0*X+0.5)
Y = apply(Xm, 1, runExperiment)
Y = matrix(Y[1,])

open3d()
spheres3d(X[,1],X[,2],Y,radius=1)
aspect3d(1,1,1)
axes3d(box=TRUE)
title3d(xlab="x1",ylab="x2",zlab="y",cex=1.5)
#snapshot3d("linReg_output")

## model
XY <- data.frame(x1=X[,1],x2=X[,2],y=Y)
mr <- lm(y ~ 1 + x1 + x2 + I((x1-x2)^2) + I((x1-x2)^3), data = XY)
# mr <- lm(y ~ 1 + I(x1-x2) + I((x1-x2)^2) + I((x1-x2)^3), data = XY)   

round(mr$coefficients,2)

summary(mr)


xg = seq(0,1,length=20)
Xnew = expand.grid(x1=xg,x2=xg)
Yp = predict(mr,Xnew)

persp3d(xg,xg,matrix(Yp,20),col='red',xlab='', ylab='', zlab='',cex=1.5,alpha=0.8)
spheres3d(X[,1],X[,2],Y,radius =1)
title3d(xlab="x1",ylab="x2",zlab="m(x)",cex=1.5)
#snapshot3d("linReg_model2.png")


surface3d(seq(0,1,0.05),seq(0,1,0.05),matrix(predmod,length(seq(0,1,0.05))),col='lightblue',xlab='x1', ylab='x2', zlab='f(x)')


points3d(X[,1],X[,2],Y,radius =1)



sigmoid <- function(x){
  1/(1+exp(-x))
}

treeSigmoid = function(x,W,B){
  n = nrow(x)
  S1 = x
  for(i in 2:length(W)){
    nbL = nrow(W[[i]])
    S2 = matrix(0,n,nbL)
    for(j in 1:nbL){
      S2[,j] = sigmoid(S1 %*% W[[i]][j,]+B[[i]][j])[,1]
    }
    S1 = S2
  }
  return(S2)
}

vec2list <- function(x,nlayer){
  j <- 1
  w = list(rep(1,nlayer[1]))
  b = list(rep(0,nlayer[1]))
  for(i in 2:length(nlayer)){
    j2 <- j-1+nlayer[i]*nlayer[i-1]
    j3 <- j2+nlayer[i]
    w = c(w,list(matrix(x[j:j2],nlayer[i])))
    b = c(b,list(matrix(x[(j2+1):j3])))
    j <- j3+1
  } 
  return(list(w,b))
}

dimParam <- function(nblayer){
  dim <- 0
  for(i in 2:length(nlayer)){
    dim <- dim + nlayer[i]*nlayer[i-1] + nlayer[i]
  } 
  return(dim)
}

RSS <- function(x){
  WB <- vec2list(x,nlayer)
  sum((treeSigmoid(Xtrain,WB[[1]],WB[[2]])-Ytrain)^2)
}

##################################################################"
### sigmoide


# tikz('neuralnet_sigmoid.tex', standAlone = TRUE, width=5, height=4)
par(mar=c(4.5,5.1,1.5,1.5),cex.axis=1.5,cex.lab=2))
plot(x,sigmoid(x), ylab='sigmoid function',col=darkBlue, type='l', lwd = 2))
# dev.off()
# tools::texi2dvi('neuralnet_sigmoid.tex',pdf=T)


##################################################################"
### example of function:

nlayer <- c(2,4,5,1)
dim <- dimParam(nlayer)
WB = vec2list(runif(dim,-5,5),nlayer)

W = WB[[1]]
B = WB[[2]]

Y = treeSigmoid(Xg,W,B)

persp3d(seq(0,1,0.05),seq(0,1,0.05),matrix(Y,length(seq(0,1,0.05))),col='red',xlab='x1', ylab='x2', zlab='f(x)')
#rgl.snapshot('neuralnetwork.png')

###############################
## fit

Xtrain <- matrix(runif(15*2),ncol=2)
Ytrain <- treeSigmoid(Xtrain,W,B)

mod <- optim(rep(0,43),RSS,control = list(maxit=15000))
WBmod <- vec2list(mod$par,nlayer)
  
predmod <- treeSigmoid(Xg,WBmod[[1]],WBmod[[2]])

persp3d(seq(0,1,0.05),seq(0,1,0.05),matrix(Y,length(seq(0,1,0.05))),col='red',xlab='x1', ylab='x2', zlab='f(x)')
surface3d(seq(0,1,0.05),seq(0,1,0.05),matrix(predmod,length(seq(0,1,0.05))),col='lightblue',xlab='x1', ylab='x2', zlab='f(x)')
spheres3d(Xtrain[,1],Xtrain[,2],Ytrain,radius =0.02)

#rgl.snapshot('neuralnetworkfit.png')

