# Lab session solutions
library(shiny)
library(DiceDesign)
library(DiceKriging)
library(DiceView)

#########################################
# Part 1 : Kriging model
#########################################

## 1.
runGitHub("shinyApps",username="NicolasDurrande",subdir="catapult") 

X0 <- lhsDesign(n = 16, dimension = 4)$design
Xopt <- maximinESE_LHS(X0, it=10)
X <- Xopt$design
colnames(X) <- c("rotation_axis", "stop", "spring_binding1", "spring_binding2")

Y <- apply(X, 1, runExperiment)[1,]

m0<- km(design=X, response=Y)

## 2.


## 3.


## 4. here's a few useful functions to asses the model quality
MSE <- function(model){
  # returns Leave one Out Mean Square error (the smaller the better)
  pred <- leaveOneOut.km(model)
  mse <- mean((model@y-pred$mean)^2)
  return(mse)
}

Q2 <- function(model){
  # returns Q2 on LOO (usually in [0,1], the closest to 1 the better)
  pred <- leaveOneOut.km(model,type="UK")
  q2 <- 1-sum((model@y-pred$mean)^2)/sum((model@y-mean(model@y))^2)
  return(q2)
}

sd_error <- function(model){
  # returns standard deviation of standardised LOO residuals (the closest to 1 the better)
  pred <- leaveOneOut.km(model,type="UK")
  standardised_residuals <- (model@y-pred$mean)/pred$sd
  return(sd(standardised_residuals))
}

Q2(m0)
sd_error(m0)

## 4. TODO: try improve the model by changing kernel and trend


## 5.
predict_wrapper <- function(x,model){
  pred <- predict(model, matrix(x,ncol=model@d), type="UK",checkNames = FALSE)
  return(pred$mean)
}

opt_model <- optim(X[which.max(Y),],predict_wrapper,model=m2,control=list(fnscale=-1))

# best inputs according to model
opt_model$par

# predicted length according to model
opt_model$value

## 5. TODO : simulator run and comparison with current optimum


#########################################
# Part 2 : Efficient global optimization
#########################################

## 1.
# First maximising f(x) <=> minimising -f(x) 
runExperimentMin <- function(x){
  -runExperiment(x)[1]
}

## 1. TODO: update Y and the kriging model accordingly
#      the new model should be named should be called "model"

## 2. Get the new location of the point maximising the EI
#      this point should be called "newx"

## 3.
newy <- runExperimentMin(newx)

cat("Previous best value:", round(-min(model@y),2),
    "\nNew computed value: ", round(-newy,2),
    "\nImprovement is:      ", round(min(Y) - newy,2), "     (Expected improvement was :", round(EGO1$value,2),")")

# The model can be updated with the new design point and observation
m <- update(model, newx, newy)

## 4.
oEGO <- EGO.nsteps(model = m, fun = runExperimentMin, nsteps = 20, 
                   lower = rep(0, 4), upper = rep(1, 4))

## look at the locations visited by the algorithm
visualizeEGO <- function(initDesign, initValues, EGOpoints, EGOvalues){
  bestIndex <- which.min(EGOvalues)
  y <- c(initValues, EGOvalues, EGOvalues[bestIndex])
  X <- rbind(initDesign, EGOpoints, EGOpoints[bestIndex, ])
  ninit <- nrow(initDesign)
  nsteps <- nrow(EGOpoints)
  pairs(cbind(y, X), 
        col = c(rep("black", ninit), rep("blue", nsteps), "red"),
        pch = c(rep(1, ninit + nsteps), 19))
}

visualizeEGO(initDesign = X, initValues = Y,
             EGOpoints = oEGO$par, EGOvalues = oEGO$val)

plot(c(Y,oEGO$value),main="convergence",xlab="evaluation number",ylab="Y values")
lines(rep(length(Y),2),range(Y,oEGO$value),lty=2,col="gray")
lines(length(Y)+0:length(oEGO$value),c(min(Y),cummin(oEGO$value)),col="red",lwd=2)

# Look at the last model around the optimum
sectionview(oEGO$lastmodel, center = oEGO$par[bestPoint, ])
par()

## 5.
bestPoint <- which.min(oEGO$value)
cat("longest shot observed:",-round(oEGO$value[bestPoint],2),
    "\ncorresponding input values:",round(oEGO$par[bestPoint,],2))

# plot best shot:
runExperiment(oEGO$par[bestPoint, ])

