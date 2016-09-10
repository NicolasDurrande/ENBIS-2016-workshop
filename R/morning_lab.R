# Lab session: file to be completed

#########################################
# Part 1 : catapult numerical simulator
#########################################

library(shiny)
runGitHub("shinyApps",username="NicolasDurrande",subdir="catapult") 

#########################################
# Part 2: Design of experiments
#########################################

## 1.

## 2.

## 3. 

## 4. plot Y against the inputs
par(mfrow=c(2,2)) # split the plotting area in 4
plot(X[,1],Y)
plot(X[,2],Y)
plot(X[,3],Y)
plot(X[,4],Y)
par(mfrow=c(1,1))

library(rgl)
plot3d(x = X[,1], y = X[,2], z = Y)

## 5. bests inputs so far

#########################################
## Kriging
#########################################

## 1. create GP model
library(DiceKriging)
?km # or help("km")

## 2. create GP model

## 3. visualization with DiceView

## 4. 