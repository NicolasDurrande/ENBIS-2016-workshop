# Lab session solutions

#########################################
# Part 1 : catapult numerical simulator
#########################################

library(shiny)
runGitHub("shinyApps",username="NicolasDurrande",subdir="catapult") 
# runApp(appDir='../../catapult') # run the app locally, the specified directory (appDir) must contain the files ui.R et server.R
# runUrl("https://sites.google.com/site/nicolasdurrandehomepage/catapult.zip")

#########################################
# Part 2: Design of experiments
#########################################

## 1.
library(DiceDesign)
X0 <- lhsDesign(n = 16, dimension = 4)$design
Xopt <- maximinESE_LHS(X0, it=10)

## you may be interested in the convergence
# plot(Xopt$critValues,type="l")

X <- Xopt$design
colnames(X) <- c("rotation_axis", "stop", "spring_binding1", "spring_binding2")
round(X,2)

## You can save X into a csv file
#write.table(X,"DoE.csv")
#Xsaved <- read.table("DoE.csv")

## 2.
pairs(X)

## 3. compute one output values
Y <- runExperiment(X[1,])

# compute directly all Y
Y <- apply(X, 1, runExperiment)
Y <- Y[1,]

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
X[which.max(Y),]

#########################################
## Kriging
#########################################

## 1. create GP model
library(DiceKriging)
?km # or help("km")

## 2. create GP model
m0<- km(design=X, response=Y)

## 3. visualization with DiceView
library(DiceView)
sectionview(m0, center = c(0.5, 0.5, 0.5, 0.5))
#sectionview3d(m0, center = c(0, 1, 0.5, 0.5))

# same thing but centred on the optimum
sectionview(m0, center = X[which.max(Y),])
#sectionview3d(m0, center = X[which.max(Y),])

## 4. The variables x3 and x4 (both are spring binding) don't seem to be very influencial