library(tikzDevice)

plot.design <- function(X, n = nrow(X), arrows.h = FALSE, arrows.v = FALSE){
  plot(X, xlab="$x_1$", ylab="$x_2$", cex.lab=2,
       xlim=c(0,1), ylim=c(0,1), 
       cex = 2, pch = 19, xaxt="n", yaxt="n")
  s <- seq(0, 1, length.out = n+1)
  abline(h = s, v = s, lty="dotted")
  if (arrows.v){
    arrows(x0 = X[, 1], y0 = X[,2], x1 = X[,1], y1 = rep(0, n), 
           col="blue", cex = 0.5, lty = "dashed")
  }
  if (arrows.h){
    arrows(x0 = X[, 1], y0 = X[,2], x1 = rep(0, n), y1 = X[,2], 
           col="blue", cex = 0.5, lty = "dashed")
  }
}


dimension <- 2
p <- 3; 
x <- ((0:(p-1))+0.5)/p
X <- expand.grid(x, x)

fileName <- 'DoE_grid.tex'
tikz(fileName, standAlone = TRUE, width=5, height=5)
plot.design(X,n=p)
dev.off()
tools::texi2dvi(fileName,pdf=T,clean=TRUE)
file.remove(fileName)

fileName <- 'DoE_gridarrow.tex'
tikz(fileName, standAlone = TRUE, width=5, height=5)
plot.design(X,n=p,arrows.v=TRUE)
dev.off()
tools::texi2dvi(fileName,pdf=T,clean=TRUE)
file.remove(fileName)


###################

library(DiceDesign)

dimension <- 2
n <- 9
X <- lhsDesign(n, dimension=2, seed=3, randomized=FALSE)$design

fileName <- 'DoE_LHS.tex'
tikz(fileName, standAlone = TRUE, width=5, height=5)
par(mar=c(4.5,4.5,.5,.5),cex.axis=1.5,cex.lab=2)
plot.design(X)
dev.off()
tools::texi2dvi(fileName,pdf=T,clean=TRUE)
file.remove(fileName)

fileName <- 'DoE_LHSarrow.tex'
tikz(fileName, standAlone = TRUE, width=5, height=5)
par(mar=c(4.5,4.5,.5,.5),cex.axis=1.5,cex.lab=2)
plot.design(X,arrows.h = TRUE, arrows.v = TRUE)
dev.off()
tools::texi2dvi(fileName,pdf=T,clean=TRUE)
file.remove(fileName)

X <- maximinESE_LHS(X, it=5)$design
fileName <- 'DoE_LHSmaximin.tex'
tikz(fileName, standAlone = TRUE, width=5, height=5)
par(mar=c(4.5,4.5,.5,.5),cex.axis=1.5,cex.lab=2)
plot.design(X)
dev.off()
tools::texi2dvi(fileName,pdf=T,clean=TRUE)
file.remove(fileName)
