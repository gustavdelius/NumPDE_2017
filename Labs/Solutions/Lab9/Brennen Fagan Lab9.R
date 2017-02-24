list.files()
source('Labs/Solutions/Lab9/doublesweep.R')

library(plot3Drgl)

#Note that this would be fairly obvious to
#write as n-dimensional. The interesting bit
#would be how the language handles it.
ADI <- function(u0, K=1, f, L1=1, N1=30, L2=1, N2=30, T=1, M=30) {
  # set up space grids
  h1 <- L1/N1
  x <- h1*(0:N1)
  h2 <- L2/N2
  y <- h2*(0:N2)
  
  # set up time grid
  tau <- T/M
  t <- tau*(0:M)
  
  gamma1 <- K*tau/h1^2
  gamma2 <- K*tau/h2^2
  # Vectors to be later used in double sweep method
  A1 <- rep(gamma1/2, N1-1)
  C1 <- rep(1+gamma1, N1-1)
  A2 <- rep(gamma2/2, N2-1)
  C2 <- rep(1+gamma2, N2-1)
  
  w <- array(0, dim=c(N1+1, N2+1, M+1))  # Array to hold the solution
  w[, , 1] <- outer(x, y, u0)  # Initial value
  # Loop over time steps
  for (n in 1:M) {
    # Matrix with contributions from inhomogeneous term $\tau/2 f^{n+1/2}$
    Fh = tau*(outer(x, y, f, t=t[n]) + outer(x, y, f, t=t[n+1]))/4
    # first half step
    wh <- matrix(0, nrow=N1+1, ncol=N2+1)  # matrix to hold w^{n+1/2}
    for (j in 2:N2) {
      F1 <- gamma2/2*(w[2:N1, j+1, n] + w[2:N1, j-1, n]) + 
        (1-gamma2)*w[2:N1, j, n] + Fh[2:N1, j]
      wh[2:N1, j] <- doublesweep(A1, A1, C1, F1, 0, 0)
    }
    # second half step
    for (k in 2:N1) {
      F2 <- gamma1/2*(wh[k+1, 2:N2] + wh[k-1, 2:N2]) + 
        (1-gamma1)*wh[k, 2:N2] + Fh[k, 2:N2]
      w[k, 2:N2, n+1] <- doublesweep(A2, A2, C2, F2, 0, 0)
    }
  }
  
  # Return a list consisting of grid and solution
  return(list(x=x, y=y, t=t, w=w))
}

#Note additionally that we have homogenous boundary conditions.
sol <- ADI(u0=function(x, y) sin(pi*x)*sin(pi*y), K=1, 
           f=function(x, y, t) 20*pi^2*sin(2*pi*x)*sin(2*pi*y),
           L1=1, L2=1, N1=40, N2=40, T=0.2, M=20)

persp3D(sol$x, sol$y, sol$w[, , 1],
        xlab="x", ylab="y", zlab="w",
        ticktype="detailed", nticks=4, phi=10, theta=90)

persp3D(sol$x, sol$y, sol$w[, , 21],
        xlab="x", ylab="y", zlab="w",
        ticktype="detailed", nticks=4, phi=10, theta=90)

for (n in 1:21) {
  persp3D(sol$x, sol$y, sol$w[, , n],
          xlab="x", ylab="y", zlab="w", zlim=c(-0.7, 1), clim=c(-0.7, 1),
          ticktype="detailed", nticks=4, phi=12, theta=90)
  title(paste("Solution: n =" , n))
}

ADI2 <- function(u0=function(x, y) sin(pi*x)*sin(pi*y),
                 f=function(x, y, t) 20*pi^2*sin(2*pi*x)*sin(2*pi*y),
                 x_at_0 = function(y,t) t*sin(2*pi*y),
                 x_at_L1 = function(y,t) t*y*(1-y),
                 y_at_0 = function (x,t) 0,
                 y_at_L2 = function (x,t) 0,
                 K=1, L1=1, N1=40, L2=1, N2=40, T=0.2, M=20) {
  # set up space grids
  #X grid
  h1 <- L1/N1
  x <- h1*(0:N1)
  #Y grid
  h2 <- L2/N2
  y <- h2*(0:N2)
  
  # set up time grid
  tau <- T/M
  t <- tau*(0:M)
  
  gamma1 <- K*tau/h1^2
  gamma2 <- K*tau/h2^2
  # Vectors to be later used in double sweep method
  A1 <- rep(gamma1/2, N1-1)
  C1 <- rep(1+gamma1, N1-1)
  A2 <- rep(gamma2/2, N2-1)
  C2 <- rep(1+gamma2, N2-1)
  
  w <- array(0, dim=c(N1+1, N2+1, M+1))  # Array to hold the solution
  w[, , 1] <- outer(x, y, u0)  # Initial value
  # Loop over time steps
  for (n in 1:M) {
    # Matrix with contributions from inhomogeneous term $\tau/2 f^{n+1/2}$
    Fh = tau*(outer(x, y, f, t=t[n]) + outer(x, y, f, t=t[n+1]))/4
    # first half step
    wh <- matrix(0, nrow=N1+1, ncol=N2+1)  # matrix to hold w^{n+1/2}
    for (j in 2:N2) {
      F1 <- gamma2/2*(w[2:N1, j+1, n] + w[2:N1, j-1, n]) + 
        (1-gamma2)*w[2:N1, j, n] + Fh[2:N1, j]
      #Primary Changes ------------------------------V------------------V
      wh[2:N1, j] <- doublesweep(A1, A1, C1, F1, x_at_0(y[j],t[n]), x_at_L1(y[j],t[n]))
      #Make sure to place the boundaries in the result.
      wh[1, j] <- x_at_0(y[j],t[n])
      wh[N1+1,j] <- x_at_L1(y[j],t[n])
    }
    # second half step
    for (k in 2:N1) {
      F2 <- gamma1/2*(wh[k+1, 2:N2] + wh[k-1, 2:N2]) + 
        (1-gamma1)*wh[k, 2:N2] + Fh[k, 2:N2]
      #Primary Changes ------------------------------V------------------V
      w[k, 2:N2, n+1] <- doublesweep(A2, A2, C2, F2, y_at_0(x[j],t[n]), y_at_L2(x[j],t[n]))
      #Make sure to place the boundaries in the result.
      w[k,1,n+1] <- y_at_0(x[j],t[n])
      w[k,N2+1,n+1] <- y_at_L2(x[j],t[n])
    }
  }
  
  # Return a list consisting of grid and solution
  return(list(x=x, y=y, t=t, w=w))
}

#Verify we have done no harm
sol <- ADI2(u0=function(x, y) sin(pi*x)*sin(pi*y), K=1, 
           f=function(x, y, t) 20*pi^2*sin(2*pi*x)*sin(2*pi*y),
           x_at_0 = function(y,t) 0,
           x_at_L1 = function(y,t) 0,
           y_at_0 = function (x,t) 0,
           y_at_L2 = function (x,t) 0,
           L1=1, L2=1, N1=40, N2=40, T=0.2, M=20)

persp3D(sol$x, sol$y, sol$w[, , 1],
        xlab="x", ylab="y", zlab="w",
        ticktype="detailed", nticks=4, phi=10, theta=90)

persp3D(sol$x, sol$y, sol$w[, , 21],
        xlab="x", ylab="y", zlab="w",
        ticktype="detailed", nticks=4, phi=10, theta=90)

#See if we have addressed the question correctly.
sol2 <- ADI2()

for (n in 1:21) {
  persp3D(sol2$x, sol2$y, sol2$w[, , n],
          xlab="x", ylab="y", zlab="w", zlim=c(-0.7, 1), clim=c(-0.7, 1),
          ticktype="detailed", nticks=4, phi=12, theta=90)
  title(paste("Exercise: n =" , n))
}

for (n in 1:21) {
  persp3D(sol2$x, sol2$y, sol2$w[, , n]-sol$w[, , n],
          xlab="x", ylab="y", zlab="w", zlim=c(-0.7, 1), clim=c(-0.7, 1),
          ticktype="detailed", nticks=4, phi=12, theta=90)
  title(paste("Difference made: n =" , n))
}

#As one can see, the difference made is _quite_ subtle.
#If we remove the bounds on zlim and clim, we see the
#the difference more clearly.
persp3D(sol2$x, sol2$y, sol2$w[, , n]-sol$w[, , n],
        xlab="x", ylab="y", zlab="w", 
        ticktype="detailed", nticks=4, phi=0, theta=105)
title(paste("Difference made: n =" , n))
plotrgl()
