#Lab 10:
#From Lab prompt.        Forcing               Initial Condition  Initial Condition'
explicitWave <- function(F=function(x, t) 0*x, f=function(x) 0*x, g=function(x) 0*x,
                         alpha=1, a=0, b=1, N=30, T=1, M=30) {
  # set up space grid
  h <- (b-a)/N
  x <- a + h*(1:(N-1))
  xLong <- c(a, x, b)  # includes the endpoints
  
  # set up time grid
  tau <- T/M
  t <- tau*(0:M)
  
  # Set up evolution matrix eq.(4.10)
  gs <- (alpha*tau/h)^2
  A <- diag(2-2*gs, N-1)
  for (k in 1:(N-2)) {
    A[k,k+1] <- gs
    A[k+1,k] <- gs
  }
  
  w <- matrix(0, N-1, M+1)  # Matrix to hold the solution
  
  # Initial conditions
  w[, 1] <- f(x)  # Initial value w_0
  fpp <- (f(x-h) -2*f(x) + f(x+h))/h^2  # Approximate derivative of f
  w[, 2] <- f(x) + tau*g(x) + tau^2/2*(alpha^2*fpp + F(x,0))  # eq.(4.14) w_1
  
  # Loop over time steps
  for (j in 2:M) {
    # Use eq.(4.9)
    w[, j+1] <- A %*% w[, j] - w[, j-1] + tau^2 * F(x, t[j])
  }
  
  # Return a list consisting of time grid, x grid and solution
  # Boundary Conditions:            v     v
  return(list(x=xLong, t=t, w=rbind(0, w, 0)))
}

#As a reminder, in order to double T, you do also need to double M.
sol <- explicitWave(f=function(x) ifelse(abs(x)<1/4, (1+cos(4*pi*x))/2, 0), 
                    a=-1, N=80, T=8, M=320)
library("plot3Drgl") 
persp3D(sol$x, sol$t, sol$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4)
title("Solution 1")

sol2 <- explicitWave(f = function(x) ifelse(abs(x)<1/4, (1+cos(4*pi*x))/2, 0),
                     g = function(x) ifelse(abs(x)<1/4, 2*pi*sin(4*pi*x), 0), 
                     a=-1, N=80, T=4, M=160)

persp3D(sol2$x, sol2$t, sol2$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4)
title("Exercise 1")


sol2 <- explicitWave(f = function(x) ifelse(abs(x)<1/4, (sin(4*pi*x))/2, 0),
                     g = function(x) ifelse(abs(x)<1/4, -2*cos(4*pi*x), 0), 
                     a=-1, N=80, T=4, M=160)

persp3D(sol2$x, sol2$t, sol2$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4,
        clim=c(-0.7, 0.7), zlim=c(-0.7, 0.7))
title("Exercise 2")

sol3 <- explicitWave(f = function(x) ifelse(abs(x)<1/4, (sin(4*pi*x))/2, 0),
                     g = function(x) ifelse(abs(x)<1/4, -2*cos(4*pi*x), 0), 
                     a=-1, N=80, T=4, M=320)

persp3D(sol3$x, sol3$t, sol3$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4,
        clim=c(-0.7, 0.7), zlim=c(-0.7, 0.7))
title("Exercise 3")
#Exercise 2 appears to be more continuous than Exercise 3.
#As discussed in practical, this is due to the error cancellation, similar
#to problem 6 in exercises 1, in which a very specific choice of tau and h
#created error cancellation. Exercise 2 is what the actual solution should
#look like, but Exercise 3 is what the approximation should look like.

sol4 <- explicitWave(f = function(x) ifelse(abs(x)<1/4, (sin(4*pi*x))/2, 0),
                     g = function(x) ifelse(abs(x)<1/4, -2*cos(4*pi*x), 0), 
                     a=-1, N=160, T=4, M=160)

persp3D(sol4$x, sol4$t, sol4$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4
        )
title("Exercise 4")
#Instabilities arise due to the explicit nature and the outpacing of the
#space positions over the time positions. Note that constraining the boundaries
#shows that it initially looks "right" before quickly disappearing off the
#boundaries.

