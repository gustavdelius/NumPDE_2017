#Lab 10:
#From Lab prompt.        Forcing               Initial Condition  Initial Condition'
#NOTE: Multiplying by x still due to the need for a return vector in the function itself.
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

#From prompt:
implicitWave <- function(F=function(x, t) 0*x, 
                         f=function(x) 0*x, g=function(x) 0*x,
                         alpha=1, a=0, b=1, N=30, T=1, M=30, sigma=0.5) {
  # set up space grid
  h <- (b-a)/N
  x <- a + h*(1:(N-1))
  xLong <- c(a, x, b)  # includes the endpoints
  
  # set up time grid
  tau <- T/M
  t <- tau*(0:M)
  
  # Set up matrices (these were not explicitly written in the notes.)
  gs <- (alpha*tau/h)^2
  A <- diag(1+2*sigma*gs, N-1)
  B <- diag(2-2*(1-2*sigma)*gs, N-1)
  C <- diag(-1-2*sigma*gs, N-1)
  for (k in 1:(N-2)) {
    A[k,k+1] <- -sigma*gs
    A[k+1,k] <- -sigma*gs
    B[k,k+1] <- (1-2*sigma)*gs
    B[k+1,k] <- (1-2*sigma)*gs
    C[k,k+1] <- sigma*gs
    C[k+1,k] <- sigma*gs
  }
  Ainv <- solve(A)
  AinvB <- Ainv %*% B
  AinvC <- Ainv %*% C
  
  w <- matrix(0, N-1, M+1)  # Matrix to hold the solution
  
  # Initial conditions
  w[, 1] <- f(x)  # Initial value
  fpp <- (f(x-h) -2*f(x) + f(x+h))/h^2  # Approximate derivative of f
  w[, 2] <- f(x) + tau*g(x) + tau^2/2*(alpha^2*fpp + F(x,0))  # eq.(4.14)
  
  # Loop over time steps
  #w_j+1 = A^-1*B*w_j + A^-1*C*w_j-1+tau^2*A^-1*F
  #A*w_j+1 = B*w_j + C*w_j-1 + tau^2*F
  #d^2(x) = w_k-1,j -2w_k,j + w_k+1,j
  #w_j+1 - (alpha*tau/h)^2*(sigma*d^2 (w_j+1)) = 2 w_j - w_j-1 + (alpha*tau/h)^2*((1-2sigma)d^2(w_j) + sigma*d^2(w_j-1)) + tau^2*F 
  #w_j+1 - 2w_j + w_j-1 - (alpha*tau/h)^2*(sigma*d^2 (w_j+1) + (1-2sigma)d^2(w_j) + sigma*d^2(w_j-1)) = tau^2*F 
  #w_j+1 - 2w_j + w_j-1 - gamma^2*(sigma*d^2 (w_j+1) + (1-2sigma)d^2(w_j) + sigma*d^2(w_j-1)) = tau^2*F 
  for (j in 2:M) {
    w[, j+1] <- AinvB %*% w[, j] +AinvC %*% w[, j-1] + tau^2 * Ainv %*% F(x, t[j])
  }
  
  # Return a list consisting of time grid, x grid and solution
  return(list(x=xLong, t=t, w=rbind(0, w, 0)))
}

#Exercise 5:
sol5 <- implicitWave(f = function(x) ifelse(abs(x)<1/4, (sin(4*pi*x))/2, 0),
                     g = function(x) ifelse(abs(x)<1/4, -2*cos(4*pi*x), 0), 
                     a=-1, N=80, T=4, M=160)

persp3D(sol5$x, sol5$t, sol5$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4
)
title("Exercise 5a")
#Much noisier than expected!

sol5 <- implicitWave(f = function(x) ifelse(abs(x)<1/4, (sin(4*pi*x))/2, 0),
                     g = function(x) ifelse(abs(x)<1/4, -2*cos(4*pi*x), 0), 
                     a=-1, N=80, T=4, M=320)

persp3D(sol5$x, sol5$t, sol5$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4
)
title("Exercise 5b")
#Much safer though!


sol5 <- implicitWave(f = function(x) ifelse(abs(x)<1/4, (sin(4*pi*x))/2, 0),
                     g = function(x) ifelse(abs(x)<1/4, -2*cos(4*pi*x), 0), 
                     a=-1, N=160, T=4, M=160)

persp3D(sol5$x, sol5$t, sol5$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4
)
title("Exercise 5c")
#Space doesn't look as important as time.


sol5 <- implicitWave(f = function(x) ifelse(abs(x)<1/4, (sin(4*pi*x))/2, 0),
                     g = function(x) ifelse(abs(x)<1/4, -2*cos(4*pi*x), 0), 
                     a=-1, N=1000, T=4, M=1000)

persp3D(sol5$x, sol5$t, sol5$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4
)
title("Exercise 5d")
#Even with significant investment, we still do not recover a solution that is as nice.
#Clearly, eliminating terms when possible can be very powerful, but, as mentioned
#in the practical session, this is unusual and not very practical.

#Exercise 6:
#Experiment with adding inhomogeneous terms to the wave equation.
#Since it says inhomogeneous terms as opposed to inhomogeneous conditions,
#we examine changes in the sense that have already been incorporated into
#implicitWave with the F forcing function.
sol6 <- implicitWave(F = function(x, t) 0*x,
                     f = function(x) sin(pi*x),
                     g = function(x) pi*cos(pi*x),
                     a = -1, N = 80, T = 8, M = 320)

persp3D(sol6$x, sol6$t, sol6$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4
)
title("Exercise 6a")
sol6 <- implicitWave(F = function(x, t) 0*x+1,
                     f = function(x) sin(pi*x),
                     g = function(x) pi*cos(pi*x),
                     a = -1, N = 80, T = 8, M = 320)

persp3D(sol6$x, sol6$t, sol6$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4
)
title("Exercise 6b")
#Can still see the wave behaviour, but the forcing is obscuring it
#by interacting with the wave differently.
sol6 <- implicitWave(F = function(x, t) -10*sin(pi*t*1/4)*x,
                     f = function(x) sin(pi*x),
                     g = function(x) pi*cos(pi*x),
                     a = -1, N = 80, T = 8, M = 320)

persp3D(sol6$x, sol6$t, sol6$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4
)
title("Exercise 6c")
#Seems like the tau^2 factor damps this more than I expected.
#Also demonstrates the difficulty of deciphering waves when they
#are layered.
sol6 <- implicitWave(F = function(x, t) -15*x,
                     f = function(x) sin(pi*x),
                     g = function(x) pi*cos(pi*x),
                     a = -1, N = 80, T = 8, M = 320)

persp3D(sol6$x, sol6$t, sol6$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4
)
title("Exercise 6d")
#Reversal of orientation of minor features.
sol6 <- implicitWave(F = function(x, t) ifelse(0*x+t<.75, 2+0*x, 0*x),
                     f = function(x) sin(pi*x),
                     g = function(x) pi*cos(pi*x),
                     a = -1, N = 80, T = 8, M = 320)

persp3D(sol6$x, sol6$t, sol6$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4
)
title("Exercise 6e")
#Temporary perturbations can have long lasting effects.
#Piecewise
sol6 <- implicitWave(F = function(x, t) ifelse(abs(x)<1/4, (abs(x)-.5)*10, cos(pi*x)),
                     f = function(x) sin(pi*x),
                     g = function(x) pi*cos(pi*x),
                     a = -1, N = 80, T = 8, M = 320)

persp3D(sol6$x, sol6$t, sol6$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4
)
title("Exercise 6f")

#Piecewise discontinuous
sol6 <- implicitWave(F = function(x, t) ifelse(abs(x)<1/4, (abs(x)-1)*7, abs(x)/x*exp(-t)),
                     f = function(x) sin(pi*x),
                     g = function(x) pi*cos(pi*x),
                     a = -1, N = 80, T = 8, M = 320)

persp3D(sol6$x, sol6$t, sol6$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4
)
title("Exercise 6g")
#It can clearly handle irregular functions.

#Exercise 7:
#Modify the methods to allow variable wave speed and experiment.
#Wave speed refers to alpha:
implicitWave2 <- function(F=function(x, t) 0*x, 
                         f=function(x) 0*x, g=function(x) 0*x,
                         alpha=function(x,t) 0*x+1, 
                         a=0, b=1, N=30, T=1, M=30, sigma=0.5) {
  # set up space grid
  h <- (b-a)/N
  x <- a + h*(1:(N-1))
  xLong <- c(a, x, b)  # includes the endpoints
  
  # set up time grid
  tau <- T/M
  t <- tau*(0:M)
  
  w <- matrix(0, N-1, M+1)  # Matrix to hold the solution
  
  # Initial conditions
  w[, 1] <- f(x)  # Initial value
  fpp <- (f(x-h) -2*f(x) + f(x+h))/h^2  # Approximate derivative of f
  w[, 2] <- f(x) + tau*g(x) + tau^2/2*(alpha(x, t[2])^2*fpp + F(x,0))  # eq.(4.14)
  
  # Loop over time steps
  #w_j+1 = A^-1*B*w_j + A^-1*C*w_j-1+tau^2*A^-1*F
  #A*w_j+1 = B*w_j + C*w_j-1 + tau^2*F
  #d^2(x) = w_k-1,j -2w_k,j + w_k+1,j
  #w_j+1 - (alpha*tau/h)^2*(sigma*d^2 (w_j+1)) = 2 w_j - w_j-1 + (alpha*tau/h)^2*((1-2sigma)d^2(w_j) + sigma*d^2(w_j-1)) + tau^2*F 
  #w_j+1 - 2w_j + w_j-1 - (alpha*tau/h)^2*(sigma*d^2 (w_j+1) + (1-2sigma)d^2(w_j) + sigma*d^2(w_j-1)) = tau^2*F 
  #w_j+1 - 2w_j + w_j-1 - gamma^2*(sigma*d^2 (w_j+1) + (1-2sigma)d^2(w_j) + sigma*d^2(w_j-1)) = tau^2*F 
  for (j in 2:M) {
    # Set up matrices (these were not explicitly written in the notes.)
    #These matrices, unfortunately, are now dependent on alpha as a function,
    #and hence will need to be evaluated at each time step.
    #The laziest implementation:
    gs <- (alpha(x, t[j])*tau/h)^2
    #Note that gs is a vector that is placed
    #along the diagonal now.
    A <- diag(1+2*sigma*gs, N-1)
    B <- diag(2-2*(1-2*sigma)*gs, N-1)
    C <- diag(-1-2*sigma*gs, N-1)
    for (k in 1:(N-2)) {
      A[k,k+1] <- -sigma*gs[k]
      A[k+1,k] <- -sigma*gs[k+1]
      B[k,k+1] <- (1-2*sigma)*gs[k]
      B[k+1,k] <- (1-2*sigma)*gs[k+1]
      C[k,k+1] <- sigma*gs[k]
      C[k+1,k] <- sigma*gs[k+1]
    }
    Ainv <- solve(A)
    AinvB <- Ainv %*% B
    AinvC <- Ainv %*% C
    
    #Could potentially use doublesweep:
    #A*w[,j+1] = B * w[,j] + C * w[, j-1] + tau^2 * F(x,t[j])
    #Use the RHS as your (A)F, take (A)A, (A)B, and (A)C in the usual way.
    w[, j+1] <- AinvB %*% w[, j] +AinvC %*% w[, j-1] + tau^2 * Ainv %*% F(x, t[j])
  }
  
  # Return a list consisting of time grid, x grid and solution
  return(list(x=xLong, t=t, w=rbind(0, w, 0)))
}

#Verify did no harm:
#Exercise 5:
sol5 <- implicitWave2(f = function(x) ifelse(abs(x)<1/4, (sin(4*pi*x))/2, 0),
                     g = function(x) ifelse(abs(x)<1/4, -2*cos(4*pi*x), 0), 
                     a=-1, N=80, T=4, M=160)

persp3D(sol5$x, sol5$t, sol5$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4
)
title("Exercise 5a")

sol7 <- implicitWave2(f = function(x) ifelse(abs(x)<1/4, (1+cos(4*pi*x))/2, 0),
                     g = function(x) ifelse(abs(x)<1/4, 2*pi*sin(4*pi*x), 0),
                     a=-1, N=80, T=4, M=160)

persp3D(sol7$x, sol7$t, sol7$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4)
title("Exercise 7a")

#Time
sol7 <- implicitWave2(f = function(x) ifelse(abs(x)<1/4, (1+cos(4*pi*x))/2, 0),
                      g = function(x) ifelse(abs(x)<1/4, 2*pi*sin(4*pi*x), 0),
                      alpha = function(x,t) ifelse(0*x + t<3, 1, 0),
                      a=-1, N=80, T=4, M=160)

persp3D(sol7$x, sol7$t, sol7$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4)
title("Exercise 7b")
#we set utt = 0 <=> ut = constant <=> u experiences linear growth with time, no change in space

sol7 <- implicitWave2(f = function(x) ifelse(abs(x)<1/4, (1+cos(4*pi*x))/2, 0),
                      g = function(x) ifelse(abs(x)<1/4, 2*pi*sin(4*pi*x), 0),
                      alpha = function(x,t) ifelse(0*x + t<2, 1, 2),
                      a=-1, N=80, T=4, M=160)

persp3D(sol7$x, sol7$t, sol7$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4)
title("Exercise 7c")
plotrgl(lighting = TRUE)
#Appears as if the sudden change in speed causes errors and opposing waves to form.

sol7 <- implicitWave2(f = function(x) ifelse(abs(x)<1/4, (1+cos(4*pi*x))/2, 0),
                      g = function(x) ifelse(abs(x)<1/4, 2*pi*sin(4*pi*x), 0),
                      alpha = function(x,t) ifelse(0*x + t<2, 1, t^-1),
                      a=-1, N=80, T=4, M=160)

persp3D(sol7$x, sol7$t, sol7$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4)
title("Exercise 7d")
plotrgl(lighting = TRUE)

#Space
#Propagate normally on right, but slowly on left
sol7 <- implicitWave2(f = function(x) ifelse(abs(x)<1/4, (1+cos(4*pi*x))/2, 0),
                      g = function(x) ifelse(abs(x)<1/4, 2*pi*sin(4*pi*x), 0),
                      alpha = function(x,t) ifelse(x>0, 1,.5),
                      a=-1, N=80, T=4, M=160)

persp3D(sol7$x, sol7$t, sol7$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4)
title("Exercise 7e")
plotrgl(lighting = TRUE)

#Propagate quickly further away from origin.
sol7 <- implicitWave2(f = function(x) ifelse(abs(x)<1/4, (1+cos(4*pi*x))/2, 0),
                      g = function(x) ifelse(abs(x)<1/4, 2*pi*sin(4*pi*x), 0),
                      alpha = function(x,t) x,
                      a=-1, N=80, T=4, M=160)

persp3D(sol7$x, sol7$t, sol7$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4)
title("Exercise 7fi")
plotrgl(lighting = TRUE)
#Can see the 0 point where nothing moves, but once outside of that point, the 0 points
#appear to induce waves

sol7 <- implicitWave2(f = function(x) ifelse(abs(x)<1/4, (1+cos(4*pi*x))/2, 0),
                      g = function(x) ifelse(abs(x)<1/4, 2*pi*sin(4*pi*x), 0),
                      alpha = function(x,t) abs(x),
                      a=-1, N=80, T=4, M=160)

persp3D(sol7$x, sol7$t, sol7$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4)
title("Exercise 7fii")
plotrgl(lighting = TRUE)
#Interesting how little difference the change in sign makes.

sol7 <- implicitWave2(f = function(x) ifelse(abs(x)<1/4, (1+cos(4*pi*x))/2, 0),
                      g = function(x) ifelse(abs(x)<1/4, 2*pi*sin(4*pi*x), 0),
                      alpha = function(x,t) -x,
                      a=-1, N=80, T=4, M=160)

persp3D(sol7$x, sol7$t, sol7$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4)
title("Exercise 7fiii")
plotrgl(lighting = TRUE)

#0 growth bands with increasing growth based on distance to origin
sol7 <- implicitWave2(f = function(x) ifelse(abs(x)<1/4, (1+cos(4*pi*x))/2, 0),
                      g = function(x) ifelse(abs(x)<1/4, 2*pi*sin(4*pi*x), 0),
                      alpha = function(x,t) ifelse(abs(x)<.75 && abs(x)>.5, 0, x+1),
                      a=-1, N=80, T=4, M=160)

persp3D(sol7$x, sol7$t, sol7$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4)
title("Exercise 7gi")
plotrgl(lighting = TRUE)
#Broke due to double condition?

sol7 <- implicitWave2(f = function(x) ifelse(abs(x)<1/4, (1+cos(4*pi*x))/2, 0),
                      g = function(x) ifelse(abs(x)<1/4, 2*pi*sin(4*pi*x), 0),
                      alpha = function(x,t) ifelse(abs(x)>.6, 0, x+1),
                      a=-1, N=80, T=4, M=160)

persp3D(sol7$x, sol7$t, sol7$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4)
title("Exercise 7gii")
plotrgl(lighting = TRUE)
#Appears so

sol7 <- implicitWave2(f = function(x) ifelse(abs(x)<1/4, (1+cos(4*pi*x))/2, 0),
                      g = function(x) ifelse(abs(x)<1/4, 2*pi*sin(4*pi*x), 0),
                      alpha = function(x,t) ifelse(abs(x)>.5, ifelse(abs(x)<.75,0,x+1), x+1),
                      a=-1, N=80, T=4, M=160)

persp3D(sol7$x, sol7$t, sol7$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4)
title("Exercise 7giii")
plotrgl(lighting = TRUE)
#Emergence of the second wave is interesting. Error or correct?

#0 space square?
sol7 <- implicitWave2(f = function(x) ifelse(abs(x)<1/4, (1+cos(4*pi*x))/2, 0),
                      g = function(x) ifelse(abs(x)<1/4, 2*pi*sin(4*pi*x), 0),
                      alpha = function(x,t) ifelse(0*x + abs(t-2)<1, ifelse(abs(x)<.5, 0, x+1),x+1),
                      a=-1, N=80, T=4, M=160)

persp3D(sol7$x, sol7$t, sol7$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4)
title("Exercise 7h")
plotrgl(lighting = TRUE)
#Can actually see the square

#Mixed space, time control. Avoiding permanent 0 growth bands.
sol7 <- implicitWave2(f = function(x) ifelse(abs(x)<1/4, (1+cos(4*pi*x))/2, 0),
                      g = function(x) ifelse(abs(x)<1/4, 2*pi*sin(4*pi*x), 0),
                      alpha = function(x,t) (cos(2*pi*x)+10^-1)*cos(t*pi),
                      a=-1, N=80, T=4, M=160)

persp3D(sol7$x, sol7$t, sol7$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4)
title("Exercise 7i")
plotrgl(lighting = TRUE)
#Very odd sharp points. Not quite sure as to explanation.