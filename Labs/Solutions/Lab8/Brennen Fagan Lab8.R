#Lab 8:
#In which I realize how bizarre it is to pass a method as
#an argument without explicitly declaring it as such, e.g. C++.

#TODO: As mentioned in the previous lab,
#replace this with a noisier version.
doublesweep <- function(A, B, C, F, a, b) {
  # Solves the equation 
  # A[i]*v[i-1] - C[i]*v[i] + B[i]*v[i+1] = F[i]
  # for v[i], i = 1,...,N-1, with boundary conditions
  # v[0]=a and v[N]=b
  
  # Check the lengths of the vectors
  N <- length(C) + 1
  if ((length(B) != N-1) || (length(A) != N-1) || (length(F) != N-1)) {
    stop("The lengths of the vector arguments need to be equal")
  }
  
  alpha <- rep(0, N)
  beta <- rep(0, N)
  beta[1] <- a
  
  #sweep up to get alphas and betas
  for (i in 1:(N-1)) {
    alpha[i+1] <- B[i] / (C[i]-alpha[i]*A[i])
    beta[i+1] <- (beta[i]*A[i] - F[i]) / (C[i] - alpha[i]*A[i])
  }
  
  v <- rep(0, N-1 )
  v[N-1] <- alpha[N]*b + beta[N]
  
  #sweep down to find v's
  for (i in (N-1):2) {
    v[i-1] <- alpha[i]*v[i] + beta[i]    
  }
  
  return(v)
}

backwardDifference2 <- function(u0, K, f, L=1, N=30, T=1, M=30) {
  #Admittedly, really want to place guards here due to
  #passing functions implicitly. It is as if a million
  #compilers cried out at once, and then were silenced!
  
  # set up space grid
  h <- L/N
  #So we do not quite care about the boundaries since 0.
  x <- h*(1:(N-1))
  #But we do need the entire grid now, in order to calculate
  #K. Notice later that we do not return the full grid.
  xFull <- h*(0:N)
  
  # set up time grid
  tau <- T/M
  t <- tau*(0:M)
  
  gamma <- tau/h^2
  
  w <- matrix(0, N-1, M+1)  # Matrix to hold the solution
  w[ , 1] <- u0(x)  # Initial value
  # Loop over time steps
  for (j in 1:M) {
    #K apparently has the function signature
    #Kv<- K(Grid, Cur'nt Time, the entire previous row of points)
    Kv <- K(xFull, t[j], c(0, w[, j], 0)) 
    #With the obvious technicality that it is an abuse to call
    #it on the previous row of points and it should only be called
    #on the current row of points.
    
    #chi here is performing some error correction by averaging
    #points due to the abuse in this method.
    chi <- (Kv[1:N] + Kv[2:(N+1)])/2
    
    #C?
    a <- 1 + gamma*(chi[2:N]+chi[1:(N-1)])
    #A?
    b <- -gamma*chi[1:(N-1)]
    #B?
    #c(b[2:(N-1)],0)
    
                #TODO: Check to see if there is an error here
                #because the function signature is 
                #doublesweep(A,B,C,F,a,b)
                #and we have (b,c,a)
    w[, j+1] <- doublesweep(b[1:(N-1)], c(b[2:(N-1)],0), -a, 
                            #Completely technically correct, 
                            #but not pretty!
                            w[, j]+tau*f(x, t[j], w[, j]), 0, 0)
  }
  
  # Return a list consisting of time grid, x grid and solution
  return(list(x=x, t=t, w=w))
}

sol <- backwardDifference2(u0=function(x) pmax(0, 1-10*abs(x-0.5)), 
                           K=function(x, t, u) u^2/2,
                           f=function(x, t, u) -u)
library("plot3Drgl") 
persp3D(sol$x, sol$t, sol$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4)

#Exercise 1:
sol <- backwardDifference2(u0=function(x) pmax(0, 1-10*abs(x-0.5)), 
                           K=function(x, t, u) x*t*u^2,
                           f=function(x, t, u) 8*x*sin(8*pi*x*t))
library("plot3Drgl") 
persp3D(sol$x, sol$t, sol$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4)
