#Lab 8:
#In which I realize how bizarre it is to pass a method as
#an argument without explicitly declaring it as such, e.g. C++.

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
  if(any(C<=0)||any(B<=0)||any(A<=0)){
    warning("Exists i: Ai or Bi or Ci <= 0")
  }
  if(any(C<A+B)){
    warning("Exists i: Ci<Ai+Bi")
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
    
    #a -> C
    C <- -(1 + gamma*(chi[2:N]+chi[1:(N-1)]))
    #b -> A
    A <- -gamma*chi[1:(N-1)]
    #B:
    #c(b[2:(N-1)],0)
    
                #Note: this gives a ton of warnings!
    w[, j+1] <- doublesweep(A[1:(N-1)], c(A[2:(N-1)],0), C, 
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

backwardDifference3 <- function(u0=function(x) pmax(0,1-10*abs(x-0.5)), 
                                K=function(x, t, u) u^2/2,
                                f=function(x, t, u) -u^2,
                                L=1, N=30, T=1, M=30,
                                max_iteration=20, tolerance=0.000001) {
  #Almost the same setup as above. 
  # set up space grid
  h <- L/N
  x <- h*(1:(N-1))
  xFull <- h*(0:N)
  
  # set up time grid
  tau <- T/M
  t <- tau*(0:M)
  
  gamma <- tau/h^2
  
  w <- matrix(0, N-1, M+1)  # Matrix to hold the solution
  iterations <- rep(0, M)  # to hold number of iterations needed at each step
  precision <- rep(0, M)  # to hold precision achieved
  
  w[ , 1] <- u0(x)  # Initial value
  # Loop over time steps
  for (j in 1:M) {
    wn <- w[, j]  # initial guess for w_{j+1}
    #max iteration guards against running away from the tolerance limit forever.
    for (s in 1:max_iteration) {
      ws <- wn
      #essentially the same, albeit safer, calculation as above.
      Kv <- K(xFull, t[j], c(0, ws, 0))
      chi <- (Kv[1:N] + Kv[2:(N+1)])/2
      a <- 1 + gamma*(chi[2:N]+chi[1:(N-1)])
      b <- -gamma*chi[1:(N-1)]
      # update guess for w_{j+1}
      wn <- doublesweep(b[1:(N-1)], c(b[2:(N-1)],0), -a, 
                        w[, j]+tau*f(x, t[j], ws), 0, 0)
      # break if tolerance limit is reached
      if (max(abs(wn-ws)) < tolerance) {
        break
      }
    }
    iterations[j] <- s
    precision[j] <- max(abs(wn-ws))
    w[, j+1] <- wn
  }
  
  # Return a list consisting of time grid, x grid and solution
  return(list(x=x, t=t, w=w, iterations=iterations, precision=precision))
}
sol <- backwardDifference3(u0=function(x) pmax(0, 1-10*abs(x-0.5)), 
                           K=function(x, t, u) u^2/2,
                           f=function(x, t, u) -u)

persp3D(sol$x, sol$t, sol$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4)

sol <- backwardDifference3(u0=function(x) pmax(0, 1-10*abs(x-0.5)), 
                           K=function(x, t, u) u^2/2,
                           f=function(x, t, u) -u,
                           max_iteration = 1)

persp3D(sol$x, sol$t, sol$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4)

sol <- backwardDifference3(u0=function(x) pmax(0, 1-10*abs(x-0.5)), 
                           K=function(x, t, u) u^2/2,
                           f=function(x, t, u) -u,
                           max_iteration = 2)

persp3D(sol$x, sol$t, sol$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4)

#Exercise 2: Part 2 applied to the simple example:
sol <- backwardDifference3(u0=function(x) pmax(0, 1-10*abs(x-0.5)), 
                           K=function(x, t, u) u^2/2,
                           f=function(x, t, u) -u,
                           max_iteration = 4)

persp3D(sol$x, sol$t, sol$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4)
cat("Precision range of 4 max_iterations:[",min(sol$precision),",",max(sol$precision),"]\n")


sol <- backwardDifference3(u0=function(x) pmax(0, 1-10*abs(x-0.5)), 
                           K=function(x, t, u) x*t*u^2,
                           f=function(x, t, u) 8*x*sin(8*pi*x*t))

persp3D(sol$x, sol$t, sol$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4)
cat("Precision range of Exercise 1: 20 max_iterations:[",min(sol$precision),",",max(sol$precision),"].\n")
cat("0 is an odd value and indicates convergence. Number of iterations:", sol$iterations, ".\n")
cat("So we probably want to neglect the first entry; its precision:", sol$precision[1],".\n")
cat("Precision range of Exercise 1: 20 max_iterations:[",min(sol$precision[0:-1]),",",max(sol$precision[0:-1]),"]\n")

print("Increasing iterations until 1 more than number used.")
sol <- backwardDifference3(u0=function(x) pmax(0, 1-10*abs(x-0.5)), 
                           K=function(x, t, u) x*t*u^2,
                           f=function(x, t, u) 8*x*sin(8*pi*x*t),
                           max_iteration = 23)

persp3D(sol$x, sol$t, sol$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4)
cat("Precision range of Exercise 1: 23 max_iterations:[",min(sol$precision[0:-1]),",",max(sol$precision[0:-1]),"]\n")
cat("Number of iterations:", sol$iterations, "\n")


print("Decreasing iterations until minimum precision is greater than 10^-6.")
sol <- backwardDifference3(u0=function(x) pmax(0, 1-10*abs(x-0.5)), 
                           K=function(x, t, u) x*t*u^2,
                           f=function(x, t, u) 8*x*sin(8*pi*x*t),
                           max_iteration = 5)

persp3D(sol$x, sol$t, sol$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4)
cat("Precision range of Exercise 1: 5 max_iterations:[",min(sol$precision[0:-1]),",",max(sol$precision[0:-1]),"]\n")
cat("Number of iterations:", sol$iterations, "\n")


print("Finally, precision for when there are only 4 iterations.")
sol <- backwardDifference3(u0=function(x) pmax(0, 1-10*abs(x-0.5)), 
                           K=function(x, t, u) x*t*u^2,
                           f=function(x, t, u) 8*x*sin(8*pi*x*t),
                           max_iteration = 4)

persp3D(sol$x, sol$t, sol$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4)
cat("Precision range of Exercise 1: 4 max_iterations:[",min(sol$precision[0:-1]),",",max(sol$precision[0:-1]),"]\n")
cat("Precision:",sol$precision,"\n")
cat("Number of iterations:", sol$iterations, "\n")
