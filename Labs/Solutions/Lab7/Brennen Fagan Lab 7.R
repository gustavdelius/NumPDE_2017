library(plot3Drgl)

#"...forwardDifference() function from lab 4:

forwardDifference <- function(f=function(x, t) 0,
                              u0=function(x) 2*sin(2*pi*x),
                              a = function(t) 0, b = function(t) 0,
                              K=1, L=1, N=30, T=0.1, M=200) {
  # set up space grid
  h <- L/N
  #Appears we are modifying the code to
  #cover the whole grid now as well!
  #Old grid h*(1:(N-1))
  x <- h*(0:N)
  
  # set up time grid
  tau <- T/M
  t <- tau*(0:M)
  
  # set up vectors with initial condition and inhomogenous term
  w <- u0(x)
  
  # Set up evolution matrix
  gamma <- K*tau/(h^2)
  # Adjust for the matrix size
  A <- diag(1-2*gamma, N+1)
  for (k in 1:(N)) {
    A[k,k+1] <- gamma
    A[k+1,k] <- gamma
  }
  
  #Adjust matrix size as well
  Temperature <- matrix(0, N+1, M+1)  # Matrix to hold the solution
  Temperature[ , 1] <- w  # Initial value
  # Loop over time steps
  for (j in 1:M) {
    #Evaluate f at current points and time
    #Assumes f can take a vector input!
    F <- f(x,t[j])
    w <- A %*% w + tau * F
    #Since we are now covering the whole grid
    #We also change our "boundaries" and force
    #them to the boundary conditions
    #Q: Confused as to the t[j+1] though.
    #Why is it here j+1, above j? Shouldn't
    #above be j+1 as well?
    #A: Looking back at the notes, we are using
    #equation 2.61, which has time of F the same
    #as the time for w when computing the next w.
    #So these make sense.
    w[1] <- a(t[j+1])
    #N-1 -> N+1
    w[N+1] <- b(t[j+1])
    Temperature[ , j+1] <- w
  }
  
  # Return a list consisting of time grid, x grid and solution
  return(list(x=x, t=t, w=Temperature))
}#..."

#Says 360 in the notes, but uses 720 as an example. Might as
#well see how much of a difference it makes!
sol <- forwardDifference(u0=function(x) x, 
                         f=function(x,t) 10*sin(10*pi*t)*exp(-10*(x-0.5)^2),
                         a=function(t) sin(20*pi*t), b=function(t) cos(20*pi*t),
                         K=1, L=1, T=0.4, N=30, M=360)
persp3D(sol$x, sol$t, sol$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4)
#Conditional Stability strikes again!

sol <- forwardDifference(u0=function(x) x, 
                         f=function(x,t) 10*sin(10*pi*t)*exp(-10*(x-0.5)^2),
                         a=function(t) sin(20*pi*t), b=function(t) cos(20*pi*t),
                         K=1, L=1, T=0.4, N=30, M=720)
persp3D(sol$x, sol$t, sol$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4)
plotrgl(smooth = TRUE)

#Recall that this version does not control for errors
#in the values of the vectors, as a warning!
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

backwardDifference <- function(u0=function(x) 2*sin(2*pi*x), 
                               K=1, L=1, N=30, T=0.1, M=30) {
  # set up space grid
  h <- L/N
  x <- h*(1:(N-1))
  
  # set up time grid
  tau <- T/M
  t <- tau*(0:M)
  
  # set up vectors with initial condition
  w <- u0(x)
  
  gamma <- K*tau/(h^2)
  
  Temperature <- matrix(0, N-1, M+1)  # Matrix to hold the solution
  Temperature[ , 1] <- w  # Initial value
  # Loop over time steps
  for (j in 1:M) {
    w <- doublesweep(rep(gamma, N-1), rep(gamma, N-1), 
                     rep(1 + 2* gamma, N-1), -w, 0, 0)
    Temperature[ , j+1] <- w
  }
  
  # Return a list consisting of time grid, x grid and solution
  return(list(x=x, t=t, w=Temperature))
}

#Clearly, this will be an extension of the mini-project.
backwardDifferenceIMPROVED_EXCLAMATION_POINT <- function(
                              #Mimic changes to function signature
                              f=function(x, t) 0,
                              u0=function(x) 2*sin(2*pi*x),
                              a = function(t) 0, b = function(t) 0,
                              K=1, L=1, N=30, T=0.1, M=30) {
  
  # set up space grid
  h <- L/N
  #Mimic the changes to the grid
  x <- h*(0:N)
  
  # set up time grid
  tau <- T/M
  t <- tau*(0:M)
  
  # set up vectors with initial condition
  w <- u0(x)
  
  gamma <- K*tau/(h^2)
  
  #Adjust matrix size as well
  Temperature <- matrix(0, N+1, M+1)  # Matrix to hold the solution
  Temperature[ , 1] <- w  # Initial value
  # Loop over time steps
  for (j in 1:M) {
    #Compute F as necessary
    #Q: As above, Cautious about t[j] vs t[j+1]
    #A: We are now using Equation 2.62. In this case
    #We should use the next time, rather than the previous.
    F <- f(x,t[j+1])
    #Temp names!
    a_ <- a(t[j+1])
    b_ <- b(t[j+1])
    w <- doublesweep(rep(gamma, N+1), rep(gamma, N+1), 
                     rep(1 + 2* gamma, N+1), -(w+tau*F), a_, b_)
    w[1] <- a_
    w[N+1] <- b_
    Temperature[ , j+1] <- w
  }
  
  # Return a list consisting of time grid, x grid and solution
  return(list(x=x, t=t, w=Temperature))
}

bd2 <- backwardDifferenceIMPROVED_EXCLAMATION_POINT
#Now we check to see if it matches for smaller values!

sol <- bd2(u0=function(x) x, 
           f=function(x,t) 10*sin(10*pi*t)*exp(-10*(x-0.5)^2),
           a=function(t) sin(20*pi*t), b=function(t) cos(20*pi*t),
           K=1, L=1, T=0.4, N=30, M=30)
persp3D(sol$x, sol$t, sol$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4)
plotrgl(smooth = TRUE)

