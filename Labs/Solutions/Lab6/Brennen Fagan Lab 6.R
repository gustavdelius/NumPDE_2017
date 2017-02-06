#Lab 6:

require("plot3Drgl")

#From Lab 5
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
  
  # Set up evolution constant
  gamma <- K*tau/(h^2)
  #We seek to solve
  #w(j-1) = A*w(j)
  #For w(j)
  
  #Now we want to use Double-sweep
  #With v = w(j), A = tridiagonal
  #and F = w(j-1).
  
  #Tridiagonal matrix doesn't change.
  Ci <- rep(1+2*gamma, N-1)
  Ai <- rep(gamma, N-1)
  Bi <- rep(gamma, N-1)
  
  Temperature <- matrix(0, N-1, M+1)  # Matrix to hold the solution
  Temperature[ , 1] <- w  # Initial value
  # Loop over time steps
  for (j in 1:M) {
    #F <- w  
    #No boundary conditions inputted
    #Will need to fix in order to add
    #boundary conditions
    
    #Note, a hilarious plot forms if you
    #forget the - in front of the w.
    #Do try and rotate the plot around.
    w <- doublesweep(Ai,Bi,Ci, -w, 0, 0)
    Temperature[ , j+1] <- w
  }
  
  # Return a list consisting of time grid, x grid and solution
  return(list(x=x, t=t, w=Temperature))
}

#Exercise 1: Copying function signature
CrankNicolson <- function(u0=function(x) 2*sin(2*pi*x), 
                          K=1, L=1, N=30, T=0.1, M=30) {
  # set up space grid
  h <- L/N
  x <- h*(1:(N-1))
  
  # set up time grid
  tau <- T/M
  t <- tau*(0:M)
  
  # set up vectors with initial condition
  w <- u0(x)
  
  # Set up evolution constant
  gamma <- K*tau/(h^2)
  #We seek to solve
  #w(j-1) = A*w(j)
  #For w(j)
  
  #Now we want to use Double-sweep
  #With v = w(j), A = tridiagonal
  #and F = w(j-1).
  
  #Tridiagonal matrix doesn't change.
  Ci <- rep(1+gamma, N-1)
  Ai <- rep(gamma/2, N-1)
  Bi <- rep(gamma/2, N-1)
  
  #Calculate transition matrix B
  B <- diag(1-gamma, N-1)
  for (k in 1:(N-2)) {
    B[k,k+1] <- gamma/2
    B[k+1,k] <- gamma/2
  }
  
  Temperature <- matrix(0, N-1, M+1)  # Matrix to hold the solution
  Temperature[ , 1] <- w  # Initial value
  # Loop over time steps
  for (j in 1:M) {
    #F <- w  
    #No boundary conditions inputted
    #Will need to fix in order to add
    #boundary conditions
    w <- doublesweep(Ai,Bi,Ci, -(B%*%w), 0, 0)
    Temperature[ , j+1] <- w
  }
  
  # Return a list consisting of time grid, x grid and solution
  return(list(x=x, t=t, w=Temperature))
}

solbd <- backwardDifference(N=300, M=300)
solCN <- CrankNicolson(N=300, M=300)
max(solbd$w - solCN$w)