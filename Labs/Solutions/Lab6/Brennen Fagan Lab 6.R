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
cat("Difference in methods: ", max(solbd$w - solCN$w))

maxError <- function(N, M, method, u0 = function(x) 2*sin(2*pi*x), omega = 2) {
  # numerical solution
  numSol <- method(M=M, N=N, u0 = u0)
  # exact solution
  x <- numSol$x
  t <- numSol$t
  xy <- mesh(x, t)
  u <- with(xy, u0(x)*exp(-(omega^2)*pi^2*y))
  
  return(max(abs(u - numSol$w)))
}

plotError <- function(N, M, method, u0 = function(x) 2*sin(2*pi*x), omega = 2) {
  #WARNING: Assumes Sine
  # numerical solution
  numSol <- method(M=M, N=N, u0 = u0)
  # exact solution
  x <- numSol$x
  t <- numSol$t
  xy <- mesh(x, t)
  u <- with(xy, u0(x)*exp(-omega^2*pi^2*y))
  
  persp3D(x,t, numSol$w,
          xlab="Length", ylab="Time", zlab="Approx. Temp.", # Provides axis labels
          ticktype="detailed", nticks=4) # Provides axis ticks
  
  persp3D(x,t, u,
          xlab="Length", ylab="Time", zlab="Exact Temp.", # Provides axis labels
          ticktype="detailed", nticks=4) # Provides axis ticks
  
  persp3D(x, t, (u - numSol$w),
          xlab="Length", ylab="Time", zlab="Errors", # Provides axis labels
          ticktype="detailed", nticks=4) # Provides axis ticks
  plotrgl()
}

N <- 15*2^(0:7) -> M
#So this is the elegant way of doing what I did in lab 5, exercise 4
#although restricted to one variable apparently.
errbd <- sapply(N, function(N) maxError(N, N, backwardDifference))
errCN <- sapply(N, function(N) maxError(N, N, CrankNicolson))
plot(errbd ~ N, type="b", log="xy", ylim=c(0.0000001, 0.1), ylab="Error")
lines(errCN ~ N, type="b", col="blue")

#Retrieve the coefficients from the plot, since I don't trust my eyes.
bd.lm <- lm(log(errbd) ~ log(N))
CN.lm <- lm(log(errCN) ~ log(N))

#Exercises 2 and 3
#We resort to the methods of Lab 5. sapply and mapply do not do quite what I am
#wanting to accomplish.

N <- 15*2^(0:8) -> M

ErrorMatbd <- matrix(0, length(N), length(M))
ErrorMatCN <- matrix(0, length(N), length(M))
for (n in 1:length(N)){
  for (m in 1:length(M)){
    ErrorMatbd[n,m] <- maxError(N[n],M[m], backwardDifference)
    ErrorMatCN[n,m] <- maxError(N[n],M[m], CrankNicolson)
  }
}

persp3D(N, M, ErrorMatbd-ErrorMatCN,
        xlab="Length Divisions", ylab="Time Divisions", zlab="Difference in Maximum Errors", # Provides axis labels
        ticktype="detailed", nticks=4) # Provides axis ticks
plotrgl()

#From which we can pull N=60 or M=60
plot(ErrorMatbd[,3] ~ N, type="b", log="xy", ylim=c(0.0000001, 0.1), ylab="Error")
lines(ErrorMatCN[,3] ~ N, type="b", col="blue")
title("Error as a function of length divisions, M=60")

plot(ErrorMatbd[3,] ~ M, type="b", log="xy", ylim=c(0.0000001, 0.1), ylab="Error")
lines(ErrorMatCN[3,] ~ M, type="b", col="blue")
title("Error as a function of time divisions, N=60")

#Exercise 4
N <- 15*2^(0:8) -> M

ErrorMatbd <- matrix(0, length(N), length(M))
ErrorMatCN <- matrix(0, length(N), length(M))
for (n in 1:length(N)){
  for (m in 1:length(M)){
    ErrorMatbd[n,m] <- maxError(N[n],M[m], backwardDifference, u0 = function(x) 2*sin(8*pi*x), omega = 8)
    ErrorMatCN[n,m] <- maxError(N[n],M[m], CrankNicolson, u0 = function(x) 2*sin(8*pi*x), omega = 8)
  }
}

persp3D(N, M, ErrorMatbd-ErrorMatCN,
        xlab="Length Divisions", ylab="Time Divisions", zlab="Difference in Maximum Errors", # Provides axis labels
        ticktype="detailed", nticks=4) # Provides axis ticks

#From which we can pull N=60 or M=60
plot(ErrorMatbd[,3] ~ N, type="b", log="xy", ylim=c(0.001, 1), ylab="Error")
lines(ErrorMatCN[,3] ~ N, type="b", col="blue")
title("Error as a function of length divisions, M=60")

plot(ErrorMatbd[3,] ~ M, type="b", log="xy", ylim=c(0.001, 1), ylab="Error")
lines(ErrorMatCN[3,] ~ M, type="b", col="blue")
title("Error as a function of time divisions, N=60")


persp3D(N, M, ErrorMatCN,
        xlab="Length Divisions", ylab="Time Divisions", zlab="Crank-Nicholson Maximum Errors", # Provides axis labels
        ticktype="detailed", nticks=4) # Provides axis ticks
plotrgl()

persp3D(N, M, ErrorMatbd,
        xlab="Length Divisions", ylab="Time Divisions", zlab="Backward Difference Maximum Errors", # Provides axis labels
        ticktype="detailed", nticks=4) # Provides axis ticks
plotrgl()

#Some conclusions and comments from Exercises 2, 3, and 4:
#   There are distinct diminishing returns in increasing the number of divisions.
#   This holds true in all exercises.
#   This partially is due to error's dependence on both length and time.
#   We can observe this from the plots of the errors in 3-D.
#     It appears that there is a minimum threshold in both length and time divisions
#     that must be passed before we reach the diminishing returns, but this threshold
#     is very small, and much smaller for Crank-Nicholson than for Backward-Difference.
#     Interestingly, Crank-Nicholson has a higher overall error initially, but a much
#     steeper drop off!
#   I am still confused as to why there is a local minimum in terms of the errors when time/space
#     has predefined divisions (M = 60 or N = 60). Perhaps it is due to terms neglected in big-O?
#       (Discussed in class)


plotError(60,60,CrankNicolson)
plotError(60,60,CrankNicolson, u0 = function(x) 2*sin(8*pi*x), omega = 8)
maxError(60,60, CrankNicolson,  u0 = function(x) 2*sin(8*pi*x), omega = 8)

