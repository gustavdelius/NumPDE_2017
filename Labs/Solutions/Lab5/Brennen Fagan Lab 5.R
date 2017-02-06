#Lab 5:

#Build Environment
N <- 30  # Number of time points after initial state
L <- 1
h <- L/N
x <- h*(1:(N-1))  # Vector of non-boundary grid points
w <- 2*sin(2*pi*x)  # Initial condition
M <- 30  # Number of time points after initial state
T <- 0.1  # Final time
tau <- T/M  # Time step size
t <- tau*(0:M)  # Vector of time steps
K <- 1  # Diffusion rate

#Build (Backward) Transform Matrix
gamma <- K*tau/(h^2)
A <- diag(1 + 2*gamma, N-1)
for (k in 1:(N-2)){
  A[k,k+1] <- -gamma
  A[k+1,k] <- -gamma
}

#Build (Forward) Transform Matrix
Ainv <- solve(A)

#Apply Transform
Temperature <- matrix(0, N-1, M+1) # Matrix to hold the solution
Temperature[ , 1] <- w
# Loop over time steps
for (j in 0:M) {
  w <- Ainv %*% w
  Temperature[ , j+1] <- w
}

#Plot
require("plot3Drgl")
persp3D(x, t, Temperature,
        xlab="x", ylab="t", zlab="Temperature", # Provides axis labels
        ticktype="detailed", nticks=4) # Provides axis ticks

#Functional version
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
  
  # Set up evolution matrix
  gamma <- K*tau/(h^2)
  A <- diag(1+2*gamma, N-1)
  for (k in 1:(N-2)) {
    A[k,k+1] <- -gamma
    A[k+1,k] <- -gamma
  }
  Ainv <- solve(A)
  
  Temperature <- matrix(0, N-1, M+1)  # Matrix to hold the solution
  Temperature[ , 1] <- w  # Initial value
  # Loop over time steps
  for (j in 1:M) {
    w <- Ainv %*% w 
    Temperature[ , j+1] <- w
  }
  
  # Return a list consisting of time grid, x grid and solution
  return(list(x=x, t=t, w=Temperature))
}

#Compare with Forward Difference from Lab 3.
function_forward_difference_heat <- function(
  Length, LengthDivisions, Time, TimeDivisions,
  HeatDispersion, InitialConditionFunction, DO_Plot = 0
){
  #Calculate h and tau
  LengthSegmentSize <- Length/LengthDivisions
  TimeSegmentSize <- Time/TimeDivisions
  #Calculate x and t
  LengthSegments <- LengthSegmentSize*(1:(LengthDivisions-1))
  TimeSegments <- TimeSegmentSize*(0:TimeDivisions)
  #Calculate gamma
  gamma <- HeatDispersion*TimeSegmentSize/(LengthSegmentSize^2)
  
  #Calculate Transition Matrix
  Transition <- diag(1-2*gamma, LengthDivisions-1)
  for (k in 1:(LengthDivisions-2)){
    Transition[k, k+1] <- gamma -> Transition[k+1,k]
  }
  
  #Calculate Initial Condition
  w <- InitialConditionFunction(LengthSegments)
  
  #Calculate propagation
  Temperature <- matrix(0, LengthDivisions-1, TimeDivisions+1) # Matrix to hold the solution
  Temperature[ , 1] <- w
  # Loop over time steps
  for (j in 0:(TimeDivisions)) {
    w <- Transition %*% w
    Temperature[ , j+1] <- w
  }
  
  if(DO_Plot){
    require("plot3Drgl")
    persp3D(LengthSegments, TimeSegments, Temperature,
            xlab="x", ylab="t", zlab="Temperature", # Provides axis labels
            ticktype="detailed", nticks=4) # Provides axis ticks
    plotrgl(smooth = TRUE, lighting = TRUE)
  }
  
  return(Temperature)
}

#Exercise 1:
#Comparing backwardDifference(M = 30) with backwardDifference(M = 500)
sol <- backwardDifference(M=500)

require("plot3Drgl")
persp3D(sol$x, sol$t, sol$w,
        xlab="x", ylab="t", zlab="Temperature", # Provides axis labels
        ticktype="detailed", nticks=4) # Provides axis ticks

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

#Exercise 2:
ptm <- proc.time()
Mat = matrix(c(-1,2,0,0,2,-1,1,0,0,3,-1,4,0,0,1,-7),nrow=4, ncol  = 4, byrow = TRUE)
Vec = c(1,2,6,-6)
ans = solve(Mat, Vec)
print(proc.time()-ptm)

ptm <- proc.time()
#alpha1 = 0, beta1 = a
#v[n] = b
#A is Diag-1, B is Diag + 1, C is Diag, F is vector
#User discretion on A1, Bn-1
A <-c(1,2,3,1)
B <-c(2,1,4,1)
C <-c(1,1,1,7)
F <-c(1,2,6,-6)
ans2 <- doublesweep(A,B,C,F,0,0)
print(proc.time()-ptm)

backwardDifference2 <- function(u0=function(x) 2*sin(2*pi*x), 
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
  
  #In the basic method, we find
  #Ainv and repeatedly apply Ainv
  #to w(j-1) to get w(j)'s.
  
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
ans1 <- backwardDifference()
ans2 <- backwardDifference2()

print(all(ans1$x-ans2$x==0))
print(all(ans1$t-ans2$t==0))

require("plot3Drgl")
#Absolute Difference:
persp3D(ans1$x, ans1$t, ans1$w-ans2$w,
        xlab="x", ylab="t", zlab="Temperature", # Provides axis labels
        ticktype="detailed", nticks=4) # Provides axis ticks
title("Absolute Difference between methods")
#8*10^-16 is small for most intents and purposes.
plotrgl()

persp3D(ans1$x, ans1$t, ans1$w,
        xlab="x", ylab="t", zlab="Temperature", # Provides axis labels
        ticktype="detailed", nticks=4) # Provides axis ticks
plotrgl()

persp3D(ans2$x, ans2$t, ans2$w,
        xlab="x", ylab="t", zlab="Temperature", # Provides axis labels
        ticktype="detailed", nticks=4) # Provides axis ticks
plotrgl()


maxError <- function(N, M) {
  # numerical solution
  numSol <- backwardDifference2(M=M, N=N)
  # exact solution
  x <- numSol$x
  t <- numSol$t
  xy <- mesh(x, t)
  u <- with(xy, 2*sin(2*pi*x)*exp(-4*pi^2*y))
  
  return(max(abs(u - numSol$w)))
}

#N, M increasing should result in decreasing errors.
#I think I still expect space to predominate over time
#if only due to the factor of h^2.
#Due to the lack of interdependence though, I imagine
#that the effects should be more similar than what I
#had last time I tried to explore the errors.

Nvals <- seq(100,1000, by = 100) -> Mvals
ErrorMat <- matrix(0, length(Nvals), length(Mvals))
for (n in 1:length(Nvals)){
  for (m in 1:length(Mvals)){
    ErrorMat[n,m] <- maxError(Nvals[n], Mvals[m])
  }
}
persp3D(Nvals, Mvals, ErrorMat,
        xlab="Length Divisions", ylab="Time Divisions", zlab="Maximum Error", # Provides axis labels
        ticktype="detailed", nticks=4) # Provides axis ticks
plotrgl()

#Well, that was unexpected.
#Space appears to suffer quite much from diminishing
#returns, and now appears to suffer from the same problem
#that time suffered under the forward difference method.
#Meanwhile, the time dimension is now more critical, and
#plays the role that the space dimension did in 
#the forward difference method's error!

#Looking back, I clearly got confused over the O(tau+h^2)
#and thought backward difference was O(tau^2+h^2) instead.
#This makes the dependence much more natural and highlights
#the aberration that is the forward difference method.

plot(Nvals, ErrorMat[,1])
plot(Mvals, ErrorMat[1,])
