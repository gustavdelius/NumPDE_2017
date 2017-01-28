#Lab 4:
forwardDifference <- function(f=function(x) 0, 
                              u0=function(x) 2*sin(2*pi*x), 
                              K=1, L=1, N=30, T=0.1, M=200) {
  # set up space grid
  h <- L/N
  x <- h*(1:(N-1))
  
  # set up time grid
  tau <- T/M
  t <- tau*(0:M)
  
  # set up vectors with initial condition and inhomogenous term
  w <- u0(x)
  F <- f(x)
  
  # Set up evolution matrix
  gamma <- K*tau/(h^2)
  A <- diag(1-2*gamma, N-1)
  for (k in 1:(N-2)) {
    A[k,k+1] <- gamma
    A[k+1,k] <- gamma
  }
  
  Temperature <- matrix(0, N-1, M+1)  # Matrix to hold the solution
  Temperature[ , 1] <- w  # Initial value
  # Loop over time steps
  for (j in 1:M) {
    w <- A %*% w + tau * F
    Temperature[ , j+1] <- w
  }
  
  # Return a list consisting of time grid, x grid and solution
  return(list(x=x, t=t, w=Temperature))
}

sol <- forwardDifference()
require("plot3Drgl")

persp3D(sol$x, sol$t, sol$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4)
plotrgl(smooth=TRUE, lighting = TRUE)


sol <- forwardDifference(f=function(x) {-25*sin(3*pi*x)},
                         u0=function(x) {-1.5*sin(2*pi*x)},
                         T=0.2, M=600
)
persp3D(sol$x, sol$t, sol$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4)
plotrgl(smooth=TRUE, lighting = TRUE)

#Exercise 1:
#Function Signature:
#forwardDifference(
#ForcingFunction,
#BoundaryCondition,
#K=1, L=1, N=30, T=0.1, M=200
#)
#Note tau <= h^2/(2*1/4), where h = L/N, tau = T/M
#     T/M <= (L/N)^2*2
#     .1/M <= (1/N)^2
#
sol <- forwardDifference(f = function(x){-16*sin(8*pi*x)},
                         u0 = function(x) {sin(pi*x)},
                         K = 1/4, L = 1, N=30, T = .2, M=200)


#Plot the solution.
persp3D(sol$x, sol$t, sol$w,
        xlab="x", ylab="t", zlab="w",
        ticktype="detailed", nticks=4)
plotrgl(smooth=TRUE, lighting = TRUE)

numSol <- forwardDifference(f=function(x) {-25*sin(3*pi*x)},
                            u0=function(x) {-3/2*sin(2*pi*x)},
                            T=0.2, M=360
)
x <- numSol$x
t <- numSol$t
xy <- mesh(x, t)
u <- with(xy, -3/2*sin(2*pi*x)*exp(-4*pi^2*y)
          -25/(9*pi^2)*sin(3*pi*x)*(1-exp(-9*pi^2*y)))
#True Solution
persp3D(x, t, u, zlab="u", ticktype="detailed", nticks=4) 
plotrgl()
#Absolute Error
persp3D(x, t, u - numSol$w, zlab="u - w", ticktype="detailed", nticks=4)
plotrgl()
#Relative Error: Issue due to how close we are to 0. 
#For strictly positive problems, rel. error is more relevant.
persp3D(x, t, (u - numSol$w)/u, zlab="(u - w)/u", ticktype="detailed", nticks=4)
plotrgl()

#Interestingly, the amount of nonsense increases as the
#distance from stability increases.
#This is expected behaviour, but nonetheless nice to
#verify.
numSol <- forwardDifference(f=function(x) {-25*sin(3*pi*x)},
                            u0=function(x) {-3/2*sin(2*pi*x)},
                            T=0.2, M=1420, N=60
)
persp3D(numSol$x, numSol$t, numSol$w, 
        xlab="x", ylab="t", zlab="w", 
        ticktype="detailed", nticks=4) 

numSol <- forwardDifference(f=function(x) {-25*sin(3*pi*x)},
                            u0=function(x) {-3/2*sin(2*pi*x)},
                            T=0.2, M=1440, N=60
)
persp3D(numSol$x, numSol$t, numSol$w, 
        xlab="x", ylab="t", zlab="w", 
        ticktype="detailed", nticks=4) 

#Selection of points
tSelect <- seq(1, 1441, by=4)
t <- numSol$t[tSelect]
w <- numSol$w[, tSelect]
x <- numSol$x
xy <- mesh(x, t)
u <- with(xy, -3/2*sin(2*pi*x)*exp(-4*pi^2*y)
          -25/(9*pi^2)*sin(3*pi*x)*(1-exp(-9*pi^2*y)))
persp3D(x, t, u - w, zlab="u - w", ticktype="detailed", nticks=4)
plotrgl()
#Clear reduction in amount of error, as expected.

#Exercise 2:
#h = L/N, L is fixed, so we set N
#h = 1/120 =>  N = 120
#Note => for 480, we need .2/M < (1/480)^2/2
#                 or .4/((1/480)^2)<M => M>92160
#To be "fair", we do this for each run
#Timer:
timer <- proc.time()
numSol120 <- forwardDifference(f=function(x) {-25*sin(3*pi*x)},
                            u0=function(x) {-3/2*sin(2*pi*x)},
                            T=0.2, M=92200, N=120
)

#Selection of points
tSelect <- seq(1, 92201, by=100)
t <- numSol120$t[tSelect]
w <- numSol120$w[, tSelect]
x <- numSol120$x
xy <- mesh(x, t)
u <- with(xy, -3/2*sin(2*pi*x)*exp(-4*pi^2*y)
          -25/(9*pi^2)*sin(3*pi*x)*(1-exp(-9*pi^2*y)))

print("Amount of Time Spent Calculating for h = 1/120")
print(proc.time()-timer)
#Calculate Maximal Error:
print("Maximal Error for h = 1/120")
print(max(abs(u - w)))


timer <- proc.time()
numSol240 <- forwardDifference(f=function(x) {-25*sin(3*pi*x)},
                               u0=function(x) {-3/2*sin(2*pi*x)},
                               T=0.2, M=92200, N=240
)

#Selection of points
tSelect <- seq(1, 92201, by=100)
t <- numSol240$t[tSelect]
w <- numSol240$w[, tSelect]
x <- numSol240$x
xy <- mesh(x, t)
u <- with(xy, -3/2*sin(2*pi*x)*exp(-4*pi^2*y)
          -25/(9*pi^2)*sin(3*pi*x)*(1-exp(-9*pi^2*y)))

print("Amount of Time Spent Calculating for h = 1/240")
print(proc.time()-timer)
#Calculate Maximal Error:
print("Maximal Error for h = 1/240")
print(max(abs(u - w)))


timer <- proc.time()
numSol480 <- forwardDifference(f=function(x) {-25*sin(3*pi*x)},
                               u0=function(x) {-3/2*sin(2*pi*x)},
                               T=0.2, M=92200, N=480
)

#Selection of points
tSelect <- seq(1, 92201, by=100)
t <- numSol480$t[tSelect]
w <- numSol480$w[, tSelect]
x <- numSol480$x
xy <- mesh(x, t)
u <- with(xy, -3/2*sin(2*pi*x)*exp(-4*pi^2*y)
          -25/(9*pi^2)*sin(3*pi*x)*(1-exp(-9*pi^2*y)))
print("Amount of Time Spent Calculating for h = 1/480")
print(proc.time()-timer)
#Calculate Maximal Error:
print("Maximal Error for h = 1/480")
print(max(abs(u - w)))

print("As expected, the amount of error decreases, although it appears to plateau. Also as expected, the amount of time taken increases quickly.")


#Exercises 3 and 4:
#Calculate Maximal Error for range of step sizes h (=>N) and number of steps M
#and plot.
#Recall Exercise 1:
# sol <- forwardDifference(f = function(x){-16*sin(8*pi*x)},
#                          u0 = function(x) {sin(pi*x)},
#                          K = 1/4, L = 1, N=30, T = .2, M=200)
# Has exact solution e^(-pi^2*t/4)*sin(pi*x)-sin(8*pi*x)/pi^2+e^(-16*pi^2*t)*sin(8*pi*x)/pi^2
#Say we stop at h = 1/100
#Then we need M>.2*2/((1/100)^2)

#First, exercise 3
m = 4000
MaxErrorsH = matrix(0,nrow = 90, ncol = 1)
for(n in seq(10,100)){
  sol <- forwardDifference(f = function(x){-16*sin(8*pi*x)},
                           u0 = function(x) {sin(pi*x)},
                           K = 1/4, L = 1, N=n, T = .2, M=m)
  tSelect <- seq(1, m+1, by=4)
  t <- sol$t[tSelect]
  w <- sol$w[, tSelect]
  x <- sol$x
  xy <- mesh(x, t) #=> y=t
  u <- with(xy, exp(-pi^2*y/4)*sin(pi*x)+sin(8*pi*x)/(pi^2)*(exp(-16*pi^2*t)-1))
  MaxErrorsH[n-9] = max(abs(u - w))
}
plot(seq(10,100), MaxErrorsH, xlab="N", ylab="Maximum Error")
title("Error as a function of the number of length divisions")

#Exercise 4: Do the same, but for M
m0 = 4000
MaxErrorsM = matrix(0,nrow = 91, ncol = 1)
for(m in seq(0,90)){
  sol <- forwardDifference(f = function(x){-16*sin(8*pi*x)},
                           u0 = function(x) {sin(pi*x)},
                           K = 1/4, L = 1, N=100, T = .2, M=m0+m*400)
  tSelect <- seq(1, m0+m*400+1, by=4)
  t <- sol$t[tSelect]
  w <- sol$w[, tSelect]
  x <- sol$x
  xy <- mesh(x, t) #=> y=t
  u <- with(xy, exp(-pi^2*y/4)*sin(pi*x)+sin(8*pi*x)/(pi^2)*(exp(-16*pi^2*t)-1))
  MaxErrorsM[m+1] = max(abs(u - w))
}
plot(seq(4000,40000, by=400), MaxErrorsM, xlab="M", ylab="Maximum Error")
title("Error as a function of the number of time divisions")