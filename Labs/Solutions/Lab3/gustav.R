# Exercises from Lab 3

#### Exercise 1 ####
# Collect all the commands from above in your own R script file and make sure 
# you can run the code from there and get the same solution.

# space grid
N <- 30
L <- 1
h <- L/N
x <- h*(1:(N-1))

# initial condition
w <- 2*sin(2*pi*x)

# time grid
T <- 0.1
M <- 200
tau <- T/M
t <- tau*(0:M) 

# parameters
K <- 1
gamma <- K*tau/(h^2)

# evolution matrix
A <- diag(1-2*gamma, N-1)
for (k in 1:(N-2)) {
    A[k,k+1] <- gamma
    A[k+1,k] <- gamma
}

# Calculate solution
Temperature <- matrix(0, N-1, M+1) # Matrix to hold the solution
Temperature[ , 1] <- w
# Loop over time steps
for (j in 0:M) {
    w <- A %*% w
    Temperature[ , j+1] <- w
}

# Plot solution
library("plot3Drgl") 
persp3D(x, t, Temperature,
        xlab="x", ylab="t", zlab="Temperature", # Provides axis labels
        ticktype="detailed", nticks=4) # Provides axis ticks
plotrgl(smooth=TRUE, lighting = TRUE)


#### Exercise 2 ####
# Evolve the heat equation from two other initial conditions of your choosing.

## First choice
w <- sin(8*pi*x)

# Calculate solution
Temperature <- matrix(0, N-1, M+1) # Matrix to hold the solution
Temperature[ , 1] <- w
# Loop over time steps
for (j in 0:M) {
    w <- A %*% w
    Temperature[ , j+1] <- w
}

# Plot solution
persp3D(x, t, Temperature,
        xlab="x", ylab="t", zlab="Temperature", # Provides axis labels
        ticktype="detailed", nticks=4) # Provides axis ticks


## Second choice
w <- 0.5-abs(x-0.5)

# Calculate solution
Temperature <- matrix(0, N-1, M+1) # Matrix to hold the solution
Temperature[ , 1] <- w
# Loop over time steps
for (j in 0:M) {
    w <- A %*% w
    Temperature[ , j+1] <- w
}
# Plot solution
persp3D(x, t, Temperature,
        xlab="x", ylab="t", zlab="Temperature", # Provides axis labels
        ticktype="detailed", nticks=4) # Provides axis ticks


#### Exercise 3 ####
# Now explore the effect of different values of $K$ on the solution. Start 
# with $K=0.5$. What do you observe? Does it agree with your expectation?

# initial condition
w <- 2*sin(2*pi*x)
# parameters
K <- 0.5
gamma <- K*tau/(h^2)

# evolution matrix
A <- diag(1-2*gamma, N-1)
for (k in 1:(N-2)) {
    A[k,k+1] <- gamma
    A[k+1,k] <- gamma
}
# Calculate solution
Temperature <- matrix(0, N-1, M+1) # Matrix to hold the solution
Temperature[ , 1] <- w
# Loop over time steps
for (j in 0:M) {
    w <- A %*% w
    Temperature[ , j+1] <- w
}
# Plot solution
persp3D(x, t, Temperature,
        xlab="x", ylab="t", zlab="Temperature", # Provides axis labels
        ticktype="detailed", nticks=4) # Provides axis ticks

# As expected, the time evolution happens more slowly than at K=1
# (by a factor of 2)


#### Exercise 4 ####
# Next try $K=1.5$. What do you observe? Does it agree with your expectation?

# initial condition
w <- 2*sin(2*pi*x)
# parameters
K <- 1.5
gamma <- K*tau/(h^2)

# evolution matrix
A <- diag(1-2*gamma, N-1)
for (k in 1:(N-2)) {
    A[k,k+1] <- gamma
    A[k+1,k] <- gamma
}
# Calculate solution
Temperature <- matrix(0, N-1, M+1) # Matrix to hold the solution
Temperature[ , 1] <- w
# Loop over time steps
for (j in 0:M) {
    w <- A %*% w
    Temperature[ , j+1] <- w
}
# Plot solution
persp3D(x, t, Temperature,
        xlab="x", ylab="t", zlab="Temperature", # Provides axis labels
        ticktype="detailed", nticks=4) # Provides axis ticks

# That does not at all look as it should. There are huge spikes towards
# the end of the time period. Stability problem.


#### Exercise 5 ####
# To avoid the instability, change the size of the time step. How small does
# the time step have to be if you want to get a solution for $K=5$?

K <- 5
# Using eq.(2.17) from the notes we find that
tau <- h^2/(2*K)
# is the smallest possible value for tau
tau
# Because $\tau=T/M$ this means we need 
M <- ceiling(T/tau)
# time steps
M
## Let's test this

t <- tau*(0:M) 
w <- 2*sin(2*pi*x)

gamma <- K*tau/(h^2)
A <- diag(1-2*gamma, N-1)
for (k in 1:(N-2)) {
    A[k,k+1] <- gamma
    A[k+1,k] <- gamma
}
Temperature <- matrix(0, N-1, M+1) # Matrix to hold the solution
Temperature[ , 1] <- w
for (j in 0:M) {
    w <- A %*% w
    Temperature[ , j+1] <- w
}
persp3D(x, t, Temperature,
        xlab="x", ylab="t", zlab="Temperature", # Provides axis labels
        ticktype="detailed", nticks=4) # Provides axis ticks

# The method is indeed stable at this stepsize


#### Exercise 6 ####
# Modify the above file appropriately to solve the following initial 
# boundary value problem

# initial condition
w <- -1.5*sin(2*pi*x)
F <- -25*sin(3*pi*x)

# time grid
T <- 0.2
M <- 400
tau <- T/M
t <- tau*(0:M) 

# parameters
K <- 1
gamma <- K*tau/(h^2)

# evolution matrix
A <- diag(1-2*gamma, N-1)
for (k in 1:(N-2)) {
    A[k,k+1] <- gamma
    A[k+1,k] <- gamma
}

# Calculate solution
Temperature <- matrix(0, N-1, M+1) # Matrix to hold the solution
Temperature[ , 1] <- w
# Loop over time steps
for (j in 0:M) {
    w <- A %*% w + tau*F
    Temperature[ , j+1] <- w
}

# Plot solution
library("plot3Drgl") 
persp3D(x, t, Temperature,
        xlab="x", ylab="t", zlab="Temperature", # Provides axis labels
        ticktype="detailed", nticks=4) # Provides axis ticks
plotrgl(smooth=TRUE, lighting = TRUE)