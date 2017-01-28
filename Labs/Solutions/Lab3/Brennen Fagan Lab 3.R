#Lab 3:
#Exercise 1:
#ROD
N <- 30            #Number of Divisions Of Rod
L <- 1             #Length of Rod
h <- L/N           #Segments Lengths of Rod
x <- h*(1:(N-1))   #Segments of Rod

function_default_inicond <- function(x){
  return (2*sin(2*pi*x))
}
w <- function_default_inicond(x) #Initial Heat Distribution at segments

#TIME
T <- .1            #Length of Time
M <- 200           #Number of Segments of Time
tau <- T/M         #Length of Time Segments
t <- tau*(0:M)     #Time Segments

#HEAT DISPERSION PARAMETER
K <- 1

#Resultant EVALUATION PARAMETER
gamma <- K*tau/(h^2)

#Time Step Matrix
A <- diag(1-2*gamma, N-1)
for (k in 1:(N-2)){
  A[k, k+1] <- gamma -> A[k+1,k]
}

#Initial Sampling
plot(x,w, type = "l")
lines(x, A%*%w, col = "blue")

#Full Sampling (would require reseting w to rerun)
Temperature <- matrix(0, N-1, M+1) # Matrix to hold the solution
Temperature[ , 1] <- w
# Loop over time steps
for (j in 0:(M)) {
  w <- A %*% w
  Temperature[ , j+1] <- w
}

require("plot3Drgl")
persp3D(x, t, Temperature,
        xlab="x", ylab="t", zlab="Temperature", # Provides axis labels
        ticktype="detailed", nticks=4) # Provides axis ticks
plotrgl(smooth = TRUE, lighting = TRUE)


#Prep for Other Exercises:
#Collecting commands into a single function
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

#Create Different Initial Conditions:
#This alters the boundary condition somewhat
function_inicond_identity <- function(x){
  return(x)
}
#This places the heat all in the middle:
function_inicond_peaked <- function(x){
  MidPoint = (length(x)+1)/2
  Output = rep(0, length(x))
  for(j in 1:MidPoint){
    Output[j] = (j-1)
    Output[length(x)-j+1] = (j-1)
  }
  return(Output)
}

#Exercise 2:
function_forward_difference_heat(L, N, T, M, K, 
                                 function_inicond_identity,
                                 DO_Plot = 1)
function_forward_difference_heat(L, N, T, M, K, 
                                 function_inicond_peaked,
                                 DO_Plot = 1)

#Exercise 3,4,5:
#We expect that for higher values of K, as
#the heat dispersion parameter, would result
#in more heat dispersed, i.e. the heat spreads
#more quickly through the rod.
lo_disp_temp <- function_forward_difference_heat(L, N, T, M, .5, 
                                 function_default_inicond,
                                 DO_Plot = 1)
#Already receiving odd errors. Need to adjust parameters.
md_disp_temp <- function_forward_difference_heat(L, N, T, M*10, 1.5, 
                                 function_default_inicond,
                                 DO_Plot = 1)
hi_disp_temp <- function_forward_difference_heat(L, N, T, M*100, 5, 
                                 function_default_inicond,
                                 DO_Plot = 1)

#Exercise 6
function_sample_nonhomogeneity <- function(x){
  return(-25*sin(3*pi*x))
}
function_inicond_ex6 <- function(x){
  return(-1.5*sin(2*pi*x))
}

function_forward_difference_heat_nonhomogenous<- function(
  Length, LengthDivisions, Time, 
  TimeDivisions, HeatDispersion, 
  InitialConditionFunction = function_inicond_ex6, 
  NonHomogenousFunction = function_sample_nonhomogeneity,
  DO_Plot = 0
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
  
  #Calculate Forcing Vector
  Forcing <- NonHomogenousFunction(LengthSegments)
  
  #Calculate Initial Condition
  w <- InitialConditionFunction(LengthSegments)
  
  #Calculate propagation
  Temperature <- matrix(-15, LengthDivisions-1, TimeDivisions+1) # Matrix to hold the solution
  Temperature[ , 1] <- w
  # Loop over time steps
  for (j in 1:(TimeDivisions)) {
    w <- Transition %*% w + TimeSegmentSize * Forcing
    Temperature[ , j+1] <- w
  }
  
  if(DO_Plot){
    require("plot3Drgl")
    persp3D(LengthSegments, TimeSegments, Temperature,
            xlab="x", ylab="t", zlab="Temperature", # Provides axis labels
            ticktype="detailed", nticks=4) # Provides axis ticks
    plotrgl(smooth = TRUE, lighting = TRUE)
  }
  
  return(list(LengthSegments, TimeSegments, Temperature))
}

function_zero <-function(x){
  dimensions <- dim(x)
  return(
          matrix(
                  0,
                  nrow = dimensions[1],
                  ncol = dimensions[2]
                )
         )
}

#So it works for the "trivial" case.
BasicTest <- function_forward_difference_heat_nonhomogenous(
  L, N, T, M, K, 
  NonHomogenousFunction = function_zero,
  DO_Plot = 0
)

#Identified N and M: tau <= h^2/2K
BasicTest2 <- function_forward_difference_heat_nonhomogenous(
  1, 30, .2, 1000, 1,
  NonHomogenousFunction = function_zero,
  DO_Plot = 0
)

Ex6Soln <- function_forward_difference_heat_nonhomogenous(
  1, 30, .2, 1000, 1, DO_Plot = 1
)