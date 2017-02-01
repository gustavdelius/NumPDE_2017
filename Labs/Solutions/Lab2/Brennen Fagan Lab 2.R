#Lab 2:
#Recreating function from previous Lab
function_circle <- function(NumberOfPoints, Radius, Colors = "default", NewPlot = 1){
  t <- (1:NumberOfPoints)/NumberOfPoints
  x <- sin(2*pi*t)
  y <- cos(2*pi*t)
  A <- list(x, y)
  par(pty="s")
  coloring <- rgb(0,0,0)
  if(Colors == "RG_Pretty"  || Colors == 1){
    coloring <- rgb(x*.5+.5,y*.5+.5,0)
  }else if(Colors == "GB_Pretty" || Colors == 2){
    coloring <- rgb(0,x*.5+.5,y*.5+.5)
  }else if(Colors == "RB_Pretty" || Colors == 3){
    coloring <- rgb(y*.5+.5,0,x*.5+.5)
  }
  A <- lapply(A,"*",Radius)
  if(NewPlot){
    plot(A[[1]], A[[2]], main="A circle", xlab="X", ylab="Y", col = coloring)
  }else{
    points(A[[1]], A[[2]], main="A circle", xlab="X", ylab="Y", col = coloring)
  }
  par(pty="m")
  return (A)
}

function_circle(100,2, "GB_Pretty")

coords <- function_circle(7,3)
par(pty="s")
plot(coords[[1]], coords[[2]], main = "A septagon", xlab = "X", ylab="Y", type = "l")
lines(c(coords[[1]][7], coords[[1]][1]), c(coords[[2]][7], coords[[2]][1]))

function_ngon <- function(n, Radius){
  coords <- function_circle(n, Radius)
  par(pty="s")
  plot(coords[[1]], coords[[2]], main = toString(paste("A", n, "-gon")), xlab = "X", ylab="Y", type = "l")
  lines(c(tail(coords[[1]],1), coords[[1]][1]), c(tail(coords[[2]],1), coords[[2]][1]))
  par(pty="m")
  return (coords)
}

function_ngon(50, 5)

function_calc_shifted_circle <- function(NumberOfPoints, Radius#, Colors = "default"
                                         ){
  t <- (1:NumberOfPoints)/NumberOfPoints
  x <- sin(2*pi*t)
  y <- cos(2*pi*t)
  A <- list(x, y)
  #par(pty="s")
  # coloring <- rgb(0,0,0)
  # if(Colors == "RG_Pretty"){
  #   coloring <- rgb(x*.5+.5,y*.5+.5,0)
  # }else if(Colors == "GB_Pretty"){
  #   coloring <- rgb(0,x*.5+.5,y*.5+.5)
  # }else if(Colors == "RB_Pretty"){
  #   coloring <- rgb(x*.5+.5,0,y*.5+.5,0)
  # }
  A <- lapply(A,"*",Radius)
  #plot(A[[1]] + SHIFT, A[[2]], main="A circle", xlab="X", ylab="Y", col = coloring)
  #par(pty="m")
  return (A <- list(x = A[[1]]+ SHIFT, y = A[[2]]))
}

SHIFT <- 1
par(pty = "s")
plot(function_calc_shifted_circle(100, 2)$x,function_calc_shifted_circle(100, 2)$y)
SHIFT <- 2
points(function_calc_shifted_circle(100, 1)$x,function_calc_shifted_circle(100, 1)$y)
par(pty="m")

function_euler_number <- function(xValue, NthExpansion){
  summand <- 1
  Taylor <- 1
  for (j in 1:NthExpansion){
    Taylor <- Taylor*xValue/j
    summand <- summand + Taylor
  }
  return (summand)
}

print(function_euler_number(1,5))

#Exercise 1:
function_p_series <- function(NthTerm, Exponent){
  #Computes 1/(1^Exponent)+...+1/(NthTerm^Exponent)
  #Request is for a loop
  summand <- 0
  for (nthTerm in 1:NthTerm){
    summand <- summand + nthTerm ^ (-Exponent)
  }
  return (summand)
}
mysum <- function(N){
  return (function_p_series(N,2))
}
print(mysum(10))
print(mysum(50)) 
print(mysum(100))

#Exercise 2:
for(j in 7:1){
  if(j==7){
    function_circle(100, j, j%%4, 1)
  }else{
    function_circle(100, j, j%%4, 0)
  }
}

#Should be able to implement a function that implements
#Taylor Series functions for general "easy" functions
#such as sin, cos, etc.
#Currently not sure how to evaluate derivative at a point.
# function_Taylor_Series <- function(ApproxValue, aPoint, EndTerm){
#   summand <- 0
#   for(nthTerm in 1:EndTerm){
#     function(f) 
#     gradient = deriv(~ f(x), "x")
#     summand <- summand + deriv(f(aPoint))*(ApproxValue - aPoint)
#   }
# }

#Exercise 3: Boring way
function_MacLaurin_Sin <- function(ApproxValue, OrderOfError){
  summand <- ApproxValue
  nthTerm <- 0
  error <- abs(ApproxValue) ^ (2*nthTerm + 3)/(factorial(2*nthTerm+3))
  while(error>10^(-OrderOfError)){
    nthTerm <- nthTerm + 1
    summand <- summand + (-1)^nthTerm * ApproxValue ^ (2*nthTerm+1)/factorial(2*nthTerm+1)
    error <- abs(ApproxValue) ^ (2*nthTerm + 3)/(factorial(2*nthTerm+3))
  }
  return(summand)
}

mysin <- function(x){
  return (function_MacLaurin_Sin(x, 4))
}

mysin(0)

require("plot3Drgl")
x <- seq(-2.5, 2.5, length.out = 100)
y <- seq(-2,2, length.out = 100)
xy <- mesh(x, y)
z <- with(xy, x*exp(-x^2 - y^2))
persp3D(x,y, z)
plotrgl(smooth = TRUE)

#Exercise 4
require("gmp")
TUPPER <- as.bigz("96093937991895888497167296212785275471500433966012930665150551927170280239526642468964284174350718121267153782770623355993237280874144307891325963941337723487857735749823926629715517173716995165232890538221612403238855866184013235585136048828693337902491454229288667081096184496091705183454067827731551705405381627380967602565625016981482083418783163849115590225610003652351370343874461848378737238198224849863465033159410054974700593138339226497249461751545728366702369745461014655997933798537483143786841806593422227898388722980000748404719")
x <- seq(0, 105)
function_TUPPER <- function(){
  Output <- matrix(0, nrow=106, ncol=17)
  x <- seq(0, 105)
  for (j in 0:16){
    y <- TUPPER + j
    for(chi in x){
      #Appears that it will not give a decimal representation in the gmp package
      #We are stuck. As a work around, we see how close performing strict modulo is.
      Output[chi+1,j+1] <-as.integer(as.bigz(
              ((y/17-((y%%17)/17))*as.bigq(2^(-17*chi-(floor(y)%%17))))
            )%%2)
    }
  }
  return(Output)
}
persp3D(x, seq(0:16), function_TUPPER())
plotrgl()
#Clearly, the strict modulo does not work in reproducing this function!