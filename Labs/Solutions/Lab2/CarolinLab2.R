#NumMethodsPDE
#Carolin Natemeyer
#Lab2

#Exercise1
mysum<-function(N){
  f=0
   for(n in 1:N){
    s = 1/n^2
  f=f+s
  }
   return(f) 
}
mysum(10)
mysum(50)
mysum(100)

#Exercise2

#function to greate circle
circle <- function(N, r) {
  t <- (1:N)/N
  A <- list(x=r*sin(2*pi*t), y=r*cos(2*pi*t))
  return(A)
}
par(pty="s")
coords <- circle(100, 7)
plot(coords$x, coords$y,  xlab="X", ylab="Y", col="green") #plots the most outstanding circle 
 for(n in 1:6){ 
   #loop plots the remaining inner circles in the same plot
coords <- circle(100, n)
lines(coords$x, coords$y,  xlab="X", ylab="Y")
}
par(pty="m")

#Exercise3
mysin<-function(x){
  #computes approximation of sin(x) with input argument x by tayler series
  #camputes summands of taylerseries in loop until the desired errow is reached
  sin = 0
  n=0
  r=1 #initial residue

while (r>=10^-4){
  s= (-1)^n *x^(2*n+1)/factorial(2*n+1)
  sin = sin+s
  r= abs(x)^(2*n+3)/factorial(2*n+3)
  n=n+1
}
    return(sin)
}

#Exercise 4
#Plot Rosenbrockfunction
x<-seq(-2,2,length.out=200)
y<-seq(-2,2,length.out=200)
xy<-mesh(x,y)
z<-with(xy,(1-x)^2 +100*(y-x^2)^2)
persp3D(x,y,z)
plotrgl