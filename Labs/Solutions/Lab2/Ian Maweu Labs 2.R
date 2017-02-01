# myList <- list (animal = c("dog", "man"), legs = c(4,2)) creates a list
#myList$legs separates out a component of the list so does myList[[2]]

circle <- function(N,r)
{
  t <- (1:N)/N     #points to be evaluated at
  A <- list(x = r*sin(2*pi*t), y = r*cos(2*pi*t))     # x and y coords
  return(A)
}

shiftedCircle <- function(N,r)
{
  t <- (1:N)/N
  A <- list (x = r*sin(2*pi*t) + w, y = r*cos(2*pi*t))
  return (A)
}


#Exercise 2.1
mysum <- function(N)
{
  c = 0.0
  for(i in 1:N)
  {
    c <- c+1/(i^2)
  }
  return(c)
}
A <- c(mysum(10),mysum(50), mysum(100))

#Exercise 2.2
plot(circle(100,7)$x,circle(100,7)$y)
for(i in 6:1)
{
  points(circle(100,i)$x,circle(100,i)$y)
}

#Exercise 2.3
fac <- function(N) #Computes the factorial
{
  if(N==0) {return (1)}
  ans = 1
  for(i in 1:N)
  {
    ans = ans*i
  }
  return (ans)
}

r_n <- function(x,n)  #Computes absolute value of R_n
{
  r <- abs(x)^(2*n+3)/fac(2*n+3)
  return (r)
}

mysin <- function(x)  #Approximates sin(x) using Taylor expansion
{
  fsin = 0;
  for(i in 0:100)
  {
    n = (2*i)+1;
    fsin = fsin + (-1)^i*(x^n)/fac(n)
    
    if (r_n(x,i) < 10^(-4))  #Checks tolerance for |R_n| is reached
        {
          #print(i)
          break
        }
  }
  return(fsin)
}

#Exercise 2.4
x <- seq(-5,5, length.out = 150)
y <- seq(-5,5, length.out = 150)
xy <- mesh(x,y)
z <- with(xy,sin(x^2*y))
persp3D(x,y,z)
plotrgl(smooth = TRUE)