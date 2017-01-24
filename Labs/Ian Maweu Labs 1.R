#Exercise 1.1
s = sin(pi/3)+log(3)
s^4/2>=exp(s)


#Exercise 1.2
t <- (1:10)
y1 <- t*cos(t)
y2 <- (t^2)*exp(-t)


#Exercise 1.3
n <- (0:20)
r <- 0.5
x <- (r^n)
S <- sum(x)
limit <- 1/(1-r)


# transpose of a matrix A is t(A)
#rbind(A, u, deparse.level = 0) adds rows and cbind adds rows
# true martix multipliction use %*%


#Exercsie 1.4
A <- t(matrix(c(1,-2,3,5,-2,1,1,1,-1,3,1,-4,1,1,1,7), nrow = 4, ncol = 4))
b <- c(-4,5,13,-4)
x <- solve(A,b)

#To obtain square plot use par(pty="s") cmd and to revert ot default use par(pty="m")

#Exercise 1.5
theta = (0:100)*2*pi/100
x <- cos(theta)
y <- sin(theta)
par(pty="s")
plot(x,y, main = "A circle", xlab = "x", ylab = "y")
plot(x,y, main = "A circle", xlab = "x", ylab = "y", col =  3)
par(pty="m")