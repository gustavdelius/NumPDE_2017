###NumMethodsPDE
###Carolin Natemeyer
###Lab1

##Exercise1
s<-sin(pi/3)+log(3)
s
(s^4)/2 >= exp(s)
#Result is TRUE

##Exercise2
t<-1:10
t
t*cos(t)
t^2*exp(-t)

#Exercise3
n<-0:10
r<-0.5
x<-r^n
sum(x) # get 1.999023
limit<-(1/(1-r))
limit #get 2 
m<-(0:20) 
y<-r^m
sum(y) #1.99999

##Exercise 4
A<-matrix(c(1,-2,3,5,-2,1,1,1,-1,3,1,-4,1,1,1,7),nrow=4, ncol=4)
b<-c(-4,5,13,-4)
x<-solve(A,b)
x

##Exercise5
theta=(0:100)*2*pi/100
x<-cos(theta)
y<-sin(theta)
par(pty="s")
plot(x,y,main = "A circle", xlab="X",ylab = "Y", col="red")# colour with 'col="'colour' "'
par(pty="m")
