D <- function(w, alpha = 2){
    return ((1-w^2)^alpha)
}
D -> D1 -> D2
alpha <- 2

P <- function(x, r){
    return(ifelse(abs(x) <= r,1,0))
}

K <- function(x, g, r) {
    Re(fft(fft(P(x, r)*x)*fft(g), inverse = TRUE))*2/length(x)
}

tau_1_1 = 1
tau_1_2 = 1
tau_2_2 = 1

M_1 = 0.95
M_2 = 0.05

w_1 <- -0.5
w_2 <- 0.5
sigma_1 <- sigma_2 <- 0.05
q_1 <- q_2 <- 0.5

g1_0 <- function(x) {
    exp(-(x/0.4)^2/2)/(sqrt(2*pi)*0.4)
}
g2_0 <- function(x) {
    q_1*exp(-((x-w_1)/sigma_1)^2/2)/(sqrt(2*pi)*sigma_1)+
         q_2* exp(-((x-w_2)/sigma_2)^2/2)/(sqrt(2*pi)*sigma_2)
}

lambda = 5*10^-3
lambda -> lambda_1_1 -> lambda_1_2 -> lambda_2_2


explicit <- function(g1_0=function(x) 0*x, g2_0=function(x) 0*x, 
                     N=400, T=10, M=10000, 
                     r1 = 0.5, r2 = 0.5, r3 = 0.5) {
    h <- 2/N
    x <- -1+(0:N)*h
    
    tau <- T/M
    M1 <- 100
    M2 <- round(M/M1)
    t <- (0:M1)*M2*tau
    
    d1 <- D1(x, alpha)^2
    d2 <- D2(x, alpha)^2
    
    w1 <- matrix(0, nrow=N+1, ncol=M1+1)
    w2 <- matrix(0, nrow=N+1, ncol=M1+1)
    
    gn1 <- g1_0(x)
    gn2 <- g2_0(x)
    w1[, 1] <- gn1
    w2[, 1] <- gn2
    
    for (j in 1:M1) {
        for (j2 in 1:M2) {
            
            g1 <- gn1
            g2 <- gn2
            
            a11 <- (K(x, g1, r1)/tau_1_1+K(x, g2, r3)/(2*tau_1_2))*g1
            a11p <- c(a11[2:(N+1)], a11[N])
            a11m <- c(a11[2], a11[1:N])
            a12 <- (lambda_1_1*M_1/(2*tau_1_1)+lambda_1_2*M_2/(4*tau_1_2))*d1*g1
            a12p <- c(a12[2:(N+1)], a12[N])
            a12m <- c(a12[2], a12[1:N])
            
            a21 <- K(x, g2, r2)/tau_2_2*g2
            a21p <- c(a21[2:(N+1)], a21[N])
            a21m <- c(a21[2], a21[1:N])
            a22 <- lambda_2_2*M_2/(2*tau_2_2)*d2*g2
            a22p <- c(a22[2:(N+1)], a22[N])
            a22m <- c(a22[2], a22[1:N])
            
            gn1 <- g1 + tau * (
                1/(2*h) * (a11p - a11m) + 1/h^2 * (a12p - 2*a12 + a12m)
            )
            
            gn2 <- g2 + tau * (
                1/(2*h) * (a21p - a21m) + 1/h^2 * (a22p - 2*a22 + a22m)
            )
        }
        w1[, j+1] <- gn1
        w2[, j+1] <- gn2
    }
    list(x=x, t=t, g1=w1, g2=w2)
}

sol <- explicit(g1_0, g2_0, N=200, T=0.02, M=1000)

persp3D(sol$x, sol$t, sol$g1,
        xlab="x", ylab="t", zlab="g1",
        ticktype="detailed", nticks=4)

persp3D(sol$x, sol$t, sol$g2,
        xlab="x", ylab="t", zlab="g2",
        ticktype="detailed", nticks=4)

plot(sol$x, sol$g2[, 101], type="l", lty="dotted", xlab="x", ylab="")
lines(sol$x, sol$g1[, 101], lty="solid")
legend("topleft", legend=c("g1", "g2"), lty=c("solid", "dotted"))
