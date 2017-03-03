library("plot3Drgl")

f <- function(P, Z) {
    beta*P*(1-P)-Z*P^2/(P^2+nu^2)
}

g <- function(P, Z) {
    gamma*Z*(P^2/(P^2+nu^2) - w)
}

D <- 0.04
nu <- 0.053
beta <- 0.43
gamma <- 0.05
w <- 0.34
hp <- 1
hm <- 1

h <- 0.02
N <- 50
x <- (1:(N-1))*h

tau <- 0.0001
M1 <- 60
M2 <- 500
t <- (0:M1)*M2*tau

kappa <- tau/h^2

p0 <- 0.045

wp <- matrix(0.035, nrow=N-1, ncol=M1+1)
wz <- matrix(0.046, nrow=N-1, ncol=M1+1)

for (j in 1:M1) {
    pn <- wp[, j]
    zn <- wz[, j]
    for (j2 in 1:M2) {
        
        p <- pn
        pp <- c(p[2:(N-1)], p[N-1])
        pm <- c(p0, p[1:(N-2)])
        
        z <- zn
        zp <- c(z[2:(N-1)], z[N-1])
        zm <- c(z[1], z[1:(N-2)])
        
        pn <- p + kappa*
            (D*(pp - 2*p + pm) + hm/2*((pp+p)*(zp-z)-(p+pm)*(z-zm))) +
            tau * f(p, z)
        
        zn <- z + kappa*
            (D*(zp - 2*z + zm) - hp/2*((zp+z)*(pp-p)-(z+zm)*(p-pm))) +
            tau * g(p, z)
    }
    wp[, j+1] <- pn
    wz[, j+1] <- zn
}

persp3D(x, t, wp,
        xlab="x", ylab="t", zlab="P",
        ticktype="detailed", nticks=4)

plotrgl(lighting = TRUE)

persp3D(x, t, wz,
        xlab="x", ylab="t", zlab="Z",
        ticktype="detailed", nticks=4)