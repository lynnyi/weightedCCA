library(ggplot2)
library(irlba)
nsamples <- 1000
W1 = matrix(c(-10, 0, 0, 0, 10, 0), 2, 3, byrow=TRUE)
W2 = matrix(c(-2, 0, -2, 0, 2 , 0), 2, 3, byrow=TRUE)
z = rnorm(nsamples*2)
z <- matrix(z, nsamples, 2)
d1 = z %*% W1
d2 = z %*% W2
noise = rnorm(nsamples*3, 0, 5)
d1 <- d1 + noise
noise = rnorm(nsamples*3, 0, 5)
d2 <- d2 + noise

t <- irlba(t(d1) %*% d2, 2)

#centroid noise
nbsamples <- 3000
C1 <- c(0)
n = rnorm(nbsamples * length(C1), C1, sd=10)
n = matrix(n, nbsamples, 3, byrow=TRUE)
d1 = rbind(d1, n)
n = rnorm(nbsamples * length(C1), C1)
n = matrix(n, nbsamples, 3, byrow=TRUE)
d2 = rbind(d2, n)

f <- irlba(t(d1) %*% d2, 2)

c11 <- d1 %*% t$u
c12 <- d2 %*% t$v
c11f <- d1 %*% f$u
c12f <- d2 %*% f$v

cor(c11[1:nsamples, 1], c12[1:nsamples, 1])
cor(c11f[1:nsamples, 1], c12f[1:nsamples, 1])
cor(c11f[1:nsamples,1], z)
cor(c11[1:nsamples,1], z)
cor(c11f[nsamples:4000,1], c12f[nsamples:4000,1])
cor(c11f, c12f)



weights <- c(rep(1, 1000), rep(.01, 3000))
d1_weighted <- d1 * weights
d2_weighted  <-  d2 * weights
new <- irlba(t(d1_weighted) %*% d2_weighted, 1)
c11new <- d1%*%new$u
c12new <- d2%*%new$v
cor(c11new, c12new)

f <- factor(c(rep(1, 1000), rep(2,3000)))
g <- ggplot(data.frame(c11new), aes(x=c11f, y=c12f, colour=f))+geom_point()
png('./sim.png')
print(g)
dev.off()


cor(d1%*%new$u, d2%*%new$v)
cor(c12new[1:1000,], z)

nbsamples <- 20000
#gaussian noise
notcor <- matrix(rnorm(nbsamples*3), nbsamples, 3)
d1 <- rbind(d1, notcor)
notcor <- matrix(rnorm(nbsamples*3), nbsamples, 3)
d2 <- rbind(d2, notcor)


