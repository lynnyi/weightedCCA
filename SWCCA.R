
#sparse weighted CCA for one dimension

#update function takes top k elements of vector z and renormalizes the vectors
update <- function(z, k)
{
	if(all(z==0))
	{
		return(z)
	}
	i <- order(-abs(z))[k+1:length(z)]
	z[i] <- 0
	z <- z / sqrt(sum(z^2))
	z
}

updatePositive <- function(z, k)
{
  if(all(z==0))
  {
    return(z)
  }
  i <- order(-(z))[k+1:length(z)]
  z[i] <- 0
  z <- z / sqrt(sum(z^2))
  z
}

iterate <- function(u, X, w, Y, v, ku, kv, kw)
{

	#4000x1
	u <- Y %*% v
	#4000x1
	u <- w * u
	#3*1
	u <- t(X) %*% matrix(u, length(u), 1)
	u <- update(u, ku)

	v <- X %*% u
	v <- w * v
	v <- t(Y) %*% matrix(v, length(v), 1)
	v <- update(v, kv)

	t1 <- X%*%u
	t2 <- Y%*%v
	w <- t1 * t2
	w <- update(w, kw)

	return(list(u, X, w, Y, v, ku, kv, kw))
}


i_old <- iterate(ulist[[1]], xlist[[1]], w, xlist[[2]], ulist[[2]], 3, 3, 1000)
i <- do.call(iterate, i_old)
while(sum((i[[3]] - i_old[[3]])^2) > 10e-16)
{
  print('loop')
  i_old <- i
  i <- do.call(iterate, i)
}

