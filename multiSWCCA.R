
#multiSWCCA, where one can take multiple views and maximize the sum of all pairwise
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

iterate <- function(xlist, klist, ulist, w, kw)
{
	K = length(xlist)
	u = list()
	for(i in 1:K)
	{
		tot = 0
		for(j in (1:K)[-i]){
			tot = tot + t(xlist[[i]]) %*% diag(as.vector(w)) %*% xlist[[j]] %*% ulist[[j]]
		}
		print(dim(tot))
		u = append(u, list(update(tot, klist[[i]])))
	}
	tot = 0
	for(i in 1:K)
	{
		for(j in 1:i)
		{
			tot= tot + (xlist[[j]]%*%u[[j]]) * (xlist[[i]]%*%u[[i]])
		}
	}
	print(dim(tot))
	w = update(tot, kw)
	return(list(xlist, klist, u, w, kw))
}

