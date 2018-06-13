
# generalized SWCCA, multiple views and multiple dimensinos
#multiSWCCA, where one can take multiple views and maximize the sum of all pairwise correlations
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

# @xlist the views, i.e. list of data matrices
# @klist number of sparse components for each x
# @ulist the initial projections
# @w the initial weighting vector
# @kw number of sparse components for weight vector w
# @ncomp the number of components to solve for
iterate <- function(xlist, klist, ulist, w, kw, ncomp)
{
	#number of views
	K = length(xlist)
	#update each component
	for(n in 1:ncomp) {
		crit=0
		crit.old=0
		while(crit.old==0 || (abs(crit-crit.old)/abs(crit.old) > 10e-2))
		{
  		#initialize starting u
  		u = list()
  		if(is.vector(ulist[[1]])){
  		  u <- ulist
  		}
      else{
        for(i in 1:K) u[[i]] <- ulist[[i]][,n]
      }
  		#optimize each ulist[[i]][,n]
  		for(i in 1:K)
  		{
  			tot = 0
  			#adding the pairwise crosscorrelations for view i with view j
  			for(j in (1:K)[-i]){
  				tot = tot + t(xlist[[i]]) %*% diag(as.vector(w)) %*% xlist[[j]] %*% u[[j]]
  				if(n>1)
  				{
  					mi = ulist[[i]][,1:n-1]
  					mj = ulist[[j]][,1:n-1]
  					if(is.vector(mi))
  					{
  						mi <- matrix(mi, length(mi), 1)
  						mj <- matrix(mj, length(mj), 1)
  					}
  					d = t(mi) %*% t(xlist[[i]]) %*% diag(w) %*% xlist[[j]] %*% mj
  					#check diagonal matrix
  					tot = tot - mi %*% d %*% t(mj) %*% u[[j]]
  				}
  			}
  			if(length(as.vector(tot)) != ncol(xlist[[1]])) print('stop! dimensionality problem!')
  			#update projection i
  			u[[i]]= update(tot, klist[[i]])
  			print(paste('updated projection:', n, i))
  			print(u[[i]])
  		}
  		if(is.vector(ulist[[1]]))
  		{
  		  ulist <- u
  		}
  		else{
  		  for(i in 1:K) ulist[[i]][,n] <- u[[i]]
  		}
  		crit.old = crit
  		crit = 0
  		for(i in 2:K)
  		{
  		  for(j in 1:(i-1))
  		  {
  		    crit <- crit + t(ulist[[i]]) %*% t(xlist[[i]]) %*% diag(w) %*% xlist[[j]] %*% ulist[[j]]
  		  }
  		}
		}
	}
	tot = 0
	for(i in 2:K)
	{
		for(j in 1:(i-1))
		{
			tot= tot + rowSums((xlist[[j]]%*%ulist[[j]]) * (xlist[[i]]%*%ulist[[i]]))
		}
	}
	if(length(tot) != nrow(xlist[[1]])) print('stop! dimensionality problem!')
	#update weights
	w = update(tot, kw)
	return(list(xlist = xlist, klist=klist, ulist = ulist, w = w, kw = kw, ncomp=ncomp))
}

SWCCA <- function(xlist, klist, ulist, w, kw, ncomp)
{
  i_old <- iterate(xlist, klist, ulist, w, kw, ncomp=2)
  i <- do.call(iterate, i_old)
  while(sum((i[[4]] - i_old[[4]])^2) > 10e-2)
  {
    print('loop')
    i_old <- i
    i <- do.call(iterate, i)
  }
  return(i)
}