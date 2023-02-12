
# normalize subjects sum to 1
normalize_rel = function(x) {
  x <- x/rowSums(x)
  return(x)
}

# calculate correlation
cal_corr = function(mat) {

  mat.centered <- t( t(mat) - colMeans(mat) )
  cov <- t(mat.centered) %*% mat.centered / nrow(mat)
  var <- diag(cov)
  sigma <- sqrt(var)
  corr <- t(1/sigma * t( 1/sigma*cov ))
  diag(corr) <- 1

  return(corr)
}


# correct correlation to avoid negative eigenvalues
correct_corr = function(mat) {
  eig <- eigen(mat)
  mat <- eig$vectors %*% tcrossprod(diag(pmax(eig$values,0)), eig$vectors)
  mat <- t((1/sqrt(diag(mat)))*t((1/sqrt(diag(mat))* mat)))
  mat <- (mat+t(mat))/2
  return(mat)
}

# calculate tetrachoric correlation
cal_tetra = function(x) {
  f11 <- as.vector( crossprod( x==1, x==1) )
  f22 <- as.vector( crossprod( x==0, x==0) )
  f12 <- as.vector( crossprod( x==1, x==0) )
  f21 <- as.vector( crossprod( x==0, x==1) )
  f <- f11 + f22 + f12 + f21 + 2
  omega <- (f11 + 0.5) * (f22 + 0.5) / (f12 + 0.5) / (f21 + 0.5)
  p1plus <- (f11 + f12)/f
  pplus1 <- (f11 + f21)/f
  pmin <- ifelse(p1plus < pplus1, p1plus, pplus1)
  c <- (1 - abs(p1plus - pplus1)/5 - (0.5-pmin)^2 )/2
  rho <- cos(pi/(1+omega^c))

  tetra.corr <- matrix(rho, nrow = ncol(x), ncol = ncol(x))
  diag(tetra.corr) <- 1

  return(tetra.corr)
}

# adjust relative abundances iteratively
adjust_iterative = function(rel.perm.f, lib.size){
  for (i in 1:dim(rel.perm.f)[1]) {
    x = rel.perm.f[i, ]
    tmp = which( x < 1/lib.size[i] & x>0 )
    while ( length(tmp)!=0 ) {
      x[-tmp] = x[-tmp] * ( 1-(1/lib.size[i]*length(tmp) - sum(x[tmp]))/sum(x[-tmp]) )
      x[tmp] = 1/lib.size[i]
      tmp = which( x < 1/lib.size[i] & x>0 )
    }
    rel.perm.f[i, ] = x
  }
  return(rel.perm.f)
}


