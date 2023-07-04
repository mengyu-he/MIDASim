
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

# calculate parameters of beta distribution for one taxon
alpha_beta = function(x) {
  x <- x[x > 0]
  xbar <- mean(x)
  s2 <- var(x)
  alpha <- xbar*(xbar*(1 - xbar)/s2 - 1)
  beta <- alpha*(1 - xbar)/xbar
  c(alpha, beta)
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

# equation solver for mu and eta

solver_mu_sigma = function(mu0, eta0, Ztotal, sample.1.prop, taxa.1.prop,
                           ids.left, n.sample, n.rm) {
  max.iter <- 100
  s1 <- sum(sample.1.prop)
  s2 <- sum(taxa.1.prop)
  for (k in 1:max.iter) {
    mu <- mu0
    eta <- eta0
    
    eta0 <- NULL
    for (i in 1:n.sample) {
      equa <- function(x) (sum(pnorm(mu - x))+n.rm)/Ztotal - sample.1.prop[i]/s1
      eta0[i] <- pracma::fzero(fun = equa , x = 0, tol = 10^-10)$x
    }
    
    mu0 <- NULL
    jj <- 1
    for (j in ids.left) {
      equa <- function(x) sum(pnorm(x - eta0))/Ztotal - taxa.1.prop[j]/s2
      mu0[jj] <- pracma::fzero(fun = equa , x = 0, tol = 10^-10)$x
      jj <- jj + 1
    }
    
    diff <- sum(c(abs(mu0-mu), abs(eta0-eta)))
    if (diff < 10^-5 ){
      break
    }
  }
  return(list(mu0 = mu0, eta0 = eta0))
}

check_taxa = function(taxa.1.prop, n.sample, Ztotal) {
  
  tmp <- -1
  exeed.id <- NULL
  
  while(tmp - length(exeed.id) != 0) {
    exeed.id <- (which(taxa.1.prop/sum(taxa.1.prop) > n.sample/Ztotal + 10^-10))
    Ztotal <- Ztotal - sum(taxa.1.prop[exeed.id]/sum(taxa.1.prop)*Ztotal - n.sample)
    
    tmp <- length(exeed.id)
  }
  
  return(list(exeed.id = exeed.id, Ztotal = Ztotal ))
  
}




