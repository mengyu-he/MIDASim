#' @importFrom stats var pnorm optimize plnorm qlnorm qgamma runif pnorm uniroot
#'
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
  rho <- cos(base::pi/(1+omega^c))

  tetra.corr <- matrix(rho, nrow = ncol(x), ncol = ncol(x))
  diag(tetra.corr) <- 1

  return(tetra.corr)
}

# equation solver for mu and eta

solver_theta_eta = function(theta0, eta0, Ztotal, sample.1.prop, taxa.1.prop,
                            ids.left, n.sample, n.rm) {
  max.iter <- 100
  Zi <- Ztotal * sample.1.prop / sum(sample.1.prop)
  Zj <- Ztotal * taxa.1.prop / sum(taxa.1.prop)
  for (k in 1:max.iter) {
    theta <- theta0
    eta <- eta0

    eta0 <- NULL
    for (i in 1:n.sample) {
      equa <- function(x) (sum(pnorm(theta + x))+n.rm) - Zi[i]
      eta0[i] <- pracma::fzero(fun = equa , x = c(-1000, 1000), tol = 10^-10)$x
    }

    theta0 <- NULL
    jj <- 1
    for (j in ids.left) {
      equa <- function(x) sum(pnorm(x + eta0)) - Zj[j]
      theta0[jj] <- pracma::fzero(fun = equa , x = c(-1000, 1000), tol = 10^-10)$x
      jj <- jj + 1
    }

    diff <- sum(c(abs(theta0-theta), abs(eta0-eta)))
    if (diff < 10^-5 ){
      break
    }
  }
  return(list(theta0 = theta0, eta0 = eta0))
}

check_taxa = function(taxa.1.prop, n.sample, Ztotal) {

  tmp <- -1
  exeed.id <- NULL
  one.id <- which(taxa.1.prop == 1)
  zero.id <- which(taxa.1.prop == 0)

  while(tmp - length(exeed.id) != 0) {
    exeed.id <- (which(taxa.1.prop/sum(taxa.1.prop) > n.sample/Ztotal + 10^-10))
    Ztotal <- Ztotal - sum(taxa.1.prop[exeed.id]/sum(taxa.1.prop)*Ztotal - n.sample)

    tmp <- length(exeed.id)
  }
  one.id <- union(one.id, exeed.id)

  return(list(one.id = one.id, zero.id = zero.id, Ztotal = Ztotal ))

}

fit_ggamma = function(cv.sq, mean.rel.abund, lib.size, mat01, m3) {

  n.taxa = length(mean.rel.abund)

  Q.est = mu.est = sigma.est = rep(0, n.taxa)
  for (j in 1:n.taxa) {

    delta = (mat01[,j] > 0)

    if (is.na(m3[j])) {
      res = suppressWarnings(optimize(f = binary.log.like,
                                      interval = c(-150,150),
                                      maximum = TRUE,
                                      m1 = mean.rel.abund[j],
                                      cv.sq = cv.sq[j],
                                      lib.size = lib.size,
                                      delta = delta))
    } else {
      res = suppressWarnings(optimize(f = obj.m3,
                                      interval = c(-150,150),
                                      maximum = TRUE,
                                      m1 = mean.rel.abund[j],
                                      cv.sq = cv.sq[j],
                                      m3 = m3[j],
                                      lib.size = lib.size,
                                      delta = delta))
    }

    Q.est[j] = res$maximum
    params = mu.sigma.of.Q( Q.est[j], mean.rel.abund[j], cv.sq[j] )
    mu.est[j] = params$mu
    sigma.est[j] = params$sigma
  }

  return(list(Q.est = Q.est, mu.est = mu.est, sigma.est = sigma.est))
}


###### gengamma #######

binary.log.like = function(Q, m1, cv.sq, lib.size, delta, eps = 10^-3,
                           p.gamma=p.gamma.nr) {
  #
  #   calculates the binary log-likelihood for Q, sigma(Q) and mu(Q, sigma(Q))
  #        given m1=mean relative abundance, cv.sq=square of empirical CV
  #
  if (abs(Q) > eps) {
    K = 1 / Q
    sigma = sigma.of.k(K, cv.sq)      # match CV2
    mu = gammln(abs(K) - sign(K) * sigma) - gammln(abs(K)) - log(m1)     # match mean
    z = sign(K) * (log(lib.size) - mu) / sigma

    if (K > 0) {
      log.like = sum( p.gamma( x = exp(z[delta == 1]), z = z[delta == 1], a = abs(K), lower.tail = TRUE, log.p = TRUE ) ) +
        sum( p.gamma( x = exp(z[delta == 0]), z = z[delta == 0], a = abs(K), lower.tail = FALSE, log.p = TRUE) )
    } else if (K < 0) {
      log.like = sum( p.gamma( x = exp(z[delta == 1]), z = z[delta == 1], a = abs(K), lower.tail = FALSE, log.p = TRUE ) ) +
        sum( p.gamma( x = exp(z[delta == 0]), z = z[delta == 0], a = abs(K), lower.tail = TRUE, log.p = TRUE) )
    }
  } else {
    sigma.sq = log(1 + cv.sq)
    mu = 0.5 * sigma.sq - log(m1)
    z = (log(lib.size) - mu) / sqrt(sigma.sq)

    log.like = sum( pnorm(z[delta == 1], lower.tail = TRUE, log.p = TRUE) ) +
      sum( pnorm(z[delta == 0], lower.tail = FALSE, log.p = TRUE) )
  }
  return(log.like)
}

#' @noRd
sigma.of.k = function(K, cv.sq) {
  #
  #   calculates sigma given a value of K (WARNING - K, not Q)
  #
  if (K > 0) {
    upper = K / 2 - 1 / (2 * gamma(K / 2) * (10 * cv.sq + 2))
    res = suppressWarnings(uniroot(f = sigma.eqn, interval=c(0, upper),
                                   K = K, cv.sq = cv.sq, tol = 10^-12))
    sigma = res$root
  }
  else if (K<0) {

    res = suppressWarnings(uniroot(f = sigma.eqn, interval = c(0, 2 * abs(K)),
                                   extendInt='upX',
                                   K = K, cv.sq = cv.sq, tol = 10^-12))
    sigma = res$root
  }
  return(sigma)
}

#' @noRd
sigma.eqn = function(x, K, cv.sq) {
  log.ratio = gammln(abs(K) - sign(K) * 2 * x) + gammln(abs(K)) - 2 * gammln(abs(K) - sign(K) * x)
  value = exp(log.ratio) - 1 - cv.sq
  return(value)
}

mu.sigma.of.Q = function(Q, m1, cv.sq) {
  #
  #	calculates mu and sigma given a value of Q, returns a list that can be used as argument for functions calling for params.
  #
  if (Q != 0) {
    K = 1/Q
    sigma = sigma.of.k(K,cv.sq)
    mu = gammln(abs(K) - sign(K) * sigma) - gammln(abs(K)) - log(m1)
  }
  else if (Q == 0) {
    sigma.sq = log(1 + cv.sq)
    mu = 0.5 * sigma.sq - log(m1)
    sigma = sqrt(sigma.sq)
  }
  res = list(mu = mu, sigma = sigma, Q = Q)
  return(res)
}

mu.of.sigma.Q = function(sigma, Q, m1) {
  K = 1/Q
  mu = gammln(abs(K) - sign(K) * sigma) - gammln(abs(K)) - log(m1)
  return(mu)
}

pr.zero = function(params, lib.size, p.gamma = p.gamma.nr) {
  #
  #	calculated the expected proportion of zero cells
  #
  mu = params[[1]]
  sigma = params[[2]]
  Q = params[[3]]
  if (Q != 0) {
    K = 1/Q
    z = sign(K) * (log(lib.size) - mu) / sigma
    use.lower.tail = ifelse(K > 0, FALSE, TRUE)
    pr.zero = p.gamma( x = exp(z), z = z, a = abs(K),
                       lower.tail = use.lower.tail, log.p=FALSE )
  } else if (Q == 0) {
    z = plnorm(lib.size, meanlog = mu, sdlog = sigma, lower.tail=FALSE)
  }
  return(pr.zero)
}

get.prob01.mat = function(mu.est, sigma.est, Q.est, lib.size, p.gamma = p.gamma.nr) {

  n.taxa = length(mu.est)
  prob01.mat = matrix(nrow = length(lib.size), ncol = n.taxa)
  for (j in 1:n.taxa) {

    params = list(mu.est[j], sigma.est[j], Q.est[j])
    prob01.mat[, j] = 1 - pr.zero(params, lib.size, p.gamma = p.gamma.nr)

  }
  return(prob01.mat)

}

p.gamma.nr = function(a, x, z, ITMAX = 1000, log.p = FALSE,
                      lower.tail = TRUE) {
  #  calculates the cdf of the gamma distribution, using code based on numerical recipes

  ind = (x < a + 1.0)
  x1 = x[ind]; z1 = z[ind]
  x2 = x[!ind]; z2 = z[!ind]

  if (lower.tail) {

    if (log.p) {
      gammp1 = gser(a, x1, z1, ITMAX)$lgamser
      tmp = gcf(a, x2, z2, ITMAX)$gammcf
      gammp2 = ifelse(tmp < 1, log(1 - tmp), - tmp)
    } else {
      gammp1 = gser(a, x1, z1, ITMAX)$gamser
      gammp2 = 1 - gcf(a, x2, z2, ITMAX)$gammcf
    }

  } else {

    if (log.p) {
      tmp = gser(a, x1, z1, ITMAX)$gamser
      gammp1 = ifelse(tmp < 1, log(1 - tmp), - tmp)
      gammp2 = gcf(a, x2, z2, ITMAX)$lgammcf
    } else {
      gammp1 = 1 - gser(a, x1, z1, ITMAX)$gamser
      gammp2 = gcf(a, x2, z2, ITMAX)$gammcf
    }

  }

  gammp = vector(length = length(x))
  gammp[ind] = gammp1; gammp[!ind] = gammp2;

  return(gammp)
}


gser = function(a, x = exp(z), z = log(x), ITMAX, EPS = 3.0e-7) {
  #USES gammln
  #Returns the incomplete gamma function P (a, x) evaluated by its series representation as
  #gamser. Also returns ln Γ(a) as gln.
  #Parameters: ITMAX is the maximum allowed number of iterations;
  #   gln=lgamma(a)
  gln = gammln(a)

  ap = a
  sum = rep(1.0 / a, length(x))
  del = sum
  count = 0
  for (n in 1:ITMAX) {
    ap = ap + 1.0
    del = del * x / ap
    sum = sum + del
    if(all(abs(del) < abs(sum) * EPS)) break
    count = count + 1
  }
  #if (count == ITMAX) print('a too large, ITMAX too small in gser')
  lgamser = log(sum) -x + a * z - gln
  gamser = sum * exp( -x + a * z - gln)
  res = list(gamser = gamser, gln = gln, sum = sum,lgamser = lgamser)
  return(res)
}



gcf = function(a, x = exp(z), z = log(x), ITMAX, EPS = 3.e-16,FPMIN = 1.e-30) {
  #USES lgamma from R
  #Returns the incomplete gamma function Q(a, x) evaluated by its continued fraction repre-
  #sentation as gammcf. Also returns ln Γ(a) as gln.
  #Parameters: ITMAX is the maximum allowed number of iterations; EPS is the relative accu-
  #racy; FPMIN is a number near the smallest representable floating-point number.
  gln = gammln(a)

  ind = (z > 709)
  z1 = z[ind]; x1 = x[ind]
  z2 = z[!ind]; x2 = x[!ind]

  lgammcf1 = -x1 + (a - 1) * z1

  lgammcf2 <- gammcf2 <- NULL
  if (length(x2) > 0) {
    b = x2 + 1.- a
    #Set up for evaluating continued fraction by modified
    #   Lentz’s method (§5.2) with b0 = 0.c=1./FPMIN
    c = 1 / FPMIN
    d = 1. / b
    h = d
    log.h = log(d)
    count = 0
    for (i in 1:ITMAX) {
      an = -i * (i - a)
      b = b + 2.
      d = an * d + b
      d[abs(d) < FPMIN] = FPMIN
      c = b + an/c
      c[abs(c) < FPMIN] = FPMIN
      d = 1. / d
      del = d * c
      h = h * del
      log.h = log.h + log(del)
      if( all(abs(del - 1.) < EPS) ) break
      count = count + 1
    }
    #if (count == ITMAX) print('a too large, ITMAX too small in gcf')
    lgammcf2 = log.h - x2 + a * z2 - gln
    gammcf2 = exp(- x2 + a * z2 - gln) * h #Put factors in front.
  }

  lgammcf  <-  gammcf <- vector(length = length(x))
  lgammcf[ind] = lgammcf1; lgammcf[!ind] = lgammcf2;
  gammcf[ind] = exp(lgammcf1); gammcf[!ind] = gammcf2;

  res = list(gammcf = gammcf, gln = gln, lgammcf = lgammcf)
  return(res)
}


gammln = function(x) {
  #Returns the value ln[Γ(x)] for x > 0.
  cof=c(76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5)
  stp=2.5066282746310005e0
  #   x=xx
  y=x
  tmp=x+5.5
  tmp=(x+0.5)*log(tmp)-tmp
  ser=1.000000000190015
  for (j in 1:6) {
    y=y+1.0
    ser=ser+cof[j]/y
  }
  gammln=tmp+log(stp*ser/x)
  return(gammln)
}


p.gen.gamma = function(t, params, p.gamma = p.gamma.nr,
                       lower.tail = TRUE, log.p = FALSE) {
  #   calculates the CDF of the generalized gamma in terms of input t.
  mu = params[[1]]             #params$mu
  sigma = params[[2]]          #params$sigma
  Q = params[[3]]              #params$Q
  if (Q != 0) {
    K = 1 / abs(Q)
    z = sign(Q) * (log(t) - mu) / sigma
    lower.tail = identical(Q > 0, lower.tail)
    cdf = rep(0,length(pi))
    cdf = p.gamma(z = z, x = exp(z), a = K, log.p = log.p, lower.tail = lower.tail)
  } else {
    z = (log(t) - mu) / sigma
    cdf = plnorm(z, meanlog = 0, sdlog = 1, lower.tail = lower.tail, log.p = log.p)
  }
  return(cdf)
}


q.gen.gamma = function(u, params, p.gamma = p.gamma.nr ) {
  #   returns the values of pi (and t) that correspond to a quantile level u
  mu = params[[1]]             #params$mu
  sigma = params[[2]]          #params$sigma
  Q = params[[3]]              #params$Q
  if (Q < 0) u = 1-u
  if (Q == 0) {
    t = qlnorm(u, meanlog = mu, sdlog = sigma)
    pi = 1 / t
  } else {
    K = 1 / abs(Q)
    tmp = qgamma(u, shape = K)
    log.q = log( tmp )
    ind = is.infinite(log.q)
    log.q[ind] = 1/K* (log(u[ind]) + gammln(K + 1))
    t = exp( mu + sign(Q) * sigma * log.q )
    pi = 1 / t
  }
  res = list(pi = pi, log.q = log.q, t = t)
  return(res)
}


r.gen.gamma = function(n, params, p.gamma = p.gamma.nr, value = 'pi') {
  #   generates a random sample from the generalized gamma distribution.
  #   returns either the value of pi or the value of t.
  u = runif(n)
  if (value == 'pi') {
    r = q.gen.gamma(u = u, params = params, p.gamma = p.gamma )$pi
  } else if (value == 't') {
    r = q.gen.gamma(u = u, params = params, p.gamma = p.gamma )$t
  } else {
    print( 'value must be either pi or t' )
    return()
  }
  return(r)
}

est.mu = function(mu, sigma, Q) {
  s = sign(Q)
  k = rep(0,length(Q))
  k[ Q!=0 ] = 1/abs( Q[Q!=0] )

  est = ifelse(Q == 0, exp( -mu + 0.5 * sigma^2 ),
               exp( -mu + gammln(k - s*sigma) - gammln(k) ) )
  return(est)
}

obj.m3 = function(Q, m1, cv.sq, m3, lib.size, delta, eps = 10^-3,
                  p.gamma = p.gamma.nr) {

  if (abs(Q) > eps) {

    K = 1 / Q
    sigma = sigma.of.k(K, cv.sq)      # match CV2
    mu = gammln(abs(K) - sign(K) * sigma) - gammln(abs(K)) - log(m1)     # match mean

    if (abs(K) - 3 * sign(K) * sigma > 0){
      log.m3 = gammln(abs(K) - 3 * sign(K) * sigma) - gammln(abs(K)) - mu * 3
      # non central moment
      f = - (m3 - exp(log.m3) )^2
    } else {
      f = -Inf
    }

  } else {
    sigma.sq = log(1 + cv.sq)
    mu = 0.5 * sigma.sq - log(m1)
    z = (log(lib.size) - mu) / sqrt(sigma.sq)

    f = - (m3 - exp(- 3 * mu + 9 * sigma.sq / 2))^2
  }
  return(f)
}
