##'Odd log-logistic Weibull distribution in Gamlss

OLLW <- function (mu.link="log", sigma.link="log", nu.link ="log") 
{
  mstats <- checklink("mu.link", "odd log-logistic W", substitute(mu.link), 
                      c("inverse", "log", "identity", "own"))
  dstats <- checklink("sigma.link", "odd log-logistic W", substitute(sigma.link), #
                      c("inverse", "log", "identity", "own"))
  vstats <- checklink("nu.link", "odd log-logistic W",substitute(nu.link), 
                      c("logshifted", "log", "identity", "own"))
  
  structure(
    list(family = c("OLLW", "odd log-logistic W"),
         parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE), 
         nopar = 3, 
         type = "Continuous",
         mu.link = as.character(substitute(mu.link)),  
         sigma.link = as.character(substitute(sigma.link)), 
         nu.link = as.character(substitute(nu.link)), 
         mu.linkfun = mstats$linkfun, 
         sigma.linkfun = dstats$linkfun, 
         nu.linkfun = vstats$linkfun,
         mu.linkinv = mstats$linkinv, 
         sigma.linkinv = dstats$linkinv,
         nu.linkinv = vstats$linkinv,
         mu.dr = mstats$mu.eta, 
         sigma.dr = dstats$mu.eta, 
         nu.dr = vstats$mu.eta,
         
         dldm = function(y,mu,sigma,nu){ #----------------------------------------------------- ok
           nd1 = gamlss:::numeric.deriv(dOLLW(y, mu, sigma, nu,log = T),"mu", delta = 1e-04)
           dldm = as.vector(attr(nd1, "gradient")) 
           dldm
         },
         d2ldm2 = function(y,mu,sigma,nu){#----------------------------------------------------- ok
           nd1 = gamlss:::numeric.deriv(dOLLW(y, mu, sigma, nu,log = T),"mu", delta = 1e-04)
           dldm = as.vector(attr(nd1, "gradient")) 
           d2ldm2 = -dldm * dldm
         },  
         dldd = function(y,mu,sigma,nu){#----------------------------------------------------- ok  
           nd1 = gamlss:::numeric.deriv(dOLLW(y, mu, sigma, nu,log =T),"sigma", delta = 1e-04)
           dldd = as.vector(attr(nd1, "gradient"))
           dldd
         } ,
         d2ldd2 = function(y,mu,sigma,nu){#----------------------------------------------------- ok
           nd1 = gamlss:::numeric.deriv(dOLLW(y, mu, sigma, nu,log =T),"sigma", delta = 1e-04)
           dldd = as.vector(attr(nd1, "gradient"))
           d2ldd2 = -dldd*dldd
           d2ldd2 
         },
         dldv = function(y,mu,sigma,nu){#----------------------------------------------------- ok
           nd1 = gamlss:::numeric.deriv(dOLLW(y, mu, sigma, nu,log =T),"nu", delta = 1e-04)
           dldv = as.vector(attr(nd1, "gradient"))
           dldv 
         },
         d2ldv2 = function(y,mu,sigma,nu){#----------------------------------------------------- ok 
           nd1 = gamlss:::numeric.deriv(dOLLW(y, mu, sigma, nu,log =T),"nu", delta = 1e-04)
           dldv = as.vector(attr(nd1, "gradient"))
           d2ldv2 = -dldv * dldv
         },
         d2ldmdd = function(y,mu,sigma,nu){#----------------------------------------------------- ok
           nd1 = gamlss:::numeric.deriv(dOLLW(y, mu, sigma, nu,log= TRUE), "mu", delta = 1e-04)
           dldm = as.vector(attr(nd1, "gradient"))
           nd1 = gamlss:::numeric.deriv(dOLLW(y, mu, sigma, nu,log=TRUE), "sigma", delta = 1e-04)
           dldd = as.vector(attr(nd1, "gradient"))           
           d2ldmdd = -dldm * dldd
           d2ldmdd               
         },
         d2ldmdv = function(y,mu,sigma,nu){#----------------------------------------------------- ok
           nd1 = gamlss:::numeric.deriv(dOLLW(y, mu, sigma, nu,log = TRUE), "mu", delta = 1e-04)
           dldm = as.vector(attr(nd1, "gradient"))
           nd1 = gamlss:::numeric.deriv(dOLLW(y, mu, sigma, nu,log = TRUE), "nu", delta = 1e-04)
           
           dldv = as.vector(attr(nd1, "gradient"))
           d2ldmdv = -dldm * dldv
           d2ldmdv			
         },
         d2ldddv = function(y,mu,sigma,nu){#----------------------------------------------------- ok
           nd1 = gamlss:::numeric.deriv(dOLLW(y, mu, sigma, nu,log = TRUE), "sigma", delta = 1e-04)
           dldm = as.vector(attr(nd1, "gradient"))
           nd1 = gamlss:::numeric.deriv(dOLLW(y, mu, sigma, nu,log = TRUE), "nu", delta = 1e-04)
           
           dldv = as.vector(attr(nd1, "gradient"))
           d2ldmdv = -dldm * dldv
           d2ldmdv	
         },
         G.dev.incr  = function(y,mu,sigma,nu,...) {
           -2*dOLLW(y,mu=mu,sigma=sigma,nu=nu,log=TRUE)
         }, 
         rqres = expression(rqres(pfun="pOLLW", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu)),
         mu.initial = expression( mu <- (y+mean(y))/2), 
         sigma.initial = expression( sigma <- rep(sd(y)/2,length(y))), 
         nu.initial = expression( nu <- rep(max(2*(mean(y)-median(y)),0.1), length(y))), # 2*max((mean(y)-median(y),0.1) # (mean(y)+sd(y))/4
         mu.valid = function(mu) all(mu > 0),  
         sigma.valid = function(sigma)  all(sigma > 0),
         nu.valid = function(nu) all(nu > 0), 
         y.valid = function(y) all(y > 0)
    ),
    class = c("gamlss.family","family"))
}

#-----------------------------------------------------------------
# Probability Density Function
dOLLW <- function(x,mu  = 0.2,sigma = 0.4, nu=1, log = FALSE){
  if (any(nu < 0)) stop(paste("nu must be positive", "\n", ""))
  if (any(mu < 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0)) stop(paste("sigma must be positive", "\n", ""))
  if (any(x < 0))  stop(paste("x must be positive", "\n", ""))
  G <- pWEI(x,mu=mu,sigma=sigma)
  g <- dWEI(x,mu=mu,sigma=sigma)
  G_ <- 1-G
  Ft <- (G^(nu))/((G^nu)+(G_^nu))
  f <- (nu*g*((G*(1-G))^(nu-1)))/(((G^nu)+(1-G)^nu)^2)
  if(log==FALSE) fx  <- f else fx <- log(f) 
  fx
}
#----------------------------------------------------------------- 
# Cumulative Density Function
pOLLW <- function(q,mu  = 0.2,sigma = 0.4, nu=1, lower.tail = TRUE, log.p = FALSE){
  if (any(nu < 0)) stop(paste("nu must be positive", "\n", ""))
  if (any(mu < 0)) stop(paste("mu must be positive", "\n", ""))
  if (any(sigma < 0)) stop(paste("sigma must be positive", "\n", ""))
  if (any(q < 0))  stop(paste("q must be positive", "\n", ""))
  G <- pWEI(q,mu=mu,sigma=sigma)
  g <- dWEI(q,mu=mu,sigma=sigma)
  G_ <- 1-G
  cdf <- (G^(nu))/((G^nu)+(G_^nu))
  if(lower.tail==TRUE) cdf  <- cdf else  cdf <- 1-cdf 
  if(log.p==FALSE) cdf  <- cdf else  cdf <- log(cdf) 
  cdf
}
#-----------------------------------------------------------------  
# Quantile Function
qOLLW <-  function(p,mu=1,sigma=1,nu=0.5, lower.tail = TRUE, log.p = FALSE){   
  if (any(mu <= 0))  stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))
  if (any(p < 0) | any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))    
  if (log.p==TRUE) p <- exp(p) else p <- p
  if (lower.tail==TRUE) p <- p else p <- 1-p
  e1 <- 1/nu
  p1 <- (p^(e1))/(((1-p)^e1)+p^e1)
  q <- qWEI(p1, mu=mu , sigma=sigma)
  q
}
#----------------------------------------------------------------- 
# Random generating function
rOLLW <- function(n, mu=1,sigma=1,nu=0.5){
  if (any(mu <= 0))  stop(paste("mu must be positive", "\n", ""))
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", ""))   
  if (any(nu <= 0))  stop(paste("nu must be positive", "\n", ""))  
  if (any(n <= 0))  stop(paste("n must be a positive integer", "\n", ""))  
  n <- ceiling(n)
  u <- runif(n,0,1)
  r <- qOLLW(u, mu =mu, sigma =sigma, nu=nu)
  r
}
