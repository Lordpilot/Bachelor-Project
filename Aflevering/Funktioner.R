library(pracma)
library(lmtest)
library(parallel)
library(MASS)
library(foreach)
library(doParallel)


#Functions
#Fourier Estimation from given periods (ex. every s=24 hours we see a period).
Fourier.Est <- function(ts, s. = c(24)){    #Use the season
  n <- length(ts)
  m <- length(s.)
  t <- seq(1,n)
  
  freq_p <- t %*% t(1/s.)
  A <- cos(2*pi*freq_p)
  B <- sin(2*pi*freq_p)
  
  factors <- B[,c(1,1)%*%t(seq(1:dim(B)[2]))]
  factors[,seq(1,dim(factors)[2],2)] <- A
  
  Lin.Mod <- lm(ts ~ factors)
  res <- list(Lin.Mod = Lin.Mod, factors = factors)
  return(res)
}

#Fourier Estimation from given frequencies (ex. the above example gives f=1/24)
Fourier.Est.Mod <- function(ts, f. = c(1/24)){
  n <- length(ts)
  m <- length(f.)
  t <- seq(1,n)
  
  freq_p <- t %*% t(s.)
  A <- cos(2*pi*freq_p)
  B <- sin(2*pi*freq_p)
  
  factors <- B[,c(1,1)%*%t(seq(1:dim(B)[2]))]
  factors[,seq(1,dim(factors)[2],2)] <- A
  
  Lin.Mod <- lm(ts ~ factors)
  res <- list(Lin.Mod = Lin.Mod, factors = factors)
  return(res)
}

#Finding all parameters through a full fourier expansion.
Fourier.Full <- function(ts){
  n <- length(ts)
  t <- seq(1,n)
  if(n %% 2 == 1){
    k <- (n-1)/2
    f_freq <- seq(1,k)/n
    coef <- NULL
    for(i in 1:k){
      A <- (2/n)*sum(ts*cos(2*pi*t*f_freq[i]))
      B <- (2/n)*sum(ts*sin(2*pi*t*f_freq[i]))
      coef <- c(coef, A, B)
    }
  }
  else{
    k <- n/2
    f_freq <- seq(1,k)/n
    coef <- NULL
    for(i in 1:(k-1)){
      A <- (2/n)*sum(ts*cos(2*pi*t*f_freq[i]))
      B <- (2/n)*sum(ts*sin(2*pi*t*f_freq[i]))
      coef <- c(coef, A, B)
    }
    Ak <- (1/n)*sum((-1)^t * ts) 
    coef <- c(coef, Ak, 0)
  }
  return(coef)
}

#Combines Fourier.Est and Fourier.Full
Fourier <- function(ts, s. = TRUE){
  s <- s.
  if(s==TRUE){
    n <- length(ts)
    t <- seq(1,n)
    if(n %% 2 == 1){
      k <- (n-1)/2
      f_freq <- seq(1,k)/n
      coef <- NULL
      for(i in 1:k){
        A <- (2/n)*sum(ts*cos(2*pi*t*f_freq[i]))
        B <- (2/n)*sum(ts*sin(2*pi*t*f_freq[i]))
        res <- c(coef, A, B)
      }
    }
    else{
      k <- n/2
      f_freq <- seq(1,k)/n
      coef <- NULL
      for(i in 1:(k-1)){
        A <- (2/n)*sum(ts*cos(2*pi*t*f_freq[i]))
        B <- (2/n)*sum(ts*sin(2*pi*t*f_freq[i]))
        coef <- c(coef, A, B)
      }
      Ak <- (1/n)*sum((-1)^t * ts) 
      res <- c(coef, Ak, 0)
    }
  }
  else{
    n <- length(ts)
    m <- length(s.)
    t <- seq(1,n)
    
    freq_p <- t %*% t(1/s.)
    A <- cos(2*pi*freq_p)
    B <- sin(2*pi*freq_p)
    
    factors <- B[,c(1,1)%*%t(seq(1:dim(B)[2]))]
    factors[,seq(1,dim(factors)[2],2)] <- A
    
    Lin.Mod <- lm(ts ~ factors)
    res <- list(Lin.Mod = Lin.Mod, factors = factors)
  }
  return(res)
}

#Scaled Periodogram
scaled.periodogram <- function(ts){
  n <- length(ts)
  if(n %% 2 == 1){
    k <- (n-1)/2
    f_freq <- seq(1,k)/n
    Is <- NULL
    coef <- Fourier(ts)
    for(i in 1:k)
      Is <- c(Is,(n/2)*(coef[2*i-1]^2 + coef[2*i]^2))
  }
  else{
    k <- n/2
    f_freq <- seq(1,k)/n
    Is <- NULL
    coef <- Fourier(ts)
    for(i in 1:(k-1)){
      Is <- c(Is,(n/2)*(coef[2*i-1]^2 + coef[2*i]^2))
    }
    Is <- c(Is, coef[k]^2)
  }
  res <- list(Frequencies = f_freq, Magnitudes = Is)
  return(res)
}

#Periodogram
my.periodogram <- function(ts){
  n <- length(ts)
  t <- seq(1,n)
  if(n %% 2 == 1){
    k <- (n-1)/2
    f_freq <- seq(1,k)/n
    ds <- NULL
    for(j in 1:k)
      ds <- c(ds, (1/n)*( as.numeric(t(ts)%*%cos(2*pi*t*f_freq[j])) ))
  }
  else{
    k <- n/2
    f_freq <- seq(1,k)/n
    ds <- NULL
    for(j in 1:(k-1)){
      ds <- c(ds, (1/n)*as.numeric(t(ts)%*%cos(2*pi*t*f_freq[j])) )
    }
    ds <- c(ds, (1/n)*as.numeric(t(ts)%*%cos(2*pi*t*f_freq[j])) )
  }
  res <- list(Frequencies = f_freq, Magnitudes = abs(ds))
  return(res)
} 

#Auto Covariance Matrix
Auto.Cov <- function(vts, h = 1){
  if(h >=0){
    X <- vts
    meanX <- apply(X, 2, mean)
    n <- dim(X)[1]
    d <- dim(X)[2]
    Xh <- X[-(1:h),]
    Xhi <- head(X,-h)
    for (i in 1:(n-h)) {
      cov <- 1/n * (Xh[i,])%*%t(Xhi[i,])
    }
  }
  else{
    h <- -h
    X <- vts
    meanX <- apply(X, 2, mean)
    n <- dim(X)[1]
    d <- dim(X)[2]
    Xh <- X[-(1:h),]
    Xhi <- head(X,-h)
    for (i in 1:(n-h)) {
      cov <- 1/n * (Xh[i,])%*%t(Xhi[i,])
    }
    cov <- t(cov)
  }
  return(cov)
}

###Test of models

##Choosing models:
#Seasonal ARMA test. (OLD)
seas.ARMA <- function(ts, p = 7, max.order = 5, test = AIC, noMA = FALSE){
  m <- max.order
  im <- m
  if (noMA == TRUE) {im <- 0}
  test.s <- NULL
  order <- NULL
  for (i in 0:m) {
    tmp.mod <- rep(0,im+1)
    tmp.mod.B <- rep(0,im+1)
    tryCatch({ 
      for (j in 0:im) {
        mod <- Arima(ts, order = c(0,0,0), seasonal = list(order = c(i,0,j), period = p), include.mean = FALSE, method = "ML")
        tmp.mod[j+1] <- AIC(mod)
        tmp.mod.B[j+1] <- BIC(mod)
        order <- rbind(order,c(i,0,j))
        cat("j =", j)
      }
    },error = function(e){})
    tmp.stat <- cbind(tmp.mod,tmp.mod.B)
    test.s <- rbind(test.s, tmp.stat)
    cat("i =", i)
  }
  bst.AIC <- which(test.s[,1] == min(test.s[test.s[,1]>0,1]))
  bst.BIC <- which(test.s[,2] == min(test.s[test.s[,2]>0,2]))
  colnames(test.s) = c("AIC","BIC")
  res <- list(Stats = test.s, model.n = rbind( c(test.s[bst.AIC,1], bst.AIC), c(test.s[bst.BIC,2], bst.BIC) ), order = rbind(order[bst.AIC,], order[bst.BIC,] ))
  return(res)
}

#Seasonal Arma Estimation for best AIC and BIC
best.SARMA <- function(ts, season = 7, max.order = 5, diffOrder = 1, sdiffOrder = 1)
  #Estimates different order of SARIMA models,
  #and it finds the best model based on AIC and BIC.
  #It also return all the tested orders.
  {
  library(astsa)
  m <- max.order
  
  NumAr <- m
  d <- diffOrder
  NumMa <- m
  NumSar <- m
  sd <- sdiffOrder
  NumSma <- m
  s <- season
  
  order <- NULL
  AICs <- NULL
  BICs <- NULL
  for (i in 0:NumAr) {
    tryCatch({ 
      for (j in 0:d) {
        for (k in 0:NumMa) {
          for (l in 0:NumSar) {
            for (o in 0:sd) {
              for (p in 0:NumSma) {
                order <- rbind(order, c(i,j,k,l,o,p,s))
                mod <- sarima(ts,i,j,k,l,o,p,s, no.constant = TRUE)
                AICs <- c(AICs,mod$AIC)
                BICs <- c(BICs,mod$BIC)
              }
            }
          }
        }
      }
    },error = function(e){})
  }
  bestAIC <- which( AICs == min( AICs[which(AICs > 0)]  ) )
  bestBIC <- which( BICs == min( BICs[which(BICs > 0)]  ) )
  res <- list(Orders = order, bestBIC)
  return(res)
}

##Remove insignificant parameters from arima and return new values
#Support function
make.AR.model <- function(ts, pointer) {
  n <- length(ts)
  m <- length(pointer)
  max <- max(pointer)
  #factors <- matrix(NULL, n-max, m)
  for (i in 1:m) {
    tmp <- rbind(ts[-c(1:pointer[i])])
  }
  factors <- matrix(c(tmp), n-max, m, byrow = TRUE)
  return(factors)
}
#The actual function
remove.insig <- function(x) {
  ts <- x$x
  p_values <- (1-pnorm(abs(x$coef)/sqrt(diag(x$var.coef))))*2
  p_values_05 <- which(p_values <= 0.05) 
  sig_coef <- x$coef[p_values_05]
  SAR_design <- make.AR.model(ts, p_values_05)
  y <- SAR_design %*% sig_coef 
  return(y)
}





#Moving Average Smoother
ma.smooth <- function(mts, level = 1, period = 12){
  #mts <- ts(mts)
  ma_filter <- rep(level,period)/period
  tsf <- filter(mts, sides = 2, filter = ma_filter)
  plot(mts, type = "l")
  lines(tsf, lwd = 2, col = 2)
  return(tsf)
}


#####Quality of Life Functions#####

#Default line plot
lplot <- function(...) plot(..., type="l")

#Col of data.frame to vector
dfCol.to.vector <- function(x) {
  as.numeric(x)
}

##QR factorisering og fra den
inner = function(v1, v2) {
  stopifnot(length(v1) == length(v2))
  return(sum(v1 * v2))
}
my_QR <- function(A) {
  dim_A <- dim(A)
  Q <- matrix(0,nrow = dim_A[1], ncol = dim_A[2])
  U <- matrix(0,nrow = dim_A[1], ncol = dim_A[2])
  U[,1] <- A[,1]
  Q[,1] <- A[,1]/sqrt(inner(U[,1],U[,1]))
  for (k in 2:dim_A[2]) {
    U[,k] <- A[,k]
    for (j in 1:(k-1)) {
      U[,k] <- U[,k] - ((inner(U[, j], U[, k])/ inner(U[,j], U[,j])) * (U[,j]))
    }
    Q[,k] <- U[,k]/sqrt(inner(U[,k],U[,k]))
  }
  R <- t(Q) %*% A
  res <- list(Q=Q,R=R,U=U)
  return(res)
}

inv_QR <- function(A) {
  QR <- my_QR(A)
  tQ <- t(QR$Q)
  R <- QR$R
  inv <- backsolve(R,tQ)
  return(inv)
}


polymult <- function(x,a,b)
  # a(B) = a_0 + a_1*B + a_2*BB2 + ... + a_p*B^p
  # b(B) = b_0 + b_1*B + b_2*BB2 + ... + a_q*B^q
  # where wlog p >= q
  # Result is a(B)b(B) x_t
{
  if (length(b) == 1){polycoeff <- b[1]*array(a)}
  else
  {alpha <- array(a)
  p <- dim(alpha) - 1
  beta <- array(b)
  q <- dim(beta) - 1
  # Assumed: p >= q
  polycoeff <- array(rep(0,p+q + 1))
  for (k in 0:(p+q))
    for (j in max(0,k-q):min(k,p))
    {polycoeff[k+1] <- polycoeff[k+1] + alpha[j+1]*beta[k-j+1]}
  }
  return(stats::filter(x,c(polycoeff),sides=1, method="convolution"))
}


polyinvers <- function (x, a, maxlag = 100)
  # a(B) = a_0 + a_1*B + a_2*BB2 + ... + a_p*B^p
  # Result is (1/a(B)) x_t
{
  phi <- array(a)
  p <- dim(phi) - 1 # assumed >= 1
  polycoeff <- array(rep(0,maxlag))
  polycoeff[1] <- 1/phi[1]
  for (k in 1:(maxlag - 1))
    for (j in (1:min(p,k)))
    {polycoeff[k+1] <- polycoeff[k+1] - phi[j+1]*polycoeff[k-j+1]/phi[1]}
  return(stats::filter(x,c(polycoeff),sides=1, method="convolution"))
}

##Plot Impulse response 
plot.irf2 <- function(mod, order = 2)
#Plots confidence intervals and IRF for class()="varirf"
  {
  conf_upper <- irf_SVARB$Upper$C
  conf_lower <- irf_SVARB$Lower$C
  m <- length(irf_SVARB$Lower$C)
  
  #Plot
  plot(1:m,irf_SVARB$irf$C, ylim = c(min(irf_SVARB$irf$C), max(irf_SVARB$irf$C)), type = "l")
  lines(1:m, conf_upper, col = "red")
  lines(1:m, conf_lower, col = "red")
}


