rm(list=ls())

# ==== commander parameter ====
args <- commandArgs(TRUE)
iloop <- as.integer(args[1])
n     <- as.integer(args[2])  

set.seed(iloop+n)
cat("Replication:", iloop, "\n")

# ==== packages ====
library(fda)

# ==== data generation function ====
a_fun = function(n,nbasis, ii,a,b,df)
{
  if(ii == 1) {a = matrix(rnorm(nbasis*n,a,b),n,nbasis)}
  else if(ii == 2) {a = matrix(rmt(nbasis*n,mean = rep(a,1),b,df=df),n,nbasis)}
  else if(ii == 3) {a = matrix(runif(nbasis*n,a,b),n,nbasis)}
  else{print("model does not exit")}
  return(a)
}

FLMR.data.generator.bsplines1 = function(n,nknots,norder,T,domain,aind,a,b,df)
{
  knots    = seq(domain[1],domain[2], length.out = nknots)
  nbasis   = nknots + norder - 2
  basis    = create.bspline.basis(knots,nbasis,norder)
  tobs = seq(domain[1],domain[2],length.out = T)
  basismat = eval.basis(tobs, basis)
  x = a_fun(n,nbasis,aind,a,b,df) %*% t(basismat)
  betaeval = sin(0.5*pi*tobs)#exp(-tobs)
  h   = (domain[2]-domain[1])/(T-1)
  cef = c(1, rep(c(4,2), (T-3)/2), 4, 1)
  y0  = rep(NA,n)
  y0  = (h/3)*x%*%diag(cef)%*%betaeval
  return(list(X = x, y0 = y0, coef = betaeval))
}

# ==== B-spline ====
compute.NV = function(X, K, d, domain)
{
  n = dim(X)[1]
  norder   = d+1
  nknots   = K+2
  knots    = seq(domain[1],domain[2], length.out = nknots)
  nbasis   = nknots + norder - 2
  basis    = create.bspline.basis(knots,nbasis,norder)
  T   = dim(X)[2]
  tobs = seq(domain[1], domain[2], length.out = T)
  basismat = eval.basis(tobs, basis)
  V = eval.penalty(basis,int2Lfd(2))
  
  h   = (domain[2]-domain[1])/(T-1)
  cef = c(1, rep(c(4,2), (T-3)/2), 4, 1)
  u   = (h/3)*X%*%diag(cef)%*%basismat
  return(list(N=u,V = V, basismat = basismat, tobs = tobs))
}

# ==== dynanmic lambda ==== 
SL <- function(i, gamma, y, z, theta_i, theta_bar, D, n){
  
  z <- matrix(z, ncol=1)
  scale <- ((n+1-i)/n)^0.5 * n^-1
  lambda_list <- 10^(seq(0, 0.95, length.out=10))  * scale
  
  best_loss <- Inf
  best_lambda <- NA
  
  for(lambda in lambda_list){
    loss <- as.numeric(y - t(z) %*% theta_i)
    grad <- -2 * loss * z + 2 * lambda * D %*% theta_i
    theta_new <- theta_i - gamma * grad
    theta_bar_new <- ((i-1)/i)*theta_bar + (1/i)*theta_new
    loss_new <- as.numeric((y - t(z) %*% theta_bar_new)^2)
    
    if(loss_new < best_loss){
      best_loss <- loss_new
      best_lambda <- lambda
    }
  }
  
  return(best_lambda)
}

# ==== IMSE ==== 
compute_IMSE <- function(beta_hat, beta_true, domain, lent){
  h <- (domain[2] - domain[1]) / (lent - 1)
  cef <- c(1, rep(c(4,2), (lent-3)/2), 4, 1)
  (h/3) * sum(cef * (beta_hat - beta_true)^2)
}

# ==== parameter setting ==== 
norder = 4
nknots1 = 20
T = 101
domain = c(0, 1)

aind1 = 1
a1 = 0
b1 = 6
df = NA

nKnot = 1
degree = 3
nbasis = nKnot + degree + 1 
n_drop <- 4000

tobs = seq(0, 1, length.out = T)
beta_1 = sin(0.5*pi*tobs)#exp(-tobs)#

# penalty
D = eval.penalty(
  create.bspline.basis(seq(0, 1, length.out = nKnot+2), nbasis, degree+1),
  int2Lfd(2)
)

# ==== data generating ==== 
data1 = FLMR.data.generator.bsplines1(n,nknots1,norder,T,domain,aind1,a1,b1,df)

Z1_list = compute.NV(data1$X, nKnot, degree, domain)
Z1 = Z1_list$N
basis1 <- Z1_list$basismat

y1 = data1$y0 
Y = y1 + rnorm(n, 0, 0.5)

# ==== main algorithm ==== 
start_time <- Sys.time()

theta_i <- matrix(0, nbasis, 1)
theta_bar <- matrix(0, nbasis, 1)
A_i <- matrix(0, nbasis, nbasis)
b_i <- matrix(0, nbasis, 1)
sum_i2 <- 0

for(i in 1:n){
  gamma <- (4/5) * i^(-0.2)
  y <- Y[i]
  z <- Z1[i,]
  
  loss <- as.numeric(y - z %*% theta_i)
  lambda <- SL(i, gamma, y, z, theta_i, theta_bar, D, n)
  #print(lambda / (((n+1-i)/n)^0.5 * n^-1))
  
  grad <- -2 * loss * matrix(z, ncol=1) + 2 * lambda * D %*% theta_i
  theta_i <- theta_i - gamma * grad
  
  if(i > n_drop){
    k <- i - n_drop
    theta_bar <- if(k==1) theta_i else ((k-1)/k)*theta_bar + (1/k)*theta_i
    
    A_i <- A_i + k^2 * (theta_bar %*% t(theta_bar))
    b_i <- b_i + k^2 * theta_bar
    sum_i2 <- sum_i2 + k^2
  }
}

m <- n - n_drop

V_mat <- (m^2)^(-1) * (
  A_i - theta_bar %*% t(b_i) - b_i %*% t(theta_bar) + sum_i2 * theta_bar %*% t(theta_bar)
)

varbeta_sgd1 <- (basis1 %*% V_mat %*% t(basis1)) / m
diag_varbeta1 <- diag(varbeta_sgd1)

beta_sgd1 <- basis1 %*% theta_bar

end_time <- Sys.time()
runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))

# ==== confidence intervals ====
Q_pt <- 6.734
lower_pt1 <- beta_sgd1 - sqrt(diag_varbeta1)*Q_pt
upper_pt1 <- beta_sgd1 + sqrt(diag_varbeta1)*Q_pt

quant_sim <- sqrt((nKnot + 4)*72.50)
lower_sim1 <- beta_sgd1 - sqrt(diag_varbeta1)*quant_sim
upper_sim1 <- beta_sgd1 + sqrt(diag_varbeta1)*quant_sim

# ==== plot ====
plot(tobs, beta_1, type="l", col="blue",
     ylim=range(c(lower_sim1, upper_sim1)),
     xlab="t", ylab="beta")
lines(tobs, beta_sgd1, col="red")
matlines(tobs, lower_pt1, col="grey", lty=2)
matlines(tobs, upper_pt1, col="grey", lty=2)
matlines(tobs, lower_sim1, col="orange", lty=2)
matlines(tobs, upper_sim1, col="orange", lty=2)

# ==== IMSE & result ====
IMSE1 <- compute_IMSE(beta_sgd1, beta_1, domain, T)

result <- list(
  IMSE1 = IMSE1,
  cover_pt1 = as.numeric(beta_1 >= lower_pt1 & beta_1 <= upper_pt1),
  cover_sim1 = as.numeric(all(beta_1 >= lower_sim1 & beta_1 <= upper_sim1)),
  length_pt1 = upper_pt1 - lower_pt1,
  length_sim1 = upper_sim1 - lower_sim1,
  beta1 = beta_sgd1,
  theta_bar = theta_bar,
  V_mat = V_mat,
  runtime = runtime
)

TargetFolder <- "./result/"
dir.create(TargetFolder, showWarnings = FALSE)
fname <- sprintf("res_n%d_J%d_rep%03d.rds", n, nbasis, iloop)
saveRDS(result, file = file.path(TargetFolder, fname))
cat("Replication", iloop, "done.", "\n")