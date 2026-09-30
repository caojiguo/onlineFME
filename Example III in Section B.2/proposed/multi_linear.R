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
  betaeval = sin(0.5*pi*tobs)
  h   = (domain[2]-domain[1])/(T-1)
  cef = c(1, rep(c(4,2), (T-3)/2), 4, 1)
  y0  = rep(NA,n)
  y0  = (h/3)*x%*%diag(cef)%*%betaeval
  return(list(X = x, y0 = y0, coef = betaeval))
}


FLMR.data.generator.bsplines2 = function(n,nknots,norder,T,domain,aind,a,b,df)
{
  knots    = seq(domain[1],domain[2], length.out = nknots)
  nbasis   = nknots + norder - 2
  basis    = create.bspline.basis(knots,nbasis,norder)
  tobs = seq(domain[1],domain[2],length.out = T)
  basismat = eval.basis(tobs, basis)
  x = a_fun(n,nbasis,aind,a,b,df) %*% t(basismat)
  betaeval = cos(0.5*pi*tobs)
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
  norder   = d+1 # order = degree +1
  nknots   = K+2 # the number of all knots = the number of inner knots + 2
  knots    = seq(domain[1],domain[2], length.out = nknots)
  nbasis   = nknots + norder - 2 # the number of basis = number of inner knots + order
  basis    = create.bspline.basis(knots,nbasis,norder)
  T   = dim(X)[2]
  tobs = seq(domain[1], domain[2], length.out = T)
  basismat = eval.basis(tobs, basis) # T x (K+d+1)
  V = eval.penalty(basis,int2Lfd(2))
  
  h   = (domain[2]-domain[1])/(T-1)
  cef = c(1, rep(c(4,2), (T-3)/2), 4, 1)
  u   = (h/3)*X%*%diag(cef)%*%basismat # n x nbasis
  return(list(N=u,V = V, basismat = basismat, tobs = tobs))
}




# ==== dynanmic lambda ==== 
SL <- function(i, gamma, y, z, theta_i, theta_bar, D, n){
  
  z <- matrix(z, ncol=1)
  scale <- ((n+1-i)/n)^0.5 * n^-1
  lambda_list <- 10^(seq(0, 1, length.out=10)) * scale
  
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
nknots1 = 36
nknots2 = 36
T = 101
domain = c(0, 1)
aind1 = 1
aind2 = 1
a1 = 0
b1 = 6
a2 = 0
b2 = 6
df = NA
nKnot = 1
degree = 3
nbasis = nKnot + degree + 1 
n_drop <- 4000

tobs = seq(0, 1, length.out = T)
beta_1 = sin(0.5*pi*tobs)
beta_2 = cos(0.5*pi*tobs)
D1 = eval.penalty(create.bspline.basis(seq(0, 1, length.out = nKnot+2), nbasis, degree+1), int2Lfd(2))
D2 = eval.penalty(create.bspline.basis(seq(0, 1, length.out = nKnot+2), nbasis, degree+1), int2Lfd(2))
D = matrix(0, nrow = 2*nbasis, ncol = 2*nbasis)
D[1:nbasis, 1:nbasis] = D1
D[(nbasis+1):(2*nbasis), (nbasis+1):(2*nbasis)] = D2
rm(D1,D2)


# ==== data generating ==== 
#n = 10000
data1 = FLMR.data.generator.bsplines1(n,nknots1,norder,T,domain,aind1,a1,b1,df)
data2 = FLMR.data.generator.bsplines2(n,nknots2,norder,T,domain,aind2,a2,b2,df)
Z1_list = compute.NV(data1$X, nKnot, degree, domain)
Z2_list = compute.NV(data2$X, nKnot, degree, domain)
Z1 = Z1_list$N
Z2 = Z2_list$N
basis1 <- Z1_list$basismat
basis2 <- Z2_list$basismat
basis <- cbind(basis1,basis2)
y1 = data1$y0 
y2 = data2$y0
Y = y1 + y2 + rnorm(n, 0, 0.5)


# ==== main algorithm ==== 
start_time <- Sys.time()


theta_i <- matrix(0, 2*nbasis, 1)
theta_bar <- matrix(0, 2*nbasis, 1)
A_i <- matrix(0, 2*nbasis, 2*nbasis)
b_i <- matrix(0, 2*nbasis, 1)
sum_i2 <- 0

for(i in 1:n){
  gamma <- (4/5) * i^(-0.2)
  y <- Y[i]
  z <- c(Z1[i,],Z2[i,])
  eta <- as.numeric(z %*% theta_i)
 
  
  loss <- as.numeric(y - z %*% theta_i)
  lambda <- SL(i, gamma, y, z, theta_i, theta_bar, D, n)#0.05*(((n+1-i)/n)^0.5 * n^-1)#
  print(lambda/(((n+1-i)/n)^0.5 * n^-1))
  grad <- -2 * loss * z + 2 * lambda * D %*% theta_i
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

V_mat <- (m^2)^(-1) * (A_i - theta_bar %*% t(b_i) - b_i %*% t(theta_bar) + sum_i2 * theta_bar %*% t(theta_bar))
varbeta_sgd1 <- (basis1 %*% V_mat[1:nbasis,1:nbasis] %*% t(basis1)) / m
varbeta_sgd2 <- (basis2 %*% V_mat[(nbasis+1):(2*nbasis),(nbasis+1):(2*nbasis)] %*% t(basis2)) / m

diag_varbeta1 <- diag(varbeta_sgd1)
diag_varbeta2 <- diag(varbeta_sgd2)

beta_sgd1 <- basis1 %*% theta_bar[1:nbasis,]
beta_sgd2 <- basis2 %*% theta_bar[(nbasis+1):(2*nbasis),]

end_time <- Sys.time()
runtime <- as.numeric(difftime(end_time, start_time, units = "secs"))

# ==== CI ====
Q_pt <- 6.734
lower_pt1 <- beta_sgd1 - sqrt(diag_varbeta1)*Q_pt
upper_pt1 <- beta_sgd1 + sqrt(diag_varbeta1)*Q_pt

lower_pt2 <- beta_sgd2 - sqrt(diag_varbeta2)*Q_pt
upper_pt2 <- beta_sgd2 + sqrt(diag_varbeta2)*Q_pt

quant_sim <- sqrt((nKnot + 4)*72.50)
lower_sim1 <- beta_sgd1 - sqrt(diag_varbeta1)*quant_sim
upper_sim1 <- beta_sgd1 + sqrt(diag_varbeta1)*quant_sim

lower_sim2 <- beta_sgd2 - sqrt(diag_varbeta2)*quant_sim
upper_sim2 <- beta_sgd2 + sqrt(diag_varbeta2)*quant_sim


# ==== plot ====
par(mfrow = c(1, 2))
plot(tobs, beta_1, type="l", col="blue", ylim=range(c(lower_sim1, upper_sim1)), 
     xlab="t", ylab="beta")
lines(tobs, beta_sgd1, col="red")
matlines(tobs, lower_pt1, col="grey", lty=2)
matlines(tobs, upper_pt1, col="grey", lty=2)
matlines(tobs, lower_sim1, col="orange", lty=2)
matlines(tobs, upper_sim1, col="orange", lty=2)





plot(tobs, beta_2, type="l", col="blue", ylim=range(c(lower_sim2, upper_sim2)), 
     xlab="t", ylab="beta")
lines(tobs, beta_sgd2, col="red")
matlines(tobs, lower_pt2, col="grey", lty=2)
matlines(tobs, upper_pt2, col="grey", lty=2)
matlines(tobs, lower_sim2, col="orange", lty=2)
matlines(tobs, upper_sim2, col="orange", lty=2)





# ==== IMSE  ====
IMSE1 <- compute_IMSE(beta_sgd1, beta_1, domain, T)
IMSE2 <- compute_IMSE(beta_sgd2, beta_2, domain, T)




result <- list(
  IMSE1 = IMSE1,
  IMSE2 = IMSE2,
  cover_pt1 = as.numeric(beta_1 >= lower_pt1 & beta_1 <= upper_pt1),  
  cover_sim1 = as.numeric(all(beta_1 >= lower_sim1 & beta_1 <= upper_sim1)),  
  cover_pt2 = as.numeric(beta_2 >= lower_pt2 & beta_2 <= upper_pt2), 
  cover_sim2 = as.numeric(all(beta_2 >= lower_sim2 & beta_2 <= upper_sim2)), 
  length_pt1 = upper_pt1 - lower_pt1, 
  length_pt2 = upper_pt2 - lower_pt2, 
  length_sim1 = upper_sim1 - lower_sim1,
  length_sim2 = upper_sim2 - lower_sim2,
  beta1 = beta_sgd1,
  beta2 = beta_sgd2,
  theta_bar = theta_bar,
  V_mat = V_mat,
  runtime = runtime
)

TargetFolder <- "./result/"
dir.create(TargetFolder, showWarnings = FALSE)
fname <- sprintf("res_n%d_J%d_rep%03d.rds", n, nbasis, iloop)
saveRDS(result, file = file.path(TargetFolder, fname))
cat("Replication", iloop, "done.", "\n")