
vcov.fn <- function(a, b, c, d, e) {
  col <- c(1, a, 0, 0, 0,
              1, b, c, d,
                 1, 0, 0, 
                    1, e,
                       1)
  mat <- xpnd(col)
  vnames <- c("x1", "x2", "x3", "z1", "z2")
  colnames(mat) <- vnames; rownames(mat) <- vnames
  return(mat)
}

rmvn.chol <- function(n, vcov.mat) {
  k <- ncol(vcov.mat)
  Q <- chol(vcov.mat)
  Z <- matrix(rnorm(k*n), nrow=n, ncol=k)
  return(Z %*% Q)
}


Wt.mat = function(G) solve( (1/n) * (t(G)%*%G) )
A = function(WW) (t(X) %*% Z) %*% WW %*% (t(Z) %*% X )
B = function(WW) (t(X) %*% Z) %*% WW %*% (t(Z) %*% y )

n <- 100000; total.reps <- 1000
vcov <- vcov.fn(0.1, 0.5, 0.4, 0.4, 0.1)
res <- rep(0,total.reps)
for (i in 1:total.reps) {
  data <- rmvn.chol(n, vcov)
  beta <- c(1, 2, -4, 1)
  X <- cbind(1, data[,1])
  Z <- cbind(data[,4], data[,5])
  eta <- rnorm(n)
  y <- cbind(1, data[,c(1,2,3)]) %*% beta + rnorm(n)
  W <- Wt.mat(Z)
  beta.hat <- solve(A(W)) %*% B(W)
  res[i] <- beta.hat[1]
}

betaHat = function(WW) {
  AA = A(WW)
  BB = B(WW)
  b = solve(AA, BB)
  return(b)
}
    
eHat = function(betaHat) y - X %*% betaHat # the residual

### 1st stage: W = (X'X)^-1
W1 = W(Z)
betaHat1 = betaHat(W1)
    
### 2nd stage: W = (G G' - g g')^-1
e = eHat(betaHat1)
G2 = Z * matrix(rep(e, 2), ncol = 2)
G2Bar = colMeans(G2)
W2 = solve( (t(G2)%*%G2)/n   -  G2Bar %*% t(G2Bar)  )
betaHat2 = betaHat(W2)

n = 10
k = 1
l = 2
d =  matrix(rnorm( (1+k+l)*n ), ncol = 1 + k + l)
colnames(d) = c("y", "x1", "z1", "z2", "z3", "z4", "z5")

d = d + 3

d[, 2:(1+k)]         = 2 + rowSums( d[, (1+k+1):(1+k+l)] )  + rnorm(n, 0, 0.2)
d[, 1] = 2 +  d[, (1 + 1):(1+k  )]   + rnorm(n, 0, 0.2)

fs = d$y ~ d$x1 
fi =     ~ d$z1 + d$z2 + d$z3 + d$z4 + d$z5

###### result from package gmm #################

### good gmm #################
source("gmm_linear.R")
gmm.linear(list(structural = fs, instruments = fi), data = d)


### big code ###
n <- 10000
x <- rnorm(n)
z1 <- 3 + 0.50*x + rnorm(n, 0, 0.1)
z2 <- 1 + 0.25*x + rnorm(n, 0, 0.1)
y <- 2 + 3*x + 4*z1 + 5*z2 + rnorm(n)

X <- cbind(1, x)
Z <- cbind(1, z1, z2)
W = function(G) solve( (1/n) * (t(G)%*%G) )
A = function(WW) ((1/n) * t(X) %*% Z) %*% WW %*% ((1/n) * t(Z) %*% X )
B = function(WW) ((1/n) * t(X) %*% Z) %*% WW %*% ((1/n) * t(Z) %*% y )
betaHat = function(WW) {
  AA = A(WW)
  BB = B(WW)
  b = solve(AA, BB)
  return(b)
}
    
eHat = function(betaHat) y - X %*% betaHat # the residual

### 1st stage: W = (X'X)^-1
W1 = W(Z)
betaHat1 = betaHat(W1)
    
### 2nd stage: W = (G G' - g g')^-1
e = eHat(betaHat1)
G2 = Z * matrix(rep(e, 3), ncol = 3)
G2Bar = colMeans(G2)
W2 = solve( (t(G2)%*%G2)/n   -  G2Bar %*% t(G2Bar)  )
betaHat2 = betaHat(W2)
    
### variance-covariance matrix
Q = crossprod(Z, X)/n
V = solve( t(Q) %*%  W2 %*% Q )
s = sqrt(diag(V)/n)
    
summ = cbind(betaHat2, s, betaHat2/s, 2*(1-pnorm(abs(betaHat2/s))) )
colnames(summ) = c("coefficient", "s.d.", "Z", "p-value")
print(summ)

data = d
formula = list(structural = fs, instruments = fi)
data1 = model.frame(formula$structural, data)
data2 = model.frame(formula$instruments, data)

##### organize the data ############

# y: endogenous variable, (n * 1) vector
# Z: explanatory variable, (n * k) matrix
# X: instruments, (n*l) matrix

y = model.response(data1, type = "numeric")
    
Z = model.matrix(formula$structural, data = data1)
X = model.matrix(formula$instrument, data = data2)
    
l = ncol(X)
k = ncol(Z)
    
if (k >l ) stop( "# of regressors must <= # of instruments" )
    
############### start gmm computation #############
    
## prepare the functions
    

    
## J stat and test ###
J = n * t(G2Bar) %*% W2 %*% G2Bar
cat("\n J-stat = ", J, 
    ". DF = ", l-k, 
    ". p-value = ", 1-pchisq(J, df = l-k), ".\n")
