# this script runs two-step gmm for linear model
# it is the very prototype function. 
# Feel free to use, and/or modify it for your own purpose.

# Note: it is NOT designed for nonliner model, though the point estimation will be correct 
# if the nonlinearity can be represented by a series expansion, but inference 
# is not applicable

# Zhentao Shi, 6/1/2011
#####################################

gmm.linear <- function(formula, data){
# syntax: gmm(formula, data, method = "twoStep")
# arguments: 
#   formula: f = list(structural, instruments)
#   data   : a matrix or data frame

# the use of the arguments are the same as the function lm()
# please refer to the companion "test.R" as an example.

######################################

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
    
    W = function(G) solve( (1/n) * (t(G)%*%G) )
    A = function(WW) ((1/n) * t(Z) %*% X) %*% WW %*% ((1/n) * t(X) %*% Z )
    B = function(WW) ((1/n) * t(Z) %*% X) %*% WW %*% ((1/n) * t(X) %*% y )
    betaHat = function(WW) {
      AA = A(WW)
      BB = B(WW)
      b = solve(AA, BB)
      return(b)
      }
    
    eHat = function(betaHat) y - Z %*% betaHat # the residual
    
    ### 1st stage: W = (X'X)^-1
    W1 = W(X)
    betaHat1 = betaHat(W1)
    
    ### 2nd stage: W = (G G' - g g')^-1
    e = eHat(betaHat1)
    G2 = X * matrix(rep(e, l), ncol = l)
    G2Bar = colMeans(G2)
    W2 = solve( (t(G2)%*%G2)/n   -  G2Bar %*% t(G2Bar)  )
    betaHat2 = betaHat(W2)
    
    ### variance-covariance matrix
    Q = crossprod(X, Z)/n
    V = solve( t(Q) %*%  W2 %*% Q )
    s = sqrt(diag(V)/n)
    
    summ = cbind(betaHat2, s, betaHat2/s, 2*(1-pnorm(abs(betaHat2/s))) )
    colnames(summ) = c("coefficient", "s.d.", "Z", "p-value")
    print(summ)
    
    ## J stat and test ###
    J = n * t(G2Bar) %*% W2 %*% G2Bar
    cat("\n J-stat = ", J, 
        ". DF = ", l-k, 
        ". p-value = ", 1-pchisq(J, df = l-k), ".\n")
        
}