# OA-CLGARCH-mode
# R codes to fit OA-CLGARCH mode
library(glmnet)
library(gregmisc)
# Function mcd_TD_lasso estimates Cholesky factors (T, D) in MCD by given order.
# data_x: data matrix;
# order: a specified order for predictors;
# return: Cholesky factor matrix T in given order;
#         the diagonal of matrix D in given order;
#         the residuals e in original order.
mcd_TD_lasso = function(data_x, order)
{
    p = ncol(data_x)
    size = nrow(data_x)
    # variable e records residuals for each p regressions
    e = matrix(data = NA, nrow = size, ncol = p)
    T = diag(rep(1,p))
    d = rep(1,p)                   # diagonal of matrix D
    new_x = data_x[, order]	       # re-order the columns of data by specified order
    e[,1] = new_x[,1]
    d[1] = var(new_x[,1])
    # fit p regressions
    for(i in c(p:2))
    {
        y = new_x[,i]
        x = new_x[,(1:(i-1))]
        # solve regression estimates
        if(i == 2)
        {
            beta = solve(t(x) %*% x) %*% t(x) %*% y	
        }
        # solve Lasso estimates
        else
        {
            lasso.fit = glmnet(x, y, family = 'gaussian', alpha = 1)
            lambda = cv.glmnet(x, y, nfolds=5, alpha = 1)$lambda.min
            beta = as.vector(coef(lasso.fit, s = lambda))[-1]
        }
        T[i,(1:(i-1))] = -beta
        e[,i] = y - x %*% beta
        d[i] = var(e[,i])
    }
    e_orig = matrix(data = NA, nrow = nrow(data_x), ncol = p)
    e_orig[, order] = e     # put residuals in the original order
    return(list(T = T, d = d, e = e_orig))
}

# Function diag_D estimates the diagonal of matrix D by the log-GARCH model.
# e: residuals obtained from the MCD
# d.sq.mw: residual variances obtained from the moving block
# return: the estimates for diagonal of matrix D
diag_D = function(e, d.sq.mw)
{
    n = nrow(e)
    p = ncol(e)
    s.sq.lag1 = log(d.sq.mw[1:(n-1),])
    s.sq.lag0 = log(d.sq.mw[2:n,])
    e.lag1 = log(e[1:(n-1),]^2+10^-6)
    I1 = I(e.lag1>=0)
    I2 = I(e.lag1<0)
    lambda = NULL
    log.d.sq = NULL
    for (j in 1:p)
    {
        # obtain the initial values from least squares
        start = lm(s.sq.lag0[,j]~I(I1[,j]*e.lag1[,j])+I(I2[,j]*e.lag1[,j])
                                    +s.sq.lag1[,j])
        s.coef = start$coef
        # fit log-GARCH model by quasi-maximum likelihood method
        fit = estimLogGARCH(s.coef[1], s.coef[2], s.coef[3], s.coef[4], e[,j], 10000)
        lambda = rbind(lambda, fit$coef)
        log.d.sq = cbind(log.d.sq, fit$log.sig2)
    }
    return(exp(log.d.sq))								
}

# Function VarAsymp is used in function estimLogGARCH.
VarAsymp = function(omega,alpha.plus,alpha.moins,beta,eps,sig2init,petit,r0)
{
    n = length(eps)
    derlogsigma2 = matrix(0, nrow=4, ncol=n)
    log.sig2 = rep(0,n)
    log.sig2[1] = log(sig2init)		
    derlogsigma2[1:4,1] = 0
    for(t in 2:n)
    {
        vec = c(1,0,0,0)
        log.sig2[t] = omega + beta*log.sig2[t-1]
        if(eps[t-1] > petit)				
        {
            log.sig2[t] = log.sig2[t] + alpha.plus * log(eps[t-1]^2)
            vec[2] = log(eps[t-1]^2)
        }
        if(eps[t-1] < -petit)				
        {
            log.sig2[t] = log.sig2[t] + alpha.moins * log(eps[t-1]^2)
            vec[3] = log(eps[t-1]^2)
        }
        vec[4] = log.sig2[t-1]
        derlogsigma2[1:4,t] = vec+beta*derlogsigma2[1:4,(t-1)]
    }
    sig2 = exp(log.sig2[(r0+1):n])
    eta = eps[(r0+1):n]/sqrt(sig2)
    eta = eta/sd(eta)
    J = derlogsigma2[1:4,(r0+1):n] %*% t(derlogsigma2[1:4,(r0+1):n])/(n-r0)
    kappa4 = mean(eta^4)
    {if(kappa(J)<1/petit) inv = solve(J) else inv = matrix(0,nrow=4,ncol=4)}
    var = (kappa4-1)*inv
    return(list(var = var, residus = eta, log.sig2 = log.sig2))
}

# Define the function which needs to be minimized in estimLogGARCH function.
objf.loggarch = function(vartheta, eps, n, sig2init, eps.plus, eps.moins, r0)
{
    omega = vartheta[1]
    alpha.plus = vartheta[2]
    alpha.moins = vartheta[3]
    beta = vartheta[4]
    log.sig2 = rep(0,n)
    log.sig2[1] = log(sig2init)
    for(t in 2:n){
    log.sig2[t] = omega + beta*log.sig2[t-1] + alpha.plus*eps.plus[t-1]
                                        + alpha.moins*eps.moins[t-1]}
    sig2 = exp(log.sig2[(r0+1):n])
    qml = mean(eps[(r0+1):n]^2/sig2 + log(sig2))
    qml
}

# Function estimLogGARCH estimates the parameters in the log-GARCh model.
estimLogGARCH = function(omega, alpha.plus, alpha.moins, beta,eps, factor,
                                petit = sqrt(.Machine$double.eps), r0=10)
{
    petit = factor*petit
    n = length(eps)
    eps.plus = rep(0,n)
    eps.moins = rep(0,n)
    for(i in 1:n)
    {
        {if(eps[i]>=0) eps.plus[i] = log(max(eps[i],petit)^2)}
        {if(eps[i]<0) eps.moins[i] = log(max(-eps[i],petit)^2)}
    }
    valinit = c(omega, alpha.plus, alpha.moins, beta)
    sig2init = var(eps[1:min(n,5)])						
    res = nlminb(valinit, objf.loggarch, lower = c(-Inf, -Inf, -Inf, -1+petit),
                upper = c(Inf,Inf,Inf,1-petit), eps=eps, n=n, sig2init=sig2init,
                eps.plus=eps.plus, eps.moins=eps.moins, r0=r0)		
	
    # record the estimates of the parameters in log-GARCH model
    omega = res$par[1]
    alpha.plus = res$par[2]
    alpha.moins = res$par[3]
    beta = res$par[4]
    var = VarAsymp(omega,alpha.plus,alpha.moins,beta,eps,sig2init,petit,r0)
    return(list(coef = res$par, residus = var$residus, var=var$var,
                            log.sig2 = var$log.sig2, loglik = res$objective))
}

# Function OA_CLGARCH provides the proposed covariance estimate.
# x: n by p data matrix
# q: moving block size
# M: number of permutations
# return: the estimated covariance matrices at each time point
OA_CLGARCH = function(x, q, M)
{
    n = nrow(x)
    p = ncol(x)
    # construct moving blocks
    s.square = NULL                 # record sample variance of predictors
    d.sq.mw = matrix(NA, n, p-1)
    for (i in 1:n)
    {
        if ((i-(q-1)/2>0) & (i+(q-1)/2<=n))
        {x.block = x[(i-(q-1)/2):(i+(q-1)/2),]}
        else if (i-(q-1)/2<=0)
        {x.block = x[1:(i+(q-1)/2),]}
        else if (i+(q-1)/2>n)
        {x.block = x[(i-(q-1)/2):n,]}
        s.square = rbind(s.square, diag(cov(x.block)))
        for (j in 2:p)
        {
           mod.mw = lm(x.block[,j] ~ x.block[,1:(j-1)]-1)
           d.sq.mw[i,j-1] = var(mod.mw$residuals)
        }
    }
    d.sq.mw = cbind(s.square[,1], d.sq.mw)
    # generate M orders of variables
    perm = rep(0,p)
    if(p < 10){perm = permutations(n = p, r = p)}
    # if p is large, randomly choose M orders of variables
    if(p >= 10)								
    {
        perm = matrix(data=NA, nrow=M, ncol=p)
        for(i in 1:M)
        {perm[i,] = sample(p, p)}
    }
    id = c(1:nrow(perm))
    if(nrow(perm) > M){id = sample(nrow(perm), M)}
    # record estimates T and D for each order of variables		
    T_l = array(data=0, dim=c(p,p,length(id)))		
    d_l = array(data=0, dim=c(n,p,length(id)))	
    # estimate matrices T and D for each order of variables
    for(i in 1:length(id))
    {
        order = perm[(id[i]),]
        # construct permutation matrix
        P_pi = matrix(data=0, nrow=p, ncol=p)		
        for(kk in 1:p)
        {P_pi[order[kk], kk] = 1}
        result_lasso = mcd_TD_lasso(x, order)				
        T_l[,,i] = P_pi %*% result_lasso$T %*% t(P_pi)
        d_l[,,i] = diag_D(result_lasso$e, d.sq.mw)
    }   	
    T_l1 = apply(T_l, 1:2, mean)	   # take average for estimates of T
    d_l1 = apply(d_l, 1:2, mean)
    # compute the covariance matrices for each time point
    S.mean_L = array(NA, dim=c(p,p,n))
    for (i in 1:n)
    {S.mean_L[,,i] = solve(T_l1) %*% diag(d_l1[i,]) %*% solve(t(T_l1))}
    return(S.mean_L)
}
