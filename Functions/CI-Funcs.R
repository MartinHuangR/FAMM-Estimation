sbwlsCI = function(X, Y, true_beta, 
                   LOOPS = 100, nOut, m = 32, largeOut = T,na.rm = F){
  n = nrow(X) 
  p = ncol(X)
  # Quantile Confidence Interval
  quantile.ci.acc = matrix(NA, nrow = LOOPS, ncol = p)
  quantile.length = matrix(NA, nrow = LOOPS, ncol = p)
  
  # BCa Confidence Interval
  bca.ci.acc = matrix(NA, nrow = LOOPS, ncol = p)
  bca.length = matrix(NA, nrow = LOOPS, ncol = p)
  for (i in 1:LOOPS){
    Y.full = Y + rnorm(n)
    n0 = sample(1:n, nOut, replace = F)
    if (largeOut == T){Y.full[n0] = runif(nOut, 1000, 2000)}
    if (largeOut == F){Y.full[n0] = runif(nOut, 30, 40)}
    Boot = 100
    rlm.full = rlm(Y.full ~ X, method='MM', maxit = 500)
    res.rank = rank(rlm.full$residuals)
    weights = rlm.full$w
    res.ind  <- rep(1,n);
    for(j1 in 1:n){
      if (res.rank[j1] <=n/8){res.ind[j1] <- res.ind[j1]+1};
      if (res.rank[j1] <=2*n/8){res.ind[j1] <- res.ind[j1]+1};
      if (res.rank[j1] <=3*n/8){res.ind[j1] <- res.ind[j1]+1};
      if (res.rank[j1] <=4*n/8){res.ind[j1] <- res.ind[j1]+1};
      if (res.rank[j1] <=5*n/8){res.ind[j1] <- res.ind[j1]+1};
      if (res.rank[j1] <=6*n/8){res.ind[j1] <- res.ind[j1]+1};
      if (res.rank[j1] <=7*n/8){res.ind[j1] <- res.ind[j1]+1};
    }
    X1 <- X[res.ind==1,]
    X2 <- X[res.ind==2,]
    X3 <- X[res.ind==3,]
    X4 <- X[res.ind==4,]
    X5 <- X[res.ind==5,]
    X6 <- X[res.ind==6,]
    X7 <- X[res.ind==7,]
    X8 <- X[res.ind==8,]
    Y1 <- Y.full[res.ind==1]
    Y2 <- Y.full[res.ind==2]
    Y3 <- Y.full[res.ind==3]
    Y4 <- Y.full[res.ind==4]
    Y5 <- Y.full[res.ind==5]
    Y6 <- Y.full[res.ind==6]
    Y7 <- Y.full[res.ind==7]
    Y8 <- Y.full[res.ind==8]
    weights1 <- weights[res.ind==1]
    weights2 <- weights[res.ind==2]
    weights3 <- weights[res.ind==3]
    weights4 <- weights[res.ind==4]
    weights5 <- weights[res.ind==5]
    weights6 <- weights[res.ind==6]
    weights7 <- weights[res.ind==7]
    weights8 <- weights[res.ind==8]
    
    beta.alpha.boot = matrix(NA, nrow = Boot, ncol = p)
    for (j in 1:Boot){
      SSAMPLE1     <- sample(1:nrow(X1), m/8, replace = T) 
      SSAMPLE2     <- sample(1:nrow(X2), m/8, replace = T)
      SSAMPLE3     <- sample(1:nrow(X3), m/8, replace = T)
      SSAMPLE4     <- sample(1:nrow(X4), m/8, replace = T)
      SSAMPLE5     <- sample(1:nrow(X5), m/8, replace = T)
      SSAMPLE6     <- sample(1:nrow(X6), m/8, replace = T)
      SSAMPLE7     <- sample(1:nrow(X7), m/8, replace = T)
      SSAMPLE8     <- sample(1:nrow(X8), m/8, replace = T)
      X.boot2 <- rbind(X1[SSAMPLE1,],X2[SSAMPLE2,],X3[SSAMPLE3,],X4[SSAMPLE4,],X5[SSAMPLE5,],X6[SSAMPLE6,],X7[SSAMPLE7,],X8[SSAMPLE8,])
      Y.boot2 <- c(Y1[SSAMPLE1],Y2[SSAMPLE2],Y3[SSAMPLE3],Y4[SSAMPLE4],Y5[SSAMPLE5],Y6[SSAMPLE6],Y7[SSAMPLE7],Y8[SSAMPLE8])
      weights.boot2 = c(weights1[SSAMPLE1],weights2[SSAMPLE2],weights3[SSAMPLE3],weights4[SSAMPLE4],weights5[SSAMPLE5],weights6[SSAMPLE6],weights7[SSAMPLE7],weights8[SSAMPLE8])
      
      beta.alpha.boot[j,] = lm(Y.boot2 ~ X.boot2, weights = weights.boot2)$coefficients[-1]; 
    }
    
    alpha = 0.05
    quantile.ci <- apply(beta.alpha.boot, 2, function(u) {
      quantile(u, prob = c(alpha/2, 1 - alpha/2), na.rm = na.rm)
    })
    
    bca.ci = apply(beta.alpha.boot, 2, function(u){
      if (na.rm == T){
        coxed::bca(na.omit(u))
      }else{
        coxed::bca(u)
      }
    })
    
    
    
    
    quantile.ci.acc[i,] = (true_beta >= quantile.ci[1,] & true_beta <= quantile.ci[2,])
    bca.ci.acc[i,] = (true_beta >= bca.ci[1,] & true_beta <= bca.ci[2,])
    
    quantile.length[i,] = quantile.ci[2,] - quantile.ci[1,]
    bca.length[i,] =  bca.ci[2,] - bca.ci[1,]
    print(i)
  }
  list("quantile.acc" = apply(quantile.ci.acc,2, function(x) {mean(x)}),
       "quantile.length" = apply(quantile.length,2, function(x) {mean(x)}),
       "bca.acc" = apply(bca.ci.acc,2, function(x) {mean(x)}),
       "bca.length" = apply(bca.length,2, function(x) {mean(x)}))
}


sbmmCI = function(X, Y, true_beta, 
                  LOOPS = 100, nOut, m = 32, verbose = T, largeOut = T, na.rm = F){
  n = nrow(X) 
  p = ncol(X)
  # Quantile Confidence Interval
  quantile.ci.acc = matrix(NA, nrow = LOOPS, ncol = p)
  quantile.length = matrix(NA, nrow = LOOPS, ncol = p)
  
  # BCa Confidence Interval
  bca.ci.acc = matrix(NA, nrow = LOOPS, ncol = p)
  bca.length = matrix(NA, nrow = LOOPS, ncol = p)
  for (i in 1:LOOPS){
    Y.full = Y + rnorm(n)
    n0 = sample(1:n, nOut, replace = F)
    if (largeOut == T){Y.full[n0] = runif(nOut, 1000, 2000)}
    if (largeOut == F){Y.full[n0] = runif(nOut, 30, 40)}
    Boot = 100
    rlm.full = rlm(Y.full ~ X, method='MM', maxit = 500)
    res.rank = rank(rlm.full$residuals)
    res.ind  <- rep(1,n);
    for(j1 in 1:n){
      if (res.rank[j1] <=n/8){res.ind[j1] <- res.ind[j1]+1};
      if (res.rank[j1] <=2*n/8){res.ind[j1] <- res.ind[j1]+1};
      if (res.rank[j1] <=3*n/8){res.ind[j1] <- res.ind[j1]+1};
      if (res.rank[j1] <=4*n/8){res.ind[j1] <- res.ind[j1]+1};
      if (res.rank[j1] <=5*n/8){res.ind[j1] <- res.ind[j1]+1};
      if (res.rank[j1] <=6*n/8){res.ind[j1] <- res.ind[j1]+1};
      if (res.rank[j1] <=7*n/8){res.ind[j1] <- res.ind[j1]+1};
    }
    X1 <- X[res.ind==1,]
    X2 <- X[res.ind==2,]
    X3 <- X[res.ind==3,]
    X4 <- X[res.ind==4,]
    X5 <- X[res.ind==5,]
    X6 <- X[res.ind==6,]
    X7 <- X[res.ind==7,]
    X8 <- X[res.ind==8,]
    Y1 <- Y.full[res.ind==1]
    Y2 <- Y.full[res.ind==2]
    Y3 <- Y.full[res.ind==3]
    Y4 <- Y.full[res.ind==4]
    Y5 <- Y.full[res.ind==5]
    Y6 <- Y.full[res.ind==6]
    Y7 <- Y.full[res.ind==7]
    Y8 <- Y.full[res.ind==8]
    
    beta.alpha.boot = matrix(NA, nrow = Boot, ncol = p)
    for (j in 1:Boot){
      while (!success) {
        tryCatch(
          {
            # Resample the data
            SSAMPLE1 <- sample(1:nrow(X1), m / 8, replace = TRUE)
            SSAMPLE2 <- sample(1:nrow(X2), m / 8, replace = TRUE)
            SSAMPLE3 <- sample(1:nrow(X3), m / 8, replace = TRUE)
            SSAMPLE4 <- sample(1:nrow(X4), m / 8, replace = TRUE)
            SSAMPLE5 <- sample(1:nrow(X5), m / 8, replace = TRUE)
            SSAMPLE6 <- sample(1:nrow(X6), m / 8, replace = TRUE)
            SSAMPLE7 <- sample(1:nrow(X7), m / 8, replace = TRUE)
            SSAMPLE8 <- sample(1:nrow(X8), m / 8, replace = TRUE)
            
            X.boot2 <- rbind(
              X1[SSAMPLE1,], X2[SSAMPLE2,], X3[SSAMPLE3,], X4[SSAMPLE4,],
              X5[SSAMPLE5,], X6[SSAMPLE6,], X7[SSAMPLE7,], X8[SSAMPLE8,]
            )
            Y.boot2 <- c(
              Y1[SSAMPLE1], Y2[SSAMPLE2], Y3[SSAMPLE3], Y4[SSAMPLE4],
              Y5[SSAMPLE5], Y6[SSAMPLE6], Y7[SSAMPLE7], Y8[SSAMPLE8]
            )
            
            beta.alpha.boot[j, ] <- rlm(Y.boot2 ~ X.boot2, method = 'MM', maxit = 500)$coefficients[-1]
            success <- TRUE  # If successful, set success to TRUE
          },
          error = function(e) {
          }
        )
      }
      
      
    }
    alpha = 0.05
    quantile.ci <- apply(beta.alpha.boot, 2, function(u) {
      quantile(u, prob = c(alpha/2, 1 - alpha/2), na.rm = na.rm)
    })
    
    bca.ci = apply(beta.alpha.boot, 2, function(u){
      if (na.rm == T){
        coxed::bca(na.omit(u))
      }else{
        coxed::bca(u)
      }
    })
    
    quantile.ci.acc[i,] = (true_beta >= quantile.ci[1,] & true_beta <= quantile.ci[2,])
    bca.ci.acc[i,] = (true_beta >= bca.ci[1,] & true_beta <= bca.ci[2,])
    
    quantile.length[i,] = quantile.ci[2,] - quantile.ci[1,]
    bca.length[i,] =  bca.ci[2,] - bca.ci[1,]
    print(i)
  }
  list("quantile.acc" = apply(quantile.ci.acc,2, function(x) {mean(x)}),
       "quantile.length" = apply(quantile.length,2, function(x) {mean(x)}),
       "bca.acc" = apply(bca.ci.acc,2, function(x) {mean(x)}),
       "bca.length" = apply(bca.length,2, function(x) {mean(x)}))
}





frbCI = function(X, Y, true_beta, 
                 LOOPS = 100, nOut, m = 32, verbose = T, largeOut = T,na.rm = F){
  n = nrow(X) 
  p = ncol(X)
  
  # Quantile Confidence Interval
  quantile.ci.acc = matrix(NA, nrow = LOOPS, ncol = p)
  quantile.length = matrix(NA, nrow = LOOPS, ncol = p)
  
  # BCa Confidence Interval
  bca.ci.acc = matrix(NA, nrow = LOOPS, ncol = p)
  bca.length = matrix(NA, nrow = LOOPS, ncol = p)
  
  
  
  for (i in 1:LOOPS){
    Y.full = Y + rnorm(n)
    n0 = sample(1:n, nOut, replace = F)
    if (largeOut == T){Y.full[n0] = runif(nOut, 1000, 2000)}
    if (largeOut == F){Y.full[n0] = runif(nOut, 30, 40)}
    Boot = 100
    lmrob.full = lmrob(Y.full ~ X, method='MM', maxit = 500, setting = "KS2014")
    res.rank = rank(lmrob.full$residuals)
    weights = lmrob.full$w
    res.ind  <- rep(1,n);
    for(j1 in 1:n){
      if (res.rank[j1] <=n/8){res.ind[j1] <- res.ind[j1]+1};
      if (res.rank[j1] <=2*n/8){res.ind[j1] <- res.ind[j1]+1};
      if (res.rank[j1] <=3*n/8){res.ind[j1] <- res.ind[j1]+1};
      if (res.rank[j1] <=4*n/8){res.ind[j1] <- res.ind[j1]+1};
      if (res.rank[j1] <=5*n/8){res.ind[j1] <- res.ind[j1]+1};
      if (res.rank[j1] <=6*n/8){res.ind[j1] <- res.ind[j1]+1};
      if (res.rank[j1] <=7*n/8){res.ind[j1] <- res.ind[j1]+1};
    }
    X1 <- X[res.ind==1,]
    X2 <- X[res.ind==2,]
    X3 <- X[res.ind==3,]
    X4 <- X[res.ind==4,]
    X5 <- X[res.ind==5,]
    X6 <- X[res.ind==6,]
    X7 <- X[res.ind==7,]
    X8 <- X[res.ind==8,]
    Y1 <- Y.full[res.ind==1]
    Y2 <- Y.full[res.ind==2]
    Y3 <- Y.full[res.ind==3]
    Y4 <- Y.full[res.ind==4]
    Y5 <- Y.full[res.ind==5]
    Y6 <- Y.full[res.ind==6]
    Y7 <- Y.full[res.ind==7]
    Y8 <- Y.full[res.ind==8]
    
    beta.alpha.boot = matrix(NA, nrow = Boot, ncol = p)
    for (j in 1:Boot){
      SSAMPLE1     <- sample(1:nrow(X1), m/8, replace = T) 
      SSAMPLE2     <- sample(1:nrow(X2), m/8, replace = T)
      SSAMPLE3     <- sample(1:nrow(X3), m/8, replace = T)
      SSAMPLE4     <- sample(1:nrow(X4), m/8, replace = T)
      SSAMPLE5     <- sample(1:nrow(X5), m/8, replace = T)
      SSAMPLE6     <- sample(1:nrow(X6), m/8, replace = T)
      SSAMPLE7     <- sample(1:nrow(X7), m/8, replace = T)
      SSAMPLE8     <- sample(1:nrow(X8), m/8, replace = T)
      X.boot2 <- rbind(X1[SSAMPLE1,],X2[SSAMPLE2,],X3[SSAMPLE3,],X4[SSAMPLE4,],X5[SSAMPLE5,],X6[SSAMPLE6,],X7[SSAMPLE7,],X8[SSAMPLE8,])
      Y.boot2 <- c(Y1[SSAMPLE1],Y2[SSAMPLE2],Y3[SSAMPLE3],Y4[SSAMPLE4],Y5[SSAMPLE5],Y6[SSAMPLE6],Y7[SSAMPLE7],Y8[SSAMPLE8])
      
      beta.alpha.boot[j,] =  lmrob(Y.boot2 ~ X.boot2, method='MM', maxit = 500, setting = "KS2014")$coefficients[-1]; 
    }
    
    # make CI (from HDCI library under ci() "quantile")
    alpha = 0.05
    quantile.ci <- apply(beta.alpha.boot, 2, function(u) {
      quantile(u, prob = c(alpha/2, 1 - alpha/2), na.rm = na.rm)
    })
    
    bca.ci = apply(beta.alpha.boot, 2, function(u){
      if (na.rm == T){
        coxed::bca(na.omit(u))
      }else{
        coxed::bca(u)
      }
    })
    
    
    
    
    quantile.ci.acc[i,] = (true_beta >= quantile.ci[1,] & true_beta <= quantile.ci[2,])
    bca.ci.acc[i,] = (true_beta >= bca.ci[1,] & true_beta <= bca.ci[2,])
    
    quantile.length[i,] = quantile.ci[2,] - quantile.ci[1,]
    bca.length[i,] =  bca.ci[2,] - bca.ci[1,]
    print(i)
  }
  list("quantile.acc" = apply(quantile.ci.acc,2, function(x) {mean(x)}),
       "quantile.length" = apply(quantile.length,2, function(x) {mean(x)}),
       "bca.acc" = apply(bca.ci.acc,2, function(x) {mean(x)}),
       "bca.length" = apply(bca.length,2, function(x) {mean(x)}))
}

