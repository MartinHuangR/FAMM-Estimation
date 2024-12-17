mmAll = function(X.full, y, ptrue, 
                 LOOPS = 100, nOut, SubM, m = 32, verbose = T, outliers,
                 largeOut = T, proteomics = F){
  
  pMax     <- dim(SubM)[1];
  n = nrow(X.full)
  M01.ind  <- rep(0,pMax)
  M02.ind  <- rep(0,pMax)
  M03.ind  <- rep(0,pMax)
  M04.ind  <- rep(0,pMax)
  M05.ind  <- rep(0,pMax)
  M06.ind  <- rep(0,pMax)
  Boot = 100
  
  for(i in 1:LOOPS){  
    y.full = y + rnorm(n)
    #ADD OUTLIERS: 32 and 16
    nO     <- sample(1:n, nOut, replace = FALSE)
    if (largeOut == T){y.full[nO] = runif(nOut, 1000, 2000)}
    if (largeOut == F){y.full[nO] = runif(nOut, 30, 40)}
    
    mad.full.ls <- mad(lm(y.full ~ X.full)$residuals);
    mad.full.ro <- mad(rlm(y.full ~ X.full, method='MM', maxit = 500)$residuals);
    
    # m out of n bootstrap
    Boot     <- 100;
    nBoot    <- rep(0,Boot);
    
    # parameters for robust model selection criterion
    b.huber  <- 2; 
    
    # robust
    MSC.1    <- rep(0,pMax)
    MSC01    <- rep(0,pMax)
    
    # ro+st
    MSC.2    <- rep(0,pMax)
    MSC02    <- rep(0,pMax)
    
    # ls
    MSC.3    <- rep(0,pMax)
    MSC03    <- rep(0,pMax)
    
    # ro+pen
    MSC.4    <- rep(0,pMax)
    MSC04    <- rep(0,pMax)
    
    # ro+st+pen
    MSC.5    <- rep(0,pMax)
    MSC05    <- rep(0,pMax)
    
    # ls+pen
    MSC.6    <- rep(0,pMax)
    MSC06    <- rep(0,pMax)
    
    
    # parameter estimation for each submodel
    beta.alpha <- list(rep(0,pMax));
    beta.alpha.ls <- list(rep(0,pMax));
    for(p in 1:pMax){ 
      X.red <- SubM[p,];
      X.red <- X.red==1;
      X.red <- X.full[,X.red];
      if (sum(X.red)<1){ beta.alpha[[p]] <- rlm(y.full ~ NULL, method='MM', maxit = 500)$coefficients;
      beta.alpha.ls[[p]] <- lm(y.full ~ NULL)$coefficients }
      if (sum(X.red)>0){ beta.alpha[[p]] <- rlm(y.full ~ X.red, method='MM', maxit = 500)$coefficients;
      beta.alpha.ls[[p]] <- lm(y.full ~ X.red)$coefficients }
    } #for end
    
    # Index sets for stratification samples
    res.rank <- rank(rlm(y.full ~ X.full, method='MM', maxit = 500)$residuals)
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
    X.full1 <- X.full[res.ind==1,]
    X.full2 <- X.full[res.ind==2,]
    X.full3 <- X.full[res.ind==3,]
    X.full4 <- X.full[res.ind==4,]
    X.full5 <- X.full[res.ind==5,]
    X.full6 <- X.full[res.ind==6,]
    X.full7 <- X.full[res.ind==7,]
    X.full8 <- X.full[res.ind==8,]
    y.full1 <- y.full[res.ind==1]
    y.full2 <- y.full[res.ind==2]
    y.full3 <- y.full[res.ind==3]
    y.full4 <- y.full[res.ind==4]
    y.full5 <- y.full[res.ind==5]
    y.full6 <- y.full[res.ind==6]
    y.full7 <- y.full[res.ind==7]
    y.full8 <- y.full[res.ind==8]
    
    
    for(j in 1:Boot){ 
      # generating bootstrap samples
      SAMPLE1      <- sample(1:n, m, replace = TRUE)
      mDim         <- length(names(table(SAMPLE1)));
      X.full.boot1 <- X.full[SAMPLE1,]
      y.full.boot1 <- y.full[SAMPLE1,]
      
      # generating stratified bootstrap samples
      SSAMPLE1     <- sample(1:nrow(X.full1), m/8, replace = TRUE)
      SSAMPLE2     <- sample(1:nrow(X.full2), m/8, replace = TRUE)
      SSAMPLE3     <- sample(1:nrow(X.full3), m/8, replace = TRUE)
      SSAMPLE4     <- sample(1:nrow(X.full4), m/8, replace = TRUE)
      SSAMPLE5     <- sample(1:nrow(X.full5), m/8, replace = TRUE)
      SSAMPLE6     <- sample(1:nrow(X.full6), m/8, replace = TRUE)
      SSAMPLE7     <- sample(1:nrow(X.full7), m/8, replace = TRUE)
      SSAMPLE8     <- sample(1:nrow(X.full8), m/8, replace = TRUE)
      X.full.boot2 <- rbind(X.full1[SSAMPLE1,],X.full2[SSAMPLE2,],X.full3[SSAMPLE3,],X.full4[SSAMPLE4,],X.full5[SSAMPLE5,],X.full6[SSAMPLE6,],X.full7[SSAMPLE7,],X.full8[SSAMPLE8,])
      y.full.boot2 <- c(y.full1[SSAMPLE1],y.full2[SSAMPLE2],y.full3[SSAMPLE3],y.full4[SSAMPLE4],y.full5[SSAMPLE5],y.full6[SSAMPLE6],y.full7[SSAMPLE7],y.full8[SSAMPLE8])
      
      
      # calculation model selection criteria for each submodel
      beta.alpha.boot1 <- list(rep(0,pMax));
      beta.alpha.boot2 <- list(rep(0,pMax));
      beta.alpha.boot3 <- list(rep(0,pMax));
      for(p in 1:pMax){ 
        X.red <- SubM[p,];
        X.red <- X.red==1;
        X.red.boot1 <- X.full.boot1[,X.red];
        X.red.boot2 <- X.full.boot2[,X.red];
        X.red.full <- cbind(rep(1,n),X.full[,X.red]);
        
        if (sum(X.red)<1){ 
          beta.alpha.boot1[[p]] <- rlm(y.full.boot1 ~ NULL, method='MM', maxit = 500)$coefficients;
          res.boot1 <- y.full - X.red.full %*%  beta.alpha.boot1[[p]];
          res.small.boot1 <- res.boot1[res.boot1<(b.huber*mad.full.ro)**2]; 
          MSC.1[p] <- MSC.1[p] + (sum(( res.small.boot1 )**2) + ((b.huber*mad.full.ro)**2)*(n-length(res.small.boot1)))/(n)
          beta.alpha.boot2[[p]] <- rlm(y.full.boot2 ~ NULL, method='MM', maxit = 500)$coefficients;
          res.boot2 <- y.full - X.red.full %*%  beta.alpha.boot2[[p]];
          res.small.boot2 <- res.boot2[res.boot2<(b.huber*mad.full.ro)**2] ;
          MSC.2[p] <- MSC.2[p] + (sum(( res.small.boot2 )**2) + ((b.huber*mad.full.ro)**2)*(n-length(res.small.boot2)))/(n);
          beta.alpha.boot3[[p]] <- lm(y.full.boot1 ~ NULL)$coefficients;
          res.boot3 <- y.full - X.red.full %*%  beta.alpha.boot3[[p]];
          res.small.boot3 <- res.boot3[res.boot3<(b.huber*mad.full.ls)**2]; 
          MSC.3[p] <- MSC.3[p] + (sum(( res.small.boot3 )**2) + ((b.huber*mad.full.ls)**2)*(n-length(res.small.boot3)))/(n)}
        
        if (sum(X.red)>0){ 
          beta.alpha.boot1[[p]] <- rlm(y.full.boot1 ~ X.red.boot1, method='MM', maxit = 500)$coefficients;
          res.boot1 <- y.full - X.red.full %*%  beta.alpha.boot1[[p]];
          res.small.boot1 <- res.boot1[res.boot1<(b.huber*mad.full.ro)**2]; 
          MSC.1[p] <- MSC.1[p] + (sum(( res.small.boot1 )**2) + ((b.huber*mad.full.ro)**2)*(n-length(res.small.boot1)))/(n)
          beta.alpha.boot2[[p]] <- rlm(y.full.boot2 ~ X.red.boot2, method='MM', maxit = 500)$coefficients;
          res.boot2 <- y.full - X.red.full %*%  beta.alpha.boot2[[p]];
          res.small.boot2 <- res.boot2[res.boot2<(b.huber*mad.full.ro)**2] ;
          MSC.2[p] <- MSC.2[p] + (sum(( res.small.boot2 )**2) + ((b.huber*mad.full.ro)**2)*(n-length(res.small.boot2)))/(n);
          beta.alpha.boot3[[p]] <- lm(y.full.boot1 ~ X.red.boot1)$coefficients;
          res.boot3 <- y.full - X.red.full %*%  beta.alpha.boot3[[p]];
          res.small.boot3 <- res.boot3[res.boot3<(b.huber*mad.full.ls)**2]; 
          MSC.3[p] <- MSC.3[p] + (sum(( res.small.boot3 )**2) + ((b.huber*mad.full.ls)**2)*(n-length(res.small.boot3)))/(n)}
        
      }#for end
      print(c(i,j));
    }#boot end
    for(p in 1:pMax){
      X.red <- SubM[p,];
      pDim  <- sum(X.red);
      X.red <- X.red==1;
      X.red.full <- cbind(rep(1,n),X.full[,X.red]);
      if (sum(X.red)<1){ 
        beta.alpha1 <- rlm(y.full ~ NULL, method='MM', maxit = 500)$coefficients;
        res         <- y.full - X.red.full %*%  beta.alpha1;
        res.small   <- res[res<(b.huber*mad.full.ro)**2]
        MSC.4[p] <- MSC.1[p]/Boot;
        MSC.5[p] <- MSC.2[p]/Boot;
        penal    <- (sum(( res.small)**2) + ((b.huber*mad.full.ro)**2)*(n-length(res.small)) +2*log(n)*pDim*(mad.full.ro**2) )/(n)
        MSC.4[p] <- penal + MSC.4[p]; 
        MSC.5[p] <- penal + MSC.5[p];
        beta.alpha1 <- lm(y.full ~ NULL)$coefficients;
        res         <- y.full - X.red.full %*%  beta.alpha1;
        res.small   <- res[res<(b.huber*mad.full.ls)**2]
        MSC.6[p]    <- MSC.3[p]/Boot + (sum(( res.small)**2) + ((b.huber*mad.full.ls)**2)*(n-length(res.small)) +2*log(n)*pDim*(mad.full.ls**2) )/(n)
      } 
      if (sum(X.red)>0){ 
        beta.alpha1 <- rlm(y.full ~ X.red.full[,-1], method='MM', maxit = 500)$coefficients;
        res         <- y.full - X.red.full %*%  beta.alpha1;
        res.small   <- res[res<(b.huber*mad.full.ro)**2]
        MSC.4[p] <- MSC.1[p]/Boot;
        MSC.5[p] <- MSC.2[p]/Boot;
        penal    <- (sum(( res.small)**2) + ((b.huber*mad.full.ro)**2)*(n-length(res.small)) +2*log(n)*pDim*(mad.full.ro**2) )/(n)
        MSC.4[p] <- penal + MSC.4[p]; 
        MSC.5[p] <- penal + MSC.5[p]; 
        beta.alpha1 <- lm(y.full ~ X.red.full[,-1])$coefficients;
        res         <- y.full - X.red.full %*%  beta.alpha1;
        res.small   <- res[res<(b.huber*mad.full.ls)**2]
        MSC.6[p]    <- MSC.3[p]/Boot + (sum(( res.small)**2) + ((b.huber*mad.full.ls)**2)*(n-length(res.small)) +2*log(n)*pDim*(mad.full.ls**2) )/(n)
      } 
    }#for end
    
    
    M01.ind[which.min(MSC.1)] = M01.ind[which.min(MSC.1)] + 1
    M02.ind[which.min(MSC.2)] = M02.ind[which.min(MSC.2)] + 1
    M03.ind[which.min(MSC.3)] = M03.ind[which.min(MSC.3)] + 1
    M04.ind[which.min(MSC.4)] = M04.ind[which.min(MSC.4)] + 1
    M05.ind[which.min(MSC.5)] = M05.ind[which.min(MSC.5)] + 1
    M06.ind[which.min(MSC.6)] = M06.ind[which.min(MSC.6)] + 1
    if ((i %% 100) == 0){
      print(paste("Outlier:",outliers))
      print(cbind(SubM,M05.ind,M04.ind,M02.ind,M01.ind,M06.ind,M03.ind));
    }
    
  } #simu end
  r = rbind(M01.ind, M02.ind, M03.ind, M04.ind, M05.ind, M06.ind) |> as.data.frame()
  if (proteomics == T){
    colnames(r) = c("True", "Model 2", "Model 3")
    r$Outliers = outliers
    r$Criterion = c("M01", "M02", "M03", "M04", "M05", "M06")
    rownames(r) <- NULL
  }else{
    colnames(r) = c("True", "Model 2", "Model 3", "Model 4")
    r$Outliers = outliers
    r$Criterion = c("M01", "M02", "M03", "M04", "M05", "M06")
    rownames(r) <- NULL
  }
  r
}


fammAll = function(X.full, y, ptrue, 
                   LOOPS = 100, nOut, SubM, m = 32, verbose = T, outliers,
                   largeOut = T, proteomics = F){
  
  pMax     <- dim(SubM)[1];
  n = nrow(X.full)
  M01.ind  <- rep(0,pMax)
  M02.ind  <- rep(0,pMax)
  M03.ind  <- rep(0,pMax)
  M04.ind  <- rep(0,pMax)
  M05.ind  <- rep(0,pMax)
  M06.ind  <- rep(0,pMax)
  Boot = 100
  
  for(i in 1:LOOPS){  
    y.full = y + rnorm(n)
    nO     <- sample(1:n, nOut, replace = FALSE)
    if (largeOut == T){y.full[nO] = runif(nOut, 1000, 2000)}
    if (largeOut == F){y.full[nO] = runif(nOut, 30, 40)}
    
    mad.full.ls <- mad(lm(y.full ~ X.full)$residuals);
    mad.full.ro <- mad(rlm(y.full ~ X.full, method='MM', maxit = 500)$residuals);
    
    # m out of n bootstrap
    Boot     <- 100;
    nBoot    <- rep(0,Boot);
    
    # parameters for robust model selection criterion
    b.huber  <- 2; 
    
    # ro
    MSC.1    <- rep(0,pMax)
    MSC01    <- rep(0,pMax)
    
    # ro+st
    MSC.2    <- rep(0,pMax)
    MSC02    <- rep(0,pMax)
    
    # ls
    MSC.3    <- rep(0,pMax)
    MSC03    <- rep(0,pMax)
    
    # ro+pen
    MSC.4    <- rep(0,pMax)
    MSC04    <- rep(0,pMax)
    
    # ro+st+pen
    MSC.5    <- rep(0,pMax)
    MSC05    <- rep(0,pMax)
    
    # ls+pen
    MSC.6    <- rep(0,pMax)
    MSC06    <- rep(0,pMax)
    
    
    # parameter estimation for each submodel
    beta.alpha <- list(rep(0,pMax));
    beta.alpha.ls <- list(rep(0,pMax));
    for(p in 1:pMax){ 
      X.red <- SubM[p,];
      X.red <- X.red==1;
      X.red <- X.full[,X.red];
      if (sum(X.red)<1){ beta.alpha[[p]] <- rlm(y.full ~ NULL, method='MM', maxit = 500)$coefficients;
      beta.alpha.ls[[p]] <- lm(y.full ~ NULL)$coefficients }
      if (sum(X.red)>0){ beta.alpha[[p]] <- rlm(y.full ~ X.red, method='MM', maxit = 500)$coefficients;
      beta.alpha.ls[[p]] <- lm(y.full ~ X.red)$coefficients }
    } #for end
    
    # Index sets for stratification samples
    rlm.full = rlm(y.full ~ X.full, method='MM', maxit = 500)
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
    X.full1 <- X.full[res.ind==1,]
    X.full2 <- X.full[res.ind==2,]
    X.full3 <- X.full[res.ind==3,]
    X.full4 <- X.full[res.ind==4,]
    X.full5 <- X.full[res.ind==5,]
    X.full6 <- X.full[res.ind==6,]
    X.full7 <- X.full[res.ind==7,]
    X.full8 <- X.full[res.ind==8,]
    y.full1 <- y.full[res.ind==1]
    y.full2 <- y.full[res.ind==2]
    y.full3 <- y.full[res.ind==3]
    y.full4 <- y.full[res.ind==4]
    y.full5 <- y.full[res.ind==5]
    y.full6 <- y.full[res.ind==6]
    y.full7 <- y.full[res.ind==7]
    y.full8 <- y.full[res.ind==8]
    weights1 <- weights[res.ind==1]
    weights2 <- weights[res.ind==2]
    weights3 <- weights[res.ind==3]
    weights4 <- weights[res.ind==4]
    weights5 <- weights[res.ind==5]
    weights6 <- weights[res.ind==6]
    weights7 <- weights[res.ind==7]
    weights8 <- weights[res.ind==8]
    
    for(j in 1:Boot){ 
      # generating bootstrap samples
      SAMPLE1      <- sample(1:n, m, replace = TRUE)
      mDim         <- length(names(table(SAMPLE1)));
      X.full.boot1 <- X.full[SAMPLE1,]
      y.full.boot1 <- y.full[SAMPLE1,]
      weights.boot1 = weights[SAMPLE1]
      # generating stratified bootstrap samples
      SSAMPLE1     <- sample(1:nrow(X.full1), m/8, replace = TRUE)
      SSAMPLE2     <- sample(1:nrow(X.full2), m/8, replace = TRUE)
      SSAMPLE3     <- sample(1:nrow(X.full3), m/8, replace = TRUE)
      SSAMPLE4     <- sample(1:nrow(X.full4), m/8, replace = TRUE)
      SSAMPLE5     <- sample(1:nrow(X.full5), m/8, replace = TRUE)
      SSAMPLE6     <- sample(1:nrow(X.full6), m/8, replace = TRUE)
      SSAMPLE7     <- sample(1:nrow(X.full7), m/8, replace = TRUE)
      SSAMPLE8     <- sample(1:nrow(X.full8), m/8, replace = TRUE)
      X.full.boot2 <- rbind(X.full1[SSAMPLE1,],X.full2[SSAMPLE2,],X.full3[SSAMPLE3,],X.full4[SSAMPLE4,],X.full5[SSAMPLE5,],X.full6[SSAMPLE6,],X.full7[SSAMPLE7,],X.full8[SSAMPLE8,])
      y.full.boot2 <- c(y.full1[SSAMPLE1],y.full2[SSAMPLE2],y.full3[SSAMPLE3],y.full4[SSAMPLE4],y.full5[SSAMPLE5],y.full6[SSAMPLE6],y.full7[SSAMPLE7],y.full8[SSAMPLE8])
      weights.boot2 = c(weights1[SSAMPLE1],weights2[SSAMPLE2],weights3[SSAMPLE3],weights4[SSAMPLE4],weights5[SSAMPLE5],weights6[SSAMPLE6],weights7[SSAMPLE7],weights8[SSAMPLE8])
      
      # calculation model selection criteria for each submodel
      beta.alpha.boot1 <- list(rep(0,pMax));
      beta.alpha.boot2 <- list(rep(0,pMax));
      beta.alpha.boot3 <- list(rep(0,pMax));
      for(p in 1:pMax){ 
        X.red <- SubM[p,];
        X.red <- X.red==1;
        X.red.boot1 <- X.full.boot1[,X.red];
        X.red.boot2 <- X.full.boot2[,X.red];
        X.red.full <- cbind(rep(1,n),X.full[,X.red]);
        
        if (sum(X.red)<1){ 
          beta.alpha.boot1[[p]] <- lm(y.full.boot1 ~ NULL, weights = weights.boot1)$coefficients;
          res.boot1 <- y.full - X.red.full %*%  beta.alpha.boot1[[p]];
          res.small.boot1 <- res.boot1[res.boot1<(b.huber*mad.full.ro)**2]; 
          MSC.1[p] <- MSC.1[p] + (sum(( res.small.boot1 )**2) + ((b.huber*mad.full.ro)**2)*(n-length(res.small.boot1)))/(n)
          beta.alpha.boot2[[p]] <- lm(y.full.boot2 ~ NULL, weights = weights.boot2)$coefficients;
          res.boot2 <- y.full - X.red.full %*%  beta.alpha.boot2[[p]];
          res.small.boot2 <- res.boot2[res.boot2<(b.huber*mad.full.ro)**2] ;
          MSC.2[p] <- MSC.2[p] + (sum(( res.small.boot2 )**2) + ((b.huber*mad.full.ro)**2)*(n-length(res.small.boot2)))/(n);
          beta.alpha.boot3[[p]] <- lm(y.full.boot1 ~ NULL)$coefficients;
          res.boot3 <- y.full - X.red.full %*%  beta.alpha.boot3[[p]];
          res.small.boot3 <- res.boot3[res.boot3<(b.huber*mad.full.ls)**2]; 
          MSC.3[p] <- MSC.3[p] + (sum(( res.small.boot3 )**2) + ((b.huber*mad.full.ls)**2)*(n-length(res.small.boot3)))/(n)}
        
        if (sum(X.red)>0){ 
          beta.alpha.boot1[[p]] <- lm(y.full.boot1 ~ X.red.boot1, weights = weights.boot1)$coefficients;
          res.boot1 <- y.full - X.red.full %*%  beta.alpha.boot1[[p]];
          res.small.boot1 <- res.boot1[res.boot1<(b.huber*mad.full.ro)**2]; 
          MSC.1[p] <- MSC.1[p] + (sum(( res.small.boot1 )**2) + ((b.huber*mad.full.ro)**2)*(n-length(res.small.boot1)))/(n)
          beta.alpha.boot2[[p]] <- lm(y.full.boot2 ~ X.red.boot2, weights = weights.boot2)$coefficients;
          res.boot2 <- y.full - X.red.full %*%  beta.alpha.boot2[[p]];
          res.small.boot2 <- res.boot2[res.boot2<(b.huber*mad.full.ro)**2] ;
          MSC.2[p] <- MSC.2[p] + (sum(( res.small.boot2 )**2) + ((b.huber*mad.full.ro)**2)*(n-length(res.small.boot2)))/(n);
          beta.alpha.boot3[[p]] <- lm(y.full.boot1 ~ X.red.boot1)$coefficients;
          res.boot3 <- y.full - X.red.full %*%  beta.alpha.boot3[[p]];
          res.small.boot3 <- res.boot3[res.boot3<(b.huber*mad.full.ls)**2]; 
          MSC.3[p] <- MSC.3[p] + (sum(( res.small.boot3 )**2) + ((b.huber*mad.full.ls)**2)*(n-length(res.small.boot3)))/(n)}
        
      }#for end
      print(c(i,j));
    }#boot end
    for(p in 1:pMax){
      X.red <- SubM[p,];
      pDim  <- sum(X.red);
      X.red <- X.red==1;
      X.red.full <- cbind(rep(1,n),X.full[,X.red]);
      if (sum(X.red)<1){ 
        beta.alpha1 <- lm(y.full ~ NULL,  weights = weights)$coefficients;
        res         <- y.full - X.red.full %*%  beta.alpha1;
        res.small   <- res[res<(b.huber*mad.full.ro)**2]
        MSC.4[p] <- MSC.1[p]/Boot;
        MSC.5[p] <- MSC.2[p]/Boot;
        penal    <- (sum(( res.small)**2) + ((b.huber*mad.full.ro)**2)*(n-length(res.small)) +2*log(n)*pDim*(mad.full.ro**2) )/(n)
        MSC.4[p] <- penal + MSC.4[p]; 
        MSC.5[p] <- penal + MSC.5[p];
        beta.alpha1 <- lm(y.full ~ NULL)$coefficients;
        res         <- y.full - X.red.full %*%  beta.alpha1;
        res.small   <- res[res<(b.huber*mad.full.ls)**2]
        MSC.6[p]    <- MSC.3[p]/Boot + (sum(( res.small)**2) + ((b.huber*mad.full.ls)**2)*(n-length(res.small)) +2*log(n)*pDim*(mad.full.ls**2) )/(n)
      } 
      if (sum(X.red)>0){ 
        beta.alpha1 <- lm(y.full ~ X.red.full[,-1], weights = weights)$coefficients;
        res         <- y.full - X.red.full %*%  beta.alpha1;
        res.small   <- res[res<(b.huber*mad.full.ro)**2]
        MSC.4[p] <- MSC.1[p]/Boot;
        MSC.5[p] <- MSC.2[p]/Boot;
        penal    <- (sum(( res.small)**2) + ((b.huber*mad.full.ro)**2)*(n-length(res.small)) +2*log(n)*pDim*(mad.full.ro**2) )/(n)
        MSC.4[p] <- penal + MSC.4[p]; 
        MSC.5[p] <- penal + MSC.5[p]; 
        beta.alpha1 <- lm(y.full ~ X.red.full[,-1])$coefficients;
        res         <- y.full - X.red.full %*%  beta.alpha1;
        res.small   <- res[res<(b.huber*mad.full.ls)**2]
        MSC.6[p]    <- MSC.3[p]/Boot + (sum(( res.small)**2) + ((b.huber*mad.full.ls)**2)*(n-length(res.small)) +2*log(n)*pDim*(mad.full.ls**2) )/(n)
      } 
    }#for end
    
    
    M01.ind[which.min(MSC.1)] = M01.ind[which.min(MSC.1)] + 1
    M02.ind[which.min(MSC.2)] = M02.ind[which.min(MSC.2)] + 1
    M03.ind[which.min(MSC.3)] = M03.ind[which.min(MSC.3)] + 1
    M04.ind[which.min(MSC.4)] = M04.ind[which.min(MSC.4)] + 1
    M05.ind[which.min(MSC.5)] = M05.ind[which.min(MSC.5)] + 1
    M06.ind[which.min(MSC.6)] = M06.ind[which.min(MSC.6)] + 1
    
    if ((i %% 100) == 0){
      print(paste("Outlier:",outliers))
      print(cbind(SubM,M05.ind,M04.ind,M02.ind,M01.ind,M06.ind,M03.ind));
    }
    
  } #simu end
  
  r = rbind(M01.ind, M02.ind, M03.ind, M04.ind, M05.ind, M06.ind) |> as.data.frame()
  if (proteomics == T){
    colnames(r) = c("True", "Model 2", "Model 3")
    r$Outliers = outliers
    r$Criterion = c("M01", "M02", "M03", "M04", "M05", "M06")
    rownames(r) <- NULL
  }else{
    colnames(r) = c("True", "Model 2", "Model 3", "Model 4")
    r$Outliers = outliers
    r$Criterion = c("M01", "M02", "M03", "M04", "M05", "M06")
    rownames(r) <- NULL
  }
  r
  
}


frbAll = function(X.full, y, ptrue, 
                  LOOPS = 100, nOut, SubM, m = 32, verbose = T, outliers,
                  largeOut = T, proteomics = F){
  
  pMax     <- dim(SubM)[1];
  n = nrow(X.full)
  M01.ind  <- rep(0,pMax)
  M02.ind  <- rep(0,pMax)
  M03.ind  <- rep(0,pMax)
  M04.ind  <- rep(0,pMax)
  M05.ind  <- rep(0,pMax)
  M06.ind  <- rep(0,pMax)
  Boot = 100
  
  for(i in 1:LOOPS){
    y.full = y + rnorm(n)
    nO     <- sample(1:n, nOut, replace = FALSE)
    if (largeOut == T){y.full[nO] = runif(nOut, 1000, 2000)}
    if (largeOut == F){y.full[nO] = runif(nOut, 30, 40)}
    
    mad.full.ls <- mad(lm(y.full ~ X.full)$residuals);
    mad.full.ro <- mad(lmrob(y.full ~ X.full, method = "MM", setting = "KS2014", maxit = 500)$residuals);
    
    # m out of n bootstrap
    Boot     <- 100;
    nBoot    <- rep(0,Boot);
    
    # parameters for robust model selection criterion
    b.huber  <- 2; 
    
    
    # robust
    MSC.1    <- rep(0,pMax)
    MSC01    <- rep(0,pMax)
    
    # ro+st
    MSC.2    <- rep(0,pMax)
    MSC02    <- rep(0,pMax)
    
    # ls
    MSC.3    <- rep(0,pMax)
    MSC03    <- rep(0,pMax)
    
    # ro+pen
    MSC.4    <- rep(0,pMax)
    MSC04    <- rep(0,pMax)
    
    # ro+st+pen
    MSC.5    <- rep(0,pMax)
    MSC05    <- rep(0,pMax)
    
    # ls+pen
    MSC.6    <- rep(0,pMax)
    MSC06    <- rep(0,pMax)
    
    
    # parameter estimation for each submodel
    beta.alpha <- list(rep(0,pMax));
    beta.alpha.ls <- list(rep(0,pMax));
    for(p in 1:pMax){ 
      X.red <- SubM[p,];
      X.red <- X.red==1;
      X.red <- X.full[,X.red];
      if (sum(X.red)<1){ beta.alpha[[p]] <- lmrob(y.full ~ NULL, method = "MM", setting = "KS2014", maxit = 500)$coefficients;
      beta.alpha.ls[[p]] <- lm(y.full ~ NULL)$coefficients }
      if (sum(X.red)>0){ beta.alpha[[p]] <- lmrob(y.full ~ X.red, method = "MM", setting = "KS2014", maxit = 500)$coefficients;
      beta.alpha.ls[[p]] <- lm(y.full ~ X.red)$coefficients }
    } #for end
    
    # Index sets for stratification samples
    res.rank <- rank(lmrob(y.full ~ X.full, method = "MM", setting = "KS2014", maxit = 500)$residuals)
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
    X.full1 <- X.full[res.ind==1,]
    X.full2 <- X.full[res.ind==2,]
    X.full3 <- X.full[res.ind==3,]
    X.full4 <- X.full[res.ind==4,]
    X.full5 <- X.full[res.ind==5,]
    X.full6 <- X.full[res.ind==6,]
    X.full7 <- X.full[res.ind==7,]
    X.full8 <- X.full[res.ind==8,]
    y.full1 <- y.full[res.ind==1]
    y.full2 <- y.full[res.ind==2]
    y.full3 <- y.full[res.ind==3]
    y.full4 <- y.full[res.ind==4]
    y.full5 <- y.full[res.ind==5]
    y.full6 <- y.full[res.ind==6]
    y.full7 <- y.full[res.ind==7]
    y.full8 <- y.full[res.ind==8]
    
    
    for(j in 1:Boot){ 
      # generating bootstrap samples
      SAMPLE1      <- sample(1:n, m, replace = TRUE)
      mDim         <- length(names(table(SAMPLE1)));
      X.full.boot1 <- X.full[SAMPLE1,]
      y.full.boot1 <- y.full[SAMPLE1,]
      
      # generating stratified bootstrap samples
      SSAMPLE1     <- sample(1:nrow(X.full1), m/8, replace = TRUE)
      SSAMPLE2     <- sample(1:nrow(X.full2), m/8, replace = TRUE)
      SSAMPLE3     <- sample(1:nrow(X.full3), m/8, replace = TRUE)
      SSAMPLE4     <- sample(1:nrow(X.full4), m/8, replace = TRUE)
      SSAMPLE5     <- sample(1:nrow(X.full5), m/8, replace = TRUE)
      SSAMPLE6     <- sample(1:nrow(X.full6), m/8, replace = TRUE)
      SSAMPLE7     <- sample(1:nrow(X.full7), m/8, replace = TRUE)
      SSAMPLE8     <- sample(1:nrow(X.full8), m/8, replace = TRUE)
      X.full.boot2 <- rbind(X.full1[SSAMPLE1,],X.full2[SSAMPLE2,],X.full3[SSAMPLE3,],X.full4[SSAMPLE4,],X.full5[SSAMPLE5,],X.full6[SSAMPLE6,],X.full7[SSAMPLE7,],X.full8[SSAMPLE8,])
      y.full.boot2 <- c(y.full1[SSAMPLE1],y.full2[SSAMPLE2],y.full3[SSAMPLE3],y.full4[SSAMPLE4],y.full5[SSAMPLE5],y.full6[SSAMPLE6],y.full7[SSAMPLE7],y.full8[SSAMPLE8])
      
      
      # calculation model selection criteria for each submodel
      beta.alpha.boot1 <- list(rep(0,pMax));
      beta.alpha.boot2 <- list(rep(0,pMax));
      beta.alpha.boot3 <- list(rep(0,pMax));
      for(p in 1:pMax){ 
        X.red <- SubM[p,];
        X.red <- X.red==1;
        X.red.boot1 <- X.full.boot1[,X.red];
        X.red.boot2 <- X.full.boot2[,X.red];
        X.red.full <- cbind(rep(1,n),X.full[,X.red]);
        
        if (sum(X.red)<1){ 
          beta.alpha.boot1[[p]] <- lmrob(y.full.boot1 ~ NULL, method = "MM", setting = "KS2014", maxit = 500)$coefficients;
          res.boot1 <- y.full - X.red.full %*%  beta.alpha.boot1[[p]];
          res.small.boot1 <- res.boot1[res.boot1<(b.huber*mad.full.ro)**2]; 
          MSC.1[p] <- MSC.1[p] + (sum(( res.small.boot1 )**2) + ((b.huber*mad.full.ro)**2)*(n-length(res.small.boot1)))/(n)
          beta.alpha.boot2[[p]] <- lmrob(y.full.boot2 ~ NULL, method = "MM", setting = "KS2014", maxit = 500)$coefficients;
          res.boot2 <- y.full - X.red.full %*%  beta.alpha.boot2[[p]];
          res.small.boot2 <- res.boot2[res.boot2<(b.huber*mad.full.ro)**2] ;
          MSC.2[p] <- MSC.2[p] + (sum(( res.small.boot2 )**2) + ((b.huber*mad.full.ro)**2)*(n-length(res.small.boot2)))/(n);
          beta.alpha.boot3[[p]] <- lm(y.full.boot1 ~ NULL)$coefficients;
          res.boot3 <- y.full - X.red.full %*%  beta.alpha.boot3[[p]];
          res.small.boot3 <- res.boot3[res.boot3<(b.huber*mad.full.ls)**2]; 
          MSC.3[p] <- MSC.3[p] + (sum(( res.small.boot3 )**2) + ((b.huber*mad.full.ls)**2)*(n-length(res.small.boot3)))/(n)}
        
        if (sum(X.red)>0){ 
          beta.alpha.boot1[[p]] <- lmrob(y.full.boot1 ~ X.red.boot1, method = "MM", setting = "KS2014", maxit = 500)$coefficients;
          res.boot1 <- y.full - X.red.full %*%  beta.alpha.boot1[[p]];
          res.small.boot1 <- res.boot1[res.boot1<(b.huber*mad.full.ro)**2]; 
          MSC.1[p] <- MSC.1[p] + (sum(( res.small.boot1 )**2) + ((b.huber*mad.full.ro)**2)*(n-length(res.small.boot1)))/(n)
          beta.alpha.boot2[[p]] <- lmrob(y.full.boot2 ~ X.red.boot2, method = "MM", setting = "KS2014", maxit = 500)$coefficients;
          res.boot2 <- y.full - X.red.full %*%  beta.alpha.boot2[[p]];
          res.small.boot2 <- res.boot2[res.boot2<(b.huber*mad.full.ro)**2] ;
          MSC.2[p] <- MSC.2[p] + (sum(( res.small.boot2 )**2) + ((b.huber*mad.full.ro)**2)*(n-length(res.small.boot2)))/(n);
          beta.alpha.boot3[[p]] <- lm(y.full.boot1 ~ X.red.boot1)$coefficients;
          res.boot3 <- y.full - X.red.full %*%  beta.alpha.boot3[[p]];
          res.small.boot3 <- res.boot3[res.boot3<(b.huber*mad.full.ls)**2]; 
          MSC.3[p] <- MSC.3[p] + (sum(( res.small.boot3 )**2) + ((b.huber*mad.full.ls)**2)*(n-length(res.small.boot3)))/(n)}
        
      }#for end
      print(c(i,j));
    }#boot end
    for(p in 1:pMax){
      X.red <- SubM[p,];
      pDim  <- sum(X.red);
      X.red <- X.red==1;
      X.red.full <- cbind(rep(1,n),X.full[,X.red]);
      if (sum(X.red)<1){ 
        beta.alpha1 <- lmrob(y.full ~ NULL, method = "MM", setting = "KS2014", maxit = 500)$coefficients;
        res         <- y.full - X.red.full %*%  beta.alpha1;
        res.small   <- res[res<(b.huber*mad.full.ro)**2]
        MSC.4[p] <- MSC.1[p]/Boot;
        MSC.5[p] <- MSC.2[p]/Boot;
        penal    <- (sum(( res.small)**2) + ((b.huber*mad.full.ro)**2)*(n-length(res.small)) +2*log(n)*pDim*(mad.full.ro**2) )/(n)
        MSC.4[p] <- penal + MSC.4[p]; 
        MSC.5[p] <- penal + MSC.5[p];
        beta.alpha1 <- lm(y.full ~ NULL)$coefficients;
        res         <- y.full - X.red.full %*%  beta.alpha1;
        res.small   <- res[res<(b.huber*mad.full.ls)**2]
        MSC.6[p]    <- MSC.3[p]/Boot + (sum(( res.small)**2) + ((b.huber*mad.full.ls)**2)*(n-length(res.small)) +2*log(n)*pDim*(mad.full.ls**2) )/(n)
      } 
      if (sum(X.red)>0){ 
        beta.alpha1 <- lmrob(y.full ~ X.red.full[,-1], method = "MM", setting = "KS2014", maxit = 500)$coefficients;
        res         <- y.full - X.red.full %*%  beta.alpha1;
        res.small   <- res[res<(b.huber*mad.full.ro)**2]
        MSC.4[p] <- MSC.1[p]/Boot;
        MSC.5[p] <- MSC.2[p]/Boot;
        penal    <- (sum(( res.small)**2) + ((b.huber*mad.full.ro)**2)*(n-length(res.small)) +2*log(n)*pDim*(mad.full.ro**2) )/(n)
        MSC.4[p] <- penal + MSC.4[p]; 
        MSC.5[p] <- penal + MSC.5[p]; 
        beta.alpha1 <- lm(y.full ~ X.red.full[,-1])$coefficients;
        res         <- y.full - X.red.full %*%  beta.alpha1;
        res.small   <- res[res<(b.huber*mad.full.ls)**2]
        MSC.6[p]    <- MSC.3[p]/Boot + (sum(( res.small)**2) + ((b.huber*mad.full.ls)**2)*(n-length(res.small)) +2*log(n)*pDim*(mad.full.ls**2) )/(n)
      } 
    }#for end
    
    
    M01.ind[which.min(MSC.1)] = M01.ind[which.min(MSC.1)] + 1
    M02.ind[which.min(MSC.2)] = M02.ind[which.min(MSC.2)] + 1
    M03.ind[which.min(MSC.3)] = M03.ind[which.min(MSC.3)] + 1
    M04.ind[which.min(MSC.4)] = M04.ind[which.min(MSC.4)] + 1
    M05.ind[which.min(MSC.5)] = M05.ind[which.min(MSC.5)] + 1
    M06.ind[which.min(MSC.6)] = M06.ind[which.min(MSC.6)] + 1
    if ((i %% 100) == 0){
      print(paste("Outlier:",outliers))
      print(cbind(SubM,M05.ind,M04.ind,M02.ind,M01.ind,M06.ind,M03.ind));
    }
    
  } #simu end
  r = rbind(M01.ind, M02.ind, M03.ind, M04.ind, M05.ind, M06.ind) |> as.data.frame()
  if (proteomics == T){
    colnames(r) = c("True", "Model 2", "Model 3")
    r$Outliers = outliers
    r$Criterion = c("M01", "M02", "M03", "M04", "M05", "M06")
    rownames(r) <- NULL
  }else{
    colnames(r) = c("True", "Model 2", "Model 3", "Model 4")
    r$Outliers = outliers
    r$Criterion = c("M01", "M02", "M03", "M04", "M05", "M06")
    rownames(r) <- NULL
  }
  r
}

getSubM=function(p){
  SubM <- list()
  for (k in 1:p) {
    SubMp <- combn(p, k)
    for (i in 1:ncol(SubMp)) {
      SubMpi <- rep(0, p)         
      SubMpi[SubMp[, i]] <- 1    
      SubM[[length(SubM) + 1]] <- SubMpi 
    }
  }
  SubM <- do.call(rbind, SubM)
  SubM
}