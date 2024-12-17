source("Functions/MS-Funcs")
set.seed(1)

# Setting 1
n = 500; p = 5; active = 2
toeshd = 0.5^abs(row(matrix(1:p, p, p)) - col(matrix(1:p, p, p)))
X.full = mvtnorm::rmvnorm(n  = n, sigma = toeshd)
beta = rep(0,p)
beta[1:active] = 1
y.full = as.matrix(X.full) %*% as.matrix(beta)
allSubM = getSubM(p)
SubM = allSubM[c(6,9,20,31),]
ptrue = 1


outliers = c(seq(0,40, by = 10), c(42,44,46,48,50))
mm1 = c()
mm1time = c()
set.seed(1)
for (i in 1:length(outliers)){
  start <- Sys.time()
  m1 = mmAll(X.full = X.full, y = y.full, ptrue = ptrue, LOOPS = 1000, SubM = SubM, m = nrow(X.full)/2,
             nOut = (nrow(X.full) * outliers[i]/100), outliers = outliers[i])
  mm1 = rbind(mm1, m1)
  mm1time[i] = print( Sys.time() - start)
  print(paste("###########################", outliers[i], "###########################"))
}



outliers = c(seq(0,40, by = 10), c(42,44,46,48,50))
famm1 = c()
famm1time = c()
set.seed(1)
for (i in 1:length(outliers)){
  start <- Sys.time()
  f1 = fammAll(X.full = X.full, y = y.full, ptrue = ptrue, LOOPS = 1000, SubM = SubM, m = nrow(X.full)/2,
               nOut = (nrow(X.full) * outliers[i]/100), outliers = outliers[i])
  famm1 = rbind(famm1, f1)
  famm1time[i] = print( Sys.time() - start)
  print(paste("###########################", outliers[i], "###########################"))
}


numCores <- 9
cl <- makeCluster(numCores)
registerDoParallel(cl)

set.seed(1)
frb1 = foreach(i = seq_along(outliers), .combine = rbind, .packages = c("MASS", "robustbase")) %dopar% {
  m1 <- frbAll(
    X.full = X.full, 
    y = y.full, 
    ptrue = ptrue, 
    LOOPS = 1000, 
    SubM = SubM, 
    m = nrow(X.full) / 2,
    nOut = (nrow(X.full) * outliers[i] / 100), 
    outliers = outliers[i]
  )
  m1
}
stopCluster(cl)


# Setting 2
set.seed(1)
n = 100; p = 10; active = 5
X.full = matrix(rnorm(n*p), ncol = p, nrow = n, byrow = T)
beta = rep(0,p)
beta[1:active] = 1
y.full = as.matrix(X.full) %*% as.matrix(beta)
allSubM = getSubM(p)
SubM = allSubM[c(386, 299, 56, 1023),]
ptrue = 1

outliers = c(seq(0,40, by = 10), c(42,44,46,48,50))
mm2 = c()
mm2time = c()
set.seed(1)
for (i in 1:length(outliers)){
  start <- Sys.time()
  m2 = mmAll(X.full = X.full, y = y.full, ptrue = ptrue, LOOPS = 1000, SubM = SubM, m = nrow(X.full)/2,
             nOut = (nrow(X.full) * outliers[i]/100), outliers = outliers[i])
  mm2 = rbind(mm2, m2)
  mm2time[i] = print( Sys.time() - start)
  print(paste("###########################", outliers[i], "###########################"))
}


outliers = c(seq(0,40, by = 10), c(42,44,46,48,50))
famm2 = c()
famm2time = c()

set.seed(1)
for (i in 1:length(outliers)){
  start <- Sys.time()
  f2 = fammAll(X.full = X.full, y = y.full, ptrue = ptrue, LOOPS = 1000, SubM = SubM, m = nrow(X.full)/2,
               nOut = (nrow(X.full) * outliers[i]/100), outliers = outliers[i])
  famm2 = rbind(famm2, f2)
  famm2time[i] = print( Sys.time() - start)
  print(paste("###########################", outliers[i], "###########################"))
}

numCores <- 9
cl <- makeCluster(numCores)
registerDoParallel(cl)

set.seed(1)
frb2 = foreach(i = seq_along(outliers), .combine = rbind, .packages = c("MASS", "robustbase")) %dopar% {
  m1 <- frbAll(
    X.full = X.full, 
    y = y.full, 
    ptrue = ptrue, 
    LOOPS = 1000, 
    SubM = SubM, 
    m = nrow(X.full) / 2,
    nOut = (nrow(X.full) * outliers[i] / 100), 
    outliers = outliers[i]
  )
  m1
}
stopCluster(cl)
