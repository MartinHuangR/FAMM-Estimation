source("Functions/CI-Funcs.R")
source("Proteomics/Proteomics-Funcs.R")
pro = read.csv("Proteomics.csv")
pro = pro |> dplyr::select(-sampleID)
X = pro |> as.matrix() |> scale()

set.seed(1)
idx = sample(1:ncol(X), 20, replace = F)
X = X[,idx]
active = 3; repeats = 1000; snr = 10; p = ncol(X)
beta = c(rep(1,active), rep(0, ncol(X) - active)) 
true = c(rep(1, active),rep(0,ncol(X) - active))
signal = sqrt(mean((as.matrix(X) %*% as.matrix(beta))^2))
sigma = as.numeric(signal/sqrt(snr))

# Compute Y with SNR
Y = as.matrix(X)%*%as.matrix(beta) + rnorm(nrow(X), 0, sd = sigma) 


outliers = c(seq(0,40, by = 10), c(42,44,46,48,50))
pwlsCIq = matrix(NA, nrow = length(outliers), ncol = p)
rownames(pwlsCIq) = outliers
pwlsCIbca = matrix(NA, nrow = length(outliers), ncol = p)
rownames(pwlsCIbca) = outliers
pwlsCIlengthq = matrix(NA, nrow = length(outliers), ncol = p)
rownames(pwlsCIlengthq) = outliers
pwlsCIlengthbca = matrix(NA, nrow = length(outliers), ncol = p)
rownames(pwlsCIlengthbca) = outliers
pwlsCItime = c()
set.seed(1)
for (i in 1:length(outliers)){
  start <- Sys.time()
  s = sbwlsCI(X,Y,beta, LOOPS = 1000, nOut = (nrow(X) * outliers[i]/100), m = nrow(X)/2, largeOut = F, na.rm = T)
  
  pwlsCIq[i,] = s[["quantile.acc"]]
  pwlsCIlengthq[i,] = s[["quantile.length"]]
  pwlsCIbca[i,] = s[["bca.acc"]]
  pwlsCIlengthbca[i,] = s[["bca.length"]]
  
  end = Sys.time()
  pwlsCItime[i] = difftime(end, start, units = "mins") |> print()
  print(paste("###########################", outliers[i], "###########################"))
}



outliers = c(seq(0,40, by = 10), c(42,44,46,48,50))
pmmCIq = matrix(NA, nrow = length(outliers), ncol = p)
rownames(pmmCIq) = outliers
pmmCIbca = matrix(NA, nrow = length(outliers), ncol = p)
rownames(pmmCIbca) = outliers
pmmCIlengthq = matrix(NA, nrow = length(outliers), ncol = p)
rownames(pmmCIlengthq) = outliers
pmmCIlengthbca = matrix(NA, nrow = length(outliers), ncol = p)
rownames(pmmCIlengthbca) = outliers
pmmCItime = c()
set.seed(1)
for (i in 1:length(outliers)){
  start <- Sys.time()
  s = sbmmCI(X,Y,beta, LOOPS = 1000, nOut = (nrow(X) * outliers[i]/100), m = nrow(X)/2, largeOut = F, na.rm = T)
  
  pmmCIq[i,] = s[["quantile.acc"]]
  pmmCIlengthq[i,] = s[["quantile.length"]]
  pmmCIbca[i,] = s[["bca.acc"]]
  pmmCIlengthbca[i,] = s[["bca.length"]]
  end = Sys.time()
  pmmCItime[i] = difftime(end, start, units = "mins") |> print()
  print(paste("###########################", outliers[i], "###########################"))
}


outliers = c(seq(0,40, by = 10), c(42,44,46,48,50))
pfrbCIq = matrix(NA, nrow = length(outliers), ncol = p)
rownames(pfrbCIq) = outliers
pfrbCIbca = matrix(NA, nrow = length(outliers), ncol = p)
rownames(pfrbCIbca) = outliers
pfrbCIlengthq = matrix(NA, nrow = length(outliers), ncol = p)
rownames(pfrbCIlengthq) = outliers
pfrbCIlengthbca = matrix(NA, nrow = length(outliers), ncol = p)
rownames(pfrbCIlengthbca) = outliers
pfrbCItime = c()
set.seed(1)
for (i in 1:length(outliers)){
  start <- Sys.time()
  s = frbCI(X,Y,beta, LOOPS = 1000, nOut = (nrow(X) * outliers[i]/100), m = nrow(X)/2, largeOut = F, na.rm = T)
  
  pfrbCIq[i,] = s[["quantile.acc"]]
  pfrbCIlengthq[i,] = s[["quantile.length"]]
  pfrbCIbca[i,] = s[["bca.acc"]]
  pfrbCIlengthbca[i,] = s[["bca.length"]]
  end = Sys.time()
  pfrbCItime[i] = difftime(end, start, units = "mins") |> print()
  print(paste("###########################", outliers[i], "###########################"))
}