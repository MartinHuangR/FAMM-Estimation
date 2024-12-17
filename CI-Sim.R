source("Functions/CI-Funcs")
set.seed(1)
n = 500; p = 5; active = 2
toeshd = 0.5^abs(row(matrix(1:p, p, p)) - col(matrix(1:p, p, p)))
X = mvtnorm::rmvnorm(n  = n, sigma = toeshd)
beta = rep(0,p)
beta[1:active] = 1
y = as.matrix(X ) %*% as.matrix(beta)



outliers = c(seq(0,40, by = 10), c(42,44,46,48,50))
wlsCIq = matrix(NA, nrow = length(outliers), ncol = p)
rownames(wlsCIq) = outliers
wlsCIbca = matrix(NA, nrow = length(outliers), ncol = p)
rownames(wlsCIbca) = outliers
wlsCIlengthq = matrix(NA, nrow = length(outliers), ncol = p)
rownames(wlsCIlengthq) = outliers
wlsCIlengthbca = matrix(NA, nrow = length(outliers), ncol = p)
rownames(wlsCIlengthbca) = outliers
wlsCItime = c()
set.seed(1)
for (i in 1:length(outliers)){
  start <- Sys.time()
  s = sbwlsCI(X,y,beta, LOOPS = 1000, nOut = (nrow(X) * outliers[i]/100), m = nrow(X)/2)
  
  wlsCIq[i,] = s[["quantile.acc"]]
  wlsCIlengthq[i,] = s[["quantile.length"]]
  wlsCIbca[i,] = s[["bca.acc"]]
  wlsCIlengthbca[i,] = s[["bca.length"]]
  
  end = Sys.time()
  wlsCItime[i] = difftime(end, start, units = "secs") |> print()
  print(paste("###########################", outliers[i], "###########################"))
}



outliers = c(seq(0,40, by = 10), c(42,44,46,48,50))
mmCIq = matrix(NA, nrow = length(outliers), ncol = p)
rownames(mmCIq) = outliers
mmCIbca = matrix(NA, nrow = length(outliers), ncol = p)
rownames(mmCIbca) = outliers
mmCIlengthq = matrix(NA, nrow = length(outliers), ncol = p)
rownames(mmCIlengthq) = outliers
mmCIlengthbca = matrix(NA, nrow = length(outliers), ncol = p)
rownames(mmCIlengthbca) = outliers
mmCItime = c()
set.seed(1)
for (i in 1:length(outliers)){
  start <- Sys.time()
  s = sbmmCI(X,y,beta, LOOPS = 1000, nOut = (nrow(X) * outliers[i]/100), m = nrow(X)/2)
  
  mmCIq[i,] = s[["quantile.acc"]]
  mmCIlengthq[i,] = s[["quantile.length"]]
  mmCIbca[i,] = s[["bca.acc"]]
  mmCIlengthbca[i,] = s[["bca.length"]]
  end = Sys.time()
  mmCItime[i] = difftime(end, start, units = "secs") |> print()
  print(paste("###########################", outliers[i], "###########################"))
}



outliers = c(seq(0,40, by = 10), c(42,44,46,48,50))
frbCIq = matrix(NA, nrow = length(outliers), ncol = p)
rownames(frbCIq) = outliers
frbCIbca = matrix(NA, nrow = length(outliers), ncol = p)
rownames(frbCIbca) = outliers
frbCIlengthq = matrix(NA, nrow = length(outliers), ncol = p)
rownames(frbCIlengthq) = outliers
frbCIlengthbca = matrix(NA, nrow = length(outliers), ncol = p)
rownames(frbCIlengthbca) = outliers
frbCItime = c()
set.seed(1)
for (i in 1:length(outliers)){
  start <- Sys.time()
  s = frbCI(X,y,beta, LOOPS = 1000, nOut = (nrow(X) * outliers[i]/100), m = nrow(X)/2)
  
  frbCIq[i,] = s[["quantile.acc"]]
  frbCIlengthq[i,] = s[["quantile.length"]]
  frbCIbca[i,] = s[["bca.acc"]]
  frbCIlengthbca[i,] = s[["bca.length"]]
  end = Sys.time()
  frbCItime[i] = difftime(end, start, units = "secs") |> print()
  print(paste("###########################", outliers[i], "###########################"))
}


##################################    SIM 2   ################################################

# Sim 2
set.seed(1)
n = 100; p = 10; active = 5
X = matrix(rnorm(n*p), ncol = p, nrow = n, byrow = T)
beta = rep(0,p)
beta[1:active] = 1
y = as.matrix(X) %*% as.matrix(beta)


outliers = c(seq(0,40, by = 10), c(42,44,46,48,50))
wlsCIq2 = matrix(NA, nrow = length(outliers), ncol = p)
rownames(wlsCIq2) = outliers
wlsCIbca2 = matrix(NA, nrow = length(outliers), ncol = p)
rownames(wlsCIbca2) = outliers
wlsCIlengthq2 = matrix(NA, nrow = length(outliers), ncol = p)
rownames(wlsCIlengthq2) = outliers
wlsCIlengthbca2 = matrix(NA, nrow = length(outliers), ncol = p)
rownames(wlsCIlengthbca2) = outliers
wlsCItime2 = c()
set.seed(1)
for (i in 1:length(outliers)){
  start <- Sys.time()
  s = sbwlsCI(X,y,beta, LOOPS = 1000, nOut = (nrow(X) * outliers[i]/100), m = nrow(X)/2)
  
  wlsCIq2[i,] = s[["quantile.acc"]]
  wlsCIlengthq2[i,] = s[["quantile.length"]]
  wlsCIbca2[i,] = s[["bca.acc"]]
  wlsCIlengthbca2[i,] = s[["bca.length"]]
  
  end = Sys.time()
  wlsCItime2[i] = difftime(end, start, units = "secs") |> print()
  print(paste("###########################", outliers[i], "###########################"))
}



mmCIq2 = matrix(NA, nrow = length(outliers), ncol = p)
rownames(mmCIq2) = outliers
mmCIbca2 = matrix(NA, nrow = length(outliers), ncol = p)
rownames(mmCIbca2) = outliers
mmCIlengthq2 = matrix(NA, nrow = length(outliers), ncol = p)
rownames(mmCIlengthq2) = outliers
mmCIlengthbca2 = matrix(NA, nrow = length(outliers), ncol = p)
rownames(mmCIlengthbca2) = outliers
mmCItime2 = c()
set.seed(1)
for (i in 1:length(outliers)){
  start <- Sys.time()
  s = sbmmCI(X,y,beta, LOOPS = 1000, nOut = (nrow(X) * outliers[i]/100), m = nrow(X)/2)
  
  mmCIq2[i,] = s[["quantile.acc"]]
  mmCIlengthq2[i,] = s[["quantile.length"]]
  mmCIbca2[i,] = s[["bca.acc"]]
  mmCIlengthbca2[i,] = s[["bca.length"]]
  end = Sys.time()
  mmCItime2[i] = difftime(end, start, units = "secs") |> print()
  print(paste("###########################", outliers[i], "###########################"))
}



frbCIq2 = matrix(NA, nrow = length(outliers), ncol = p)
rownames(frbCIq2) = outliers
frbCIbca2 = matrix(NA, nrow = length(outliers), ncol = p)
rownames(frbCIbca2) = outliers
frbCIlengthq2 = matrix(NA, nrow = length(outliers), ncol = p)
rownames(frbCIlengthq2) = outliers
frbCIlengthbca2 = matrix(NA, nrow = length(outliers), ncol = p)
rownames(frbCIlengthbca2) = outliers
frbCItime2 = c()
set.seed(1)
for (i in 1:length(outliers)){
  start <- Sys.time()
  s = frbCI(X,y,beta, LOOPS = 1000, nOut = (nrow(X) * outliers[i]/100), m = nrow(X)/2)
  
  frbCIq2[i,] = s[["quantile.acc"]]
  frbCIlengthq2[i,] = s[["quantile.length"]]
  frbCIbca2[i,] = s[["bca.acc"]]
  frbCIlengthbca2[i,] = s[["bca.length"]]
  end = Sys.time()
  frbCItime2[i] = difftime(end, start, units = "secs") |> print()
  print(paste("###########################", outliers[i], "###########################"))
}



