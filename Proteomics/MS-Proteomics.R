source("Functions/MS-Funcs.R")
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

# eats_sel = eats(X,Y) 
eats_sel = c("IL32", "TANK", "CCL24")

# lasso_sel = lasso(X,Y)
lasso_sel = c("IL32", "TANK", "CCL24", "IKBKG", "ARG2", "PLA2G4A", "TEK", "MAZ", "PLAUR")

# ko_sel = ko(X,Y)
ko_sel = c("IL32", "TANK", "CCL24", "POLB", "ARG2", 
           "COL5A1", "MLN", "PLA2G4A", "TEK", "PSMC3", 
           "HSP90AB1", "TPMT", "MAZ", "PIK3R1", "PLAUR", 
           "HNMT")

SubM = matrix(0, nrow = 3, ncol = ncol(X))
SubM[1, which(colnames(X) %in% eats_sel)] = 1 # EATS
SubM[2, which(colnames(X) %in% lasso_sel)] = 1 # LASSO
SubM[3, which(colnames(X) %in% ko_sel)] = 1 # Knockoff
ptrue = 1


outliers = c(seq(0,40, by = 10), c(42,44,46,48,50))
pfamm = c()
pfammtime = c()

set.seed(1)
for (i in 1:length(outliers)){
  start <- Sys.time()
  pf = fammAll(X.full = X, y = Y, ptrue = ptrue, LOOPS = 1000, SubM = SubM, m = nrow(X)/2,
               nOut = (nrow(X) * outliers[i]/100), outliers = outliers[i], largeOut = F,
               proteomics = T)
  pfamm = rbind(pfamm, pf)
  end = Sys.time()
  pfammtime[i] = difftime(end, start, units = "mins") |> print()
  print(paste("###########################", outliers[i], "###########################"))
}


outliers = c(seq(0,40, by = 10), c(42,44,46,48,50))
pmm = c()
pmmtime = c()

set.seed(1)
for (i in 1:length(outliers)){
  start <- Sys.time()
  pm = mmAll(X.full = X, y = Y, ptrue = ptrue, LOOPS = 1000, SubM = SubM, m = nrow(X)/2,
             nOut = (nrow(X) * outliers[i]/100), outliers = outliers[i], largeOut = F,
             proteomics = T)
  pmm = rbind(pmm, pm)
  end = Sys.time()
  pmmtime[i] = difftime(end, start, units = "mins") |> print()
  print(paste("###########################", outliers[i], "###########################"))
}


outliers = c(seq(0,40, by = 10), c(42,44,46,48,50))
pfrb = c()
pfrbtime = c()

set.seed(1)
for (i in 1:length(outliers)){
  start <- Sys.time()
  pfr = frbAll(X.full = X, y = Y, ptrue = ptrue, LOOPS = 1000, SubM = SubM, m = nrow(X)/2,
               nOut = (nrow(X) * outliers[i]/100), outliers = outliers[i], largeOut = F,
               proteomics = T)
  pfrb = rbind(pfrb, pfr)
  end = Sys.time()
  pfrbtime[i] = difftime(end, start, units = "mins") |> print()
  print(paste("###########################", outliers[i], "###########################"))
}
