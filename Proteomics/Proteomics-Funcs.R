library(Rcpp)
library(stabs)
library(glmnet)
library(ncvreg)
eats = function(X, Y){
  
  # Stability Selection
  s = stabs::stabsel(x = X, y = Y, B = 100,
                     fitfun = stabs::lars.lasso, PFER = 5, cutoff = 0.75,
                     sampling.type = "MB")
  
  idx = sample(1:nrow(X), replace = F)
  rX = X[idx,]
  idxPushed = c(tail(idx, 1), head(idx, -1))
  rY = Y[idxPushed] |> as.matrix(ncol = 1)
  
  sMix = stabs::stabsel(x = rX, y = rY, B = 100,
                        fitfun = stabs::lars.lasso, PFER = 5, cutoff = 0.75,
                        sampling.type = "MB")
  
  sMix_prob = sort(sMix$max, decreasing = T)
  mix_exclusion = quantile(sMix_prob, 0.95)
  
  # Exclusion ATS  
  EATS =  convert(s)[1:(length(convert(s)))][convert(s)[1:(length(convert(s)))] >= 100*mix_exclusion] 
  if (length(EATS) == 2){EATS = c(EATS, EATS[1])}
  EATS = EATS |> getR()
  EATS_selected = sort(s$max, decreasing = T)[1:EATS] |> names()
  
  EATS_selected
}


convert = function(s){
  return((as.vector(s$max) |> sort(decreasing = T))*100)
}

Rcpp::cppFunction('
int getR(const NumericVector& d) {
  int p = d.size();
  NumericVector lq(p, 0.0);
  NumericVector sigma2(p);
  for (int q = 0; q < p; q++) {
    NumericVector d1 = head(d, q + 1);
    NumericVector d2 = tail(d, p - (q + 1));
    double mu1 = mean(d1);
    double mu2 = mean(d2);
    sigma2[q] = (sum(pow(d1 - mu1, 2)) + sum(pow(d2 - mu2, 2))) / (p - 2);
    lq[q] = sum(dnorm(d1, mu1, sqrt(sigma2[q]), true)) +
      sum(dnorm(d2, mu2, sqrt(sigma2[q]), true));
  }
  return which_max(lq) + 1;
}
')

lasso = function(X,Y){
  l = cv.glmnet(X, Y)
  c1se = which(coef(l, lambda = l$lambda.1se)[-1] != 0)
  c1se_selected = colnames(X)[c1se]
  c1se_selected
}

ko = function(X,Y){
  kf = knockoff::knockoff.filter(X, Y, fdr = 5)
  ko_selected = names(kf$selected)
  ko_selected
}
