# install.packages("SMME")
library("SMME")

ssize <- c(10000, 40000, 60000)
p <- 50000
zeta = c(5, 10)
## read genotype files
vcf_afr <- file.path("simu_afr_chr4_0.vcf.gz")
vcf_eas <- file.path("simu_eas_chr4_0.vcf.gz")
vcf_eur <- file.path("simu_eur_chr4_0.vcf.gz")
pop <- c("afr", "eas", "eur")
vcf_test <- paste0("simu_test_", pop, "_chr4_0")

vcf_vec <- c(vcf_afr, vcf_eas, vcf_eur)

x <- list()
for (i in 1:length(pop)) {
  vcf <- sim1000G::readVCF(vcf_vec[i], maxNumberOfVariants = p, min_maf = 0 ,max_maf = NA)
  x[[pop[i]]] = t(vcf$gt1 + vcf$gt2)
  remove(vcf)
}

y <- list()
for (i in pop) {
  y[[i]] = read.table(paste0("simu_", i, "_chr4_0.txt"), quote="\"", comment.char="")$V3
}

##fit model for range of lambda and zeta
system.time(fit <- softmaximin(x, y, zeta = zeta, penalty = "lasso", alg = "npg", nlambda = 5, nthreads = 10))
# Multithreading enabled using 10 threads
# 用户    系统    流逝 
# 1454.09 1217.09 1636.80 
save(fit, file = "fit.RData")
save(x, file = "x.RData")
save(y, file = "y.RData")
betahat <- fit$coef
save(betahat, file = "betahat.RData")


eval <- function(pred, ytest, family){
  if (family == 'binomial') {
    ## type.measure="class"
    y = as.factor(ytest)
    ntab = table(y)
    nc = as.integer(length(ntab))
    y = diag(nc)[as.numeric(y), , drop=FALSE]
    predmat=1 / (1+exp(-pred))
    class = y[, 1] * (predmat > 0.5) + y[, 2] * (predmat <= 0.5)
    return(sum(class)/length(class))
  }
  if (family == 'gaussian') {
    ## type.measure="mse"
    mse = mean((ytest - pred)^2)
    rsq = (cor(ytest, pred))^2
    return(list(mse = mse, nrsq = -rsq))
  }
}

cross_validation <- function(x, y, betahat, family = "gaussian"){
  metric = list()
  met_all = c()
  lam_idx = c()
  if(typeof(betahat) != "list"){
    tmpLis = list()
    tmpLis[[1]] = betahat
    betahat = tmpLis
  }
  for (z_idx in 1:length(betahat)) {
    temp = c()
    for (lam in 1:dim(betahat[[z_idx]])[2]) {
      vg = c()
      for (i in 1:length(y)) {
        beta <- betahat[[z_idx]][, lam]
        if(sum(beta)==0){
          vg = c(vg, NA)
        }else{
          pred <- x[[i]] %*% beta
          vg = c(vg, eval(pred, y[[i]], family = family))$nrsq
        }
      }
      vg = ifelse(is.na(vg), 0, vg)
      vg = min(vg, na.rm = T)
      temp = c(temp, vg)
    }
    met_all = c(met_all, min(temp, na.rm = T))
    lam_idx = c(lam_idx, which(temp==met_all[z_idx])[1])
    metric[[z_idx]] = temp
  }
  z_idx =  which(met_all == min(met_all, na.rm = T))
  best_m = -metric[[z_idx]][lam_idx[z_idx]]
  return(list(z_idx = z_idx, lam_idx = lam_idx, metric = metric, best_m = best_m))
}

res = cross_validation(x, y, betahat)
z = fit$zeta[res$z_idx]
if(typeof(fit$lambda) == "list"){
  lamda_cv = fit$lambda[[res$z_idx]][res$lam_idx[res$z_idx]]
}else{
  lamda_cv = fit$lambda[res$lam_idx[res$z_idx]]
}
if(typeof(fit$coef) == "list"){
  beta = fit$coef[[res$z_idx]][, res$lam_idx[res$z_idx]]
}else{
  beta = fit$coef[, res$lam_idx[res$z_idx]]
}

write.table(beta, paste0("maximin_beta.txt"), row.names = F, col.names = F, quote = F)
#######################################################################
### evaluation
## metrics for selected beta
vg = c()
for (i in 1:length(y)) {
  pred <- x[[i]] %*% beta
  vg = c(vg, (cor(pred, y[[i]])^2))
  vg = ifelse(is.na(vg), 0, vg)
}
write.table(vg, paste0("rsq_ins.txt"), append = T, col.names = F)

## metrics for all beta
for (dx in 1:length(res[["metric"]][[res$z_idx]])) {
  vg = c()
  beta = fit$coef[[res$z_idx]][, dx]
  if(sum(beta, na.rm = T)==0){
    print("all zero, skip")
  }else{
    for (i in 1:length(y)) {
      pred <- x[[i]] %*% beta
      vg = c(vg, (cor(pred, y[[i]])^2))
      vg = ifelse(is.na(vg), 0, vg)    
  }
  print(vg)
  }
}

## metrics in testing afr data
outs_metric = c()
for(vcf_nm in vcf_test){
  vcf <- sim1000G::readVCF(paste0(vcf_nm, ".vcf.gz"), maxNumberOfVariants = p, min_maf = 0 ,max_maf = NA)
  x_test = t(vcf$gt1 + vcf$gt2)
  remove(vcf)
  y_test = read.table(paste0(vcf_nm, ".txt"), quote="\"", comment.char="")$V3
  if(typeof(fit$coef) == "list"){
    pred <- x_test %*% fit$coef[[res$z_idx]][, res$lam_idx[res$z_idx]]
  }else{
    pred <- x_test %*% fit$coef[, res$lam_idx[res$z_idx]]
  }
  outs_metric = c(outs_metric, cor(pred, y_test)^2)
  print(cor(pred, y_test)^2)
}
write.table(outs_metric, paste0("rsq_outs.txt"), append = T, col.names = F)

