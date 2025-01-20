setwd("/users/mhe/midas")
args=(commandArgs(TRUE))

library(RcppArmadillo)
library(Rcpp)
# sourceCpp("LOCOM.cpp")
# source("LOCOM_fun.R")
library(LOCOM)
library(phyloseq)
library(survival)
library(permute)
library(parallel)

alln = 200 # total sample size
allb = c(0, 0.5, 1, 1.5, 2)  # effect size
rnt = c(10, 20) # number of causal taxa

eg = expand.grid(n = alln, beta = allb, rnt = rnt, 
                 stringsAsFactors = F)

args = (commandArgs(trailingOnly=TRUE))
iscen = as.numeric(args[1])  # which setup from 'eg'
sim.batch = as.numeric(args[2]) # every 100 replicates is a batch
print(iscen)
thresh = as.numeric(args[3]) # filter out taxa that exist in less than thresh(eg.thresh=0.2 means 20%) samples

condition = eg[iscen, ]

condition = eg[iscen,]

n.data = condition$n
beta1 = condition$beta # effect of the variable of interest
rnt = condition$rnt

if (!dir.exists("result_locom")) dir.create("result_locom", recursive = TRUE)
fn =sprintf("result_locom/n_%d_b_%s_rnt_%d.txt",
            n.data, beta1, rnt) 

start.index = 100 * (sim.batch - 1) + 1
end.index = 100 * sim.batch

for (i in start.index:end.index){
  
  data.fn = sprintf("MIDASim_generatedData/simulated.ibd.midas_n_%d_b_%s_rnt_%d_sim%d.rds", 
                    n.data, beta1, rnt, i)
  
  if (file.access(data.fn)!= 0 ){
    next()
  }
  
  load(data.fn)
  
  if (class(midas.data) != "list"|length(midas.data) != 6){
    next()
  }
  trueDA = midas.data[[1]]
  data.all = midas.data[[2]]
  cens.prob = midas.data[[3]]
  x.true = midas.data[[4]]
  x2 = midas.data[[5]]
  n.taxa = ncol(data.all)
  ra.all <- data.all/rowSums(data.all)
  
  data.all = ifelse(is.na(data.all), 0, data.all)
  
  if (any(rowSums(data.all) == 0)){
    flag = which(rowSums(data.all) == 0)
    data.all = data.all[-flag, ]
    x.true = x.true[-flag]
    x2 = x2[-flag]
  }
  
  lib.size.use = rowSums(data.all)
  
  res <- locom(otu.table = data.all, Y = x.true, C = x2,
               fdr.nominal = 0.2, seed = 1, n.cores = 4,
               filter.thresh = thresh)
  
  set.seed(NULL)
  
  out = data.frame(effect = drop(res$effect.size), p = drop(res$p.otu),
                   q = drop(res$q.otu))
  
  keep.otu = which(colMeans(data.all > 0) >= thresh)
  
  allout = as.data.frame(matrix(NA,  n.taxa, 4))
  allout[keep.otu, 1:3] = out
  allout$trueDA = trueDA
  colnames(allout)[1:3] = colnames(out)
  allout$skip.rare = rep(1, n.taxa)
  allout$skip.rare[keep.otu] = 0
  
  
  null.taxa = which(trueDA == 0)
  alt.taxa = which(trueDA == 1)
  fdr = allout$q
  pos.true = sum(fdr < 0.05 & allout$trueDA %in% c(1, -1), na.rm =T)
  pos.false = sum(fdr < 0.05 & allout$trueDA == 0, na.rm =T)
  pos.tot = sum(fdr < 0.05, na.rm =T)
  trueDA.tot = sum(allout$trueDA != 0 & allout$skip.rare == 0, na.rm= T)
  
  pos.true.20 = sum(fdr < 0.2 & allout$trueDA %in% c(1, -1), na.rm =T)
  pos.false.20 = sum(fdr < 0.2 & allout$trueDA == 0, na.rm =T)
  pos.tot.20 = sum(fdr < 0.2, na.rm =T)
  
  mse = mean(c((allout$effect[alt.taxa] - beta1)^2, (allout$effect[null.taxa])^2), na.rm =T)
  if (pos.tot == 0){
    fdr = 0
  }else{
    fdr = pos.false/pos.tot
  }
  
  if (pos.tot.20 == 0){
    fdr.20 = 0
  }else{
    fdr.20 = pos.false.20/pos.tot.20
  }
  
  allp = allout$p[which(allout$trueDA == 0& allout$skip.rare == 0)]
  typeI = mean(allp < 0.05, na.rm = T)
  
  uni.test = suppressWarnings(ks.test(x = allp, y = "punif")$p)
  n.taxa.included = sum(allout$skip.rare == 0)
  
  out2 = c(pos.tot, pos.true, pos.false, trueDA.tot, fdr, mse,
           pos.tot.20, pos.true.20, pos.false.20, fdr.20, n.taxa.included)
  
  sink(fn, append = T)
  cat(out2,sep = "\t")
  cat("\n")
  sink()
  
}




