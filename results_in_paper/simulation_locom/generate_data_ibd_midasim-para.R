rm(list = ls() )

setwd("")
library(MIDASim)
# ml r/4.2.0
library(nleqslv)
library(survival)
load("count.ibd.setup.RData") # fitted midasim model: count.ibd.setup <- MIDASim.setup(count.ibd, mode = 'parametric')

alln = 200 # total sample size
allb = c(0, 0.5, 1, 1.5, 2)  # effect size
rnt = c(10, 20) # number of causal taxa

eg = expand.grid(n = alln, beta = allb, rnt = rnt, 
                 stringsAsFactors = F)

args = (commandArgs(trailingOnly=TRUE))
iscen = as.numeric(args[1])  # which setup from 'eg'
sim.batch = as.numeric(args[2]) # every 100 replicates is a batch

condition = eg[iscen,]

n.data = condition$n
beta1 = condition$beta # effect of the variable of interest
rnt = condition$rnt

n.taxa = count.ibd.setup$n.taxa
n.data = count.ibd.setup$n.sample

start.index = 100 * (sim.batch - 1) + 1
end.index = 100 * sim.batch

m = 50
top.index = order(count.ibd.setup$mean.rel.abund, decreasing =T)[1:m]
x.index = sample(top.index, rnt, replace =F) # get causal taxa set: randomly select rnt=10/20 taxa from the most abundant 50 taxa

# There are five of them that are the same and the rest are not the same
x2.index = c(sample(x.index, 5, replace =F),
             sample(setdiff(top.index, x.index), 5, replace =F )) # the causal set of the covariate has overlap with the causal set of the variable of interest

set.seed(NULL)
x.true = c(rep(0, n.data/2), rep(1, n.data/2))
x2 =  c(rep(0, n.data/4), rep(1, n.data/4), rep(0, n.data/4), rep(1, n.data/4))

x.true = x.true - mean(x.true)
x2 = x2 - mean(x2)


if (!dir.exists("MIDASim_generatedData")) dir.create("MIDASim_generatedData", recursive = TRUE)
for (sim in start.index:end.index){
  
  fn = sprintf("MIDASim_generatedData/simulated.ibd.midas_n_%d_b_%s_rnt_%d_sim%d.rds",
               n.data, beta1, rnt, sim)
  
  data.all = ra.all = matrix(NA, n.data, n.taxa)
  
  effect1 = x.true * beta1
  effect2 = x2 * 1  # effect of the covariate is always 1
  new.mean.rel.abund = count.ibd.setup$mean.rel.abund
  new.mean.rel.abund = outer(rep(1,n.data), new.mean.rel.abund)
  new.mean.rel.abund[, x.index] = new.mean.rel.abund[, x.index] * outer(exp(effect1), rep(1, length(x.index)))
  new.mean.rel.abund[, x2.index] = new.mean.rel.abund[, x2.index] * outer(exp(effect2), rep(1, length(x2.index)))
  new.mean.rel.abund = new.mean.rel.abund/rowSums(new.mean.rel.abund)
  
  count.ibd.modified.2 =  try(MIDASim.modify(count.ibd.setup,
                                             rel.abund = new.mean.rel.abund,
                                             x.cov = cbind(x.true, x2)))
  
  if ("try-error" %in% class(count.ibd.modified.2)){
    # out.list = list(count.ibd.setup, lib.size , new.mean.rel.abund)
    # save(out.list, file = "Midas_debug.RData")
    break()
  }else{
    temp = MIDASim(count.ibd.modified.2)
    data.all <- temp$sim_count
    ra.all<- temp$sim_rel
  }
  
  if (any(is.na(rowSums(data.all)))){
    next()
  }
  
  cens.prob = colMeans(data.all == 0)

  trueDA = rep(0, n.taxa)
  trueDA[x.index] = 1
  
  trueDA.x2 = rep(0, n.taxa)
  trueDA.x2[x2.index] = 1
  colnames(data.all) = paste0("otu", 1:ncol(data.all))
  
  midas.data = list(trueDA, data.all, cens.prob, x.true, x2, trueDA.x2)
  
  
  save(midas.data, file = fn)
}








