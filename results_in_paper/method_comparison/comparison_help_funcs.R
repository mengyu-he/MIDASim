normalize.rel = function(x) {
  lib.size = rowSums(x)
  lib.size[lib.size == 0] = 1 # to prevent Inf, won't change data
  x <- x/lib.size
  return(x)
}

draw_pcoa_compare = function(dat1, dat2, color, shape, title, xy = c(1, 1), distance){
  tmp = rbind(dat1, dat2)
  rownames(tmp) = 1:nrow(tmp)
  OTU <- otu_table( tmp, taxa_are_rows = F)
  SAM = sample_data( data.frame(
    name = rownames(tmp),
    type = rep(c("Real","Sim"), 
               each = nrow(tmp)/2) ) )
  sample_names(SAM) <- SAM$name
  dt = phyloseq(OTU,SAM)
  
  if (distance == "jaccard") {
    ord1.jc = ordinate(dt, "PCoA", "jaccard", binary = T)
  } else {
    ord1.jc = ordinate(dt, "PCoA", "bray")
  }
  
  g = plot_ordination(dt, ord1.jc, type="samples",
                      color = "type", shape = "type")+
    scale_color_manual(values=c(color[1],color[2]) )+
    scale_shape_manual(values=c(shape[1], shape[2]))+
    xlab("PCoA 1")+
    ylab("PCoA 2")+
    theme_minimal()+ 
    theme(legend.position = "none")+
    ggtitle(title)
  
  if(xy[1] == 0) {
    g = g + theme(axis.title.x = element_blank())
  } 
  if(xy[2] == 0) {
    g = g + theme(axis.title.y = element_blank())
  } 
  
  return(g)
}


draw_tsne_compare = function(dat1, dat2, color, shape, title, xy = c(1, 1), distance){
  tmp <- rbind(dat1, dat2)
  rownames(tmp) <- 1:nrow(tmp)
  OTU <- otu_table(tmp, taxa_are_rows = F)
  SAM <- sample_data( data.frame(
    name = rownames(tmp),
    type = rep(c("Real","Sim"), 
               each = nrow(tmp)/2) ) )
  sample_names(SAM) <- SAM$name
  dt <- phyloseq(OTU,SAM)
  
  if (distance == "jaccard") {
    dist.jac <- phyloseq::distance(dt, method = "jaccard", binary = T)
  } else {
    dist.jac <- phyloseq::distance(dt, method = "bray")
  }
  
  tsne.res <- Rtsne(as.matrix(dist.jac), is_distance = T)
  
  df_tsne <- data.frame(
    Method = SAM$type,
    X = tsne.res$Y[, 1],
    Y = tsne.res$Y[, 2]
  )
  
  g <- ggplot(data = df_tsne, aes(x = X, y = Y))+
    geom_point(aes(color = Method, shape = Method))+
    scale_color_manual(values=c(color[1],color[2]) )+
    scale_shape_manual(values=c(shape[1], shape[2]))+
    xlab("t-SNE 1")+
    ylab("t-SNE 2")+
    theme_minimal()+ 
    theme(legend.position = "none")+
    ggtitle(title)
  
  if(xy[1] == 0) {
    g <- g + theme(axis.title.x = element_blank())
  } 
  if(xy[2] == 0) {
    g <- g + theme(axis.title.y = element_blank())
  } 
  
  return(g)
}


draw_umap_compare = function(dat1, dat2, color, shape, title, xy = c(1, 1), distance){
  tmp = rbind(dat1, dat2)
  rownames(tmp) = 1:nrow(tmp)
  OTU <- otu_table( tmp, taxa_are_rows = F)
  SAM = sample_data( data.frame(
    name = rownames(tmp),
    type = rep(c("Real","Sim"), 
               each = nrow(tmp)/2) ) )
  sample_names(SAM) <- SAM$name
  dt = phyloseq(OTU,SAM)
  
  if (distance == "jaccard") {
    dist.jac = phyloseq::distance(dt, method = "jaccard", binary = T)
  } else {
    dist.jac = phyloseq::distance(dt, method = "bray")
  }
  
  umap.res <- umap(as.matrix(dist.jac))
  
  df_umap <- data.frame(
    Method = SAM$type,
    X = umap.res$layout[, 1],
    Y = umap.res$layout[, 2]
  )
  
  g = ggplot(data = df_umap, aes(x = X, y = Y))+
    geom_point(aes(color = Method, shape = Method))+
    scale_color_manual(values=c(color[1],color[2]) )+
    scale_shape_manual(values=c(shape[1], shape[2]))+
    xlab("UMAP 1")+
    ylab("UMAP 2")+
    theme_minimal()+ 
    theme(legend.position = "none")+
    ggtitle(title)
  
  if(xy[1] == 0) {
    g = g + theme(axis.title.x = element_blank())
  } 
  if(xy[2] == 0) {
    g = g + theme(axis.title.y = element_blank())
  } 
  
  return(g)
}


beta_alpha_tests = function(otu.tab, simulated) {
  
  rel.tab = normalize.rel(otu.tab)
  
  p.jaccard.ibd = NULL
  p.bc.ibd = NULL
  p.t.observed = p.ks.observed = p.t.shannon = p.ks.shannon = NULL
  for (i in 1:(nrow(simulated)/nrow(otu.tab))) {
    print(paste0("sim: ", i))
    
    ind = ((i-1)*nrow(otu.tab)+1): (i*nrow(otu.tab))
    tmp = rbind(rel.tab, simulated[ind,  ])
    rownames(tmp) = 1:dim(tmp)[1]
    OTU <- otu_table( tmp, taxa_are_rows = F)
    SAM = sample_data( data.frame(name = rownames(tmp), 
                                  type = c(rep("Original", nrow(rel.tab)), 
                                           rep("Simulated", nrow(rel.tab))) ))
    sample_names(SAM) <- SAM$name
    dt = phyloseq(OTU,SAM)
    
    
    D.jaccard = vegdist(dt@otu_table@.Data, method="jaccard", binary = T)
    p.jaccard.ibd = c(p.jaccard.ibd, 
                      adonis2( as.matrix(D.jaccard) ~ dt@sam_data$type , permutations=1000)$`Pr(>F)`[1] )
    
    
    D.bc = vegdist(dt@otu_table@.Data, method="bray")
    p.bc.ibd = c(p.bc.ibd, 
                 adonis2( as.matrix(D.bc) ~ dt@sam_data$type , permutations=1000)$`Pr(>F)`[1] )
    
    lib.size = rowSums(otu.tab)
    tmp = rbind(otu.tab, lib.size * simulated[ind,  ])
    rownames(tmp) = 1:dim(tmp)[1]
    OTU <- otu_table( tmp, taxa_are_rows = F)
    SAM = sample_data( data.frame(name = rownames(tmp), 
                                  type = c(rep("Original", nrow(rel.tab)), 
                                           rep("Simulated", nrow(rel.tab))) ))
    sample_names(SAM) <- SAM$name
    dt = phyloseq(OTU,SAM)
    alpha.sub = microbiome::alpha(dt, index = c("observed", "diversity_shannon"))
    alpha.sub$Method = factor(dt@sam_data$type, 
                              levels = c("Original","MIDAS"))
    
    p.t.observed = c(p.t.observed, as.numeric(t.test( alpha.sub[1:nrow(alpha.sub)/2, "observed"], 
                                                      alpha.sub[(nrow(alpha.sub)/2+1):nrow(alpha.sub), "observed"])$p.value))
    p.ks.observed = c(p.ks.observed, as.numeric(ks.boot( alpha.sub[1:(nrow(alpha.sub)/2), "observed"], 
                                                         alpha.sub[(nrow(alpha.sub)/2+1):nrow(alpha.sub), "observed"], 
                                                         nboots = 10^4)$ks.boot.pvalue))
    
    p.t.shannon = c(p.t.shannon, as.numeric(t.test( alpha.sub[1:nrow(alpha.sub)/2, "diversity_shannon"], 
                                                    alpha.sub[(nrow(alpha.sub)/2+1):nrow(alpha.sub), "diversity_shannon"])$p.value))
    p.ks.shannon = c(p.ks.shannon, as.numeric(ks.boot( alpha.sub[1:(nrow(alpha.sub)/2), "diversity_shannon"], 
                                                       alpha.sub[(nrow(alpha.sub)/2+1):nrow(alpha.sub), "diversity_shannon"], 
                                                       nboots = 10^4)$ks.boot.pvalue))
    
  }
  
  print("p.jaccard: ")
  print(mean(p.jaccard.ibd))
  
  print("p.bc: ")
  print(mean(p.bc.ibd))
  
  print("p.t.observed: ")
  print(mean(p.t.observed))
  
  print("p.ks.observed: ")
  print(mean(p.ks.observed))
  
  print("p.t.shannon: ")
  print(mean(p.t.shannon))
  
  print("p.ks.shannon: ")
  print(mean(p.ks.shannon))
  
}


get_ks_p = function(otu.tab, simulated) {
  
  rel.tab <- normalize.rel(otu.tab)
  n.sam <- nrow(otu.tab)
  
  p.ks.jac <- p.ks.bc <- NULL
  for (i in 1:(nrow(simulated)/n.sam) ) {
    ind <- ((i-1)*nrow(otu.tab)+1): (i*nrow(otu.tab))
    tmp <- rbind(rel.tab, simulated[ind,  ])
    rownames(tmp) = 1:dim(tmp)[1]
    OTU <- otu_table( tmp, taxa_are_rows = F)
    SAM <- sample_data( data.frame(name = rownames(tmp), 
                                   type = c(rep("Original", nrow(rel.tab)), 
                                            rep("Simulated", nrow(rel.tab))) ))
    sample_names(SAM) <- SAM$name
    dt <- phyloseq(OTU,SAM)
    D.jaccard.ibd <- vegdist(dt@otu_table@.Data, 
                             method="jaccard", binary = T)
    D.bc.ibd <- vegdist(dt@otu_table@.Data, method="bray")
    bd.jaccard.ibd <- betadisper(D.jaccard.ibd,
                                 as.factor(dt@sam_data$type), 
                                 type = c("centroid"), 
                                 bias.adjust = FALSE)
    bd.bc.ibd <- betadisper(D.bc.ibd, 
                            as.factor(dt@sam_data$type),
                            type = c("centroid"), 
                            bias.adjust = FALSE)
    
    group1 <- 1:n.sam
    group2 <- (n.sam+1):(2*n.sam)
    p.ks.jac <- c(p.ks.jac, ks.test(bd.jaccard.ibd$distances[group1],
                                    bd.jaccard.ibd$distances[group2])$p.value)
    p.ks.bc <- c(p.ks.bc, ks.test(bd.bc.ibd$distances[group1],
                                  bd.bc.ibd$distances[group2])$p.value)
  }
  
  return(list(p.ks.jac = p.ks.jac,
              p.ks.bc = p.ks.bc))
}


# This function from MetaSPARSim is debugged
estimate_variability <-
  function (data,
            variability_func = "dispersion",
            perc_not_zeros = 0.2) {
    variability <- array(NA, dim = nrow(data))
    
    # variability as edgeR dispersion of RAW counts in data
    if (variability_func == "dispersion") {
      # compute percentage of entry > 0 for each OTU
      perc_not_zeros_per_OTU <- rowSums(data > 0) / ncol(data)
      
      # Indices of OTUs passing the filter
      ind_pass_filter <- (perc_not_zeros_per_OTU > perc_not_zeros)
      
      # number of OTUs passing the filter
      N_OTU_pass_filter <- sum(ind_pass_filter)
      message(
        N_OTU_pass_filter,
        " OTUs have values >0 in more than ",
        perc_not_zeros * 100,
        "% of samples"
      )
      
      if (N_OTU_pass_filter > 1) {
        f <- edgeR::DGEList(as.matrix(data[ind_pass_filter, ]))
        f <- edgeR::estimateGLMCommonDisp(f)
        f <- edgeR::estimateGLMTrendedDisp(f)
        f$trended.dispersion = rep(f$common.dispersion,length(f$AveLogCPM)) # (debugged!)
        f <- edgeR::estimateGLMTagwiseDisp(f)
        
        variability[ind_pass_filter] <- f$tagwise.dispersion
      }
    }
    
    # variability as variance of NORMALIZED counts in data
    if (variability_func == "variance") {
      variability<- apply(data, 1, var)
    }
    return(variability)
    
  }



# estimate_parameter_from_data() remains unchanged, the update is in estimate_variability()
estimate_parameter_from_data <-
  function (raw_data,
            norm_data,
            conditions,
            intensity_func = "mean",
            keep_zeros = TRUE,
            perc_not_zeros = 0.2) {
    # Number of samples group
    N_cond <- length(conditions)
    
    # Number of OTUs
    N_feature <- nrow(raw_data)
    
    dataset_parameter <- list()
    if (any(lengths(conditions) <= 1)) {
      stop("The conditions should have at least two samples per group.")
    }
    
    cond_i <- 1
    for (cond in conditions) {
      cond_parameter <- list()
      # estimate intesity
      cond_parameter$intensity <-
        estimate_intensity (data = norm_data[, cond, drop = FALSE],
                            aggregate_func = intensity_func,
                            keep_zeros = keep_zeros)
      
      # estimate variability
      cond_parameter$variability <-
        estimate_variability (
          data = raw_data[, cond, drop = FALSE],
          variability_func = "dispersion",
          perc_not_zeros = perc_not_zeros
        )
      
      # estimate library size
      cond_parameter$lib_size <-
        estimate_library_size(data = raw_data[, cond, drop = FALSE])
      
      dataset_parameter[[cond_i]] <- cond_parameter
      cond_i <- cond_i + 1
      
    }
    
    
    names(dataset_parameter) <- names(conditions)
    
    return(dataset_parameter)
  }




