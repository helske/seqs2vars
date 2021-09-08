library(TraMineR)
library(cluster)
library(foreach)
library(doParallel)

ncores <- 32
nsim <- 10000
cl <- makeCluster(ncores)
registerDoParallel(cl)


simulation <- function(
  sim_seq,
  y,
  n = 1000L, nsim = 1000L, k = 4,
  seed = sample(1:.Machine$integer.max,size=1)) {
  
  # run FANNY
  cluster_with_fanny <- function(d, memb.exp, k = 4) {
    out <- suppressWarnings(try(fanny(k = k, diss = TRUE, x = d, memb.exp = memb.exp), silent = TRUE))
    if (!inherits(out, "try-error") &&
        abs(out$coeff["normalized"]) > 1e-06 &&
        out$convergence[2] == 1 &&
        length(unique(out$clustering)) == k &&
        out$k.crisp == k &&
        qr(out$membership)$rank == k &&
        qr(head(out$membership, nrow(d) / 2))$rank == k &&
        qr(tail(out$membership, nrow(d) / 2))$rank == k) {
      return(list(out = out))
    } else "error"
  }
  
  set.seed(seed)
  
  # Compute RMSEs, repeat nsim times
  rmse <- foreach(j = 1:nsim, .combine = cbind, 
    .packages = c("cluster", "TraMineR")) %dopar% {
      
      
      rmse <- numeric(4)
      names(rmse) <-  c("Hard classification",
        "Representativeness",
        "Soft classification",
        "Pseudoclass")
      
      # take a subsample of the full data
      idx <- sample(1:nrow(sim_seq), size = n)
      seq_sample <- sim_seq[idx,]
      y_sample <- y[idx]
      # Sequence analysis (using OMspell dissimilarity measure which is sensitive to sequencing) and PAM clustering
      sm <- suppressMessages(seqcost(seq_sample, method = "CONSTANT"))
      dist_OMsp <- suppressMessages(seqdist(seq_sample, method = "OMspell", sm = sm$sm, tpow = 1.5))
      # Dissimilarities from OMspell
      d <- as.matrix(dist_OMsp)
      
      # PAM
      cluster_pam <- pam(k = k, diss = TRUE, x = dist_OMsp)
      
      # representativeness to cluster medoids (based on PAM)
      
      d_to_medoid <- 1 - d[, cluster_pam$id.med] / max(d)
      
      # FANNY, with memb.exp 1.4 (can fail)
      cluster_fanny <- cluster_with_fanny(d, 1.4, k = k)
      
      ## Modelling alternatives ##
      
      # Estimate model using first half of the sample, and predict latter half
      # note that the X is built using full sample
      
      # Traditional prediction using hard classification
      cl_max <- factor(cluster_pam$clustering)
      df <- data.frame(y = y_sample, cl_max = cl_max) # create data
      fit <- lm(y ~ ., data = head(df, n/2))
      if (fit$rank == k) {
        # Predict rest
        rmse["Hard classification"] <- sqrt(mean((tail(df$y, n/2) -
            predict(fit, newdata = tail(df, n/2)))^2))
      }
      
      # Using representativeness
      df <- data.frame(y = y_sample, x = d_to_medoid)
      fit <- lm(y ~  ., data = head(df, n/2))
      if (fit$rank == (k+1)) {
        rmse["Representativeness"] <- sqrt(mean((tail(df$y, n/2) -
            predict(fit, newdata = tail(df, n/2)))^2))
      }
      
      if (is.list(cluster_fanny)) {
        
        x <- cluster_fanny$out$membership
        df <- data.frame(y = y_sample, x = x)
        fit <- lm(y ~ -1 + .,data = head(df, n/2))
        if (fit$rank == k) {
          rmse["Soft classification"] <- sqrt(mean((tail(df$y, n/2) -
              predict(fit, newdata = tail(df, n/2)))^2))
        }
        
        # Averaging predictions 
        # This is simpler than averaging coefficients but gives equal results in the linear model case
        pseudo_preds <- matrix(NA, n/2, 50)
        df <- data.frame(y = y_sample, x = factor(NA, levels = 1:k))
        for(jj in 1:50) {
          df$x <- factor(apply(x, 1, function(p) sample(1:length(p), size = 1, prob = p)))
          if(length(unique(df$x)) == k){
            f <- lm(y ~ -1 + x, data = head(df, n/2))
            pseudo_preds[, jj] <- predict(f, newdata = tail(df, n/2))
          }
        }
        rmse["Pseudoclass"] <- sqrt(mean((tail(y_sample, n/2) - rowMeans(pseudo_preds, na.rm = TRUE))^2))
      }
      rmse
    }
  rmse
}




# Created with create_data.R
load("fuzzy_data.rds")
load("crisp_data.rds")
load("super_crisp_data.rds")

res <- simulation(sim_seq_fuzzy, y_repr_fuzzy, n = 1000, nsim = nsim, seed = 1)
saveRDS(res, file = "repr_fuzzy.rds")
res <- simulation(sim_seq_fuzzy, y_class_fuzzy, n = 1000, nsim = nsim, seed = 1)
saveRDS(res, file = "class_fuzzy.rds")
res <- simulation(sim_seq_crisp, y_repr_crisp, n = 1000, nsim = nsim, seed = 1)
saveRDS(res, file = "repr_crisp.rds")
res <- simulation(sim_seq_crisp, y_class_crisp, n = 1000, nsim = nsim, seed = 1)
saveRDS(res, file = "class_crisp.rds")
res <- simulation(sim_seq_super_crisp, y_repr_super_crisp, n = 1000, nsim = nsim, seed = 1)
saveRDS(res, file = "repr_super_crisp.rds")
res <- simulation(sim_seq_super_crisp, y_class_super_crisp, n = 1000, nsim = nsim, seed = 1)
saveRDS(res, file = "class_super_crisp.rds")

stopCluster(cl)