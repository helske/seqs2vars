library(TraMineR)
library(cluster)
library(foreach)
library(doParallel)

ncores <- 32
cl <- makeCluster(ncores)
registerDoParallel(cl)


simulation <- function(
  sim_seq,
  y,
  n = 1000L, 
  nsim = 1000L, 
  k = 4,
  seed = sample.int(.Machine$integer.max, size = 1)) {
  
  # run FANNY
  cluster_with_fanny <- function(d, memb.exp, k = 4) {
    out <- suppressWarnings(
      try(fanny(k = k, diss = TRUE, x = d, memb.exp = memb.exp), silent = TRUE))
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
      
      
      rmse <- rep(NA, 6)
      names(rmse) <-  c("Hard classification (PAM)",
        "Hard classification (FANNY)",
        "Soft classification",
        "Representativeness",
        "Gravity center",
        "Pseudoclass")
      
      # take a subsample of the full data
      idx <- sample.int(nrow(sim_seq), size = n)
      seq_sample <- sim_seq[idx,]
      y_sample <- y[idx]
      # Sequence analysis 
      # (using OMspell dissimilarity measure which is sensitive to sequencing) 
      sm <- suppressMessages(seqcost(seq_sample, method = "CONSTANT"))
      dist_OMsp <- suppressMessages(
        seqdist(seq_sample, method = "OMspell", sm = sm$sm, tpow = 1.5))
      # Dissimilarities from OMspell
      d <- as.matrix(dist_OMsp)
      
      # PAM
      cluster_pam <- pam(k = k, diss = TRUE, x = dist_OMsp)
      
      # representativeness to cluster medoids (based on PAM)
      
      # scaled distance to medoid
      d_to_medoid <- 1 - d[, cluster_pam$id.med] / max(d)
      
      # distance to cluster's gravity center (can be negative!)
      d_to_gc <- 
        as.matrix(disscenter(d, group = cluster_pam$clustering, allcenter = TRUE))
      
      # FANNY, with memb.exp 1.4 (can fail with large values)
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
        rmse["Hard classification (PAM)"] <- sqrt(mean((tail(df$y, n/2) -
            predict(fit, newdata = tail(df, n/2)))^2))
      }
      
      # Using representativeness
      df <- data.frame(y = y_sample, x = d_to_medoid)
      fit <- lm(y ~  ., data = head(df, n/2))
      if (fit$rank == (k+1)) {
        rmse["Representativeness"] <- sqrt(mean((tail(df$y, n/2) -
            predict(fit, newdata = tail(df, n/2)))^2))
      }
      
      # Using gravity centers
      df <- data.frame(y = y_sample, x = d_to_gc)
      fit <- lm(y ~  ., data = head(df, n/2))
      if (fit$rank == (k+1)) {
        rmse["Gravity center"] <- sqrt(mean((tail(df$y, n/2) -
            predict(fit, newdata = tail(df, n/2)))^2))
      }
      # FANNY based approaches
      if (is.list(cluster_fanny)) {
        
        # Membership probabilities
        x <- cluster_fanny$out$membership
        df <- data.frame(y = y_sample, x = x)
        fit <- lm(y ~ -1 + .,data = head(df, n/2))
        if (fit$rank == k) {
          rmse["Soft classification"] <- sqrt(mean((tail(df$y, n/2) -
              predict(fit, newdata = tail(df, n/2)))^2))
        }
        
        # Hard classification using FANNY
        cl_max <- factor(cluster_fanny$out$clustering)
        df <- data.frame(y = y_sample, cl_max = cl_max)
        fit <- lm(y ~ ., data = head(df, n/2))
        if (fit$rank == k) {
          # Predict rest
          rmse["Hard classification (FANNY)"] <- sqrt(mean((tail(df$y, n/2) -
              predict(fit, newdata = tail(df, n/2)))^2))
        }
        # Pseudoclass
        # Averaging predictions, this is simpler than averaging coefficients
        # but gives equal results in the linear model case
        pseudo_preds <- matrix(NA, n/2, 50)
        df <- data.frame(y = y_sample, x = factor(NA, levels = 1:k))
        for(jj in 1:50) {
          df$x <- factor(apply(x, 1, function(p) sample(1:length(p), size = 1, prob = p)))
          if(length(unique(df$x)) == k){
            f <- lm(y ~ -1 + x, data = head(df, n/2))
            pseudo_preds[, jj] <- predict(f, newdata = tail(df, n/2))
          }
        }
        rmse["Pseudoclass"] <- 
          sqrt(mean((tail(y_sample, n/2) - rowMeans(pseudo_preds, na.rm = TRUE))^2))
      }
      
      rmse
    }
  rmse
}



nsim <- 10000
set.seed(1)
# use the same noise for all cases
noise <- 0.25 * rnorm(1e4)
n <- c(200, 1000) # sample size
y_type <- c("class_pam", "class_fanny", "prob", "repr", "gc")
cluster_strength <- c("fuzzy", "crisp", "super_crisp")

opts <- expand.grid(n = n, y = y_type, strength = cluster_strength)
for (i in seq_len(nrow(opts))) {
  # data created with create_data.R
  load(paste0(opts$strength[i], "_data.rds"))
  y <- get(paste0("y_", opts$y[i])) + noise
  res <- simulation(sim_seq, y, n = opts$n[i], nsim = nsim, seed = 1)
  saveRDS(res, 
    paste0(opts$y[i], "_", opts$strength[i], "_", opts$n[i], ".rds"))
}

## Same with "Null-model", i.e. response variable y does not depend on sequences

y <- rnorm(1e4)

# Created with create_data.R
load("fuzzy_data.rds")
res <- simulation(sim_seq, y, n = 1000, nsim = nsim, seed = 1)
saveRDS(res, file = "fuzzy_1000_noeffect.rds")

load("crisp_data.rds")
res <- simulation(sim_seq, y, n = 1000, nsim = nsim, seed = 1)
saveRDS(res, file = "crisp_1000_noeffect.rds")

load("super_crisp_data.rds")
res <- simulation(sim_seq, y, n = 1000, nsim = nsim, seed = 1)
saveRDS(res, file = "super_crisp_1000_noeffect.rds")

stopCluster(cl)