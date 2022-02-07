library(seqHMM)
library(TraMineR)
library(cluster)

# Simulate "population" data
# Observations y do not contain noise at this point, add that later

set.seed(1)
# create MMM
# Everyone starts from state 1
init1 <- c(1, 0, 0, 0)
init <- list(init1, init1, init1, init1)
emiss1 <- diag(1, nrow = 4)
emiss <- list(emiss1, emiss1, emiss1, emiss1)

# regression coefficients

beta <- c(0, 1, 1, -1)
# Fuzzy case

trans1 <- matrix(c(0.7, 0.2, 0.1, 0,
  0, 0.93, 0.04, 0.03,
  0, 0.1, 0.87, 0.03,
  0, 0.15, 0.05, 0.8), nrow = 4, byrow = TRUE)
trans2 <- matrix(c(0.4, 0.3, 0.25, 0.05,
  0, 0.7, 0.17, 0.13,
  0, 0.05, 0.85, 0.1,
  0, 0.1, 0.05, 0.85), nrow = 4, byrow = TRUE)
trans3 <- matrix(c(0.7, 0.3, 0, 0,
  0, 0.83, 0.04, 0.13,
  0, 0, 0.9, 0.1,
  0, 0.07, 0.03, 0.9), nrow = 4, byrow = TRUE)
trans4 <- matrix(c(0.85, 0.1, 0.05, 0,
  0, 0.9, 0.03, 0.07,
  0, 0.03, 0.9, 0.07,
  0, 0.01, 0.01, 0.98), nrow = 4, byrow = TRUE)
transitions <- list(trans1, trans2, trans3, trans4)

sim_seq <- suppressMessages(simulate_mhmm(n_sequences = 1e4, 
  initial_probs = init, 
  transition_probs = transitions, 
  emission_probs = emiss,
  sequence_length = 20)$observations)

sm <- suppressMessages(seqcost(sim_seq, method = "CONSTANT"))
# Dissimilarities from OMspell
dist_OMsp <- suppressMessages(seqdist(sim_seq, method = "OMspell", sm = sm$sm, 
  tpow = 1.5))

cluster_pam <- pam(k = 4, diss = TRUE,  x = dist_OMsp)
cluster_fanny <- fanny(k = 4, diss = TRUE,  x = dist_OMsp, memb.exp =  1.4)

# PAM is used for representativeness
d <- as.matrix(dist_OMsp)
# Dissimilarity to cluster medoids
d_to_medoid <- d[, cluster_pam$id.med]

X_repr <- 1 - d_to_medoid /  max(d)
X_class_pam <- model.matrix(~ -1 + factor(cluster_pam$clustering))
X_class_fanny <- model.matrix(~ -1 + factor(cluster_fanny$clustering))
X_prob <- cluster_fanny$membership
X_gc <- as.matrix(disscenter(d, group = cluster_pam$clustering, allcenter = TRUE))

# scale to make the effects and RMSEs comparable
X_gc <- X_gc / max(abs(X_gc))

y_repr <- X_repr %*% beta 
y_class_pam <- X_class_pam %*% beta 
y_class_fanny <- X_class_fanny %*% beta 
y_prob <- X_prob %*% beta 
y_gc <- X_gc %*% beta 

save(sim_seq, y_repr, y_class_pam, y_class_fanny, y_prob, y_gc, 
  file = "fuzzy_data.rds")

####################################################################

set.seed(1)
# Crisp case

trans1 <- matrix(c(0.8, 0.2, 0, 0,
  0, 0.95, 0.04, 0.01,
  0, 0, 0.995, 0.005,
  0, 0.15, 0.05, 0.8), nrow = 4, byrow = TRUE)
trans2 <- matrix(c(0.5, 0.45, 0.05, 0,
  0, 0.8, 0.19, 0.01,
  0, 0.00, 0.99, 0.01,
  0, 0.01, 0.01, 0.98), nrow = 4, byrow = TRUE)
trans3 <- matrix(c(0.6, 0.35, 0.05, 0,
  0, 0.75, 0.07, 0.18,
  0, 0.00, 0.82, 0.18,
  0, 0.02, 0.01, 0.97), nrow = 4, byrow = TRUE)
trans4 <- matrix(c(0.999, 0.0005, 0.0005, 0,
  0, 0.9, 0.03, 0.07,
  0, 0.03, 0.9, 0.07,
  0, 0.01, 0.01, 0.98), nrow = 4, byrow = TRUE)
transitions <- list(trans1, trans2, trans3, trans4)

sim_seq <- suppressMessages(simulate_mhmm(n_sequences = 1e4, 
  initial_probs = init, 
  transition_probs = transitions, 
  emission_probs = emiss,
  sequence_length = 20)$observations)

sm <- suppressMessages(seqcost(sim_seq, method = "CONSTANT"))
# Dissimilarities from OMspell
dist_OMsp <- suppressMessages(seqdist(sim_seq, method = "OMspell", sm = sm$sm, 
  tpow = 1.5))

cluster_pam <- pam(k = 4, diss = TRUE,  x = dist_OMsp)
cluster_fanny <- fanny(k = 4, diss = TRUE,  x = dist_OMsp, memb.exp =  1.4)

# PAM is used for representativeness
d <- as.matrix(dist_OMsp)
# Dissimilarity to cluster medoids
d_to_medoid <- d[, cluster_pam$id.med]

X_repr <- 1 - d_to_medoid /  max(d)
X_class_pam <- model.matrix(~ -1 + factor(cluster_pam$clustering))
X_class_fanny <- model.matrix(~ -1 + factor(cluster_fanny$clustering))
X_prob <- cluster_fanny$membership
X_gc <- as.matrix(disscenter(d, group = cluster_pam$clustering, allcenter = TRUE))

# scale to make the effects and RMSEs comparable
X_gc <- X_gc / max(abs(X_gc))

y_repr <- X_repr %*% beta 
y_class_pam <- X_class_pam %*% beta 
y_class_fanny <- X_class_fanny %*% beta 
y_prob <- X_prob %*% beta 
y_gc <- X_gc %*% beta

save(sim_seq, y_repr, y_class_pam, y_class_fanny, y_prob, y_gc, 
  file = "crisp_data.rds")

###########################################################

set.seed(1)
# Super crisp case
trans1 <- matrix(c(0.72, 0.25, 0.03, 0,
  0, 1, 0, 0,
  0, 0, 1, 0,
  0, 0, 0, 1), nrow = 4, byrow = TRUE)

trans2 <- matrix(c(0.5, 0.45, 0.05, 0,
  0, 0.7, 0.27, 0.03,
  0, 0, 1, 0,
  0, 0, 0, 1), nrow = 4, byrow = TRUE)

trans3 <- matrix(c(0.6, 0.4, 0.00, 0,
  0, 0.65, 0.05, 0.3,
  0, 0.00, 0.75, 0.25,
  0, 0, 0, 1), nrow = 4, byrow = TRUE)

trans4 <- matrix(c(1, 0, 0, 0,
  0, 1, 0, 0,
  0, 0, 1, 0,
  0, 0, 0, 1), 
  nrow = 4, byrow = TRUE)
transitions <- list(trans1, trans2, trans3, trans4)

sim_seq <- suppressMessages(simulate_mhmm(n_sequences = 1e4, 
  initial_probs = init, 
  transition_probs = transitions, 
  emission_probs = emiss,
  sequence_length = 20)$observations)

sm <- suppressMessages(seqcost(sim_seq, method = "CONSTANT"))
# Dissimilarities from OMspell
dist_OMsp <- suppressMessages(seqdist(sim_seq, method = "OMspell", sm = sm$sm, tpow = 1.5))

cluster_pam <- pam(k = 4, diss = TRUE,  x = dist_OMsp)
cluster_fanny <- fanny(k = 4, diss = TRUE,  x = dist_OMsp, memb.exp =  1.4)

# PAM is used for representativeness
d <- as.matrix(dist_OMsp)
# Dissimilarity to cluster medoids
d_to_medoid <- d[, cluster_pam$id.med]

X_repr <- 1 - d_to_medoid /  max(d)
X_class_pam <- model.matrix(~ -1 + factor(cluster_pam$clustering))
X_class_fanny <- model.matrix(~ -1 + factor(cluster_fanny$clustering))
X_prob <- cluster_fanny$membership
X_gc <- as.matrix(disscenter(d, group = cluster_pam$clustering, allcenter = TRUE))

# scale to make the effects and RMSEs comparable
X_gc <- X_gc / max(abs(X_gc))

y_repr <- X_repr %*% beta 
y_class_pam <- X_class_pam %*% beta 
y_class_fanny <- X_class_fanny %*% beta 
y_prob <- X_prob %*% beta 
y_gc <- X_gc %*% beta 

save(sim_seq, y_repr, y_class_pam, y_class_fanny, y_prob, y_gc, 
  file = "super_crisp_data.rds")

