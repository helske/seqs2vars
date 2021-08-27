library(seqHMM)
library(TraMineR)
library(forcats)
library(cluster)

# Simulate "population" data:

set.seed(1)
# create MMM
init1 <- c(1, 0, 0, 0)
init <- list(init1, init1, init1, init1)
emiss1 <- diag(1, nrow = 4)
emiss <- list(emiss1, emiss1, emiss1, emiss1)

# Fuzzy
trans1 <- matrix(c(0.7, 0.2, 0.1, 0,
  0, 0.93, 0.04, 0.03,
  0, 0.1, 0.87, 0.03,
  0, 0.15, 0.05, 0.8), nrow = 4, byrow = TRUE)
trans2 <- matrix(c(0.4, 0.3, 0.25, 0,
  0,   0.7, 0.17, 0.13,
  0,   0.05, 0.85, 0.1,
  0,   0.1, 0.05, 0.85), nrow = 4, byrow = TRUE)
trans3 <- matrix(c(0.7, 0.3, 0, 0,
  0,   0.83, 0.04, 0.13,
  0,   0, 0.9, 0.1,
  0,   0.07, 0.03, 0.9), nrow = 4, byrow = TRUE)
trans4 <- matrix(c(0.85, 0.1, 0.05, 0,
  0,    0.9, 0.03, 0.07,
  0,    0.03, 0.9, 0.07,
  0,    0.01, 0.01, 0.98), nrow = 4, byrow = TRUE)
trans_fuzzy <- list(trans1, trans2, trans3, trans4)



sim_seq_fuzzy <- suppressMessages(simulate_mhmm(n_sequences = 10000, 
  initial_probs = init, 
  transition_probs = trans_fuzzy, 
  emission_probs = emiss,
  sequence_length = 20)$observations)

sm <- suppressMessages(seqcost(sim_seq_fuzzy, method = "CONSTANT"))
dist_OMsp <- suppressMessages(seqdist(sim_seq_fuzzy, method = "OMspell", sm = sm$sm, tpow = 1.5))
# Dissimilarities from OMspell
cluster_pam_fuzzy <- pam(k = 4, diss = TRUE,  x = dist_OMsp)
# PAM is used for representativeness
d_fuzzy <- as.matrix(dist_OMsp)
# Dissimilarity to cluster medoids
d_to_medoid <- d_fuzzy[, cluster_pam_fuzzy$id.med]
scale_fuzzy <-  max(d_fuzzy)
X_repr_fuzzy <- 1 - d_to_medoid / scale_fuzzy
X_class_fuzzy <- matrix(model.matrix(~ -1 + factor(cluster_pam_fuzzy$clustering)), ncol = 4)

y_repr_fuzzy <- X_repr_fuzzy %*% c(0, 1, 1, -1) + rnorm(1e4, sd = 0.25)
y_class_fuzzy <- X_class_fuzzy %*% c(0, 1, 1, -1) + rnorm(1e4, sd = 0.25)

save(sim_seq_fuzzy, y_repr_fuzzy, y_class_fuzzy, 
  file = "fuzzy_data.rds")

####################################################################

set.seed(1)
# Crisp

trans1 <- matrix(c(0.8, 0.2, 0, 0,
  0, 0.95, 0.04, 0.01,
  0, 0, 0.995, 0.005,
  0, 0.15, 0.05, 0.8), nrow = 4, byrow = TRUE)
trans2 <- matrix(c(0.5, 0.45, 0.05, 0,
  0,   0.8, 0.19, 0.01,
  0,   0.00, 0.99, 0.01,
  0,   0.01, 0.01, 0.98), nrow = 4, byrow = TRUE)
trans3 <- matrix(c(0.6, 0.35, 0.05, 0,
  0,   0.75, 0.07, 0.18,
  0,   0.00, 0.82, 0.18,
  0,   0.02, 0.01, 0.97), nrow = 4, byrow = TRUE)
trans4 <- matrix(c(0.999, 0.0005, 0.0005, 0,
  0,    0.9, 0.03, 0.07,
  0,    0.03, 0.9, 0.07,
  0,    0.01, 0.01, 0.98), nrow = 4, byrow = TRUE)
trans_crisp <- list(trans1, trans2, trans3, trans4)

sim_seq_crisp <- suppressMessages(simulate_mhmm(n_sequences = 10000, 
  initial_probs = init, 
  transition_probs = trans_crisp, 
  emission_probs = emiss,
  sequence_length = 20)$observations)

sm <- suppressMessages(seqcost(sim_seq_crisp, method = "CONSTANT"))
dist_OMsp <- suppressMessages(seqdist(sim_seq_crisp, method = "OMspell", sm = sm$sm, tpow = 1.5))
# Dissimilarities from OMspell
cluster_pam_crisp <- pam(k = 4, diss = TRUE,  x = dist_OMsp)
# PAM is used for representativeness
d_crisp <- as.matrix(dist_OMsp)
# Dissimilarity to cluster medoids
d_to_medoid <- d_crisp[, cluster_pam_crisp$id.med]
scale_crisp <-  max(d_crisp)
X_repr_crisp <- 1 - d_to_medoid / scale_crisp
X_class_crisp <- matrix(model.matrix(~ -1 + factor(cluster_pam_crisp$clustering)), ncol = 4)

y_repr_crisp <- X_repr_crisp %*% c(0, 1, 1, -1) + rnorm(1e4, sd = 0.25)
y_class_crisp <- X_class_crisp %*% c(0, 1, 1, -1) + rnorm(1e4, sd = 0.25)

save(sim_seq_crisp, y_repr_crisp, y_class_crisp, 
  file = "crisp_data.rds")

###########################################################

set.seed(1)
# Super crisp
# Green
trans1 <- matrix(c(0.72, 0.25, 0, 0,
  0, 1, 0, 0,
  0, 0, 0, 0,
  0, 0, 0, 0), nrow = 4, byrow = TRUE)
# Red
trans2 <- matrix(c(0.5, 0.45, 0.05, 0,
  0,   0.7, 0.27, 0,
  0,   0, 1, 0,
  0,   0, 0, 1), nrow = 4, byrow = TRUE)
# Grey
trans3 <- matrix(c(0.6, 0.4, 0.00, 0,
  0,   0.65, 0.05, 0.3,
  0,   0.00, 0.75, 0.25,
  0,   0, 0, 1), nrow = 4, byrow = TRUE)
# Purple
trans4 <- matrix(c(1, 0, 0, 0,
  0,    1, 0, 0,
  0,    0, 1, 0,
  0,    0, 0, 0), nrow = 4, byrow = TRUE)
trans_super_crisp <- list(trans1, trans2, trans3, trans4)


sim_seq_super_crisp <- suppressMessages(simulate_mhmm(n_sequences = 10000, 
  initial_probs = init, 
  transition_probs = trans_super_crisp, 
  emission_probs = emiss,
  sequence_length = 20)$observations)

sm <- suppressMessages(seqcost(sim_seq_super_crisp, method = "CONSTANT"))
dist_OMsp <- suppressMessages(seqdist(sim_seq_super_crisp, method = "OMspell", sm = sm$sm, tpow = 1.5))
# Dissimilarities from OMspell
cluster_pam_super_crisp <- pam(k = 4, diss = TRUE,  x = dist_OMsp)
# PAM is used for representativeness
d_super_crisp <- as.matrix(dist_OMsp)
# Dissimilarity to cluster medoids
d_to_medoid <- d_super_crisp[, cluster_pam_super_crisp$id.med]
scale_super_crisp <-  max(d_super_crisp)
X_repr_super_crisp <- 1 - d_to_medoid / scale_super_crisp
X_class_super_crisp <- matrix(model.matrix(~ -1 + factor(cluster_pam_super_crisp$clustering)), ncol = 4)

y_repr_super_crisp <- X_repr_super_crisp %*% c(0, 1, 1, -1) + rnorm(1e4, sd = 0.25)
y_class_super_crisp <- X_class_super_crisp %*% c(0, 1, 1, -1) + rnorm(1e4, sd = 0.25)

save(sim_seq_super_crisp, y_repr_super_crisp, y_class_super_crisp, 
  file = "super_crisp_data.rds")


