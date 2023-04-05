###### Calculate representativeness values ######

# Note: here seq_sample is the sequence data created with the seqdef function from TraMineR.

# Representativeness values are calculated from sequence dissimilarities. 
# Here is an example for using cluster medoids as representative sequences.
# Note: You can also opt for any other relevant representative sequences (see below)

# Compute sequence dissimilarities using OMspell dissimilarity measure
# Note: See the TraMineR manual for more information
# Costs for substitutions (here: constant costs)
sm <- seqcost(seq_sample, method = "CONSTANT")
# Calculate dissimilarities 
dist_OMsp <- seqdist(seq_sample, method = "OMspell", sm = sm$sm, tpow = 1.5)
# Dissimilarities from OMspell collected in a matrix
d <- as.matrix(dist_OMsp)

# Using cluster medoids for representative sequences
# Clustering with PAM
cluster_pam <- pam(k = k, diss = TRUE, x = dist_OMsp)
      
# Calculate representativeness to cluster medoids (based on PAM)
      
# Scaled distance to medoid
d_to_medoid <- 1 - d[, cluster_pam$id.med] / max(d)

# If you want to use other sequences as representative sequences (other than medoids)

# In this case there's no need to calculate the full dissimilarity matrix or perform 
# clustering (these steps take a lot of time and memory).
# Instead, use the refseq argument in the seqdist fuction to calculate dissimilarities to
# each representative sequence of your choice. This is faster and takes less memory.
# See the seqdist function for more information.
