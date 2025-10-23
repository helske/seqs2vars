###### Calculate representativeness values ######

# Note: here seq_sample is the sequence data created with the seqdef function from TraMineR.
# See below for reproducible example


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

# use d_to_medoid in regression model etc.

# If you want to use other sequences as representative sequences (other than medoids):
# In this case there's no need to calculate the full dissimilarity matrix or perform 
# clustering (these steps take a lot of time and memory).
# Instead, use the refseq argument in the seqdist fuction to calculate dissimilarities to
# each representative sequence of your choice. This is faster and takes less memory.
# See the seqdist function for more information.


### Fully working example using TraMineR's biofam data ###

# Load packages
library(dplyr)
library(TraMineR)
library(ggplot2)
library(cluster)
library(marginaleffects) # for average marginal effects

#### Preparations for the Biofam data ####

# Load the Biofam data
data("biofam")

# Create the state sequence object
biofamseq <- 
  seqdef(
    biofam,
    var = 10:25, # Sequence data from the 15th to the 25th variable
    alphabet = c(0, 4, 2, 1, 5, 3, 6, 7), # Change the order of states
    labels = c("Parent", "Child", "Married", "Left", "Left+Child", "Left+Marr", 
               "Left+Marr+Child", "Divorced"), # Labels
    start = 15, # Starting from age 15,
    cpal = c("lightblue", "pink", "palegreen1", "darkblue", "brown1", 
             "palegreen3", "darkgreen", "orange") # Colour palette
  )


# OM with constant costs
dist_const_OM <- seqdist(
  biofamseq, 
  method = "OM", 
  sm = "CONSTANT"
)

# Using cluster medoids for representative sequences
# Clustering with PAM
cluster_pam <- pam(k = 5, diss = TRUE, x = dist_const_OM)
# Scaled distance to medoid
d_to_medoid <- 1 - dist_const_OM[, cluster_pam$id.med] / max(dist_const_OM)

# Name variables with representativeness scores to representative sequences
d_to_medoid <- data.frame(d_to_medoid)
names(d_to_medoid) <- c(
  "Married_child", 
  "Single", 
  "Delayed_leaver", 
  "Parents_married", 
  "Married_no_child"
)
# Combine with covariates
d <- cbind(
  biofam |> select(
    birthyr, 
    sex, 
    nationality = nat_1_02
  ), 
  d_to_medoid
) |> 
  na.exclude() |> 
  droplevels()


# OLS with representative sequences
# There is no sum-to-one constraint as in membership degrees, 
# so you can use all representativeness variables
# Predict birth_year (model is substantively meaningless; for illustration only)
fit <- lm(
  birthyr ~ sex + nationality + 
    Married_child + Single + Delayed_leaver +
    Parents_married + Married_no_child,
  data = d
)
summary(fit)

# Note: Interpreting the coefficients directly is difficult. Instead, we use
# average marginal predictions.

amps <- vector("list", length = 5) # list of length 5 (number of clusters)
names(amps) <- names(d_to_medoid) # name elements with medoid names
for (variable in names(d_to_medoid)) { # loop over medoids
  tmp <- d # new data
  # replace all values of current variable with 1 (perfectly representative)
  tmp[[variable]] <- 1
  # average marginal prediction over the data
  amps[[variable]] <- avg_predictions(fit, newdata = tmp)
}
# combine predictions and convert to data frame
amps <- bind_rows(amps, .id = "variable") |> data.frame()

amps |> 
  ggplot(aes(variable, estimate)) +
  geom_pointrange(aes(ymin = conf.low, ymax = conf.high)) +
  ylab("Birth year") +
  xlab("Group") + 
  theme_minimal()

