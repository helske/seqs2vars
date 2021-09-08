library(dplyr)
library(ggplot2)

# read files
result_files <- list.files(pattern = "^(class|prob|repr)")
# collect RMSEs into a data frame

rmse <- do.call("rbind", lapply(result_files, function(i) {
  x <- readRDS(i)
  s <- unlist(strsplit(i, "(_|\\.)"))
  data.frame(cluster_quality = s[2],
    variable_type = s[1],
    estimation_method = rownames(x),
    mean_rmse = as.numeric(rowMeans(x, na.rm = TRUE)),
    lwr = apply(x, 1, quantile, 0.025, na.rm = TRUE),
    upr = apply(x, 1, quantile, 0.975, na.rm = TRUE),
    n_na = rowSums(is.na(x)))
}))

methods <- c("Hard classification", "Soft classification", "Pseudoclass", "Representativeness")

p <- rmse %>%
  mutate(cluster_quality =
      factor(cluster_quality,
        levels = c("fuzzy", "crisp", "super"),
        labels = c("Weak", "Moderate", "Strong"), ordered = TRUE),
    estimation_method = factor(estimation_method,
      levels = methods, ordered = TRUE)) %>%
  ggplot(aes(x = cluster_quality, y = mean_rmse,
    shape = estimation_method)) +
  geom_pointrange(aes(ymin = lwr, ymax = upr), position = position_dodge(0.6)) +
  theme_bw() + scale_y_continuous("RMSE") +
  scale_x_discrete("Clustering tendency", labels = c("Weak",
    "Moderate",
    "Strong"))+
  scale_shape_discrete("Estimation method", labels = c(
    "Hard classification", "Soft classification", "Pseudoclass", "Representativeness")) +
  facet_grid(~variable_type,
    labeller = as_labeller(c("class" = "Response based on classification",
      "repr" = "Response based on representativeness"))) +
  theme(legend.position = "bottom")
ggsave(p, file = "simulation.png")
