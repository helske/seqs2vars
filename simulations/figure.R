library(dplyr)
library(ggplot2)
# read files
setwd("simulations")
files <- list.files()
result_files <- files[grep(files, 
  pattern = "^(class|prob|repr|gc).*\\.rds", perl = TRUE)]

# collect RMSEs into a data frame

rmse <- do.call("rbind", lapply(result_files, function(i) {
  x <- readRDS(i)
  i <- sub(".rds", "", i)
  s <- unlist(strsplit(i, "_"))
  type <- s[1]
  if(type == "class") {
    type <- paste0(s[1], "_", s[2])
    strength <- s[3]
  } else {
    strength <- s[2]
  }
  n <- as.numeric(s[length(s)])
  ranks <- apply(x, 2, rank)
  x[x==0] <- NA
  data.frame(
    variable_type = factor(type, 
      levels = c("class_pam", "class_fanny", "prob", "repr", "gc"),
      labels = c("hard classification (PAM)", "hard classification (FANNY)",
        "soft classification", "representativeness", "gravity center")),
    strength = factor(strength, levels = c("fuzzy", "crisp", "super"),
      labels = c("Weak", "Moderate", "Strong"), ordered = TRUE),
    n = factor(n),
    estimation_method = factor(rownames(x), 
      levels = c("Hard classification (PAM)", "Hard classification (FANNY)", 
        "Soft classification", "Pseudoclass", "Gravity center",
        "Representativeness"), ordered = TRUE),
    mean_rmse = as.numeric(rowMeans(x, na.rm = TRUE)),
    lwr = apply(x, 1, quantile, 0.025, na.rm = TRUE),
    upr = apply(x, 1, quantile, 0.975, na.rm = TRUE),
    avg_rank = rowMeans(ranks),
    min_rank = apply(ranks, 1, quantile, 0.25, na.rm = TRUE),
    max_rank = apply(ranks, 1, quantile, 0.75, na.rm = TRUE),
    n_na = rowSums(is.na(x)))
}))

p <- rmse %>%
  filter(n == 200) %>%
  ggplot(aes(x = strength, y = mean_rmse,
    shape = estimation_method, colour = estimation_method)) +
  geom_pointrange(aes(ymin = lwr, ymax = upr), 
    position = position_dodge(0.6)) +
  theme_bw() + scale_y_continuous("RMSE") + 
  scale_x_discrete("Clustering tendency", labels = c("Weak",
    "Moderate",
    "Strong"))+
  scale_shape_discrete("Estimation method") +
  scale_colour_discrete("Estimation method") +
  facet_grid(~ variable_type,
    labeller = as_labeller(function(x) paste0("Outcome based on ", x))) +
  theme(legend.position = "bottom")
p

ggsave(p, file = "simulation_n200.png")


p <- rmse %>%
  filter(n == 1000) %>%
  ggplot(aes(x = strength, y = mean_rmse,
    shape = estimation_method, colour = estimation_method)) +
  geom_pointrange(aes(ymin = lwr, ymax = upr), 
    position = position_dodge(0.6)) +
  theme_bw() + scale_y_continuous("RMSE") + 
  scale_x_discrete("Clustering tendency", labels = c("Weak",
    "Moderate",
    "Strong"))+
  scale_shape_discrete("Estimation method") +
  scale_colour_discrete("Estimation method") +
  facet_grid(~ variable_type,
    labeller = as_labeller(function(x) paste0("Outcome based on ", x))) +
  theme(legend.position = "bottom")
p

ggsave(p, file = "simulation_n1000.png")


# for the paper
rmse_subset <- rmse %>% filter(n == 1000) %>%
  filter(estimation_method %in% 
      c("Hard classification (PAM)", "Soft classification", 
        "Pseudoclass", "Representativeness")) %>%
  filter(variable_type %in% c("hard classification (PAM)",
    "soft classification", "representativeness"))
rmse_subset$estimation_method <- droplevels(rmse_subset$estimation_method)
levels(rmse_subset$estimation_method) <- 
  c("Hard classification", "Soft classification", "Pseudoclass", 
    "Representativeness")
rmse_subset$variable_type <- droplevels(rmse_subset$variable_type)
levels(rmse_subset$variable_type) <- 
  c("hard classification", "soft classification", "representativeness")

p <- rmse_subset %>%
  ggplot(aes(x = strength, y = mean_rmse,
    shape = estimation_method)) +
  geom_pointrange(aes(ymin = lwr, ymax = upr), 
    position = position_dodge(0.6)) +
  theme_bw() + scale_y_continuous("RMSE") + 
  scale_x_discrete("Clustering tendency", labels = c("Weak",
    "Moderate",
    "Strong"))+
  scale_shape_discrete("Estimation method") +
  facet_grid(~ variable_type,
    labeller = as_labeller(function(x) paste0("Outcome based on ", x))) +
  theme(legend.position = "bottom", text = element_text(size = 19))
p
# adjust size for the paper (see also element_text above)
ggsave(p, file = "simulation.png", width = 12, height = 8)

# Null-model case
# 

# read files
result_files <- list.files(pattern = "noeffect.rds")
# collect RMSEs into a data frame

rmse <- do.call("rbind", lapply(result_files, function(i) {
  x <- readRDS(i)
  s <- unlist(strsplit(i, "(_|\\.)"))
  data.frame(strength = s[1],
    estimation_method = factor(rownames(x), 
      levels = c("Hard classification (PAM)", "Hard classification (FANNY)", 
        "Soft classification", "Pseudoclass", "Gravity center",
        "Representativeness"), ordered = TRUE),
    mean_rmse = as.numeric(rowMeans(x, na.rm = TRUE)),
    lwr = apply(x, 1, quantile, 0.025, na.rm = TRUE),
    upr = apply(x, 1, quantile, 0.975, na.rm = TRUE),
    avg_rank = rowMeans(ranks),
    min_rank = apply(ranks, 1, quantile, 0.05, na.rm = TRUE),
    max_rank = apply(ranks, 1, quantile, 0.95, na.rm = TRUE),
    n_na = rowSums(is.na(x)))
}))

p <- rmse %>%
  ggplot(aes(x = strength, y = mean_rmse,
    shape = estimation_method, colour = estimation_method)) +
  geom_pointrange(aes(ymin = lwr, ymax = upr), 
    position = position_dodge(0.6)) +
  theme_bw() + scale_y_continuous("RMSE") + 
  scale_x_discrete("Clustering tendency", labels = c("Weak",
    "Moderate",
    "Strong"))+
  scale_shape_discrete("Estimation method") +
  scale_colour_discrete("Estimation method") +
  theme(legend.position = "bottom")
p
ggsave(p, file = "simulation_noeffect.png")

