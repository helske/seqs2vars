# read files
files <- list.files()
result_files <- files[grep(files, 
  pattern = "^(class|prob|repr|gc).*\\.rds", perl = TRUE)]

# collect RMSEs into a data frame

bic <- do.call("rbind", lapply(result_files, function(i) {
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
  data.frame(
    variable_type = factor(type, 
      levels = c("class_pam", "class_fanny", "prob", "repr", "gc"),
      labels = c("hard classification (PAM)", "hard classification (FANNY)",
        "membership degree", "representativeness", "gravity center")),
    strength = factor(strength, levels = c("fuzzy", "crisp", "super"),
      labels = c("Weak", "Moderate", "Strong"), ordered = TRUE),
    estimation_method = factor(rownames(x), 
      levels = c("Hard classification (PAM)", "Hard classification (FANNY)", 
        "Soft classification", "Pseudoclass", "Gravity center",
        "Representativeness"), ordered = TRUE),
    mean_bic = as.numeric(rowMeans(x, na.rm = TRUE)),
    lwr = apply(x, 1, quantile, 0.025, na.rm = TRUE),
    upr = apply(x, 1, quantile, 0.975, na.rm = TRUE),
    n_na = rowSums(is.na(x)))
}))


p <- bic %>% 
  filter(estimation_method != "Pseudoclass") %>%
  ggplot(aes(x = strength, y = mean_bic,
    shape = estimation_method)) +
  geom_pointrange(aes(ymin = lwr, ymax = upr), 
    position = position_dodge(0.6)) +
  theme_bw() + scale_y_continuous("BIC") + 
  scale_x_discrete("Clustering tendency", labels = c("Weak",
    "Moderate",
    "Strong"))+
  scale_shape_discrete("Estimation method") +
  facet_wrap(~ variable_type, ncol = 1, scales = "free",
    labeller = as_labeller(function(x) paste0("Outcome based on ", x))) +
  theme(legend.position = "right") 
p

ggsave(p, file = "simulation_BIC.png")
