

#### set-up ####
# set workdir
setwd("C:\\Users\\camay\\Dropbox\\Postdocs_research\\Edwards-lab\\bolivia\\pgls_and_reduced_envPCA\\")

# bring in data
data <- read.csv("env_and_morpho_matching_with_tree.csv", header = T, sep = ",", na.strings = "")

subset <- data[, 16:ncol(data)]
colnames(subset)


df_sorted <- subset[, order(colnames(subset))]


# load package
library(ggstatsplot)

# correlogram
ggstatsplot::ggcorrmat(
  data = subset,
  type = "parametric", # parametric for Pearson, nonparametric for Spearman's correlation
  # pch = "*"
  # colors = c("darkred", "white", "steelblue"), # change default colors
)

# 
# 
# # correlogram
# ggstatsplot::ggcorrmat(
#   data = subset,
#   type = "nonparametric", # parametric for Pearson, nonparametric for Spearman's correlation
# )
# 
# # correlogram
# ggstatsplot::ggcorrmat(
#   data = subset,
#   type = "robust",
# )
# 
# 
# # correlogram
# ggstatsplot::ggcorrmat(
#   data = subset,
#   type = "bayes", 
# )
