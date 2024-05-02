#The idea here is to see if there is a correlation between the morphology (represented by the NMDS axis) and the environmental PCs after a pca





library(phylolm)

library(ape)
# library(nlme)
library(geiger)



#load data
data <-read.table("./envPCsNohighCorrelated_and_morpho_matching_with_tree.csv", sep = ",", header = TRUE, row.names=1)  # env pca without high correlated vars
# data <-read.table("./envPCs_and_morpho_matching_with_tree.csv", sep = ",", header = TRUE, row.names=1) #env pca with all chelsa vars
tree <-read.tree("./RAxML_bipartitions.10-bolivia-initial_mcov0.25_rcov0")


#check correspondence
obj<-name.check(tree, data)
obj


#prune tree 
pruned_tree <- drop.tip(tree, obj$tree_not_data)

#check again correspondence
name.check(pruned_tree, data)


# Extract the row names as a separate column called ids
data$ids <- row.names(data)



## Try all response and predictors individually


response_columns <- colnames(data)[4:13]
response_columns
predictor_columns <- colnames(data)[17:ncol(data)-1]
predictor_columns

# Load the 'gtools' package
library(gtools)

# Sort the predictor_columns
predictor_columns <- mixedsort(predictor_columns)
predictor_columns

#create empty matrices
results <- matrix(NA, nrow = length(response_columns), ncol = length(predictor_columns),
                  dimnames = list(response_columns, predictor_columns))

pvalues <- matrix(NA, nrow = length(response_columns), ncol = length(predictor_columns),
                  dimnames = list(response_columns, predictor_columns))


aics_bm <- matrix(NA, nrow = length(response_columns), ncol = length(predictor_columns),
                  dimnames = list(response_columns, predictor_columns))


# Try using brownian motion
for (i in 1:length(response_columns)) {
  for (j in 1:length(predictor_columns)) {
    formula_str <- paste(response_columns[i], "~", predictor_columns[j])
    
    # Attempt to fit the model
    tryCatch({
      # model <- gls(as.formula(formula_str), data = data, correlation = bm)
      # summary_table <- summary(model)$tTable
      # results[i, j] <- summary_table[predictor_columns[j], "Value"]
      # pvalues[i, j] <- summary_table[predictor_columns[j], "p-value"]
      # 
      # print(formula_str)
      # print(summary_table)
      
      #new phyloml model
      model <- phylolm(as.formula(formula_str), data = data, pruned_tree, model = "BM")
      results[i, j] <- summary(model)$coefficients[predictor_columns[j], "Estimate"]
      pvalues[i, j] <- summary(model)$coefficients[predictor_columns[j], "p.value"]
      aics_bm[i, j] <- summary(model)$aic
      
    }, error = function(e) {
      # If an error occurs, populate the results with NAs
      results[i, j] <- NA
      pvalues[i, j] <- NA
    })
  }
}





# Create a new dataframe with significant coefficient estimates
masked_results <- results
masked_results[pvalues >= 0.055] <- NA
masked_results[pvalues == "NaN"] <- NA

masked_dataframe <- as.data.frame(masked_results)
results_dataframe <- as.data.frame(results)


write.csv(masked_dataframe, "1e-masked_coefficients_results_OnlySignificants_BM.csv", row.names = FALSE)
write.csv(results_dataframe, "1e-coefficients_results_dataframe_BM.csv", row.names = FALSE)



library(ggplot2)
library(reshape2)

# Convert results to a data frame
df_results <- melt(masked_results, varnames = c("Response", "Predictor"), value.name = "Coefficient")


# Create the heatmap plot using ggplot2 with log scale
heatmap_plot <- ggplot(df_results, aes(x = Predictor, y = Response, fill = Coefficient)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "#fafafa") + #, limits = c(-1, 1)) +
  labs(title = "Coefficient Estimates Heatmap (model: BM)", x = "Predictor Variables", y = "Response Variables") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_vline(xintercept = seq(0.5, nlevels(df_results$Predictor) + 0.5, by = 1), color = "gray", linetype = "dashed") +
  geom_hline(yintercept = seq(0.5, nlevels(df_results$Response) + 0.5, by = 1), color = "gray", linetype = "dashed") +
  geom_text(aes(label = ifelse(is.na(Coefficient), "", sprintf("%.3f", Coefficient))), color = "black", size = 3)
 

heatmap_plot







# Save the plot with fixed height
ggsave("1e-heatmap_morphoVSenvPCs_BM.svg", plot = heatmap_plot, width = 6, height = 10, units = "in")
ggsave("1e-heatmap_morphoVSenvPCs_BM.png", plot = heatmap_plot, width = 6, height = 10, units = "in")



### Conclusion
# Results are very different to the nmds but still everything is kinda weak


############## now using OU model

# empty matrices
results <- matrix(NA, nrow = length(response_columns), ncol = length(predictor_columns),
                  dimnames = list(response_columns, predictor_columns))

pvalues <- matrix(NA, nrow = length(response_columns), ncol = length(predictor_columns),
                  dimnames = list(response_columns, predictor_columns))


aics_ou <- matrix(NA, nrow = length(response_columns), ncol = length(predictor_columns),
               dimnames = list(response_columns, predictor_columns))




# Try using brownian motion
for (i in 1:length(response_columns)) {
  for (j in 1:length(predictor_columns)) {
    formula_str <- paste(response_columns[i], "~", predictor_columns[j])
    
    # Attempt to fit the model
    tryCatch({
      # model <- gls(as.formula(formula_str), data = data, correlation = bm)
      # summary_table <- summary(model)$tTable
      # results[i, j] <- summary_table[predictor_columns[j], "Value"]
      # pvalues[i, j] <- summary_table[predictor_columns[j], "p-value"]
      # 
      # print(formula_str)
      # print(summary_table)
      
      #new phyloml model
      model <- phylolm(as.formula(formula_str), data = data, pruned_tree, model = "OUrandomRoot")
      results[i, j] <- summary(model)$coefficients[predictor_columns[j], "Estimate"]
      pvalues[i, j] <- summary(model)$coefficients[predictor_columns[j], "p.value"]
      aics_ou[i, j] <- summary(model)$aic
      
    }, error = function(e) {
      # If an error occurs, populate the results with NAs
      results[i, j] <- NA
      pvalues[i, j] <- NA
    })
  }
}





# Create a new dataframe with significant coefficient estimates
masked_results <- results
masked_results[pvalues >= 0.055] <- NA
masked_results[pvalues == "NaN"] <- NA

masked_dataframe <- as.data.frame(masked_results)
results_dataframe <- as.data.frame(results)


write.csv(masked_dataframe, "1e-masked_coefficients_results_OnlySignificants_OU.csv", row.names = FALSE)
write.csv(results_dataframe, "1e-coefficients_results_dataframe_OU.csv", row.names = FALSE)



library(ggplot2)
library(reshape2)

# Convert results to a data frame
df_results <- melt(masked_results, varnames = c("Response", "Predictor"), value.name = "Coefficient")



# Create the heatmap plot using ggplot2 with log scale
heatmap_plot <- ggplot(df_results, aes(x = Predictor, y = Response, fill = Coefficient)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", na.value = "#fafafa") + #, limits = c(-1, 1)) +
  labs(title = "Coefficient Estimates Heatmap (model: OU)", x = "Predictor Variables", y = "Response Variables") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_vline(xintercept = seq(0.5, nlevels(df_results$Predictor) + 0.5, by = 1), color = "gray", linetype = "dashed") +
  geom_hline(yintercept = seq(0.5, nlevels(df_results$Response) + 0.5, by = 1), color = "gray", linetype = "dashed") +
  geom_text(aes(label = ifelse(is.na(Coefficient), "", sprintf("%.3f", Coefficient))), color = "black", size = 3)


heatmap_plot



# Save the plot with fixed height
ggsave("1e-heatmap_morphoVSenvPCs_OU.svg", plot = heatmap_plot, width = 6, height = 10, units = "in")
ggsave("1e-heatmap_morphoVSenvPCs_OU.png", plot = heatmap_plot, width = 6, height = 10, units = "in")



## Comparing AICs of both models


### Using simple means: BM is better
#remove problematic stellate_hairs_stalked categorical var
#this is removed by the model but the aic is still reported
# Replace -Inf with NA in aics_bm
aics_bm[aics_bm == -Inf] <- NA
aics_ou[aics_ou == -Inf] <- NA

# Calculate mean AIC values for each matrix
mean_aic_bm <- mean(aics_bm, na.rm = TRUE)
mean_aic_ou <- mean(aics_ou, na.rm = TRUE)

print(paste("bm: ",mean_aic_bm, " ou: ", mean_aic_ou))


# Calculate mean AIC values for each column in 'aics_bm'
mean_aics_bm <- apply(aics_bm, 2, mean, na.rm = TRUE)

# Calculate mean AIC values for each column in 'aics_ou'
mean_aics_ou <- apply(aics_ou, 2, mean, na.rm = TRUE)

# Print mean AIC values
print(mean_aics_bm)
print(mean_aics_ou)


# Calculate mean AIC values for each row in 'aics_bm'
mean_aics_bm_r <- apply(aics_bm, 1, mean, na.rm = TRUE)

# Calculate mean AIC values for each row in 'aics_ou'
mean_aics_ou_r <- apply(aics_ou, 1, mean, na.rm = TRUE)

# Print mean AIC values
print(mean_aics_bm_r)
print(mean_aics_ou_r)



# Combine mean AIC values into a data frame
mean_aics <- data.frame(
  Variable = row.names(aics_bm),
  BM = mean_aics_bm_r,
  OU = mean_aics_ou_r
)

# Reshape data for plotting
library(tidyr)
mean_aics_long <- pivot_longer(mean_aics, cols = c(BM, OU), names_to = "Model", values_to = "Mean_AIC")

# Plot side-by-side bar plot
library(ggplot2)
ggplot(mean_aics_long, aes(x = Variable, y = Mean_AIC, fill = Model)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Mean AIC Values for Each Variable",
       x = "Variable",
       y = "Mean AIC Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




# Perform a paired t-test to compare mean AIC values
t.test(aics_bm, aics_ou, paired = TRUE)
# This is showing that it is significant the difference between both matrices
#t = -3.1393, df = 39, p-value = 0.003223


# Perform Wilcoxon signed-rank test
wilcox.test(as.vector(aics_bm), as.vector(aics_ou), paired = TRUE)
