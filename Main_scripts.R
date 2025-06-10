## Analysis1 : Relationship between plastome size and structural variables ####

### OLS model ######################################################
# Load required package
library(readxl)

# Load data
data <- read.csv("plastid_structure.csv")
data[1:13] <- lapply(data[1:13], as.numeric)

# Define variables
traits_to_fit <- c("LSC_length", "SSC_Length", "Total_IR_Length", 
                   "CDS_Length", "tRNA_Length", "rRNA_Length", 
                   "Speacer_Length", "Dispeard_repeat_Length", 
                   "Tandem_repeat_Length", "ShannonEntropy", "GC","KC_ratio")
trait_to_response <- c("All_Length")

# Prepare result dataframe
ols_results <- data.frame(
  Variable = character(),
  Estimate = numeric(),
  Std_Error = numeric(),
  t_value = numeric(),
  P_value = numeric(),
  R2adj = numeric(),
  Residual_SE = numeric(),
  F_statistic = numeric(),
  stringsAsFactors = FALSE
)

# Loop over predictor traits
for (predictor in traits_to_fit) {
  for (response in trait_to_response) {
    # Remove NA
    df <- data[, c(response, predictor)]
    df <- df[complete.cases(df), ]
    
    if (nrow(df) < 5) next  # skip if not enough data
    
    formula <- as.formula(paste(response, "~", predictor))
    fit <- lm(formula, data = df)
    summary_fit <- summary(fit)
    
    estimate <- coef(summary_fit)[2, 1]
    std_error <- coef(summary_fit)[2, 2]
    t_val <- coef(summary_fit)[2, 3]
    p_val <- coef(summary_fit)[2, 4]
    r2adj <- summary_fit$adj.r.squared
    resid_se <- summary_fit$sigma
    f_stat <- summary_fit$fstatistic[1]
    
    ols_results <- rbind(ols_results, data.frame(
      Variable = predictor,
      Estimate = estimate,
      Std_Error = std_error,
      t_value = t_val,
      P_value = p_val,
      R2adj = r2adj,
      Residual_SE = resid_se,
      F_statistic = f_stat
    ))
  }
}

# Save results
write.csv(ols_results, file = "OLS_results_summary_PP.csv", row.names = FALSE)

### PGLS model ######################################################
# Load required package
library(ape)
library(phytools)
library(phylolm)
library(ggplot2)
library(dplyr)
library(future.apply)

# Load phylogenetic trees
phy_tree <- read.nexus("vascular_plant_tree.randrom100tres.v1.nex")

# Load data
data <- read.csv("plastid_structure.csv")
data[1:13] <- lapply(data[1:13], as.numeric)

# Define response and predictor traits
traits_to_fit <- c("LSC_length", "SSC_Length", "Total_IR_Length", 
                   "CDS_Length", "tRNA_Length", "rRNA_Length", 
                   "Speacer_Length", "Dispeard_repeat_Length", 
                   "Tandem_repeat_Length", "ShannonEntropy", "GC", "KC_ratio")
trait_to_response <- c("All_Length")

# Define phylogenetic models to test
models <- c("BM", "OUfixedRoot", "lambda", "kappa", "EB")

# Initialize results list
results <- list()

# Run PGLS for each tree and trait pair
for (tree_index in 1:length(phy_tree)) {
  current_tree <- phy_tree[[tree_index]]
  tree_species <- current_tree$tip.label
  
  for (predictor_trait in traits_to_fit) {
    for (dependent_trait in trait_to_response) {
      
      complete_data <- traits[!is.na(traits[[predictor_trait]]) & !is.na(traits[[dependent_trait]]), ]
      if (nrow(complete_data) < 3) {
        cat("Warning: Tree", tree_index, ": insufficient data for", predictor_trait, "and", dependent_trait, "\n")
        next
      }
      
      matched_species <- intersect(complete_data$wcvp_name, tree_species)
      matched_data <- complete_data[complete_data$wcvp_name %in% matched_species, ]
      pruned_tree <- drop.tip(current_tree, setdiff(tree_species, matched_data$wcvp_name))
      if (length(matched_species) < 3) next
      rownames(matched_data) <- matched_data$wcvp_name
      
      for (model in models) {
        cat("Tree:", tree_index, ", Model:", model, ", Predictor:", predictor_trait, ", Response:", dependent_trait, "\n")
        formula <- as.formula(paste(dependent_trait, "~", predictor_trait))
        
        fit <- phylolm(formula, data = matched_data, phy = pruned_tree, model = model)
        
        results <- append(results, list(list(
          tree_index = tree_index,
          model = model,
          predictor_trait = predictor_trait,
          dependent_trait = dependent_trait,
          AIC = fit$aic,
          fit = fit
        )))
      }
    }
  }
}

saveRDS(results, file = "PGLS_results_PP.rds")

# Extract AIC and logLik
result_table <- data.frame(
  tree_number = integer(),
  dependent_trait = character(),
  predictor_trait = character(),
  Model = character(),
  logLik = numeric(),
  AICc = numeric(),
  stringsAsFactors = FALSE
)

for (res in results) {
  result_table <- rbind(result_table, data.frame(
    tree_number = res$tree_index,
    dependent_trait = res$dependent_trait,
    predictor_trait = res$predictor_trait,
    Model = res$model,
    logLik = res$fit$logLik,
    AICc = res$fit$aic,
    stringsAsFactors = FALSE
  ))
}

# Summarize model AICc across trees
summary_table <- result_table %>%
  group_by(dependent_trait, predictor_trait, Model) %>%
  summarise(
    logLik = mean(logLik, na.rm = TRUE),
    AICc = mean(AICc, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  group_by(dependent_trait, predictor_trait) %>%
  mutate(
    min_AICc = min(AICc, na.rm = TRUE),
    delta_AICc = AICc - min_AICc,
    AICc_weight = exp(-0.5 * delta_AICc) / sum(exp(-0.5 * delta_AICc))
  ) %>%
  ungroup()

final_table <- summary_table %>%
  select(dependent_trait, predictor_trait, Model, logLik, AICc, delta_AICc, AICc_weight)

colnames(final_table) <- c("dependent_trait", "predictor_trait", "Model", "logLik", "AICc", "△AICc", "AICc weight")
write.csv(final_table, file = "AICc_weight_summary_PP.csv", row.names = FALSE)

library(dplyr)

# Initialize results dataframe
lambda_results <- data.frame(
  dependent_trait = character(),
  predictor_trait = character(),
  Estimate = numeric(),
  Std_Error = numeric(),
  t_value = numeric(),
  P_value = numeric(),
  R_squared_adj = numeric(),
  Residual_SE = numeric(),
  F_statistic = numeric(),
  lambda = numeric(),
  Sample_size = integer(),
  stringsAsFactors = FALSE
)

# Loop through model results
for (tree_index in seq_along(results)) {
  tree_results <- results[[tree_index]]
  
  for (model_result in tree_results) {
    if (model_result$model == "lambda") {
      fit <- model_result$fit
      
      # Extract traits
      dependent_trait <- model_result$dependent_trait
      predictor_trait <- model_result$predictor_trait
      
      # Estimate statistics
      beta <- fit$coefficients[2]
      se <- sqrt(diag(fit$vcov))[2]
      t_val <- beta / se
      df <- fit$n - length(fit$coefficients)
      p_val <- 2 * pt(-abs(t_val), df)
      
      # Additional stats
      r2_adj <- fit$adj.r.squared
      resid_se <- sqrt(fit$sigma2)
      f_stat <- t_val^2  # For single predictor
      lambda_val <- fit$optpar
      n <- fit$n
      
      # Store in results
      lambda_results <- rbind(lambda_results, data.frame(
        dependent_trait = dependent_trait,
        predictor_trait = predictor_trait,
        Estimate = beta,
        Std_Error = se,
        t_value = t_val,
        P_value = p_val,
        R_squared_adj = r2_adj,
        Residual_SE = resid_se,
        F_statistic = f_stat,
        lambda = lambda_val,
        Sample_size = n,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Summarize by trait pair
lambda_summary <- lambda_results %>%
  group_by(dependent_trait, predictor_trait) %>%
  summarise(
    Estimate = mean(Estimate, na.rm = TRUE),
    Std_Error = mean(Std_Error, na.rm = TRUE),
    t_value = mean(t_value, na.rm = TRUE),
    P_value = mean(P_value, na.rm = TRUE),
    R_squared_adj = mean(R_squared_adj, na.rm = TRUE),
    Residual_SE = mean(Residual_SE, na.rm = TRUE),
    F_statistic = mean(F_statistic, na.rm = TRUE),
    lambda = mean(lambda, na.rm = TRUE),
    Sample_size = mean(Sample_size, na.rm = TRUE),
    .groups = 'drop'
  )

# View summary
write.csv(lambda_summary, file = "lambda_model_summary_PP.csv", row.names = FALSE)

## Analysis2: Relationship between plastid and nuclear genome sizes (n = 1,291) ####

### OLS model ######################################################
# Load required package
library(readxl)

# Load data
data <- read_excel("10kp_structure_detial.xlsx", sheet = 13)

# Define regression formulas
formula_raw <- Chloroplast_genome_size ~ Nuclear_genome_size
formula_log <- log_Chloroplast_genome_size ~ log_Nuclear_genome_size

# Fit OLS models
lm_model_raw <- lm(formula_raw, data = data)
lm_model_log <- lm(formula_log, data = data)

# Summarize models
summary_lm_raw <- summary(lm_model_raw)
summary_lm_log <- summary(lm_model_log)

# Print summaries
cat("\n### Raw Genome Size Regression ###\n")
print(summary_lm_raw)

cat("\n### Log-Transformed Genome Size Regression ###\n")
print(summary_lm_log)

### PGLS model ######################################################
# Load required packages
library(ape)
library(phytools)
library(phylolm)
library(readxl)
library(dplyr)

# Read the phylogenetic trees in Newick format
phy_tree <- read.tree("pruned_100_plastid&nuclear_trees.nwk")

# Read the trait data
traits <- read_excel("10kp_structure_detial.xlsx", sheet = 21)
traits$wcvp_name <- as.character(traits$wcvp_name)
traits[, 2:5] <- lapply(traits[, 2:5], as.numeric)

# Define model pairs (raw and log-transformed)
model_pairs <- list(
  list(response = "Chloroplast_genome_size", predictor = "Nuclear_genome_size"),
  list(response = "log_Chloroplast_genome_size", predictor = "log_Nuclear_genome_size")
)

# Define phylogenetic models to test
models <- c("BM", "OUfixedRoot", "lambda", "kappa", "EB")

# Initialize result list
results <- list()

# PGLS across trees and models
for (tree_index in seq_along(phy_tree)) {
  current_tree <- phy_tree[[tree_index]]
  tree_species <- current_tree$tip.label
  
  for (pair in model_pairs) {
    response <- pair$response
    predictor <- pair$predictor
    
    data_complete <- traits[!is.na(traits[[response]]) & !is.na(traits[[predictor]]), ]
    matched_species <- intersect(data_complete$wcvp_name, tree_species)
    if (length(matched_species) < 3) next
    
    data_matched <- data_complete[data_complete$wcvp_name %in% matched_species, ]
    pruned_tree <- drop.tip(current_tree, setdiff(tree_species, matched_species))
    rownames(data_matched) <- data_matched$wcvp_name
    
    for (model in models) {
      cat("Tree:", tree_index, "Model:", model, "Predictor:", predictor, "Response:", response, "\n")
      
      formula <- as.formula(paste(response, "~", predictor))
      fit <- phylolm(formula, data = data_matched, phy = pruned_tree, model = model)
      
      results <- append(results, list(list(
        tree_index = tree_index,
        model = model,
        predictor_trait = predictor,
        dependent_trait = response,
        AIC = fit$aic,
        fit = fit
      )))
    }
  }
}

# Save results
saveRDS(results, file = "PGLS_results_PN.rds")

# Extract AIC and logLik
result_table <- data.frame(
  tree_number = integer(),
  dependent_trait = character(),
  predictor_trait = character(),
  Model = character(),
  logLik = numeric(),
  AICc = numeric(),
  stringsAsFactors = FALSE
)

for (res in results) {
  result_table <- rbind(result_table, data.frame(
    tree_number = res$tree_index,
    dependent_trait = res$dependent_trait,
    predictor_trait = res$predictor_trait,
    Model = res$model,
    logLik = res$fit$logLik,
    AICc = res$fit$aic,
    stringsAsFactors = FALSE
  ))
}

# Calculate ΔAICc and AICc Weight
summary_table <- result_table %>%
  group_by(dependent_trait, predictor_trait, Model) %>%
  summarise(
    logLik = mean(logLik, na.rm = TRUE),
    AICc = mean(AICc, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  group_by(dependent_trait, predictor_trait) %>%
  mutate(
    min_AICc = min(AICc, na.rm = TRUE),
    delta_AICc = AICc - min_AICc,
    AICc_weight = exp(-0.5 * delta_AICc) / sum(exp(-0.5 * delta_AICc))
  ) %>%
  ungroup()

# Rename and export
colnames(summary_table) <- c("dependent_trait", "predictor_trait", "Model", "logLik", "AICc", "min_AICc", "△AICc", "AICc weight")
write.csv(summary_table, "AICc_weight_summary.csv", row.names = FALSE)

# According to the AICc_weight_summary.csv results, the model with the biggest AIC wight was selected for statistics
best_models <- summary_table %>%
  group_by(dependent_trait, predictor_trait) %>%
  filter(`AICc weight` == max(`AICc weight`, na.rm = TRUE)) %>%
  select(dependent_trait, predictor_trait, Model)

print(best_models)

# Initialize list to store final results
final_results <- list()

for (i in seq_len(nrow(best_models))) {
  target_dep <- best_models$dependent_trait[i]
  target_pred <- best_models$predictor_trait[i]
  target_model <- best_models$Model[i]
  
  # Match relevant results from 'results'
  subset_results <- lapply(results, function(res) {
    if (res$dependent_trait == target_dep &&
        res$predictor_trait == target_pred &&
        res$model == target_model) {
      
      fit <- res$fit
      
      # Ensure coefficient and vcov exist
      if (length(fit$coefficients) < 2 || is.null(fit$vcov) || nrow(fit$vcov) < 2) return(NULL)
      
      beta <- fit$coefficients[2]
      se <- sqrt(diag(fit$vcov))[2]
      t_val <- beta / se
      df <- fit$n - length(fit$coefficients)
      p_val <- 2 * pt(-abs(t_val), df)
      r2_adj <- fit$adj.r.squared
      resid_se <- sqrt(fit$sigma2)  # or use: sqrt(fit$sigma2_error)
      fstat <- tryCatch(t_val^2, error = function(e) NA)  # F = t² for single predictor
      
      data.frame(
        tree_index = res$tree_index,
        dependent_trait = target_dep,
        predictor_trait = target_pred,
        model = target_model,
        Estimate = beta,
        Std_Error = se,
        t_value = t_val,
        P_value = p_val,
        R2 = r2_adj,
        Residual_SE = resid_se,
        F_statistic = fstat
      )
    } else {
      NULL
    }
  })
  
  # Filter out NULLs and bind
  filtered_results <- do.call(rbind, subset_results[!sapply(subset_results, is.null)])
  
  # Average summary
  averaged <- filtered_results %>%
    summarise(
      across(c(Estimate, Std_Error, t_value, P_value, R2, Residual_SE, F_statistic), mean, na.rm = TRUE)
    )
  
  # Add trait/model info
  averaged$dependent_trait <- target_dep
  averaged$predictor_trait <- target_pred
  averaged$model <- target_model
  
  final_results[[i]] <- averaged
}

# Combine all to final table
final_df <- do.call(rbind, final_results)
write.csv(final_df, "Final_Best_Model_Averaged_Stats_PN.csv", row.names = FALSE)

## Analysis3: Association between functional traits and plastid genome size / structural variables ####

### Eta (η) ######################################################
# Load required package
library(readxl)
library(dplyr)
library(DescTools)

# Load data
data <- read_excel("10kp_structure_detial.xlsx", sheet = 3)

# Define continuous variables
continuous_vars <- c("All_Length", "LSC_Length", "Total_IR_Length", "CDS_Length", "Speacer_Length")

# Convert continuous variables to numeric
data[, continuous_vars] <- lapply(data[, continuous_vars], function(x) as.numeric(as.character(x)))

# Define categorical (discrete) variables
discrete_vars <- c("Lifecycle", "Deciduousness", "Growth_form", "Photosynthetic_pathway",
                   "Lifecycle_imput", "Deciduousness_imput", "Growth_form_imput", "Photosynthetic_pathway_imput")

# Initialize result table
result_df <- data.frame(
  Continuous = character(),
  Discrete = character(),
  Eta_squared = numeric(),
  Eta = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

# Loop over all variable combinations
for (cont_var in continuous_vars) {
  for (disc_var in discrete_vars) {
    # Filter non-missing values
    clean_data <- data[!is.na(data[[cont_var]]) & !is.na(data[[disc_var]]), ]
    clean_data[[disc_var]] <- as.factor(clean_data[[disc_var]])
    
    # Fit ANOVA and compute eta²
    formula <- as.formula(paste(cont_var, "~", disc_var))
    model <- aov(formula, data = clean_data)
    eta_sq_result <- EtaSq(model, type = 1)
    
    eta2 <- as.numeric(eta_sq_result[1, "eta.sq"])
    eta <- sqrt(eta2)
    p_value <- summary(model)[[1]]$`Pr(>F)`[1]
    
    # Save results
    result_df <- rbind(result_df, data.frame(
      Continuous = cont_var,
      Discrete = disc_var,
      Eta_squared = eta2,
      Eta = eta,
      P_value = p_value,
      stringsAsFactors = FALSE
    ))
  }
}

# Output results
write.csv(result_df, file = "Eta_result_TP.csv", row.names = FALSE)

### OLS model ######################################################
# Load required package
library(readxl)

# Load data
data <- read_excel("10kp_structure_detial.xlsx", sheet = 3)

# Define variables
traits_to_fit <- c("N_mean","C_mean", "leaf_size",
                   "N_mean_imput","C_mean_imput", "leaf_size_imput")
trait_to_response <- c("All_Length", "LSC_Length", "Total_IR_Length", "CDS_Length", "Speacer_Length")

# Prepare result dataframe
ols_results <- data.frame(
  Predictor = character(),
  Response = character(),
  Estimate = numeric(),
  Std_Error = numeric(),
  t_value = numeric(),
  P_value = numeric(),
  R2adj = numeric(),
  Residual_SE = numeric(),
  F_statistic = numeric(),
  stringsAsFactors = FALSE
)

# Loop over predictor traits
for (predictor in traits_to_fit) {
  for (response in trait_to_response) {
    # Extract and clean data
    df <- data[, c(response, predictor)]
    df <- df[complete.cases(df), ]
    
    if (nrow(df) < 5) next  # skip if insufficient data
    
    # Apply log-transformation to predictor (add a small constant to avoid log(0) if needed)
    df[[predictor]] <- log(df[[predictor]])
    
    # Fit OLS model
    formula <- as.formula(paste(response, "~", predictor))
    fit <- lm(formula, data = df)
    summary_fit <- summary(fit)
    
    estimate <- coef(summary_fit)[2, 1]
    std_error <- coef(summary_fit)[2, 2]
    t_val <- coef(summary_fit)[2, 3]
    p_val <- coef(summary_fit)[2, 4]
    r2adj <- summary_fit$adj.r.squared
    resid_se <- summary_fit$sigma
    f_stat <- summary_fit$fstatistic[1]
    
    ols_results <- rbind(ols_results, data.frame(
      Predictor = predictor,
      Response = response,
      Estimate = estimate,
      Std_Error = std_error,
      t_value = t_val,
      P_value = p_val,
      R2adj = r2adj,
      Residual_SE = resid_se,
      F_statistic = f_stat
    ))
  }
}

# Save results
write.csv(ols_results, file = "OLS_results_summary_TP.csv", row.names = FALSE)

## Analysis4: Association between plastid genome size / structural variables and latitude####

### First- to third-degree polynomial regression models ######################################################
# Load required packages
library(readxl)
library(dplyr)

# Load data
data <- read_excel("10kp_structure_detial.xlsx", sheet = 15)

# Convert all columns except the first to numeric
data[, -1] <- lapply(data[, -1], function(x) as.numeric(as.character(x)))

# Define dependent variables (genome structure features)
variables <- c("All_Length", "LSC_Length", "Total_IR_Length", "CDS_Length", "Speacer_Length")

# Initialize result storage
result_df <- data.frame(
  Variable = character(),
  Degree = integer(),
  P_Value = numeric(),
  F_Statistic = numeric(),
  R2_Adj = numeric(),
  AIC = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each response variable
for (var in variables) {
  # Fit polynomial models of degree 1 to 3
  for (degree in 1:3) {
    # Remove rows with missing values in response or predictor
    clean <- data[!is.na(data[[var]]) & !is.na(data$Latitude), ]
    
    # Construct formula for polynomial regression
    formula <- as.formula(paste(var, "~ poly(Latitude,", degree, ")"))
    
    # Fit linear model
    lm_model <- lm(formula, data = clean)
    
    # Extract model summary and ANOVA
    model_summary <- summary(lm_model)
    anova_result <- anova(lm_model)
    
    # Extract statistics
    p_value <- anova_result$`Pr(>F)`[1]
    f_statistic <- anova_result$`F value`[1]
    r2_adj <- model_summary$adj.r.squared
    aic_value <- AIC(lm_model)
    
    # Save results
    result_df <- rbind(result_df, data.frame(
      Variable = var,
      Degree = degree,
      P_Value = p_value,
      F_Statistic = f_statistic,
      R2_Adj = r2_adj,
      AIC = aic_value,
      stringsAsFactors = FALSE
    ))
  }
}

# Output results
write.csv(result_df, file = "Latitude_pattern_result.csv", row.names = FALSE)

## Analysis5: Association between plastid genome size and key climatic factors ####

### Variance Inflation Factor ######################################################
# Load required packages
library(car)
library(tidyr)
library(readxl)
library(dplyr)

# Load data
data <- read_excel("10kp_structure_detial.xlsx", sheet = 20)  
data[, 1:11] <- lapply(data[, 1:11], as.numeric)
clean_data <- data %>% drop_na()  

# Extract weights and predictors separately
Weight <- clean_data[1]  # First column assumed to be weights
VIF_data <- clean_data[, -1]  # Remove weight column, keep predictors

# Initial full model with weighted least squares
full_model <- lm(All_Length ~ ., data = VIF_data, weights = Weight$Weight)

# Compute initial VIF values
VIF_values <- vif(full_model)
print(VIF_values)

# Iteratively remove high-VIF variables and recalculate
# Step 1: Remove "BIO1"
VIF_data <- VIF_data %>% select(-"BIO1")
model_step1 <- lm(All_Length ~ ., data = VIF_data, weights = Weight$Weight)
VIF_values <- vif(model_step1)
print(VIF_values)

# Step 2: Remove "Mean.UVB1"
VIF_data <- VIF_data %>% select(-"Mean.UVB1")
model_step2 <- lm(All_Length ~ ., data = VIF_data, weights = Weight$Weight)
VIF_values <- vif(model_step2)
print(VIF_values)

# Step 3: Remove "BIO3"
VIF_data <- VIF_data %>% select(-"BIO3")
final_model <- lm(All_Length ~ ., data = VIF_data, weights = Weight$Weight)
VIF_values <- vif(final_model)
print(VIF_values)

# Save final VIF plot
png("VIF_result_env.png", width = 20, height = 12, units = "cm", res = 300)
par(mar = c(2, 10, 2, 2))  # Set margins: bottom, left, top, right

barplot(
  VIF_values,
  main = "VIF Values of Climatic Predictors",
  horiz = TRUE,
  col = "#7bbabe",
  las = 1,
  beside = TRUE,
  space = 0.5,
  names.arg = names(VIF_values),
  cex.names = 0.8,
  xlim = c(0, 5),
  cex.axis = 0.8
)

abline(v = 5, col = "darkred", lty = 2)  # Threshold line
dev.off()

### Backward step wise regression ######################################################
# Load required packages
library(MASS)
library(readxl)

# Load data
data <- read_excel("10kp_structure_detial.xlsx", sheet = 21)
data[, -1] <- lapply(data[, -1], as.numeric)
data_clean <- na.omit(data)

# Define response variable and candidate predictors
response_variable <- "All_Length"
predictor_variables <- c("GST", "GSL", "BIO12", "GDD0", "BIO7", "BIO7_2", "PAR", "PAR_2")

# Construct the full model formula
full_formula <- as.formula(paste(response_variable, "~", paste(predictor_variables, collapse = " + ")))

# Fit the full model using a Gamma distribution with log link
full_model <- glm(full_formula, data = data_clean, weights = data_clean$Weight, family = Gamma(link = "log"))

# Perform backward stepwise selection based on AIC
final_model <- stepAIC(full_model, direction = "backward", trace = TRUE)

# Output final model summary
summary(final_model)

## Analysis6: the relative importance of structural, functional trait, climatic, and phylogenetic factors in plastome size variation ####

### Variance Inflation Factor (VIF) Analysis ######################################################
# Load required libraries
library(readxl)
library(car)
library(dplyr)

# Load data
data <- read_excel("10kp_structure_detial.xlsx", sheet = 22)
data[, -1] <- lapply(data[, -1], as.numeric)
clean_data <- na.omit(data)

# Extract weights and predictors separately
Weight <- clean_data[1]  # First column assumed to be weights
VIF_data <- clean_data[, -1]  # Remove weight column, keep predictors

# Fit initial Gamma GLM model with log link
full_model <- glm(All_Length ~ ., data = VIF_data, weights = Weight$Weight, family = Gamma(link = "log"))

# Compute VIF values
VIF_values <- vif(full_model)
print(VIF_values)

# Stepwise manual removal of high-VIF variables
# Step 1: Remove "C3"
VIF_data <- VIF_data %>% select(-"C3")

# Fit reduced linear model (OLS used here for VIF calculation)
reduced_model <- lm(All_Length ~ ., data = VIF_data, weights = Weight$Weight)

# Recompute VIF
VIF_values <- vif(reduced_model)
print(VIF_values)

# Save VIF values plot
png("VIF_result_all.png", width = 20, height = 12, units = "cm", res = 300)
par(mar = c(2, 10, 2, 2))  # Set plot margins: bottom, left, top, right

barplot(
  VIF_values,
  main = "VIF Values of All Selected Predictors",
  horiz = TRUE,
  col = "#7bbabe",
  las = 1,
  beside = TRUE,
  space = 0.5,
  names.arg = names(VIF_values),
  cex.names = 0.8,
  xlim = c(0, 5),
  cex.axis = 0.8
)

# Add threshold line at VIF = 5
abline(v = 5, col = "darkred", lty = 2)

dev.off()

### GLM Regression and Variance Partitioning ######################################################
# Load required packages
library(readxl)
library(glmm.hp)
library(dplyr)
library(broom)
library(gridExtra)

# Load data
data <- read_excel("10kp_structure_detial.xlsx", sheet = 23)
data[, -1] <- lapply(data[, -1], as.numeric)
clean_data <- na.omit(data)

# Define a list of GLM models for different plastome components
models <- list(
  All_Length  = glm(All_Length ~ GST + BIO12 + PAR + BIO7 + C4 + PAR2 + CAM + HET + phy,
                    data = clean_data, weights = clean_data$Weight, family = Gamma(link = "log")),
  
  LSC_Length  = glm(LSC_Length ~ GST + BIO12 + PAR + BIO7 + C4 + PAR2 + CAM + HET + phy,
                    data = clean_data, weights = clean_data$Weight, family = Gamma(link = "log")),
  
  IR_Length   = glm(IR_Length ~ GST + BIO12 + PAR + BIO7 + C4 + PAR2 + CAM + HET + phy,
                    data = clean_data, weights = clean_data$Weight, family = Gamma(link = "log")),
  
  CDS_Length  = glm(CDS_Length ~ GST + BIO12 + PAR + BIO7 + C4 + PAR2 + CAM + HET + phy,
                    data = clean_data, weights = clean_data$Weight, family = Gamma(link = "log")),
  
  Speacer_Length = glm(Speacer_Length ~ GST + BIO12 + PAR + BIO7 + C4 + PAR2 + CAM + HET + phy,
                       data = clean_data, weights = clean_data$Weight, family = Gamma(link = "log"))
)

# Define independent variable groupings
iv <- list(
  c4 = c("C4"),
  HET = c("HET"),
  CAM = c("CAM"),
  Phy = c("phy"),
  GST = c("GST"),
  BIO12 = c("BIO12"),
  BIO7 = c("BIO7"),
  PAR = c("PAR", "PAR2")
)

iv_1 <- list(
  Functional_trait = c("C4", "HET", "CAM"),
  Phy = c("phy"),
  Climate = c("GST", "BIO12", "PAR", "BIO7", "PAR2")
)

# Perform hierarchical partitioning using glmm.hp
glmm_results_raw <- lapply(models, function(m) glmm.hp(m, iv))
glmm_results_grouped <- lapply(models, function(m) glmm.hp(m, iv_1))

# Extract model coefficients and significance
result_tables <- list()

for (model_name in names(models)) {
  model <- models[[model_name]]
  
  coef_df <- tidy(model) %>%
    mutate(
      Model = model_name,
      Significant = ifelse(p.value < 0.05, "Yes", "No")
    ) %>%
    dplyr::select(Model, term, estimate, std.error, statistic, p.value, Significant) %>%
    rename(
      Variable = term,
      Estimate = estimate,
      Std.Error = std.error,
      Statistic = statistic,
      P.Value = p.value
    )
  
  result_tables[[model_name]] <- coef_df
}

# Combine all coefficient results into one table
final_table <- bind_rows(result_tables)

# Output coefficient table
print(final_table)
write.csv(final_table, file = "Glmm_all_factors.csv", row.names = FALSE)

### Structural equation model ######################################################
# Load required packages
library(readxl)
library(lavaan)
library(semPlot)

# Load data
data <- read_excel("10kp_structure_detial.xlsx", sheet = 24)
data[, 1:14] <- lapply(data[, 1:14], as.numeric)

# Remove rows with missing values
clean_data <- na.omit(data)

# Log-transform and standardize all variables except the first column (Weight)
clean_data_log <- clean_data
clean_data_log[, -1] <- log(clean_data_log[, -1] + 1)

clean_data_scaled <- clean_data_log
clean_data_scaled[, -1] <- scale(clean_data_scaled[, -1])

# Define SEM model structure
model <- '
All_Length ~ CDS_Length + Speacer_Length + IR_Length + LSC_Length +
             BIO7 + GST + PAR + BIO12 + C4 + HET + CAM + phy
CDS_Length ~ LSC_Length + IR_Length BIO7 + GST + PAR + BIO12 + C4 + HET + CAM + phy
Speacer_Length ~ LSC_Length + IR_Length BIO7 + GST + PAR + BIO12 + C4 + HET + CAM + phy
IR_Length ~ BIO7 + GST + PAR + BIO12 + C4 + HET + CAM + phy
LSC_Length ~ BIO7 + GST + PAR + BIO12 + C4 + HET + CAM + phy
HET ~ BIO7 + GST + PAR + BIO12 + phy
C4 ~ BIO7 + GST + PAR + BIO12 + phy
CAM ~ BIO7 + GST + PAR + BIO12 + phy
phy ~ BIO7 + GST + PAR + BIO12

GST ~~ BIO12
GST ~~ BIO7
GST ~~ PAR
BIO12 ~~ BIO7
BIO12 ~~ PAR
BIO7 ~~ PAR
LSC_Length ~~ IR_Length
CDS_Length ~~ Speacer_Length
C4 ~~ HET
C4 ~~ CAM
HET ~~ CAM

'

# Fit SEM model with weighted least squares and sampling weights
fit <- sem(
  model,
  data = clean_data_scaled,
  sampling.weights = "Weight",
  estimator = "WLS"
)

# Summarize model fit and standardized estimates
summary(fit, fit.measures = TRUE, standardized = TRUE)
summary(fit, standardized = TRUE, rsquare = TRUE)

# Modification indices: identify possible model improvements
modificationIndices(fit, sort = TRUE, minimum.value = 3.84)

# Plot standardized path diagram
semPaths(
  fit,
  what = "std",               # Show standardized path coefficients
  layout = "spring",          # Automatic spring layout
  residuals = FALSE,          # Hide residual variances
  edge.label.cex = 0.9,       # Edge label size
  nCharNodes = 5,             # Node label abbreviation
  curvePivot = TRUE,          # Smooth curves to reduce overlap
  as.edges = "all"            # Draw all edges
)

## Analysis7: Sampling biases ####

### Bootstrap Resampling ######################################################
# Load necessary libraries
library(dplyr)
library(tidyr)
library(readxl)

# Load data
data <- read_excel("10kp_structure_detial.xlsx", sheet = 24)
data[, -1] <- lapply(data[, -1], as.numeric)
data <- na.omit(data)

# Calculate observed (non-zero) regional means
observed_summary <- data %>%
  group_by(TDWG_level_3_code) %>%
  summarise(across(1:13, ~ mean(.x[.x != 0], na.rm = TRUE), .names = "mean_{col}"))

# Define columns to include in resampling
cols_to_analyze <- colnames(data)[2:14]  # Exclude region code

# Function for bootstrap resampling
bootstrap <- function(df, ratio) {
  df %>%
    group_by(TDWG_level_3_code) %>%
    summarise(across(all_of(cols_to_analyze), list(
      b_mean = ~ mean(sample(.x[.x != 0], size = floor(length(.x) * ratio), replace = TRUE), na.rm = TRUE)
    ), .names = "{col}_{fn}"))
}

# Set resampling parameters
set.seed(1234)
sample_ratio <- 0.5

# Perform bootstrap resampling with 1000, 5000, 10000 replicates
bootstrap_1000 <- replicate(1000, bootstrap(data, sample_ratio), simplify = FALSE)
bootstrap_5000 <- replicate(5000, bootstrap(data, sample_ratio), simplify = FALSE)
bootstrap_10000 <- replicate(10000, bootstrap(data, sample_ratio), simplify = FALSE)

# Function to summarize results across replicates
compute_means_by_group <- function(bootstrap_list) {
  bind_rows(bootstrap_list) %>%
    group_by(TDWG_level_3_code) %>%
    summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE), .names = "{col}_avg"))
}

# Summarize bootstrap results
summary_1000 <- compute_means_by_group(bootstrap_1000)
summary_5000 <- compute_means_by_group(bootstrap_5000)
summary_10000 <- compute_means_by_group(bootstrap_10000)

# Combine summaries
bootstrap_summary <- bind_rows(
  "bootstrap_1000" = summary_1000,
  "bootstrap_5000" = summary_5000,
  "bootstrap_10000" = summary_10000,
  .id = "Bootstrap_Replicate"
)

# Save summary table
write.csv(bootstrap_summary, "bootstrap_summary.csv", row.names = FALSE)

# Compute mean of bootstrapped means per region
compute_means_and_medians_by_group <- function(bootstrap_results) {
  combined_results <- bind_rows(bootstrap_results)
  
  summary_stats <- combined_results %>%
    group_by(TDWG_level_3_code) %>%
    summarise(across(
      where(is.numeric),
      list(mean_of = ~ mean(.x, na.rm = TRUE)),
      .names = "{col}_{fn}"
    ))
  
  return(summary_stats)
}

# Compute summaries for each bootstrap replicate set
summary_1000 <- compute_means_and_medians_by_group(bootstrap_1000)
summary_5000 <- compute_means_and_medians_by_group(bootstrap_5000)
summary_10000 <- compute_means_and_medians_by_group(bootstrap_10000)

# Combine summaries and add identifier
total_summary <- bind_rows(
  "bootstrap_1000" = summary_1000,
  "bootstrap_5000" = summary_5000,
  "bootstrap_10000" = summary_10000,
  .id = "bootstrap_set"
)

# Save combined bootstrap summary
write.csv(total_summary, "bootstrap_grouped_summary.csv", row.names = FALSE)

# Relative Deviation Calculation
bootstrap_summary <- read.csv("bootstrap_grouped_summary.csv")
actual_results <- read.csv("structure_mean_by_region.csv")

# Prepare column name alignment
colnames(actual_results)[colnames(actual_results) == "TDWG_level_3_code"] <- "TDWG_level_3_code"

# Merge for comparison
merged_data <- merge(bootstrap_summary, actual_results, by = "TDWG_level_3_code")

# Define trait names (assumes matching column structure)
traits <- c("All_Length", "LSC_length", "Total_IR_Length", "CDS_Length", "Speacer_Length", "Tandem_repeat_Length")

# Calculate relative deviations from bootstrap mean_of_means
for (trait in traits) {
  bootstrap_col <- paste0(trait, "_b_mean_mean_of")
  actual_col <- trait
  rel_dev_col <- paste0(trait, "_rel_dev")
  
  if (bootstrap_col %in% colnames(merged_data)) {
    merged_data[[rel_dev_col]] <- abs(merged_data[[bootstrap_col]] - merged_data[[actual_col]]) / 
      merged_data[[actual_col]] * 100
  }
}

# Save result
write.csv(merged_data, file = "absolute_deviation_results.csv", row.names = FALSE)

### Null Model Testing ##########################################################
# Load required package
library(dplyr)
library(tidyr)
library(ggplot2)

# Generate and Save Null Model Functions for Random Assignment
save_null_model_script <- function() {
  dir.create("D:/10KP", showWarnings = FALSE)
  file_path <- file.path("D:/10KP", "run_zero_model.R")
  
  script_content <- c(
    "library(dplyr)",
    "library(tidyr)",
    "",
    "# Simulate null distributions via random grid assignment",
    "calculate_random_distribution <- function(sample_data, grid_data, n_repeats = 999) {",
    "  sample_names <- sample_data$SampleName",
    "  grid_ids <- grid_data$GridID",
    "",
    "  if (length(sample_names) == 0) stop(\"No sample names provided.\")",
    "  if (length(grid_ids) == 0) stop(\"No grid IDs provided.\")",
    "",
    "  random_results <- vector(\"list\", n_repeats)",
    "",
    "  for (i in 1:n_repeats) {",
    "    random_assignment <- sample(sample_names, length(grid_ids), replace = TRUE)",
    "    assigned <- data.frame(GridID = grid_ids, SampleName = random_assignment)",
    "    merged <- merge(assigned, sample_data, by = 'SampleName')",
    "    random_results[[i]] <- merged %>%",
    "      group_by(GridID) %>%",
    "      summarise(across(2:14, mean, na.rm = TRUE))",
    "  }",
    "",
    "  bind_rows(random_results, .id = 'Simulation') %>%",
    "    group_by(GridID) %>%",
    "    summarise(across(2:14, list(mean = ~ mean(.x, na.rm = TRUE), sd = ~ sd(.x, na.rm = TRUE)),",
    "                     .names = \"{col}_{fn}\"))",
    "}",
    "",
    "# Calculate standardized effect sizes (SES)",
    "calculate_ses <- function(observed_data, random_summary) {",
    "  observed_data %>%",
    "    left_join(random_summary, by = 'GridID') %>%",
    "    rowwise() %>%",
    "    mutate(across(2:14, ~ (. - get(paste0(cur_column(), '_mean'))) / get(paste0(cur_column(), '_sd')),",
    "                  .names = 'SES_{col}')) %>%",
    "    ungroup()",
    "}",
    "",
    "# Main wrapper function",
    "run_zero_model <- function(sample_file, grid_file, observed_file, n_repeats = 999) {",
    "  sample_data <- read.csv(sample_file)",
    "  grid_data <- read.csv(grid_file)",
    "  observed_data <- read.csv(observed_file)",
    "  random_summary <- calculate_random_distribution(sample_data, grid_data, n_repeats)",
    "  ses_results <- calculate_ses(observed_data, random_summary)",
    "  return(ses_results)",
    "}"
  )
  
  writeLines(script_content, con = file_path)
}

# Save the script to file
save_null_model_script()

# Load and export sample-by-trait matrix
data <- read_excel("10kp_structure_detial.xlsx", sheet = 25)
samples <- data[c(5,12,15,21,25:29,31:35)]
colnames(samples)[colnames(samples) == "Species"] <- "SampleName"
write.csv(samples, "Null_model_1.5wamples_structure.csv", row.names = FALSE)

# Export grid ID file
grid <- read.csv("D:/10KP/grid_file.csv")
write.csv(grid, "Null_model_grid.csv", row.names = FALSE)

# Load and format observed region-level summary
observed <- read.csv("structure_mean_by_region.csv")
observed <- observed[c(2:15)]
colnames(observed)[colnames(observed) == "TDWG_level_3_code"] <- "GridID"
write.csv(observed, "Null_model_observed.csv", row.names = FALSE)

# Run the null model
source("D:/10KP/run_zero_model.R")

ses_results <- run_zero_model(
  sample_file = "Null_model_1.5wamples_structure.csv",
  grid_file = "Null_model_grid.csv",
  observed_file = "Null_model_observed.csv",
  n_repeats = 999
)

# Save SES results
write.csv(ses_results, file = "Null_ses_result.csv", row.names = FALSE)
