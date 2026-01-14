# GeMLR: Generative Mixture of Logistic Regression<a/></a>

**GeMLR** is an R package for predictive clustering, particularly well-suited for small datasets common in vaccine studies. It simultaneously delivers strong predictive performance and interpretability by identifying latent subgroups with shared biomarker-outcome relationships.

The model integrates a Gaussian Mixture Model (GMM) for clustering and logistic regression within each cluster, allowing discovery of subgroup-specific predictive biomarkers.

## Installation
------------------------------------------------------------------------

You can use the **remotes** package, a lightweight replacement of the install_* functions in devtools, to install the package as follows:
```r
install.packages("remotes")
remotes::install_github("llin-lab/GeMLR")
```

## Example Usage
------------------------------------------------------------------------
This example walks through using the VAST dataset (included in the \data folder) which contains 18 immune features.

### Step 1. Load the package
```r
library(GeMLR)
```

### Step 2. Read the data
```r
# Replace with full path to your data file with header (column names)
result <- read_data(
  dat_path = "data/VASTd0_Indi.txt",
  ycol = 20,           # Column index for outcome
  Indi_col = 1         # Column index for indicator
)
```

**Parameters:**

- `dat_path`: Full path to your data file. Supports `.txt`, `.csv`, `.tsv`, `.xlsx`, `.rds`, `.sav`, `.sas7bdat`, `.dta` with column headers.
- `ycol`: Column index (or name) for the binary outcome (e.g., infection). If `NULL`, uses the last column.
- `Indi_col`: Column index (or name) for indicator variable(s) (e.g., vaccine group). Can be a single value (`Indi_col = 1`), a vector (`Indi_col = c(1, 5)`), or `NULL` to auto-detect. **Recommended: specify explicitly.**

```r
# Example data format
> head_data <- head(read.table("data/VASTd0_Indi.txt"))
> head_data
  V1   V2   V3   V4   V5   V6   V7   V8   V9  V10  V11  V12  V13  V14 V15  V16  V17  V18  V19 V20
   0 5.33 1.08 4.11 5.60 1.30 5.32 0.85 5.43 0.85 5.72 1.46 5.72 5.52   1 5.84 4.21 4.05 1.20   0
   1 6.35 1.28 3.89 6.52 1.59 5.90 1.28 6.36 1.26 5.18 1.18 5.35 3.64 0.3 6.23 4.10 3.37 1.18   1
   1 5.52 0.78 3.50 5.47 1.04 5.58 0.78 3.48 0.30 4.67 1.20 5.54 3.95 0.7 4.53 4.33 3.95 0.30   1
   0 5.72 1.26 5.53 6.35 1.41 6.12 1.11 5.13 0.70 6.02 1.71 7.03 5.37 1.1 6.44 4.56 3.82 1.32   1
   0 5.96 1.66 5.12 6.02 1.60 6.30 0.90 5.57 0.48 6.33 2.00 7.18 5.87 1.6 6.96 4.69 4.18 1.41   1
   1 3.59 0.70 3.87 4.39 1.04 4.01 0.78 2.90 0.70 3.52 1.45 5.41 4.18   1 4.13 4.49 3.49 0.00   0
... 
```

### Step 3. Extract Model Inputs

After running `read_data()`, the function returns a **list** of components that serve as inputs to the GeMLR model.
**You must unpack this list manually** to access the relevant variables for clustering, modeling, and visualization.
```r
# Unpack necessary components from the result list
dim     <- result$dim      # Number of immune features (excluding indicator and outcome). For VAST data, dim = 18
numdata <- result$numdata  # Number of samples (rows) in the dataset
rawdat  <- result$rawdat   # Full raw data matrix (with all columns)

# Predictor matrix
X       <- result$X        # Predictor features used in model (excluding indicator and outcome)
Xs      <- result$Xs       # Standardized version of X (continuous variables standardized, binary kept as 0/1)

# Outcome and indicators
Y       <- result$Y        # Binary outcome (default: last column)
Indi    <- result$Indi     # Indicator variable(s) (default: first column). For VAST data, means vaccinated or not.
```

---

## Quick Start: Direct Model Fitting
------------------------------------------------------------------------

If you already know the optimal number of clusters, use `fit_model()` for quick fitting:

```r
# Quick model fitting without cross-validation
fit <- fit_model(
  X = X,               # Raw features from read_data()
  Xs = Xs,             # Standardized features from read_data()
  Y = Y,
  Indi = Indi,
  K = 3,               
  vargmm = NULL,       
  VS = 5,              
  vlasso = NULL,       
  nseeds = 10,         
  alphaLasso = 0.8,    
  verbose = 1          
)

# View results
print(fit$metrics)     # Accuracy, AUC, logloss
print(fit$beta)        # Cluster-specific coefficients

# Visualize
plot_beta_heatmap(fit$beta)
```

**Parameters:**

- `K`: Number of clusters (e.g., 2, 3, or 4)
- `vargmm`: Features for GMM clustering. `NULL` = use all features (recommended for initial analysis)
- `VS`: Variable selection. `NA` = use all features in `vargmm` pool; integer (e.g., `5`) = select top 5 features by variance
- `vlasso`: Lasso penalty strength. `NULL` = auto-compute via cross-validation (recommended)
- `nseeds`: Number of random initializations (10-20 recommended for stability)
- `alphaLasso`: Elastic net mixing parameter. `1` = Lasso, `0` = Ridge, `0.5` = equal mix
- `verbose`: Verbosity level. `0` = quiet, `1` = show progress

---

## Full Workflow: Cross-Validation for Optimal Number of Clusters
------------------------------------------------------------------------

If you don't know the optimal number of clusters, use cross-validation:

### Step 4. Run Cross-Validation

```r

# Perform cross-validation
result_cv <- runCV(
  k = 5,                   
  ncmp = c(2, 3, 4),       
  nseeds = 20,             
  rangeSeed = 30,          
  vargmm = 1:ncol(X),        # Use all features for GMM 
  vlasso = NULL,           
  Y = Y,                   
  X = X,                   
  Indi = Indi,             
  alphaLasso = 0.8,        
  verbose = 1              
)

# Extract results
cvAUC <- result_cv$cvAUCfinal
cv_mean <- apply(cvAUC, 2, mean)

```

**Parameters:**

- `k`: Number of cross-validation folds (typically 5 or 10)
- `ncmp`: Vector of cluster numbers to evaluate (e.g., `c(2, 3, 4)`)
- `rangeSeed`: Maximum seed range for sampling

`result_cv$cvAUCfinal` is a matrix where:
- Each **row** = one cross-validation fold
- Each **column** = one cluster number from `ncmp`

```r
# Example output
> result_cv$cvAUCfinal
         cluster=2 cluster=3 cluster=4
1 fold      0.5111     0.7778     0.8889
2 fold      0.6667     0.8333     0.8889
3 fold      0.6400     0.7600     0.7000
4 fold      0.6667     0.6667     0.7556
5 fold      0.5111     0.7333     0.4889


```

The column with the highest mean AUC indicates the optimal K.

### Step 5. Fit the Final Model

```r
# Fit the final model using the best cluster setting
result_final <- finalModel(
  cvAUCfinal = cvAUC,      
  ncmp = c(2, 3, 4),       
  nseeds = 20,             
  rangeSeed = 30,          
  vargmm = vargmm,         
  vlasso = NULL,           
  Y = Y,                   
  Xs = Xs,                 
  X = X,                   
  Indi = Indi,             
  alphaLasso = 0.8,        
  verbose = 1              
)
```

**Note:** `finalModel()` automatically selects the best K based on highest mean AUC from `cvAUCfinal`.

### Step 6. Visualize Cluster-Specific Coefficients

You can inspect the model's coefficients for each cluster by visualizing the heatmap of β coefficients using the `plot_beta_heatmap` function:
```r
# Show the picture in the sidebar
plot_beta_heatmap(result_final$beta)

# Save the picture
plot_beta_heatmap(result_final$beta, output_file = "beta_heatmap.png")
```

This function will generate a heatmap where:

* Rows correspond to variables (features + Indi)
* Columns correspond to clusters
* Cell color indicates magnitude and sign of coefficients

Then you will see a plot similar to this:
![](https://github.com/llin-lab/GeMLR/blob/main/example.png "Example Image")

---

## Example: CVIA078 Malaria Vaccine Study
------------------------------------------------------------------------

This example demonstrates how GeMLR can be applied to different analysis tasks using the CVIA078 dataset, which includes 117 subjects, 66 antibody features (33 at V5/1 month, 33 at V20/6 months), and demographics (AGE, SEX, SICKLE).

### Data Preprocessing

```r
# Load and preprocess CVIA078 data
source("CVIA078_preprocess_corrected.R")

# After preprocessing, you will have:
# Y: Binary infection outcome (117 subjects)
# X_full_raw: Raw antibody data (117 × 66)
# X_full_scaled: Standardized antibody data (117 × 66)
# Indi_full: Demographics matrix (AGE, SEX, SICKLE)
```

### Task 1: Single Feature with Demographics

**Purpose:** Run analysis using a single antibody feature with AGE, sickle cell trait, and Sex as covariates.

```r
library(glmnet)
library(pROC)
library(caret)
library(GeMLR)

# Prepare data
X_full <- cbind(V5_features, V20_features)
X_full_raw <- X_full
X_full_scaled <- as.data.frame(scale(X_full))
colnames(X_full_scaled) <- colnames(X_full)

# Task 1: Use first 2 V5 features (minimum for GMM stability)
vargmm_idx_1 <- c(1, 2)
varreg_idx_1 <- c(1, 2)

# Fit initial model
fit_task1 <- fit_model(
  X = X_full_raw,
  Xs = X_full_scaled,
  Y = Y,
  Indi = as.matrix(Indi_full),
  K = 2,
  vargmm = vargmm_idx_1,
  VS = NA,
  varreg = varreg_idx_1,
  nseeds = 3,
  alphaLasso = 0.8,
  verbose = 1
)

vargmm_1 <- fit_task1$vargmm
varreg_1 <- fit_task1$varreg
vlasso_1 <- fit_task1$vlasso

# Cross-validation
result_cv_1 <- runCV(
  k = 5,
  ncmp = c(2, 3, 4),
  nseeds = 20,
  rangeSeed = 30,
  vargmm = vargmm_1,
  vlasso = vlasso_1,
  Y = Y,
  X = X_full_raw,
  Indi = as.matrix(Indi_full),
  alphaLasso = 0.8
)

print(result_cv_1$cvAUCfinal)

# Final model
result_final_1 <- finalModel(
  cvAUCfinal = result_cv_1$cvAUCfinal,
  ncmp = c(2, 3, 4),
  nseeds = 20,
  rangeSeed = 30,
  vargmm = vargmm_1,
  vlasso = vlasso_1,
  Y = Y,
  Xs = X_full_scaled,
  X = X_full_raw,
  Indi = as.matrix(Indi_full),
  varreg = varreg_1,
  alphaLasso = 0.8,
  verbose = 1
)

# Visualize
png("CVIA078_task1_heatmap.png", width = 10, height = 8, units = "in", res = 300)
plot_beta_heatmap(result_final_1$beta)
dev.off()
```

**Cross-validation results:**
```r
> result_cv_1$cvAUCfinal
         cluster=2 cluster=3 cluster=4
1 fold   0.5929    0.4500    0.5000
2 fold   0.6875    0.5714    0.7143
3 fold   0.7778    0.7847    0.7014
4 fold   0.5500    0.5833    0.6417
5 fold   0.6515    0.7121    0.6667
```

### Task 2: V5 for GMM, V20 for Logistic Regression

**Purpose:** Use V5 antibodies (1 month post-vaccination) for GMM clustering and V20 antibodies (6 months) for logistic regression prediction.

```r
# Define pools
v5_indices <- 1:ncol(V5_features)
v20_indices <- (ncol(V5_features) + 1):ncol(X_full)

# Fit model with feature selection
fit_task2 <- fit_model(
  X = X_full_raw,
  Xs = X_full_scaled,
  Y = Y,
  Indi = as.matrix(Indi_full),
  K = 2,
  vargmm = v5_indices,
  VS = 5,
  varreg = v20_indices,
  nseeds = 3,
  alphaLasso = 0.8,
  verbose = 1
)

vargmm_2 <- fit_task2$vargmm
varreg_2 <- fit_task2$varreg
vlasso_2 <- fit_task2$vlasso

# Cross-validation
result_cv_2 <- runCV(
  k = 5,
  ncmp = c(2, 3, 4),
  nseeds = 20,
  rangeSeed = 30,
  vargmm = vargmm_2,
  vlasso = vlasso_2,
  Y = Y,
  X = X_full_raw,
  Indi = as.matrix(Indi_full),
  alphaLasso = 0.8
)

print(result_cv_2$cvAUCfinal)

# Final model
result_final_2 <- finalModel(
  cvAUCfinal = result_cv_2$cvAUCfinal,
  ncmp = c(2, 3, 4),
  nseeds = 20,
  rangeSeed = 30,
  vargmm = vargmm_2,
  vlasso = vlasso_2,
  Y = Y,
  Xs = X_full_scaled,
  X = X_full_raw,
  Indi = as.matrix(Indi_full),
  varreg = varreg_2,
  alphaLasso = 0.8,
  verbose = 1
)

# Visualize
png("CVIA078_task2_heatmap.png", width = 10, height = 10, units = "in", res = 300)
plot_beta_heatmap(result_final_2$beta)
dev.off()
```

**Cross-validation results:**
```r
> result_cv_2$cvAUCfinal
         cluster=2 cluster=3 cluster=4
1 fold   0.5857    0.6000    0.6571
2 fold   0.5982    0.5625    0.6696
3 fold   0.5833    0.5694    0.5833
4 fold   0.6167    0.5250    0.6250
5 fold   0.5758    0.4621    0.7348
```

### Task 3: Top 5 Variable Features for GMM, Remaining for LR

**Purpose:** Select top 5 highly variable features among all V5 and V20 antibodies for GMM clustering, and use the remaining features for logistic regression.

```r
# Define pool (all features)
all_indices <- 1:ncol(X_full)

# Fit model with feature selection
fit_task3 <- fit_model(
  X = X_full_raw,
  Xs = X_full_scaled,
  Y = Y,
  Indi = as.matrix(Indi_full),
  K = 2,
  vargmm = all_indices,
  VS = 5,
  varreg = all_indices,
  nseeds = 3,
  alphaLasso = 0.8,
  verbose = 1
)

vargmm_3 <- fit_task3$vargmm
varreg_3 <- fit_task3$varreg
vlasso_3 <- fit_task3$vlasso

# Cross-validation
result_cv_3 <- runCV(
  k = 5,
  ncmp = c(2, 3, 4),
  nseeds = 20,
  rangeSeed = 30,
  vargmm = vargmm_3,
  vlasso = vlasso_3,
  Y = Y,
  X = X_full_raw,
  Indi = as.matrix(Indi_full),
  alphaLasso = 0.8
)

print(result_cv_3$cvAUCfinal)

# Final model
result_final_3 <- finalModel(
  cvAUCfinal = result_cv_3$cvAUCfinal,
  ncmp = c(2, 3, 4),
  nseeds = 20,
  rangeSeed = 30,
  vargmm = vargmm_3,
  vlasso = vlasso_3,
  Y = Y,
  Xs = X_full_scaled,
  X = X_full_raw,
  Indi = as.matrix(Indi_full),
  varreg = varreg_3,
  alphaLasso = 0.8,
  verbose = 1
)

# Visualize
png("CVIA078_task3_heatmap.png", width = 10, height = 10, units = "in", res = 300)
plot_beta_heatmap(result_final_3$beta)
dev.off()
```

**Cross-validation results:**
```r
> result_cv_3$cvAUCfinal
         cluster=2 cluster=3 cluster=4
1 fold   0.5857    0.6000    0.6571
2 fold   0.5982    0.5625    0.6696
3 fold   0.5833    0.5694    0.5833
4 fold   0.6167    0.5250    0.6250
5 fold   0.5758    0.4621    0.7348
```

---

## Citation
------------------------------------------------------------------------
If you use this package, please cite:

[1] Lin, Lin, et al. "GeM-LR: Discovering predictive biomarkers for small datasets in vaccine studies." PLoS computational biology 20.11 (2024): e1012581.
