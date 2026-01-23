# GeMLR: Generative Mixture of Logistic Regression

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
library(glmnet)
```

### Step 2. Read the data
#### Option A: Read from text file (e.g., .txt, .csv)
```r
# Replace with full path to your data file with header (column names)
result <- read_data(
  dat_path = "data/VASTd0_Indi.txt",
  ycol = 20,           # Column index for outcome
  Indi_col = 1         # Column index for indicator
)
```
#### Option B: Read from RData file (data frame)
```r
# For demonstration only - preprocess example RData
# If your data is already clean, skip to read_data() directly
load("data/VASTd0_log.RData") 
dfd0$Vaccine <- ifelse(dfd0$Vaccine == "Vi-TT", 0, 1)
dfd0$Diagnosis <- NULL

result <- read_data(
  dat = dfd0,          
  ycol = "Y",
  Indi_col = "Vaccine"
)

```
**Parameters:**

- `dat_path`: Full path to your data file. Supports `.txt`, `.csv`, `.tsv`, `.xlsx`, `.rds`, `.rda`, `.RData`, `.sav`, `.sas7bdat`, `.dta`, `.dat` with column headers.
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
  varreg = NULL,              
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
- `varreg`: Features for logistic regression. `NULL` = use all features
- `vlasso`: Lasso penalty strength. `NULL` = auto-compute via cross-validation (recommended)
- `nseeds`: Number of random initializations (10-20 recommended for stability)
- `alphaLasso`: Elastic net mixing parameter. `1` = Lasso, `0` = Ridge, `0.5` = equal mix
- `verbose`: Verbosity level. `0` = quiet, `1` = show progress

---

## Full Workflow: Cross-Validation for Optimal Number of Clusters
------------------------------------------------------------------------

If you don't know the optimal number of clusters, use cross-validation. 

### Step 4. Initial Model Fitting

First, fit an initial model to determine which features to use:

```r
# Fit initial model with variable selection
fit <- fit_model(
  X = X,
  Xs = Xs,
  Y = Y,
  Indi = Indi,
  K = 3,                
  vargmm = NULL,        
  VS = 5,               
  varreg = NULL,        
  vlasso = NULL,
  nseeds = 10,
  alphaLasso = 0.8,
  verbose = 1
)
```

### Step 5. Run Cross-Validation

Use the selected features from `fit` to evaluate different cluster numbers:

```r
# Perform cross-validation
result_cv <- runCV(
  k = 5,                   
  ncmp = c(2, 3, 4),       
  nseeds = 20,             
  rangeSeed = 30,          
  vargmm = fit$vargmm,  
  varreg = fit$varreg,  
  vlasso = fit$vlasso,            
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

### Step 6. Fit the Final Model

```r
# Fit the final model using the best cluster setting
result_final <- finalModel(
  cvAUCfinal = cvAUC,      
  ncmp = c(2, 3, 4),       
  nseeds = 20,             
  rangeSeed = 30,          
  vargmm = fit$vargmm,  
  varreg = fit$varreg,  
  vlasso = fit$vlasso,            
  Y = Y,                   
  Xs = Xs,                 
  X = X,                   
  Indi = Indi,             
  alphaLasso = 0.8,        
  verbose = 1              
)
```

**Note:** `finalModel()` automatically selects the best K based on highest mean AUC from `cvAUCfinal`.

### Step 7. Visualize Cluster-Specific Coefficients

You can inspect the model's coefficients for each cluster by visualizing the heatmap of β coefficients using the `plot_beta_heatmap` function:
```r
# Show the picture in the sidebar
plot_beta_heatmap(result_final$beta)

# Save the picture
plot_beta_heatmap(result_final$beta, output_file = "beta_heatmap.png")
```

This function will generate a heatmap where:

* Rows correspond to variables 
* Columns correspond to clusters
* Cell color indicates magnitude and sign of coefficients

Then you will see a plot similar to this:
![](https://github.com/llin-lab/GeMLR/blob/main/example.png "Example Image")

## Example: VAST Vaccine Study - Multiple Analysis Tasks
------------------------------------------------------------------------

This example demonstrates three different analysis strategies using the VAST dataset (VASTd0_log.RData), which contains subjects with 18 immune features plus a vaccine indicator.

### Data Preparation

```r

# Load and preprocess data
load("data/VASTd0_log.RData")
dfd0$Vaccine <- ifelse(dfd0$Vaccine == "Vi-TT", 0, 1)
dfd0$Diagnosis <- NULL

# Read data using GeMLR
result <- read_data(
  dat = dfd0,
  ycol = "Y",
  Indi_col = "Vaccine"
)

# Extract components
X <- result$X
Xs <- result$Xs
Y <- result$Y
Indi <- result$Indi
```

### Task 1: Single Feature Analysis with Vaccine Indicator

**Purpose:** Analyze using only the first two features with vaccine indicator.

```r
# Use first two features (columns 1-2 in X)
vargmm_1 <- c(1, 2)
varreg_1 <- c(1, 2)

# Fit initial model
fit_task1 <- fit_model(
  X = X,
  Xs = Xs,
  Y = Y,
  Indi = Indi,
  K = 2,
  vargmm = vargmm_1,
  VS = NA,
  varreg = varreg_1,
  nseeds = 3,
  alphaLasso = 0.8,
  verbose = 1
)

# Cross-validation
result_cv_1 <- runCV(
  k = 5,
  ncmp = c(2, 3, 4),
  nseeds = 20,
  rangeSeed = 30,
  vargmm = fit_task1$vargmm,
  varreg = fit_task1$varreg,
  vlasso = fit_task1$vlasso,
  Y = Y,
  X = X,
  Indi = Indi,
  alphaLasso = 0.8
)

print(result_cv_1$cvAUCfinal)

# Final model
result_final_1 <- finalModel(
  cvAUCfinal = result_cv_1$cvAUCfinal,
  ncmp = c(2, 3, 4),
  nseeds = 20,
  rangeSeed = 30,
  vargmm = fit_task1$vargmm,
  varreg = fit_task1$varreg,
  vlasso = fit_task1$vlasso,
  Y = Y,
  Xs = Xs,
  X = X,
  Indi = Indi,
  alphaLasso = 0.8,
  verbose = 1
)

# Visualize
png("VAST_task1_heatmap.png", width = 10, height = 8, units = "in", res = 300)
plot_beta_heatmap(result_final_1$beta)
dev.off()
```

**Cross-validation results:**
```r
> result_cv_1$cvAUCfinal
         cluster=2 cluster=3 cluster=4
1 fold    0.5795    0.5000    0.6591
2 fold    0.5222    0.5333    0.5111
3 fold    0.6786    0.6161    0.5893
4 fold    0.9615    0.9231    0.9615
5 fold    0.5313    0.4792    0.5417
```

### Task 2: Feature Grouping Strategy

**Purpose:** Use first half of features for GMM clustering and second half for logistic regression.

```r
# Define feature groups
group1_features <- 1:9   # First half features (columns 1-9 in X)
group2_features <- 10:18  # Second half features (columns 10-18 in X)

# Fit model
fit_task2 <- fit_model(
  X = X,
  Xs = Xs,
  Y = Y,
  Indi = Indi,
  K = 2,
  vargmm = group1_features,
  VS = NA,
  varreg = group2_features,
  nseeds = 3,
  alphaLasso = 0.8,
  verbose = 1
)

# Cross-validation
result_cv_2 <- runCV(
  k = 5,
  ncmp = c(2, 3, 4),
  nseeds = 20,
  rangeSeed = 30,
  vargmm = fit_task2$vargmm,
  varreg = fit_task2$varreg,
  vlasso = fit_task2$vlasso,
  Y = Y,
  X = X,
  Indi = Indi,
  alphaLasso = 0.8
)

print(result_cv_2$cvAUCfinal)

# Final model
result_final_2 <- finalModel(
  cvAUCfinal = result_cv_2$cvAUCfinal,
  ncmp = c(2, 3, 4),
  nseeds = 20,
  rangeSeed = 30,
  vargmm = fit_task2$vargmm,
  varreg = fit_task2$varreg,
  vlasso = fit_task2$vlasso,
  Y = Y,
  Xs = Xs,
  X = X,
  Indi = Indi,
  alphaLasso = 0.8,
  verbose = 1
)

# Visualize
png("VAST_task2_heatmap.png", width = 10, height = 10, units = "in", res = 300)
plot_beta_heatmap(result_final_2$beta)
dev.off()
```

**Cross-validation results:**
```r
> result_cv_2$cvAUCfinal
         cluster=2 cluster=3 cluster=4
1 fold    0.6705    0.6364    0.6023
2 fold    0.4556    0.5556    0.5556
3 fold    0.6786    0.5536    0.5625
4 fold    0.8846    0.8077    0.8462
5 fold    0.5521    0.6042    0.6562
```

### Task 3: Top Variable Features

**Purpose:** Select top 5 most variable features from all available features for GMM, use all features for logistic regression.

```r
# Use all features as pool, select top 5 by variance
all_features <- 1:18

# Fit model with feature selection
fit_task3 <- fit_model(
  X = X,
  Xs = Xs,
  Y = Y,
  Indi = Indi,
  K = 2,
  vargmm = all_features,
  VS = 5,
  varreg = all_features,
  nseeds = 3,
  alphaLasso = 0.8,
  verbose = 1
)

# Cross-validation
result_cv_3 <- runCV(
  k = 5,
  ncmp = c(2, 3, 4),
  nseeds = 20,
  rangeSeed = 30,
  vargmm = fit_task3$vargmm,
  varreg = fit_task3$varreg,
  vlasso = fit_task3$vlasso,
  Y = Y,
  X = X,
  Indi = Indi,
  alphaLasso = 0.8
)

print(result_cv_3$cvAUCfinal)

# Final model
result_final_3 <- finalModel(
  cvAUCfinal = result_cv_3$cvAUCfinal,
  ncmp = c(2, 3, 4),
  nseeds = 20,
  rangeSeed = 30,
  vargmm = fit_task3$vargmm,
  varreg = fit_task3$varreg,
  vlasso = fit_task3$vlasso,
  Y = Y,
  Xs = Xs,
  X = X,
  Indi = Indi,
  alphaLasso = 0.8,
  verbose = 1
)

# Visualize
png("VAST_task3_heatmap.png", width = 10, height = 10, units = "in", res = 300)
plot_beta_heatmap(result_final_3$beta)
dev.off()
```

**Cross-validation results:**
```r
> result_cv_3$cvAUCfinal
         cluster=2 cluster=3 cluster=4
1 fold    0.9091    0.6932    0.4432
2 fold    0.6222    0.6111    0.5778
3 fold    0.5536    0.5714    0.5893
4 fold    0.8462    0.7692    0.8462
5 fold    0.5625    0.7083    0.5000
```

---

## Citation
------------------------------------------------------------------------
If you use this package, please cite:

[1] Lin, Lin, et al. "GeM-LR: Discovering predictive biomarkers for small datasets in vaccine studies." PLoS computational biology 20.11 (2024): e1012581.
