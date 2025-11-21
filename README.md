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
  ycol = 20,
  Indi_col = 1
)
```

**Important:**

- `dat_path` must be the full path to a data file. Supports `.txt`, `.csv`, `.tsv`, `.xlsx`, `.rds`, `.sav`, `.sas7bdat`, `.dta` with column headers.
- `ycol` is the column index (or name) for the binary outcome (e.g. infection). If NULL, uses the last column.
- `Indi_col` is the column index (or name) for indicator variable(s) (e.g. vaccine group). Can be:
  - A single value: `Indi_col = 1` or `Indi_col = "vaccination"`
  - A vector: `Indi_col = c(1, 5)` for multiple indicators
  - `NULL`: Auto-detects all binary variables (excluding Y) as indicators. **Warning**: This may include unwanted variables (e.g., gender, age groups). **Recommended to specify explicitly!**
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


### Step 4. Prepare Model Parameters

```r
# Compute vargmm and vlasso
params <- compute_gemlr_params(
  X = X,
  Y = Y,
  Indi = Indi,   
  num_gmm = 5
)

# Extract parameters
vargmm <- params$vargmm
vlasso <- params$vlasso
```
- `num_gmm` defines how many top-variable features (by variance) to use in the Gaussian Mixture Model (GMM)-based clustering. 
  
  - If you set `num_gmm` to a positive integer (recommended is 5), the function will select the top `num_gmm` variables with the highest variance for GMM input.
  - If you set `num_gmm` = 0, you must provide a separate `gmm_var` argument to specify which variables to use. `gmm_var` can be a vector of column names or column indices.
  - If you omit `num_gmm` entirely, all available variables (excluding Indi and Y) will be used in GMM by default.

```r
# Initialize model parameters
MLMoption <- init_MLMoption(
  alphaLasso = 0.8, 
  vlasso = vlasso,       
  numcmp = 1,            
  stopratio = 1.0e-5,
  verbose = 1, 
  minloop = 3, 
  maxloop = 50,         
  constrain = "DIAS",
  diagshrink = 0.9, 
  kmseed = 0, 
  algorithm = 1, 
  kappa = -1,
  AUC = 1, 
  DISTR = "binomial", 
  NOEM = 0, 
  Yalpha = 1.0
)
```

All model settings are specified within `MLMoption`. You may modify these values based on your study requirements, but we recommend doing so carefully and with reference to the original documentation.

**Key parameters explained:**

* `vargmm`: Column indices of features used for GMM clustering 
* `vlasso`: Lasso penalty strength λ 
* `alphaLasso`: Elastic net mixing parameter (1=Lasso, 0=Ridge, 0.5=equal mix)
* `stopratio`: Convergence threshold controlling the number of EM iterations (default: `1.0e-5`)
* `kappa`: Controls whether sample weights are used. Default is `-1`, which disables weighting
* `verbose`: Verbosity flag (default: `1` = show messages)
* `minloop`: Minimum number of EM iterations (default: `3`, must be ≥2)
* `maxloop`: Maximum number of EM iterations
* `constrain`: Covariance structure for GMMs (default: `'DIAS'`). Options include `'N'` (no constraint), `'EI'`, `'VI'`, `'EEE'`, `'VVV'`, `'DIA'`, `'DIAE'`, `'DIAS'`, `'EEV'`, `'VEV'`
* `diagshrink`: Shrinkage toward diagonal in constrained models (default: `0.9`, only used for `'DIAS'`)
* `algorithm`: Model fitting method (1 = Lasso-regularized logistic regression; 0 = logistic without variable selection)
* `numcmp`: Number of clusters (default: `2`)
* `AUC`: Selection metric for best seed (1 = AUC, 0 = accuracy)
* `DISTR`: Distribution (default: `'binomial'` for classification; 'normal' for regression)
* `NOEM`: Whether to use EM algorithm. If set to `1`, disables EM updates and only runs initialization

At this point, all the raw materials needed to build the model are ready.

**Alternative Quick Start:** If you already know the optimal number of clusters K (from prior analysis or domain knowledge), you can skip Steps 5-6 and directly use the `fit_model()` function for quick fitting:
```r
# Quick model fitting without cross-validation

fit <- fit_model(
  X = Xs,              # Standardized features from read_data()
  Y = Y,
  Indi = Indi,
  K = 3,               # Specify number of clusters
  vargmm = vargmm,
  vlasso = vlasso,       # Automatically estimates lambda (step 4.2); or specify your own
  nseeds = 10,         # Number of random initializations
  alphaLasso = 0.8,
  verbose = 0          # Set to 1 to see progress
)

# View results
print(fit$metrics)     # Accuracy, AUC, logloss
print(fit$beta)        # Cluster-specific coefficients

# Visualize
plot_beta_heatmap(fit$beta)

```

Otherwise, continue with **Step 5** to use cross-validation to determine the optimal number of clusters.


***

### Step 5. Select Optimal Model via Cross-Validation

To determine the optimal number of clusters, you can use the built-in cross-validation function `runCV()`.
```r
# Perform cross-validation to select the optimal number of clusters
result2 <- runCV(
  k = 5,                   # Number of folds used in cross-validation (user defined)
  ncmp = c(2, 3, 4),       # Number of clusters to evaluate (user defined)
  nseeds = 20,             # Number of random seeds used in k-means (user defined)
  rangeSeed = 30,          # Maximum range of seeds to draw from (user defined)
  vargmm = vargmm,           
  Y = Y, 
  X = X,                 
  Indi = Indi, 
  MLMoption = MLMoption
)

# Extract and summarize CV results
cvAUC   <- result2$cvAUCfinal
cv_mean <- apply(cvAUC, 2, mean)
```

**Details**:

* `k`: Number of folds used in cross-validation.
* `ncmp`: A vector of possible cluster counts to consider (e.g., `c(2, 3, 4)`). Users may define this based on prior knowledge or modeling goals.
* `nseeds`: Number of k-means initializations per candidate model.
* `rangeSeed`: The upper bound of the random seed range used to draw `nseeds`.


`runCV()` returns a list that includes `cvAUCfinal`, a matrix of AUC values:
  
- Each **column** corresponds to a different number of clusters (`ncmp`).
- Each **row** represents one fold of cross-validation.
```r
# Example output for ncmp = c(2, 3, 4)
> result2$cvAUCfinal
         cluster=2 cluster=3 cluster=4
1 fold      0.5111     0.7778     0.8889
2 fold      0.6667     0.8333     0.8889
3 fold      0.6400     0.7600     0.7000
4 fold      0.6667     0.6667     0.7556
5 fold      0.5111     0.7333     0.4889
```

***

### Step 6. Fit the Final Model

Once you have identified the preferred number of clusters (based on average AUC or interpretability), you can fit the final model:
```r
# Fit the final model using the best cluster setting
result3 <- finalModel(
  cvAUC = cvAUC, 
  ncmp = c(2, 3, 4), 
  nseeds = 20, 
  rangeSeed = 30, 
  vargmm = vargmm,       
  Y = Y, 
  Xs = Xs, 
  X = X, 
  Indi = Indi, 
  MLMoption = MLMoption
)
```

### Step 7. Visualize Cluster-Specific Coefficients

You can inspect the model's coefficients for each cluster by visualizing the heatmap of β coefficients using the `plot_beta_heatmap` function:
```r
# Show the picture in the sidebar
plot_beta_heatmap(result3$beta)

# Save the picture
plot_beta_heatmap(result3$beta, output_file = "beta_heatmap.png")
```

This function will generate a heatmap where:

* Rows correspond to variables (features + Indi)
* Columns correspond to clusters
* Cell color indicates magnitude and sign of coefficients

Then you will see a plot similar to this:
![](https://github.com/llin-lab/GeMLR/blob/main/example.png "Example Image")

## Citation
------------------------------------------------------------------------
The content of this package is sourced from the following article. If you use it, please quote:

[1] Lin, Lin, et al. "GeM-LR: Discovering predictive biomarkers for small datasets in vaccine studies." PLoS computational biology 20.11 (2024): e1012581.
