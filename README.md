
# GeM-LR<a/></a>

**GeM-LR(Generative Mixture of Logistic Regression)** is a package for predictive clustering. 
It works well with small data sets, and at the same time performs well in both predictive results and interpretability.

## Installation
------------------------------------------------------------------------

You can use the **remote** packages, a lightweight replacement of the install_* functions in devtools, to install our package as follows:
```r
install.packages("remote")
remotes::install_github("llin-lab/GeMLR")
```


## Example
------------------------------------------------------------------------
This example uses the VAST dataset in the \data folder, which contains 18 immune features.
```r
library(GeM-LR)

result = read_data(dat_road = "data\\VASTd0_Indi.txt");
# load some necessary variables of the dataset
dim = result$dim;
MLMoption = result$MLMoption;
ncmp = result$ncmp;
nseeds = result$nseeds;
numcmp = result$numcmp;
numdata = result$numdata;
rangeSeed = result$rangeSeed;
rawdat = result$rawdat;
vargmm = result$vargmm;
vlasso = result$vlasso;
X1 = result$X1;
X1s = result$X1s;
Y1 = result$Y1
kkk = result$kkk
Indi = result$Indi

# use cross-validation to choose the seed with best performance
result2 = runCV(kkk, ncmp,nseeds,rangeSeed,vargmm, Y1,X1, Indi, MLMoption)
cv_mean = apply(result2$cvAUCfinal,2,mean)
cvAUC = result2$cvAUCfinal

# find the final model for every cluster
result3 = finalModel(cvAUC,ncmp,nseeds,rangeSeed,vargmm, Y1,X1s, Indi, MLMoption)
```


## Citation
------------------------------------------------------------------------

