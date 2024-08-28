
# GeMLR<a/></a>

**GeMLR(Generative Mixture of Logistic Regression)** is a package for predictive clustering. 
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
library(GeMLR)
```

```r
result = read_data(dat_road = "data\\VASTd0_Indi.txt");
```
Make sure that your input data (for example, VAST there) have the right format.It means that the first column in your data is vaccinated or not, and the last column means infected or not. The rest columns are some features that you want to use in your GMM and logitstic model.

```r
# There are some requires about the data format in read_data().
# The first column indicates vaccination, and the last column indicates infection.
> head_data = head(read.table( "data\\VASTd0_Indi.txt"))
> head_data
"V1"	"V2"	"V3"	"V4"	"V5"	"V6"	"V7"	"V8"	"V9"	"V10"	"V11"	"V12"	"V13"	"V14"	"V15"	"V16"	"V17"	"V18"	"V19"	"V20"
0	5.33	1.08	4.11	5.6	1.3	5.32	0.85	5.43	0.85	5.72	1.46	5.72	5.52	1	5.84	4.21	4.05	1.2	0
1	6.35	1.28	3.89	6.52	1.59	5.9	1.28	6.36	1.26	5.18	1.18	5.35	3.64	0.3	6.23	4.1	3.37	1.18	1
1	5.52	0.78	3.5	5.47	1.04	5.58	0.78	3.48	0.3	4.67	1.2	5.54	3.95	0.7	4.53	4.33	3.95	0.3	1
0	5.72	1.26	5.53	6.35	1.41	6.12	1.11	5.13	0.7	6.02	1.71	7.03	5.37	1.11	6.44	4.56	3.82	1.32	1
0	5.96	1.66	5.12	6.02	1.6	6.3	0.9	5.57	0.48	6.33	2	7.18	5.87	1.58	6.96	4.69	4.18	1.41	1
1	3.59	0.7	3.87	4.39	1.04	4.01	0.78	2.9	0.7	3.52	1.45	5.41	4.18	1	4.13	4.49	3.49	0	0
```
After run the read_data(), you will acquire some ingredients for your model in the 'list' format (which means you need to unpack them by yourself) .

```r
# load some necessary variables of the dataset
dim = result$dim; # The number of all features except vaccination and infection.in VAST data,dim = 18.
ncmp = result$ncmp; # The number of clusters available for selection.In VAST data, ncmp = [2,3,4]
nseeds = result$nseeds; # The number of random seeds for k-means initiation.
numdata = result$numdata; # The number of samples in your data.
rangeSeed = result$rangeSeed; # The maximum value of random seeds.
rawdat = result$rawdat; # The dataframe that have all samples and all features.
vargmm = result$vargmm; # The indexs of the top 5 highly variance features among all features for GMM.
vlasso = result$vlasso; # The value of lambda that come with minimum bias in lasso.
X1 = result$X1; # The columns(except first and last) in rawdat, means all features you want to use in the model.
X1s = result$X1s; # standardized X1.
Y1 = result$Y1 # The last column in rawdat, means infected or not.
kkk = result$kkk # The number of folds in cross validation.
Indi = result$Indi # The first column in rawdat, means vaccinated or not.
```

Apart from the data ingredients above, you also need to initiate some parameters (in GeMLR, we call it 'MLMoption') for your algorithm. All parameters are packed in a list called MLMoption.
```r
# read necessary parameters for model in MLMopton
MLMoption = init_MLMoption()
```
In our package, we set all parameters in MLMoption. If you want to change it for your need, you should read the annotations and change the values carefully.



```r
# use cross-validation to choose the seed with best performance
result2 = runCV(kkk, ncmp,nseeds,rangeSeed,vargmm, Y1,X1, Indi, MLMoption)
cv_mean = apply(result2$cvAUCfinal,2,mean)
cvAUC = result2$cvAUCfinal

```r
# find the final model for every cluster
result3 = finalModel(cvAUC,ncmp,nseeds,rangeSeed,vargmm, Y1,X1s, Indi, MLMoption)
```


## Citation
------------------------------------------------------------------------

