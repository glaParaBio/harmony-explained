<!-- vim-markdown-toc GFM -->

* [What and why integrating](#what-and-why-integrating)
* [Method overview](#method-overview)
* [Setup](#setup)
    * [Helpers](#helpers)
* [Data setup](#data-setup)
    * [Initial clustering](#initial-clustering)
* [Intermezzo](#intermezzo)
    * [kmeans](#kmeans)
    * [Regress-out](#regress-out)
    * [Evaluating initial cluster diversity](#evaluating-initial-cluster-diversity)
* [Maximum-diversity soft-clustering](#maximum-diversity-soft-clustering)
    * [Built-in clustering](#built-in-clustering)
    * [By hand](#by-hand)
    * [Update cluster centroids](#update-cluster-centroids)
* [Correction](#correction)
    * [Mixture of Experts model](#mixture-of-experts-model)
    * [Cell specific corrections](#cell-specific-corrections)
    * [Correction](#correction-1)

<!-- vim-markdown-toc -->

This is a rehash of the excellent [detailed
walk-through](https://htmlpreview.github.io/?https://github.com/immunogenomics/harmony/blob/master/doc/detailedWalkthrough.html)
of the harmony integration method. Mostly for my own (dariober) understanding. The original paper is [Fast, sensitive and accurate integration of single-cell data with Harmony](https://www.nature.com/articles/s41592-019-0619-0).

# What and why integrating

Different datasets from the same or similar biological replicates may have
systematic differences of technical nature. That is, there may be clusters that
are the biologically the same across datasets, but technical differences
sepratae them. We would like to merge these datasets to remove this technical
variation while keeping intact the biological differences.

* **Input**
    * Cells from multiple dataset embedded in lower dimensional space of gene
    expression. I.e. a matrix of cells (rows) and top N PCs of gene expression
    (columns). NB: All cells from all dataset are concatenated and embeded as if they were a single dataset
    * Metadata assigning each cell to its respective dataset and possibly to other covariates

* **Main tuning parameters**: 
    * **sigma**: How much weight we want to give to distance of a cell from cluster centroids
    * **theta**: How much weight we want to give to cluster diversity
    * **nclust**: Number of clusters to work with (should be more or less the number of real bioloigical clusters)
    * **lambda**: Shrinkage of dataset coefficients (how heavily we want to correct)

* **Output**: Same cell embedding as input with cell coordinates corrected to remove batch effect(s)

# Method overview

* Cluster cells
* To each cell give a score for each cluster. This score depends on: How distant is the cell from the cluster centroid and how diverse is the cluster
* For each dimension (e.g. each principal component) and for each cluster regress the cell coordinate position against a common biological term (intercept) and a dataset specific term. The dataset term has the probability of a cell to belong to that cluster
* Correct the cell position 

# Setup 

```
mamba create -n harmony-explained
mamba activate harmony-explained
mamba install 'r-harmony=1.2.1' r-data.table r-tidyverse r-ggthemes r-ggrepel r-patchwork r-tidyr
```

## Helpers

```
source('helpers.R')
```

# Data setup

```
library(harmony)

data(cell_lines)
meta_data <- cell_lines$meta_data
V <- data.frame(cell_lines$scaled_pcs)
rownames(V) <- sprintf('cell%04d', 1:nrow(V))
colnames(V) <- sprintf('PC%02d', 1:ncol(V))
dim(V)
V[1:5,]
            PC01     PC02    PC03        PC19     PC20
cell0001  0.0028 -0.00145 -0.0064 ... 1.1e-03  0.00026
cell0002 -0.0117  0.00088  0.0009     4.1e-04 -0.00171
cell0003  0.0093 -0.00697 -0.0026     7.4e-05 -0.00066
cell0004  0.0063 -0.00252 -0.0044     9.3e-04  0.00037
cell0005  0.0085  0.00709 -0.0023 ... 1.8e-04  0.00048
```

Map cell to spherical space of radius 1

```
V_cos <- cosine_normalize(V, 1)
stopifnot(sqrt(sum(V_cos[1,]^2)) == 1) # L2 norm is just Pythagorean theorem sqrt(x^2 + y^2 + ...)


```

To get a hang of cosine distance vs euclidean distance see
(https://stats.stackexchange.com/a/508675/31142)

```
do_scatter(V, meta_data, 'dataset', no_guides = TRUE, do_labels = TRUE) + 
    labs(title = 'Colored by dataset', x = 'PC1', y = 'PC2') +
do_scatter(V, meta_data, 'cell_type', no_guides = TRUE, do_labels = TRUE) + 
    labs(title = 'Colored by cell type', x = 'PC1', y = 'PC2') +
NULL
```

## Initial clustering

```
set.seed(1)
hobj <- RunHarmony(
    data_mat = V,          # PCA embedding matrix of cells
    meta_data = meta_data, # dataframe with cell labels
    sigma = 0.1,           # Soft clustering: higher -> softer
    theta = 1,             # cluster diversity enforcement
    vars_use = 'dataset',  # variable to integrate out
    nclust = 5,            # number of clusters
    max_iter = 0,          # stop after initialization
    lambda=1,              # Penalty term for ridge regression. 0: No penalty.
    return_object = TRUE
)
```

If you want to check R kmeans matches harmony's:

```
set.seed(1)
km <- kmeans(t(hobj$Z_cos), 5)
cluster_centroids <- t(km$centers)
cluster_centroids <- cosine_normalize(cluster_centroids, 2)
cor(hobj$Y, cluster_centroids)
          1      2      3      4      5
[1,]  1.000 -0.801 -0.895 -0.653  0.800
[2,] -0.884  0.558  0.998  0.850 -0.895
[3,] -0.640  0.399  0.810  0.999 -0.906
[4,]  0.806 -0.704 -0.893 -0.907  1.000
[5,] -0.802  1.000  0.583  0.394 -0.702
```

Probability of each cell to belong to cluster *k*. Get Distance of each cell
from cluster k and rescale according to tuning param *sigma*. Sigma is set in
`RunHarmony`.

```
distance_matrix <- 2 * (1 - t(hobj$Y) %*% hobj$Z_cos)
colnames(distance_matrix) <- sprintf('cell_%04d', 1:ncol(distance_matrix))
rownames(distance_matrix) <- sprintf('k_%s', 1:nrow(distance_matrix))
cosine_normalize(distance_matrix[, 1:10], 2)
    cell_0001 cell_0002 cell_0003 cell_0004 cell_0005 cell_0006 cell_0007 cell_0008 cell_0009 cell_0010
k_1      0.54     0.054     0.698     0.655      0.59     0.071     0.694     0.022     0.043     0.689
k_2      0.41     0.594     0.181     0.217      0.17     0.605     0.219     0.607     0.618     0.236
k_3      0.33     0.574     0.254     0.224      0.05     0.610     0.324     0.522     0.579     0.295
k_4      0.63     0.042     0.643     0.683      0.71     0.029     0.600     0.085     0.051     0.616
k_5      0.18     0.559     0.052     0.085      0.33     0.506     0.076     0.593     0.527     0.052
```

Distance to score using sigma tuning param, then rescale to probability:

```
sigma <- 0.1
distance_score <- exp(-distance_matrix / sigma)
R <- apply(distance_score, 2, function(x) x/sum(x))
rownames(R) <- sprintf('k_%s', 1:nrow(R))
colnames(R) <- sprintf('cell%04d', 1:ncol(R))

options(scipen=99)
R[, 1:5] # same:
hobj$R[, 1:5]
options(scipen=0)
```

Try to use `sigma = 1` and see how `R` becomes less black-or-white.

## Evaluating initial cluster diversity

Count cells in each cluster separating by dataset. We don't really count, we
instead sum the weights of each cell according to R. 

`Phi` is a matrix of 0/1 indicating if a cell belongs to a dataset

```
hobj$Phi[, 1:10]
[1,] 1 . 1 1 . 1 1 . 1 1
[2,] . . . . 1 . . . . .
[3,] . 1 . . . . . 1 . .

observed_counts <- hobj$R %*% t(as.matrix(hobj$Phi))
colnames(observed_counts) <- sprintf('dataset_%s', 1:ncol(observed_counts))
rownames(observed_counts) <- sprintf('k_%s', 1:nrow(observed_counts))

round(observed_counts, 1)
    dataset_1 dataset_2 dataset_3
k_1     158.5       0.0       295
k_2       2.8     419.4         0
k_3       7.8     398.5         0
k_4     247.7       0.0       405
k_5     429.3       6.1         0

# Same as 
round(hobj$O, 1)
```

Expected counts: Expected number of cells in each cluster as if there is no
batch effect at all.

```
dataset_size <- rowSums(as.matrix(hobj$Phi))
cluster_props <- rowSums(hobj$O) / sum(hobj$O)
expected_counts <- t(matrix(dataset_size, ncol=1) %*% cluster_props)
colnames(expected_counts) <- sprintf('dataset_%s', 1:ncol(expected_counts))
rownames(expected_counts) <- sprintf('k_%s', 1:nrow(expected_counts))

# same as
hobj$E
```

# Maximum-diversity soft-clustering

## Built-in clustering

```
hobj$max_iter_kmeans <- 10
hobj$cluster_cpp()
```

## By hand

* Diversity score

*theta*: Diversity clustering penalty parameter. theta=0 does not encourage any
diversity. Larger values of theta result in more diverse clusters. *Theta*
decides how much weight to give diversity versus distance to cluster centroids.

The diversity matrix

```
theta <- 1 #  Also Set by user in RunHarmony
diversity <- (hobj$E / hobj$O) ^ theta 

colnames(diversity) <- sprintf('dataset_%s', 1:ncol(diversity))
rownames(diversity) <- sprintf('k_%s', 1:nrow(diversity))
diversity
    dataset_1  dataset_2     dataset_3
k_1      1.02 3742177.95          0.45
k_2     53.30       0.35 1450794326.65
k_3     18.59       0.35   33512710.10
k_4      0.94 2091631.98          0.48
k_5      0.36      24.93 2872367558.04

round(hobj$O, 1)
     [,1]  [,2] [,3]
[1,] 158.5   0.0  295
[2,]   2.8 419.4    0
[3,]   7.8 398.5    0
[4,] 247.7   0.0  405
[5,] 429.3   6.1    0
```

----

MEMO: Phi is the indicator matrix assigning cells to datasets.

Since *Phi* has 1 when a cell belongs to a dataset and 0 otherwise (1-hot
encoding), this multiplication gives to each cell the column in diversity
corresponoding to thta dataset. E.g. if a cell is from dataset 1, assign to
that cell the diversity in `diversity[,1]`, if cell is from dataset 2 then
assign `diversity[, 2]`, etc.

```
diversity_score <- diversity %*% hobj$Phi

colnames(diversity_score) <- sprintf('cell_%04d', 1:ncol(diversity_score))
rownames(diversity_score) <- sprintf('k_%s', 1:nrow(diversity_score))
diversity_score[, 1:5]
    cell_0001     cell_0002 cell_0003 cell_0004  cell_0005
k_1      1.02          0.45      1.02      1.02 3742177.95
k_2     53.30 1450794326.65     53.30     53.30       0.35
k_3     18.59   33512710.10     18.59     18.59       0.35
k_4      0.94          0.48      0.94      0.94 2091631.98
k_5      0.36 2872367558.04      0.36      0.36      24.93

Phi <- hobj$Phi
rownames(Phi) <- colnames(diversity)
colnames(Phi) <- colnames(diversity_score)
Phi[, 1:5]
          cell_0001 cell_0002 cell_0003 cell_0004 cell_0005
dataset_1         1         .         1         1         .
dataset_2         .         .         .         .         1
dataset_3         .         1         .         .         .

round(hobj$O, 1)
     [,1]  [,2] [,3]
[1,] 158.5   0.0  295
[2,]   2.8 419.4    0
[3,]   7.8 398.5    0
[4,] 247.7   0.0  405
[5,] 429.3   6.1    0
```

* new assignments are based on **distance** and **diversity**

MEMO: The diversity score comes from `exp(-distance)`, note minus sign, so
higher distance means lower score. I.e.: Higher score means closer and the cell is encouraged to belong to that cluster.

Same for diversity: higher score makes a cluster attractive.

```
R_new <- distance_score * diversity_score 
R_new <- apply(R_new, 2, function(x) x / sum(x))

# Example from this cell:
cell_id <- 2
round(cbind(dist_score=distance_score[, cell_id], clst_diversity=diversity_score[,cell_id], R_new_scaled=R_new[,cell_id]), 6)
    dist_score clst_diversity R_new_scaled
k_1      0.037           0.45 0.312708
k_2      0.000  1450794326.65 0.000005
k_3      0.000    33512710.10 0.000000
k_4      0.077           0.48 0.687201
k_5      0.000  2872367558.04 0.000085
```

## Update cluster centroids

Because we change the cell-to-cluster assignments, we need to update the coordinates of the cluster centroids:

For example, to obtain the new PC1 coordinate of cluster 1 centroid, the matrix multiplication says:

* Take the current PC1 coord of *all* cells
* For each cell take probability of a cell to belong to cluster 1
* For each cell multiply its PC1 coord by its prob of belonging to cluster 1. That is, if prob(k=1) is high, this cell contributes more to redifine the cluster centroid
* Sum all these products to get the new coordinate of cluster 1 as a weighted contribution of each cell. Weight given by the probability of a cell top belong to cluster 1

```
Y_unscaled <- hobj$Z_cos %*% t(hobj$R)
# PC1 coord of cluster 1 is Y_unscaled[1,1] and comes from
sum(hobj$Z_cos[1,] * hobj$R[1,])
# PC2 coord of cluster 1 is Y_unscaled[1,2] and comes from
sum(hobj$Z_cos[2,] * hobj$R[1,])
# Etc...
```

Rescale the new coordinates to length 1 (i.e. get back cosine dist):

```
Y_new <- cosine_normalize(Y_unscaled, 2)
rownames(Y_new) <- sprintf('PC%02d', 1:nrow(Y_new))
colnames(Y_new) <- sprintf('k%01d', 1:ncol(Y_new))
Y_new
          k1      k2       k3      k4      k5
PC01 -0.9362  0.9563  0.88525 -0.9737  0.7020
PC02  0.3017  0.1052  0.44127 -0.1291 -0.5944
PC03 -0.1526  0.2458 -0.02968  0.1717 -0.3855
...
PC18  0.0118  0.0025 -0.00631  0.0065  0.0190
PC19 -0.0044  0.0036 -0.00488  0.0039 -0.0105
PC20 -0.0283 -0.0231  0.01598  0.0195 -0.0013

# Compared to hobj$Y from initialised object
```

# Correction

In the previous section, we performed clustering in order to identify shared
groups of cells between batches. Now we make use of these groups in order to
correct the data in a sensitive way.

For each cell, we estimate how much its batch identity contributes to its PCA
scores. We then subtract this contribution from that cell's PCA scores.

## Mixture of Experts model

The theory behind this algorithm is based on the Mixture of Experts model. This
is a natural extension of linear modeling, in which each cluster is deemed an
expert and is assigned its own linear model.

We model each PC coordinate with a combination of linear factors.

Note that we work with `Z_orig`, that is the original PC coordinates, not with
the normalised, corrected `Z_cos`


```
W <- list()
Phi.moe <- as.matrix(hobj$Phi_moe)
lambda <- diag(c(hobj$lambda))
## Get beta coeeficients for all the clusters
for (k in 1:hobj$K) {
    design <- Phi.moe %*% diag(hobj$R[k, ]) # 4x2370 * 2370x2370 = 4x2370
    coeffs <- solve((design %*% t(Phi.moe)) + lambda) %*% design %*% t(hobj$Z_orig)

    rownames(coeffs) <- c('Biol_effect', sprintf('dataset_%s', 1:(nrow(coeffs)-1)))
    colnames(coeffs) <- sprintf('PC%02d', 1:ncol(coeffs))
    W[[k]] <- coeffs
}
names(W) <- sprintf('k_%s', 1:length(W))
```

The design matrix has non-zero coefficient for the biological effect and for the dataset it belongs to:

```
mat <- t(design)
colnames(mat) <- c('Biol_effect', sprintf('dataset_%s', 1:(ncol(mat)-1)))
rownames(mat) <- sprintf('cell_%04d', 1:nrow(mat))
mat[,1] <- 1
mat[1:10,]
mat[1:10,]
          Biol_effect dataset_1 dataset_2 dataset_3
cell_0001           1   1.0e+00   0.0e+00   0.0e+00
cell_0002           1   0.0e+00   0.0e+00   3.4e-14
cell_0003           1   9.9e-01   0.0e+00   0.0e+00
cell_0004           1   9.9e-01   0.0e+00   0.0e+00
cell_0005           1   0.0e+00   6.8e-06   0.0e+00
cell_0006           1   8.7e-14   0.0e+00   0.0e+00
cell_0007           1   1.0e+00   0.0e+00   0.0e+00
cell_0008           1   0.0e+00   0.0e+00   5.2e-15
cell_0009           1   3.4e-14   0.0e+00   0.0e+00
cell_0010           1   1.0e+00   0.0e+00   0.0e+00
```

So the model for say PC1 should be something like (code runs but is not right):

```
fit <- lm(hobj$Z_orig[1,] ~ 0 + mat)
```

Each item in list `W` has a table of coefficient for each cluster (so `length(W) == 5`)
Each table has 20 columns, one for each PC dimension and 4 rows: Intercept,
dataset 1, 2, and 3.

```
W[[1]][, 1:3]
            [,1]        [,2]        [,3] ... [,20]
[1,] -0.01089167  0.00350216 -0.00205706 ...       < Cell type effect (intercept) to be preserved
[2,]  0.00018787 -0.00051940 -0.00074868 ...       < Effect for cells in dataset 1 to be corrected
[3,]  0.00000067  0.00000022 -0.00000014 ...       < Effect for cells in dataset 2 to be corrected
[4,] -0.00018854  0.00051919  0.00074882 ...       < Effect for cells in dataset 3 to be corrected
W[[5]][, 1:3]
               [,1]           [,2]           [,3] ... [,20]
[1,]  0.00753574696 -0.00466306588 -0.00170242230 ... 
[2,] -0.00058120795 -0.00141331825 -0.00189154145 ... 
[3,]  0.00058120859  0.00141331809  0.00189154131 ... 
[4,] -0.00000000064  0.00000000016  0.00000000014 ... 
```

## Cell specific corrections

Now we want to correct the position of cells according to the batch specific
effect and cluster. We do this correction for each of the 20 PC dimensions.
Say we want to correct cell #5:

```
Z_i <- hobj$Z_orig[, 5] # NB: we correct the original coords
length(Z_1) # 20
```

Here `hobj$Phi_moe[, 5]` is and indicator vector of length 4 pertaining to cell
5.  It is `1 0 1 0` meaning that cell #5 has cell type effect (like all cells,
first 1) and cell #5 belongs to dataset 2 (second 1)

Each cluster contributes and estimate of the cell type and batch (dataset)
specific effect. We cumulate all these estimate by iterating for each cluster
and extracting the coefficient for that cluster. The contribution of each
cluster is weighted by the probability that cell #5 belongs to that cluster.

```
cell_id <- 5
cluster_id <- 1:hobj$K

# Initialize the matrix of coefficients. One row per coefficient, one column per dimension
pred <- matrix(0, nrow=nrow(hobj$Phi_moe), ncol=nrow(hobj$Z_orig))
for (k in cluster_id) { 
    prob_cluster_k <- hobj$R[k, cell_id]
    cat(sprintf('Cluster %s; prob: %.4f\n', k, prob_cluster_k))
    pred <- pred + W[[k]] * hobj$Phi_moe[, cell_id] * prob_cluster_k 
}
```

Since cell #5 belongs to dataset 2, rows 2 and 4 are all zeros since the effect
of dataset 1 and 3 don't apply to cell #5. The cumulative correction
(biological + batch) is the sum of these two coefficients. The result is a
vector of length 20 (the number of dimensions to correct):

```
pred <- colSums(pred)
```

## Correction

To get the amount to correct in each dimension we do the same as above except
we omit the intercept since we don't want to correct the biological effect.

```
cell_id <- 5
cluster_id <- 1:hobj$K

delta <- matrix(0, nrow=nrow(hobj$Phi), ncol=nrow(hobj$Z_orig))
for (k in cluster_id) { 
    prob_cluster_k <- hobj$R[k, cell_id]
    coeffs <- W[[k]]
    coeffs <- coeffs[2:nrow(coeffs),] # Omit intercept (do not correct biol effect)
    delta <- delta + coeffs * hobj$Phi[, cell_id] * prob_cluster_k 
}
delta <- colSums(delta)
```

Finally correct (to be done for all cells): 

```
Z_corrected <- hobj$Z_orig[, cell_id] - delta
```
