<!-- vim-markdown-toc GFM -->

* [Setup](#setup)
    * [Helpers](#helpers)
* [Data setup](#data-setup)
    * [Initial clustering](#initial-clustering)
    * [Evaluating initial cluster diversity](#evaluating-initial-cluster-diversity)
* [Maximum-diversity soft-clustering](#maximum-diversity-soft-clustering)
    * [Built-in clustering](#built-in-clustering)
    * [By hand](#by-hand)
    * [Update cluster centroids](#update-cluster-centroids)
* [Correction](#correction)
    * [Mixture of Experts model](#mixture-of-experts-model)
    * [Cell specific corrections](#cell-specific-corrections)
* [5.  It is `1 0 1 0` meaning that cell #5 has cell type effect (like all cells,](#5--it-is-1-0-1-0-meaning-that-cell-5-has-cell-type-effect-like-all-cells)
    * [Correction](#correction-1)

<!-- vim-markdown-toc -->

This is a rehash of the excellent [detailed
walk-through](https://htmlpreview.github.io/?https://github.com/immunogenomics/harmony/blob/master/doc/detailedWalkthrough.html)
of the harmony integration method. Mostly for my own (dariober) understanding.

# Setup 

```
mamba create -n harmony-explained
mamba activate harmony-explained
mamba install 'r-harmony=1.2.1' r-data.table r-tidyverse r-ggthemes r-ggrepel r-patchwork r-tidyr
```

## Helpers

```
cosine_normalize <- function(X, margin) {
    if (margin == 1) {
        res <- sweep(as.matrix(X), 1, sqrt(rowSums(X ^ 2)), '/')
        row.names(res) <- row.names(X)
        colnames(res) <- colnames(X)
    } else {
        res <- sweep(as.matrix(X), 2, sqrt(colSums(X ^ 2)), '/')
        row.names(res) <- row.names(X)
        colnames(res) <- colnames(X)
    }
    return(res)
}

onehot <- function(vals) {
    t(model.matrix(~0 + as.factor(vals)))
}
```

# Data setup

```
library(harmony)

data(cell_lines)
V <- cell_lines$scaled_pcs
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

meta_data <- cell_lines$meta_data
```

To get a hang of cosine distance vs euclidean distance see
(https://stats.stackexchange.com/a/508675/31142)

## Initial clustering

```
set.seed(1)
hobj <- RunHarmony(
    data_mat = V, ## PCA embedding matrix of cells
    meta_data = meta_data, ## dataframe with cell labels
    sigma = 0.1,
    theta = 1, ## cluster diversity enforcement
    vars_use = 'dataset', ## variable to integrate out
    nclust = 5, ## number of clusters in Harmony model
    max_iter = 0, ## stop after initialization
    lambda=1, # Penalty term for ridge regression. 0: No penalty.
    return_object = TRUE ## return the full Harmony model object
)
```

Using `base::kmeans` we get similar results, but not exactly, depending on the
seed. Cluster names are arbitrary so they differ.

MEMO: In K-means you first decide the number of clusters you want, then you
find the cluster centroids that minimise the within cluster sum of squares.
That is, you assign cells to clusters in such way that the sum of distances
between cells and the cluster centroid they belong to is smallest.

```
set.seed(1)
km <- kmeans(t(hobj$Z_cos), 5)
cluster_centroids <- t(km$centers)
cor(hobj$Y, cluster_centroids)
          1      2      3      4      5
[1,]  0.800 -0.646 -0.889  1.000 -0.801
[2,] -0.895  0.841  0.999 -0.884  0.558
[3,] -0.906  0.999  0.820 -0.640  0.399
[4,]  1.000 -0.906 -0.895  0.806 -0.704
[5,] -0.702  0.395  0.572 -0.802  1.000
```

Probability of each cell to belong to cluster *k*. Get Distance of each cell
from cluster k and rescale according to tuning param *sigma*. Sigma is set in
`RunHarmony`.

```
distance_matrix <- 2 * (1 - t(hobj$Y) %*% hobj$Z_cos)
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
t(matrix(dataset_size, ncol=1) %*% cluster_props)
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

* Distance of each cell from each cluster centroid

MEMO: Y: coords of centroids; Z_cos: coords of cell in L2 normalised space

```
distance_matrix <- 2 * (1 - t(hobj$Y) %*% hobj$Z_cos)
colnames(distance_matrix) <- sprintf('cell_%04d', 1:ncol(distance_matrix))
rownames(distance_matrix) <- sprintf('k_%s', 1:nrow(distance_matrix))

distance_matrix[, 1:10]
    cell_0001 cell_0002 cell_0003 cell_0004 cell_0005 cell_0006 cell_0007 cell_0008 cell_0009 cell_0010
k_1       2.4      0.33      3.52      3.22      2.84      0.44      3.44      0.13      0.26      3.41
k_2       1.8      3.62      0.91      1.07      0.82      3.71      1.09      3.61      3.77      1.17
k_3       1.5      3.49      1.28      1.10      0.24      3.74      1.61      3.10      3.53      1.46
k_4       2.8      0.26      3.24      3.35      3.44      0.18      2.98      0.51      0.31      3.05
k_5       0.8      3.41      0.26      0.42      1.59      3.10      0.38      3.52      3.21      0.26
```

* Convert the distance to a score using the tuning param *sigma*

```
sigma <- 0.1 # also stored in hobj$sigma where this is a 1-column matrix
distance_score <- exp(-distance_matrix / sigma)
```

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

> Recall that in cluster k=1, batches 1 and 3 were well represented. Below, note
> that in that cluster (k=1), the penalties for batches 1 and 3 are relatively
> low (0.98 and 0.47). On the other hand, batch 2 gets a penalty score of 30914.
> This means that cells from batches 1 and 3 will be encouraged to move into
> cluster k=1. On the other hand, cluster k=2 is replete with batch 2. The
> penalty for batch 2 in cluster k=2 is relatively low, and noticeably smaller
> than the penalty score for batch 2 in cluster k=1. Thus, cells from batch 2
> will be discouraged from moving into cluster k=1, as this cluster has a higher
> penalty score for cells from batch 2 compared to other clusters (such as k=1).

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
round(cbind(dist_score=distance_score[, cell_id], clst_diversity=diversity_score[,cell_id], R_new=R_new[,cell_id]), 6)
    dist_score clst_diversity    R_new
k_1      0.037           0.45 0.312708
k_2      0.000  1450794326.65 0.000005
k_3      0.000    33512710.10 0.000000
k_4      0.077           0.48 0.687201
k_5      0.000  2872367558.04 0.000085
```

## Update cluster centroids

Because we change the cell-to-cluster assignments (rather, their assignment
probabilities), we need to update the coordinates of the cluster centroids:

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
```

# Correction

In the previous section, we performed clustering in order to identify shared
groups of cells between batches. Now we make use of these groups in order to
correct the data in a sensitive way.

For each cell, we estimate how much its batch identity contributes to its PCA
scores. We then subtract this contribution from that cell’s PCA scores.

## Mixture of Experts model

The theory behind this algorithm is based on the Mixture of Experts model. This
is a natural extension of linear modeling, in which each cluster is deemed an
expert and is assigned its own linear model.

We model each PC coordinate with a combination of linear factors.

MEMO: Linear regression minimizes the residual sum of squares (RSS). In ridge regression instead we minimize:

```
RSS + lambda * foreach(beta^2)
```

The solution via OLS and for ridge, respectively are:

```
beta_hat = solve(t(X) * X) * t(X) * y
beta_hat = solve(t(X) * X + lambda * I) * t(X) * y
```

The penalty *lambda* may be assigned for individual coefficient, it doesn't
have to be the same for all.

Note that we work with `Z_orig`, that is the original PC coordinates, not with
the normalised, corrected `Z_cos`


```
W <- list()
Phi.moe <- as.matrix(hobj$Phi_moe)
lambda <- diag(c(hobj$lambda))
## Get beta coeeficients for all the clusters
for (k in 1:hobj$K) {
    design <- Phi.moe %*% diag(hobj$R[k, ]) 
    W[[k]] <- solve((design %*% t(Phi.moe)) + lambda) %*% design %*% t(hobj$Z_orig)
}

# Each item in list W has a table of coefficient for each cluster (so length(W) == 5)
# Each table has 20 columns, one for each PC dimension and 4 rows: Intercept,
dataset 1, 2, and 3.

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
#5.  It is `1 0 1 0` meaning that cell #5 has cell type effect (like all cells,
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
(biological + batch) is the su of these two coefficients. The result ois a
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
