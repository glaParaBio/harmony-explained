# Intermezzo

## kmeans

```
set.seed(1234)
x <- rbind(matrix(rnorm(100, sd = 0.3), ncol = 2),
         matrix(rnorm(100, mean = 1, sd = 0.3), ncol = 2))
colnames(x) <- c("x", "y")
cl <- kmeans(x, 2)
plot(x, col=cl$cluster, pch=19, xlab='gene 1', ylab='gene 2')
mtext(side=3, adj=0,
    text=sprintf('Within cluster sum of squares: %s', paste(round(cl$withinss, 1), collapse=', ')))
points(cl$centers, col = 1:2, pch = 8, cex = 3)
```

Using `base::kmeans` we get similar results, but not exactly, depending on the
seed. Cluster names are arbitrary so they differ.

MEMO: In K-means you first decide the number of clusters you want, then you
find the cluster centroids that minimise the within cluster sum of squares.
That is, you assign cells to clusters in such way that the sum of distances
between cells and the cluster centroid they belong to is smallest.

## Regress-out

```
set.seed(1234)
N <- 10
batch <- rep(LETTERS[1:3], each=N)
gex <- 10 + as.numeric(as.factor(batch)) + rnorm(length(batch))
dat <- data.table(batch, gex)
gg <- ggplot(data=dat, aes(x=batch, y=gex, colour=batch)) +
    geom_jitter(width=0.2) +
    ylab('Gene expression')
```

Regress gene expression (or whatever) on batch:

```
fit <- lm(gex ~ batch, data=dat)
Coefficients:
            Estimate 
(Intercept)   10.617 
batchB         1.265 
batchC         1.995 

Corrected B = gex - batchB
Corrected C = gex - batchC
```

The corrected values are the residuals from the model, optionally added to the
intercept to have them aligned with the baseline (batch A in this case)

```
dat[, corr_gex := fit$residuals + fit$coefficients[1]]

gg <- ggplot(data=dat, aes(x=batch, y=corr_gex, colour=batch)) +
    geom_jitter(width=0.2)
```

Equivalently, remove from each raw value the coefficient estimated for that
batch. Leave batch unchanged if you want to have that as the baseline. This is
closed to what Harmony does.

```
dat[, corr_gex2 := ifelse(batch == 'A', gex,
    ifelse(batch == 'B', gex - fit$coefficients[2],
        ifelse(batch == 'C', gex - fit$coefficients[3], NA)))]

gg <- ggplot(data=dat, aes(x=batch, y=corr_gex2, colour=batch)) +
    geom_jitter(width=0.2) +
    ylab('Gene expression')
```
