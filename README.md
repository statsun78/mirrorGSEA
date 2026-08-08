# mirrorGSEA

An R package/script for Gene Set Enrichment Analysis (GSEA) using mirror statistics.

Please download `function.R` to compute mirror statistics from gene expression data.

---

## Functions

### 1. `mirror_gsea`

Performs Gene Set Enrichment Analysis (GSEA) after computing the mirror statistics of individual genes from gene expression data.

#### Arguments

* `x` (*matrix*): An $n \times p$ matrix containing expression levels for $n$ samples and $p$ genes.
* `y` (*vector*): A binary outcome vector (`0` for controls, `1` for cases).
* `pathway.list` (*list*): A list of pathways containing sets of gene indices/names.
* `lam.opt` (*numeric*, default: `NULL`): The optimal $\lambda$ regularization parameter controlling model sparsity. If `NULL`, 5-fold cross-validation is automatically performed to select the optimal value.
* `type` (*character*, default: `"ridge"`): The regularization method used to compute mirror statistics. Options:
  * `"ridge"`: Ridge regression
  * `"lasso"`: Lasso regression
  * `"elastic"`: Elastic-net regression
* `mirror.type` (*integer*, default: `"1"`): Choice among three different symmetric functions to compute mirror statistics. Options : `"1"`, `"2"` or `"3"`. See the manuscript for details.
* `alpha` (*numeric*): The elastic-net mixing parameter ($0 \le \alpha \le 1$). Automatically set to `0` for `"ridge"` and `1` for `"lasso"`.
* `num.split` (*integer*, default: `100`): The total number of resamplings to compute mirror statistics.

#### Returns

* Returns a data frame identical to the output format of the `fgsea` R package.

---

### 2. `sim.fun`

Generates synthetic gene expression data with an AR(1) correlation structure as described in the manuscript.

#### Arguments

* `n` (*integer*): Sample size.
* `p` (*integer*): Total number of genes.
* `path.size` (*integer*, default: `NULL`): Pathway size. Options are `10`, `50`, `100`, or `NULL` (mixed sizes combining 10, 50, and 100).
* `global.seed` (*integer*): Random seed for reproducibility.
* `cov.rho` (*numeric*): Correlation coefficient ($\rho$) for the AR(1) covariance structure.

#### Returns

A list with the following components:
* `x`: Synthetic gene expression matrix ($n \times p$).
* `y`: Binary outcome vector (`0` for controls, `1` for cases).
* `sig_path`: List of pathways genuinely associated with the outcome.
* `sig_index`: List of individual genes associated with the outcome.

---

## Example Usage

```r
# Source required functions
source("function.R")

# 1. Initialize parameters
n <- 100
p <- 1000
size <- 50
rho <- 0.5

# 2. Simulate gene expression data
data <- sim.fun(n = n, p = p, path.size = size, global.seed = 123, cov.rho = rho)
x <- data$x
y <- data$y

# 3. Define pathways for GSEA
pathways <- split(1:p, rep(1:(p / size), each = size))
names(pathways) <- paste0("Pathway_", size, "_", 1:(p / size))

# 4. Perform GSEA using mirror statistics
g <- mirror_gsea(x, y, pathways)
head(g[order(g$pval), ])
```

---

## Reference
Kim, Y. and Sun, H. (2026) Mirror statistics-based gene set enrichment analysis for correlated gene expression data, *submitted*.
