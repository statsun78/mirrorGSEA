# mirrorGSEA

R code for GSEA (Gene Set Enrichment Analysis) using mirror statistics.

Please download `function.R` to compute mirror statistics from gene expression data.

---

## Functions

### 1. `mirror_gsea`

Performs gene set enrichment analysis after computing the mirror statistics of individual genes from gene expression data.

#### Arguments:

* `x`: An $n \times p$ matrix with $n$ samples and $p$ genes.
* `y`: A vector of binary outcome values consisting of `0` for controls and `1` for cases.
* `pathway.list`: A list of pathways annotated to the corresponding genes.
* `lam.opt`: The optimal lambda value for the regularization parameter to control model sparsity. If `NULL`, 5-fold cross-validation is performed to select the optimal value. (Default: `NULL`)
* `type`: The regularization method used to compute mirror statistics. Options are `"ridge"` for ridge regression, `"lasso"` for lasso regression, and `"elastic"` for elastic-net. (Default: `"ridge"`)
* `mirror.type`: Three different symmetric functions to compute mirror statistics. See the manuscript for details. (Default: `1`)
* `alpha`: The tuning parameter for elastic-net regularization. It must be between 0 and 1. If `type` is `"ridge"` or `"lasso"`, `alpha` is automatically set to `0` or `1`, respectively.
* `num.split`: The total number of resamplings to compute mirror statistics. (Default: `100`)

---

### 2. `sim.fun`

Generates pseudo gene expression data with an AR(1) correlation structure described in the manuscript.

#### Arguments:

* `n`: The sample size.
* `p`: The total number of genes.
* `path.size`: The pathway size. Options are `10`, `50`, `100`, or `NULL`, which correspond to small, moderate, large, and mixed pathways, respectively. For mixed pathways, 10, 50, and 100 size pathways are generated. Please see the manuscript for details.
* `global.seed`: A seed number to generate random values.
* `cov.rho`: The correlation value used for the AR(1) covariance matrix.

#### Value:

* `x`: Pseudo gene expression data represented as an $n \times p$ matrix.
* `y`: A vector of binary outcome values consisting of `0` for controls and `1` for cases.
* `sig_path`: A list of pathways associated with the outcome values.
* `sig_index`: A list of genes associated with the outcome values.

---

## Example
