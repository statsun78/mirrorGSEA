library(mnormt) 
library(glmnet)
library(fgsea) 

mirror_gsea <- function(x, y, pathway.list, lam.opt = NULL,
                        type = c("ridge", "lasso", "elastic"), 
                        mirror.type = c("1", "2", "3"), 
                        alpha = 0.1, num.split = 100) {
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  type <- match.arg(type)
  mtype <- match.arg(mirror.type)
  
  var_names <- colnames(x)
  Mstat <- vector("list", num.split)
  
  if (type == "ridge") {
    alpha <- 0 
  } else if (type == "lasso") {
    alpha <- 1
  }

  if (is.null(lam.opt)) {
    cvfit <- cv.glmnet(x, y, family = "binomial", alpha = alpha, 
                       nfolds = 5, nlambda = 100, intercept = FALSE, standardize = FALSE)
    lam.opt <- cvfit$lambda.min / sqrt(2) 
  }

  for (k in 1:num.split) {
    sample_index1 <- sample(x = 1:n, size = floor(0.5 * n), replace = FALSE)
    sample_index2 <- setdiff(1:n, sample_index1)
    
    X1 <- x[sample_index1, ]
    y1 <- y[sample_index1]
    X2 <- x[sample_index2, ]
    y2 <- y[sample_index2]

    beta1 <- as.vector(glmnet(X1, y1, family = "binomial", alpha = alpha, lambda = lam.opt, intercept = FALSE, standardize = FALSE)$beta)
    beta2 <- as.vector(glmnet(X2, y2, family = "binomial", alpha = alpha, lambda = lam.opt, intercept = FALSE, standardize = FALSE)$beta)
    
    signs <- sign(beta1 * beta2)
    abs_b1 <- abs(beta1)
    abs_b2 <- abs(beta2)
    
    if (mtype == "1") {
      Mstat[[k]] <- signs * 2 * pmin(abs_b1, abs_b2) 
    } else if (mtype == "2") {
      Mstat[[k]] <- signs * (abs_b1 * abs_b2)  
    } else if (mtype == "3") {
      Mstat[[k]] <- signs * (abs_b1 + abs_b2) 
    }
  } 
  
  ms <- do.call(rbind, Mstat) 
  mirror <- abs(colMeans(ms))
  names(mirror) <- var_names
  mirror <- sort(mirror, decreasing = TRUE)
  
  res <- fgsea(pathways = pathway.list, stats = mirror, scoreType = "pos", eps = 0)
return(res)
}

sim.fun <- function(n, p, path.size, global.seed, cov.rho) {
  path <- path.generate(p, size=path.size, seed=global.seed)
  cov <- cov.generate(p, cov.rho, pathway=path$pathway)
  data <- data.generate(n, p, signal_index=path$signal_index, seed=global.seed, covariance=cov)
  x <- data$x
  y <- data$y
return(list(x = x, y = y, sig_path = path$sig_path, sig_index = path$signal_index))
}


## 1. Generate pathway structure ----
## Input 
## - p : total number of genes / - size : pathway size(10 or 50 or 100 or NULL(mix))
## - seed : random seed        / - sparse : proportion of significant pathways (used only when size is not NULL)
path.generate <- function(p, size, seed, sparse=0.1){
  
  ## mixed pathway 
  if (is.null(size)) { 
    
    num1 = 55; num2 = 61; num3 = 64  #num1=size10, num2=size50, num3=size100
    
    # generate size10 pathways 
    start.10 <- seq(1, 10*num1, 10) 
    end.10 <- seq(10, 10*num1, 10) 
    size10 <- mapply(function(start, end) start:end, start = start.10, end = end.10, SIMPLIFY = FALSE)
    names(size10) <- paste0("Pathway_", 10, "_", 1:num1)
    
    # assign signal genes
    set.seed(seed)  #3개 pathway에는 signal gene이 50% 존재
    sig50_1 <- sort(sapply(0:2, function(i) sample((1+(i*10)):(10+(10*i)), 5, replace = F)))
    set.seed(seed*2)  #3개 pathway에는 signal gene이 20% 존재
    sig20_1 <- sort(sapply(0:2, function(i) sample((31+(i*10)):(40+(10*i)), 2, replace = F)))
    
    # generate size50 pathways 
    start.50 <- seq(10*num1+1, 10*num1+50*num2, 50)   
    end.50 <- seq(10*num1+50, 10*num1+50*num2, 50)
    size50 <- mapply(function(start, end) start:end, start = start.50, end = end.50, SIMPLIFY = FALSE)
    names(size50) <- paste0("Pathway_", 50, "_", 1:num2)
    
    # assign signal genes
    set.seed(seed)
    sig50_2 <- sort(sapply(0:2, function(i) sample((1+(i*50)+10*num1):(50+(50*i)+10*num1), 25, replace = F)))
    set.seed(seed*2)
    sig20_2 <- sort(sapply(0:2, function(i) sample((151+(i*50)+10*num1):(200+(50*i)+10*num1), 10, replace = F)))
    
    # generate size100 pathways
    start.100 <- seq(10*num1+50*num2+1, p, 100)   
    end.100 <- seq(10*num1+50*num2+100, p, 100)
    size100 <- mapply(function(start, end) start:end, start = start.100, end = end.100, SIMPLIFY = FALSE)
    names(size100) <- paste0("Pathway_", 100, "_", 1:num3)
    
    # assign signal genes
    set.seed(seed)
    sig50_3 <- sort(sapply(0:2, function(i) sample((1+(i*100)+10*num1+50*num2):(100+(100*i)+10*num1+50*num2), 50, replace = F)))
    set.seed(seed*2)
    sig20_3 <- sort(sapply(0:2, function(i) sample((301+(i*100)+10*num1+50*num2):(400+(100*i)+10*num1+50*num2), 20, replace = F)))
    
    pathway <- c(size10, size50, size100)   #mix인 경우 총 180개 pathway 생성
    
    # Construct signal and non-signal gene indices
    signal_index <- c(sig50_1, sig20_1, sig50_2, sig20_2, sig50_3, sig20_3)   
    no_sig_index <- setdiff(1:p, signal_index)
    
    # Identify significant and non-significant pathways
    sig_gene_sum <- sapply(pathway, function(x) sum(x %in% signal_index))
    sig_path <- pathway[which(sig_gene_sum > 0)]   #significant pathway list
    no_sig_path <- setdiff(names(pathway), names(sig_path))
    
  } else {  
    
    ## fixed pathway size(10 or 50 or 100)
    start_point <- seq(1, p-size+1, size)
    end_point <- seq(size, p, size)
    pathway <- mapply(function(start, end) start:end, start = start_point, end = end_point, SIMPLIFY = FALSE)
    
    # determine the numbers of total and significant pathways
    total_path_num <- length(pathway)
    sig_path_num <- total_path_num * sparse
    sig_50_num <- sig_20_num <- sig_path_num/2  #
    names(pathway) <- paste0("Pathway_", size, "_", 1:total_path_num)
    
    # assign signal genes
    set.seed(seed)   #signal gene 50% 포함
    sig_50_idx <- do.call(c, lapply(0:(sig_50_num-1), function(i) sort(sample((1 + i*size):(size + i*size), size*0.5, replace = F))))
    set.seed(seed*2)  #signal gene 20% 포함
    point <- size * sig_50_num  
    sig_20_idx <- do.call(c, lapply(0:(sig_20_num-1), function(i) sort(sample(((point+1) + i*size):((point+size) + i*size), size*0.2, replace = F))))
    
    # Construct signal and non-signal gene indices
    signal_index <- c(sig_50_idx, sig_20_idx) #signal gene : 350개
    no_sig_index <- setdiff(1:p, signal_index)
    
    # Identify significant and non-significant pathways
    sig_gene_sum <- sapply(pathway, function(x) sum(x %in% signal_index)) 
    sig_path <- pathway[which(sig_gene_sum > 0)]   
    no_sig_path <- setdiff(names(pathway), names(sig_path))
  }
  
  return(list(pathway = pathway, sig_path = sig_path, no_sig_path = no_sig_path,
              signal_index = signal_index, no_sig_index = no_sig_index))
}

## 2. Generate signal effects ----
## Input
## - n : number of signal genes
signal.fun <- function(n) {  
  signs <- sample(c(-1, 1), n, replace = TRUE)
  runif(n, min = 1, max = 3) * signs
}


## 3.Generate covariance matrix (Identity or AR1) ----
## Input
## - p : total number of genes / - rho, cov.rho : correlation parameter for AR(1)
## - pathway.list, pathway : pathway information
BlockDiag.cov <- function(p, rho, pathway.list) {
  
  cov.mat <- adj.mat <-  matrix(0, p, p)
  link <- list()
  
  for (k in seq_along(pathway.list)) {
    
    block <- pathway.list[[k]]
    adj.mat[block,block] <- 1
    
    for (i in block) {
      for (j in block) {
        cov.mat[i,j] <- rho^(abs(i-j))  
      } 
    }
    
    link[[k]] <- adj.mat[block, block]
  }
  
  return(list(cov.mat = cov.mat, adj.mat = adj.mat, link = link))
}

cov.generate <- function(p, cov.rho, pathway) {
  
  if (is.null(cov.rho)) {   #identity matrix
    
    cov.mat <- diag(p)
    
  } else {    #AR1 covariance matrix
    
    cov.mat <- BlockDiag.cov(p, rho = cov.rho, pathway.list = pathway)$cov.mat
  }
  
  return(cov.mat)
}


## 4. Generate simulation data(x,y) ----
## Input
## - n : sample size    / - p : total number of genes   / - signal_index : indices of signal genes
## - seed : random seed / - covariance : covariance matrix
data.generate <- function(n, p, signal_index, seed, covariance){
  
  # generate true beta coefficients
  beta <-  numeric(p)
  set.seed(seed)
  beta[signal_index] <- signal.fun(length(signal_index))
  
  # generate design matrix X and continuous response
  set.seed(seed)
  x2 = rmnorm(2*n, rep(0,p), varcov = covariance)
  mu = x2 %*% beta
  y0 = mu + rnorm(2*n)
  
  # convert continuous response to binary response using quartile
  qq <- quantile(y0, c(0.25, 0.75))
  by <- rep(NA, 2*n)  
  by[y0 <= qq[1]] <- 0
  by[y0 > qq[2]] <- 1
  y <- by[!is.na(by)]
  x <- x2[!is.na(by),]
  colnames(x) <- as.character(1:p) 
  
  return(list(x = x, y = y))
}
