#dependencies
library(LiblineaR)
library(parallel)
library(preprocessCore)

#main function
CIBERSORT <- function(sig_matrix, mixture_file, QN=TRUE){
  X = sig_matrix
  Y = mixture_file
  
  #read in data
  if (length(as.numeric(Y[,dim(Y)[2]])[!is.na(as.numeric(Y[,dim(Y)[2]]))])==0){
    Y = Y[,1:(dim(Y)[2]-1)];
  }
  Y = Y[match(row.names(X),row.names(Y))[!is.na(match(row.names(X),row.names(Y)))],];
  
  X <- data.matrix(X)
  Y <- data.matrix(Y)
  
  #quantile normalization of mixture file
  if(QN == TRUE){
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }
  
  #intersect genes
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX,]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY,]
  
  #standardize sig matrix
  X <- (X - mean(X)) / sd(as.vector(X))
  
  #print(nulldist)
  
  header <- c('Mixture',colnames(X))
  #print(header)
  
  output <- matrix()
  
  itor <- 1
  mixtures <- dim(Y)[2]
  
  C = heuristicC(X)
  
  #iterate through mixtures
  
  while(itor <= mixtures){

    y <- Y[,itor]

    #standardize mixture
    y <- (y - mean(y)) / sd(y)

    #run SVR core algorithm
    w = (LiblineaR(data = X, target = y, type = 11, cost = C)$W)[1:dim(X)[2]]

    #print output
    out <- c(colnames(Y)[itor],w)
    if(itor == 1) {output <- out}
    else {output <- rbind(output, out)}

    itor <- itor + 1
  }
  
  #return matrix object containing all results
  obj <- rbind(header,output)
  obj <- obj[,-1]
  obj <- obj[-1,]
  #obj = t(matrix(obj))
  obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
  rownames(obj) <- colnames(Y)
  colnames(obj) <- c(colnames(X))
  obj
}