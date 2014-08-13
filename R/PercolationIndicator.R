PercolationIndicator <- 
  function(data, adjacency, separate = FALSE, ...){
    Lag0Nodes <- colnames(data)
    adj <- adjacency
    col <- ncol(adj)
    
    # if adjacency matrix is not lag0-lag1, then adjust the adjacency matrix
    if (col != 2 * ncol(data)){
      adjlag01 <- matrix(0,2*col,2*col)
      adjlag01[1:col,1:col] <- adjlag01[(col+1):(2*col),(col+1):(2*col)] <- adj
      diag(adjlag01[(col+1):(2*col),(col+1):(2*col)]) <- 1
      adj <- adjlag01
    }
    
    # Discard nodes that are not connected
    NodesToAnalyze <- rowSums(adj[1:ncol(adj),1:ncol(adj)]) != 0
    if (!any(NodesToAnalyze)) stop("There is no network")
    
    x0x1 <- cbind(data[-1,],data[-nrow(data),]) # data set with lag1 variables
    x0x1 <- x0x1[,NodesToAnalyze,drop=FALSE] # without unconnected nodes
    data <- data[,NodesToAnalyze[1:(length(NodesToAnalyze)/2)],drop=FALSE]
    adj <- adj[NodesToAnalyze,NodesToAnalyze,drop=FALSE] # without unconnected nodes
    
    # counting the infected neighbors of healthy (not infected) nodes
    k <- x0x1%*%adj
    k[x0x1==1] <- 0
    A <- sum(k)
    A_x <- colSums(k[,1:ncol(data)])
    
    # counting upward jumps (from 0 to 1, but only when at least one neighbor was infected)
    d <- rep(0,ncol(data))
    for(j in 1:(nrow(data) - 1)) for(i in 1:ncol(data))
      if(data[j,i]==0 & data[j + 1,i]==1) d[i] <- d[i] + 1
    U <- sum(d)
    U_x <- d
    
    # counting infected sites
    B <- sum(data)
    B_x <- colSums(data)
    
    # counting downward jumps (1 to 0)
    e <- rep(0,ncol(data))
    for(j in 1:(nrow(data) - 1)) for(i in 1:ncol(data))
      if(data[j,i]==1 & data[j+1,i]==0) e[i] <- e[i] + 1
    D <- sum(e)
    D_x <- e
    
    # computing parameters
    inf.rate <- U/A
    inf.rate_x <- U_x/A_x
    inf.rate_x[is.nan(inf.rate_x)]  <- 0
    inf.rate_x[is.infinite(inf.rate_x)] <- 0
    recov.rate <- D/B
    recov.rate_x <- D_x/B_x
    recov.rate_x[is.nan(recov.rate_x)] <- 0
    recov.rate_x[is.infinite(recov.rate_x)] <- 0
    perc.ind <- inf.rate/recov.rate
    perc.ind_x <- inf.rate_x/recov.rate_x
    
    # results
    res <- list(inf.rate = inf.rate, recov.rate = recov.rate, perc.ind = perc.ind)
    res_x <- list(inf.rate = inf.rate_x, recov.rate = recov.rate_x, perc.ind = perc.ind_x)
    class(res) <- class(res_x) <- "PercolationIndicator"
    ifelse (separate==TRUE, result <- res_x, result <- res)
    return(result)
  }


## Methods:
print.PercolationIndicator <- function(x, ...)
{
  cat("perc.ind:\n")
  
  print(x$perc.ind)
  
  cat("\n\n inf.rate:\n")
  
  print(x$inf.rate)  
  
  cat("\n\n recov.rate:\n")
  
  print(x$recov.rate)  
}

