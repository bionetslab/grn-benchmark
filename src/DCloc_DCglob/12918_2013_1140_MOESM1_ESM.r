library(igraph)

#########################################
##########  Helper functions  ###########
#########################################

### Significance test for Fisher-transformed correlation coefficients
# n: number of samples
# x: first correlation to compare
# y: second correlation
#
c.test <- function(n, x, y) {
    stan.dev <- sqrt(2/(n - 3))
    dist <- x-y
    z <- abs(dist/stan.dev)
    p <- 2*(1 - pnorm(z))
    return(c(p,z))
}

### Fisher transformation
fisherr2z <- function(r) {
  z <- log((1 + r)/(1 - r)) / 2
  return(z)
}

### Inverse Fisher transformation
fisherz2r <- function(z) {
  r <- (exp(2 * z) - 1)/(1 + exp(2 * z))
  return(r)
}

#########################################
##########  DCglob() ####################
#########################################
# DESCRIPTION
# This functions calculates the differential correlation between disease
# condition A and B by comparison of the global topology of correlation
# networks.
#
# ARGUMENTS
# Mat.A: gene expression matrix for disease condition A.
# Mat.B: gene expression matrix for disease condition B.
# Mat.A and Mat.B should have the same dimension.
# Rows = genes, columns = samples.
# n.supp: number of supporting points between r.min and r.max.
# r.min: minimal correlation to be considered.
# z.max: maximal correlation to be considered.
# min.clus: minimum cluster size to considered.
#
# VALUE
# The function returns a matrix with the differential correlation of each gene
# and the condition where it shows higher correlation.
# 1 = Higher correlation in disease condition A
# 2 = Higher correlation in disease condition B
#
DCglob <- function(Mat.A, Mat.B, n.supp=200, r.min=0, r.max=fisherz2r(2.5), min.clus=3) {

    cat("Computing correlation matrices ...\n")
    ngenes <- nrow(Mat.A)
    cor.A <- cor(t(Mat.A))
    cor.B <- cor(t(Mat.B))
    for(i in 1:ngenes){
        cor.A[i,i]<-0
        cor.B[i,i]<-0
    }
    cor.A <- fisherr2z(cor.A)
    cor.B <- fisherr2z(cor.B)
    z.min <- fisherr2z(r.min)
    z.max <- fisherr2z(r.max)
    z <- seq(z.max, z.min, length.out=n.supp)
    result <- matrix(nrow=n.supp, ncol=ngenes, 0)
    rownames(result) <- z
    colnames(result) <- colnames(cor.A)

    # Computing the differentially correlated genes for the series of correlation thresholds.
    for (j in 1:n.supp) {
        z.thres <- z[j]
        cat(j, ". threshold r >= ", fisherz2r(z.thres), "\n", sep="")
        network.A <- (cor.A >= z.thres)
        network.B <- (cor.B >= z.thres)
        graph.A <- graph.adjacency(network.A, mode=c("undirected"))
        graph.B <- graph.adjacency(network.B, mode=c("undirected"))
        clus.A <- clusters(graph.A)
        clus.B <- clusters(graph.B)
        if (0 %in% clus.A$membership) clus.A$membership <- clus.A$membership + 1
        if (0 %in% clus.B$membership) clus.B$membership <- clus.B$membership + 1

        comp.A <- NULL
        comp.B <- NULL
        index.clus.A <-which(clus.A$csize >= min.clus)
        index.clus.B <-which(clus.B$csize >= min.clus)
        if (length(index.clus.A) > 0)
         for(i in 1:length(index.clus.A))
          comp.A[[i]] <- which(clus.A$membership == index.clus.A[i])
        if (length(index.clus.B) > 0)
         for(i in 1:length(index.clus.B))
          comp.B[[i]] <- which(clus.B$membership == index.clus.B[i])
        select.A <- NULL
        select.B <- NULL
        if (is.null(comp.A)) select.B <- unlist(comp.B)
        if (is.null(comp.B)) select.A <- unlist(comp.A)

        if (!is.null(comp.A) && !is.null(comp.B)) {
            index.gene.A <- setdiff(unlist(comp.A),unlist(comp.B))
            index.gene.B <- setdiff(unlist(comp.B),unlist(comp.A))
            if (length(index.gene.A) >= min.clus) {
                cor.index.A <- network.A[index.gene.A,index.gene.A]
                graph.index.A <- graph.adjacency(cor.index.A)
                subclus.A <- clusters(graph.index.A)
                if (0 %in% subclus.A$membership) subclus.A$membership <- subclus.A$membership + 1
                subclus.index.A <- which(subclus.A[[2]] >= min.clus)
                for(i in subclus.index.A) select.A <- union(select.A, index.gene.A[which(subclus.A[[1]] == i)])
            }
            if (length(index.gene.B) >= min.clus) {
                cor.index.B <- network.B[index.gene.B, index.gene.B]
                graph.index.B <- graph.adjacency(cor.index.B)
                subclus.B <- clusters(graph.index.B)
                if (0 %in% subclus.B$membership) subclus.B$membership <- subclus.B$membership + 1
                subclus.index.B <- which(subclus.B[[2]] >= min.clus)
                for(i in subclus.index.B) select.B <- union(select.B, index.gene.B[which(subclus.B[[1]] == i)])
            }
        }
        result[j, select.A] <- 1
        result[j, select.B] <- -2
    }

    # Assessing the significance of differential correlation
    cat("Averaging the results ... \n")

    clus.par <- 50
    n.samples <- length(Mat.A[1, ])
    result <- rbind(result, matrix(nrow=7 + 2*clus.par, ncol=ngenes,0))
    for(i in 1:ngenes){result[n.supp+2, i] <- sum(result[1:n.supp, i])}
    res.A <- which(result[n.supp+2,] > 1)
    res.B <- which(result[n.supp+2,] < -2)
    rownames(result)[(n.supp + 1):(n.supp + 2*clus.par + 7)] <- (n.supp + 1):(n.supp + 2*clus.par + 7)
    rownames(result)[n.supp + 2] <- "sum"
    rownames(result)[n.supp + 3] <- "p-value"
    result[n.supp+3, ] <- 1

    for(k in res.A) {
        par <- n.supp+6
        for(i in 1:n.supp){
            ind<-result[i+1,k]-result[i,k]
            if(ind== 1){
                result[par,k]<-as.numeric(rownames(result)[i+1])
                par<-par+1
            }
            if(ind== -1){
                result[par,k]<-as.numeric(rownames(result)[i])
                par<-par+1
            }
        }
    }
    for(k in res.B) {
        par <- n.supp + 6
        for(i in 1:n.supp){
            ind<-result[i+1,k]-result[i,k]
            if(ind== -2){
                result[par,k]<-as.numeric(rownames(result)[i+1])
                par<-par+1
            }
            if(ind== 2){
                result[par,k]<- as.numeric(rownames(result)[i])
                par<-par+1
            }
        }
    }

    for(i in c(res.A, res.B)) {
        index<-n.supp+8
        while(result[index, i] != 0) index<-index+1
        result[n.supp+4,i]<-index
    }
    for(i in c(res.A, res.B)) result[n.supp + 5, i] <- (result[n.supp+4, i]-(n.supp+6))/2
    for(i in c(res.A, res.B)){
        index<-0
        while(index<(result[(n.supp+5),i])){
            result[n.supp+7+clus.par+index,i] <- c.test(n.samples, result[n.supp+6+2*index,i], result[n.supp+6+2*index+1,i])[1]
            index<-index+1
        }
        result[n.supp+3,i]<-min(result[(n.supp+7+clus.par):(n.supp+7+clus.par+index-1),i])
    }

    ret <- matrix(ncol=2, nrow=ngenes, 0)
    rownames(ret) <- colnames(cor.A)
    colnames(ret) <- c("p-value","High. Corr.")
    ret[, 1] <- result[n.supp+3, ]
    ret[res.A, 2] <- 1
    ret[res.B, 2] <- 2
    return(ret)
}

#########################################
##########  DCloc()  ####################
#########################################
# DESCRIPTION
# This functions calculates the differential correlation between disease
# condition A and B by comparison of the local topology of correlation
# networks.
#
# ARGUMENTS
# Mat.A: gene expression matrix for disease condition A.
# Mat.B: gene expression matrix for disease condition B.
# Mat.A and Mat.B should have the same dimension.
# Rows = genes, columns = samples.
# n.supp: number of supporting points between r.min and r.max.
# r.min: minimal correlation to be considered.
# r.max: maximal correlation to be considered.
# min.neigh: minimum number of neighbors of a gene to be considered as
# differentially correlated.
#
# VALUE
# The function returns a matrix with the differential correlation of each gene.
# The columns of the matrix contain the absolute topological dissimilarity,
# the topological dissimilarity and the mean number of network neighbors for each
# of the disease conditions.
#
DCloc  <- function(Mat.A, Mat.B, n.supp=100, r.min=0, r.max=fisherz2r(2.5), min.neigh=3) {

    cat("Computing correlation matrices ...\n")
    ngenes <- nrow(Mat.A)
    cor.A <- cor(t(Mat.A))
    cor.B <- cor(t(Mat.B))
    for(i in 1:ngenes){
        cor.A[i,i]<-0
        cor.B[i,i]<-0
    }
    cor.A <- fisherr2z(cor.A)
    cor.B <- fisherr2z(cor.B)
    z.min <- fisherr2z(r.min)
    z.max <- fisherr2z(r.max)
    z <- seq(z.min, z.max, length.out=n.supp)

    #Computing the topological dissimilarity
    result <- list()
    for (x in c("nA", "nB", "one", "both", "d")) {
        result[[x]] <- matrix(nrow=ngenes, ncol=n.supp, 0)
        rownames(result[[x]]) <- colnames(cor.A)
        colnames(result[[x]]) <- z
    }

    for (j in 1:n.supp) {
        z.thres <- z[j]
        cat(j, ". threshold r >= ", fisherz2r(z.thres), "\n", sep="")
        network.A <- (cor.A >= z.thres)
        network.B <- (cor.B >= z.thres)
        for (i in 1:ngenes) {
            index.A <- which(network.A[i, ]==1)
            index.B <- which(network.B[i, ]==1)
            result$nA[i, j] <- length(index.A)
            result$nB[i, j] <- length(index.B)
            result$both[i, j] <- length(intersect(index.A, index.B))
            result$one[i, j] <- length(union(index.A, index.B))

            if( result$one[i, j] >= min.neigh ) {
                 result$d[i, j] <-(( 1 - result$both[i, j]/result$one[i, j] )* sign(result$nA[i, j] - result$nB[i, j]))
             }

        }
    }

    cat("Averaging of the results... \n")
    ret <- matrix(ncol=4,nrow=ngenes,0)
    colnames(ret) <- c("abs.top.dissim", "top.dissim", "mean.neigh.A", "mean.neigh.B")
    rownames(ret) <- colnames(cor.A)
    ret[ ,1] <- apply(result$d, 1, function(x){mean(abs(x))})
    ret[ ,2] <- apply(result$d, 1, mean)
    ret[ ,3] <- apply(result$nA, 1, mean)
    ret[ ,4] <- apply(result$nB, 1, mean)
    return(ret)
}
