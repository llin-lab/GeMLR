#' The algorithmic prototype of Wasserstein Barycenter using Bregman ADMM
#'
#' @param instanceW the weights of each sample
#' @param supp mu and sigma informations in GMM models
#' @param w the weighhts of each cluster
#' @param c0 the initialization
#' @param options MLMoptions
#'
#' @return the final GMM result
#' @export
#'
centroid_sphBregman_GMM <- function(stride, instanceW, supp, w, c0, options){

  d <- floor(sqrt(dim(supp)[1]))
  n <- length(stride)
  m <- length(w)
  posvec <- c(1,cumsum(stride)+1)

  if (is.null(c0)) {
    stop('Please give a GMM barycenter initialization,not empty.')
    #c <- centroid_init(stride, supp, w, options)
  }
  else {
    c <- c0
  }
  support_size <- length(c$w)

  X <- matrix(0, nrow = support_size, ncol = m)
  Y <- matrix(0, nrow = dim(X)[1], ncol = dim(X)[2])
  Z <- X
  spIDX_rows <- matrix(0, nrow = support_size * m, ncol = 1)
  spIDX_cols <- matrix(0, nrow = support_size * m, ncol = 1)
  for (i in 1:n) {
    temp <- pracma::meshgrid((i-1)*support_size + (1:support_size), posvec[i]:(posvec[i+1]-1))
    xx <- temp$X
    yy <- temp$Y
    ii <- support_size*(posvec[i]-1) + (1:(support_size*stride[i]))
    spIDX_rows[ii,1] <- t(xx)
    spIDX_cols[ii,1] = t(yy)
  }
  spIDX <- pracma::repmat(diag(x = 1, nrow = support_size, ncol = support_size), 1, n)

  # initialization
  suppW <- matrix(0, nrow = 1, ncol = m)
  for (i in 1:n) {
    suppW[posvec[i]:(posvec[i+1]-1)] <- instanceW[i]
    Z[,posvec[i]:(posvec[i+1]-1)] <- 1 / (support_size * stride[i])
  }
  if (d > 1) {
    C <- pracma::pdist2(t(c$supp[1:d,]), t(supp[1:d,]))^2
  }
  else {
    C <- pracma::pdist2(t(t(c$supp[1:d,])), t(t(supp[1:d,])))^2
  }
  C <- C + gaussian_wd(c$supp[(d+1):nrow(c$supp),], supp[(d+1):nrow(supp),])
  nIter <- 2000
  if (!is.null(options$badmm_max_iters)) {
    nIter <- options$badmm_max_iters
  }

  if (!is.null(options$badmm_rho)) {
    rho <- options$badmm_rho*mean(apply(C,2,median))
  }
  else {
    rho <- 10 * mean(apply(C,2,median))
  }

  C <- C * as.vector(suppW)
  rho <- rho * median(instanceW)

  if (!is.null(options$tau)) {
    tau <- options$tau
  }

  else {
    tau <- 10
  }


  if (!is.null(options$badmm_tol)) {
    badmm_tol <- options$badmm_tol
  }
  else {
    badmm_tol <- 1e-4
  }


  for (iter in 1:nIter) {
    # update X
    X <- Z * exp((C+Y)/(-rho)) + pracma::eps(x = 1)
    X <-  t(t(X) * (as.vector(w) / colSums(X)))

    # update Z
    Z0 <- Z
    Z <- X * exp(Y/rho) + pracma::eps(x = 1)
    spZ <- Matrix::sparseMatrix(i = as.vector(spIDX_rows), j = as.vector(spIDX_cols), x = as.vector(Z), dims = c(support_size * n, m))
    tmp <- rowSums(as.matrix(spZ))
    tmp <- matrix(tmp, nrow = support_size, ncol = n)
    dg <-  1/tmp * as.vector(c$w)
    dg <- Matrix::sparseMatrix(i = 1:(support_size*n), j = 1:(support_size*n), x = as.vector(dg))
    Z <- as.matrix(spIDX %*% dg %*% spZ)

    # update Y
    Y = Y + rho * (X - Z);

    # update c.w
    tmp <- t(t(tmp) * 1/colSums(tmp))
    sumW <- rowSums(sqrt(tmp))^2 # (R2)
    #sumW = sum(tmp,2)' # (R1)
    c$w <- sumW / sum(sumW)
    #c.w = Fisher_Rao_center(tmp')

    # update c.supp and compute C (lazy)
    if (iter %% tau == 0 & is.null(options$support_points)) {
      tmpX <- t(t(X) * as.vector(suppW))
      c$supp[1:d,] <- (supp[1:d,] %*% t(tmpX)) * 1 / pracma::repmat(rowSums(tmpX), d, 1)
      c$supp[(d+1):nrow(c$supp),] <- gaussian_mean(supp[(d+1):nrow(supp),], tmpX, c$supp[(d+1):nrow(c$supp),])
      if (d > 1) {
        C <- pracma::pdist2(t(c$supp[1:d,]), t(supp[1:d,]))^2
      }
      else {
        C <- pracma::pdist2(t(t(c$supp[1:d,])), t(t(supp[1:d,])))^2
      }
      C <- C + gaussian_wd(c$supp[(d+1):nrow(c$supp),], supp[(d+1):nrow(supp),])
      C = t(t(C) * as.vector(suppW))
    }



    # The constraint X=Z are not necessarily strongly enforced
    # during the update of w, which makes it suitable to reset
    # lagrangian multipler after a few iterations
    # if (mod(iter, 100) == 0)
    #          Y(:,:) = 0;
    #           if primres > 10*dualres
    #             rho = 2 * rho;
    #             fprintf(' *2');
    #           elseif 10*primres < dualres
    #             rho = rho / 2;
    #             fprintf(' /2');
    #           end
    #end

    # output
    if (iter %% 100 == 0 | iter == 10) {
      primres <- norm(X - Z, type = "F") / norm(Z, type = "F")
      dualres <- norm(Z - Z0, type = "F") / norm(Z, type = "F")
      cat(sprintf('\t %d %f %f %f \n', iter, sum(t(t(C * X) * as.vector(suppW))) / sum(instanceW),
                  primres, dualres))
      if (sqrt(dualres * primres) < badmm_tol) {
        break
      }
    }

  }
  return(list(c = c, X = X))
}
