#' Compute the pairwise squared Wasserstein distance between two GMM model
#'
#' @param dim the dimension of the data
#' @param supp1 mu and sigma in model1
#' @param supp2 mu and sigma in model2
#'
#' @return The distance matrix
#' @export
#'
GaussWasserstein <- function(dim, supp1, supp2) {
  library(pracma)

  numcmp1 <- dim(supp1)[2]
  numcmp2 <- dim(supp2)[2]
  pairdist <- matrix(0, nrow = numcmp1,ncol = numcmp2)

  for (ii in 1:numcmp1) {
    for (jj in 1:numcmp2) {
      sigma1 <- pracma::Reshape(supp1[(dim+1):(dim+dim*dim),ii],dim,dim)
      sigma2 <- pracma::Reshape(supp2[(dim+1):(dim+dim*dim),jj],dim,dim)

      # b2=sqrtm_eig(sigma2);
      b1 <- sqrtm_eig(sigma1)
      b2 <- sqrtm_eig(sigma2)

      mudif <- supp1[1:dim,ii] - supp2[1:dim,jj]
      pairdist[ii,jj] = sum(mudif*mudif) + sum((b1-b2)*(b1-b2))

      # a=sigma1+sigma2-2*sqrtm(b1*sigma2*b1);
      # pairdist(ii,jj)=sum(mudif.*mudif)+sum(diag(a));
    }
  }
  return(pairdist)
}
