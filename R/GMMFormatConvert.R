#' GMM format convert
#' @description convert the compressed variables to the Gaussian mixtured model format.
#'
#' @param dim dimensions of Gaussian distribution for each cluster
#' @param c a list that store the weight, mean and sigma of each cluster
#'
#' @return a list of informations of each cluster
#' @export
#'
GMMFormatConvert <- function(dim,c){
  numcmp=length(c$w);
  a=c$w;
  mu=c$supp[1:dim,];
  start_ind = dim+1;
  end_ind = dim+dim*dim
  sigma = as.vector(unlist(c$supp[start_ind:end_ind,]))
  dim(sigma) = c(dim,dim,numcmp)
  return(list(numcmp = numcmp, a = a, mu = mu, sigma = sigma))
}
