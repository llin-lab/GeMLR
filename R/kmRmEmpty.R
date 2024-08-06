#' kmRmEmpty
#' @description
#' to remove the empty clusters of k-means.
#'
#' @param dat the input data with every sample stored in one row
#' @param nkm the number of clusters prespecified
#' @param cdbkinit the randomly initiated centers of each cluster
#'
#' @return the center of each cluster and the identity of each sample
#' @export
kmRmEmpty <- function(dat,nkm,cdbkinit){

  len = nrow(dat)
  dim = ncol(dat)
  if (nkm > len){
    stop('Abort: number of clusters=%d is larger than data size %d',nkm,len)
  }
  ptdist = matrix(0,nrow=1,ncol=len)
  ind = matrix(0,nrow=1,ncol=len)
  cdbk = cdbkinit
  aa <- rep(0, nkm)

  for(i in 1:len){
    for (j in 1:nkm){
      aa[j] = sum((dat[i,]-cdbk[j,])^2)
    }
    ptdist[i] = min(aa)
    ind[i] = which.min(aa)
  }

  ndatpercls = matrix(0,nrow = 1,ncol = nkm)
  for (i in 1:len){
    ndatpercls[ind[i]] = ndatpercls[ind[i]]+1
  }

  numempty = sum(ndatpercls==0)
  if (numempty == 0){
    return(invisible(NA))
  }

  ptdist = sort(ptdist,decreasing = TRUE)
  id2 = order(ptdist,decreasing = TRUE)

  ii=1
  for (i in 1:nkm){
    if (!is.na(numempty) && ndatpercls[i]==0){
      cdbk[i,] = dat[id2[ii],]
      ii = ii+1
    }
  }

  for (i in 1:len){
    for (j in 1:nkm){
      aa[j] = sum((dat[i,]-cdbk[j,])^2)
    }
    ptdist[i] = min(aa)
    ind[i] = which.min(aa)
  }

  ndatpercls = matrix(0,nrow=1,ncol=nkm)
  for (i in 1:len){
    ndatpercls[ind[i]] = ndatpercls[ind[i]]+1
  }

  numempty = sum(ndatpercls==0)
  if (numempty>0){
    stop('Abort: Cannot get rid of empty cluster via kmRmEmpty, data size too small or too many identical data entries')
    return(invisible(NA))
  }

  return(list(cdbk = cdbk,ind=ind))
}
