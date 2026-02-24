# Function for checking what is estimable -----------------------------------------------------------------------
#' @param designMat an X.matrix from netmeta
#' @param v a numeric vector representing a contrast or treatment of interest
#' @param verbose logical indicating if information should be printed
#' @returns logical; TRUE if v is estimable from designMat, FALSE otherwise
checkIdentifiable <- function(designMat, v, verbose = T) {
  
  # check lengths
  if(ncol(designMat) != length(v)) {
    
    stop("v must be a vector of length ncol(designMat)")
    
  }
  
  rankX <- Matrix::rankMatrix(designMat)[1]
  
  rankAug <- Matrix::rankMatrix(rbind(designMat, v))[1]
  
  estimable <- rankX == rankAug
  
  if(verbose) {
    
    print(paste0("The rank of the design matrix is ", rankX))
    print(paste0("The rank of the augmented matrix is ", rankAug))
    
    print(paste0('The contrast represented by v is ', ifelse(estimable, "", "not "), "estimable."))
    
  }
  
  return(estimable)
  
}

#' Internal function compatible with cnmaRank
#' @param M a design matrix M with column names giving component names. Equivalently, an X.matrix from netmeta object
#' @param set a character vector of components/treatments
#' @param verbose logical indicating if information should be printed
#' @returns logical; TRUE if v is estimable from designMat, FALSE otherwise
identified <- function(M, set, verbose = T, refix = 1) {
  
  ss <- length(set)
  ncomp <- ncol(M)
  
  # create V matrix
  v <- matrix(0, nrow = ss, ncol = ncomp, dimnames = list(NULL, colnames(M)))
  if(refix>0) {
    
    refcomps <- unlist(str_split(set[refix], "[+]"))
    
    for(comp in refcomps) {
      
      v[,which(str_detect(colnames(M), comp))] <- -1
      
    }
    
    setix <- seq(1, ss)[-c(refix)]
    
  } else {
    
    setix <- 1:length(set)
    
  }
  
  
  for(i in setix) {
    
    compsi <- unlist(str_split(set[i], "[+]"))
    
    for(c in 1:length(compsi)) {
      
      compix <- which(str_detect(colnames(M), compsi[c]))
      v[i,compix] <- v[i,compix] + 1 
      
    }
    
  }
  
  if(is.null(ncol(v))) {
    
    v <- matrix(v, nrow = 1)
    
  }
  
  rankX <- Matrix::rankMatrix(M)[1]
  
  estimable <- logical(nrow(v))
  
  for(i in 1:nrow(v)) {
    
    rankAug <- Matrix::rankMatrix(rbind(M, v[i,]))[1]
    
    estimable[i] <- rankX == rankAug
    
    if(verbose) {
      
      # print(paste0("The rank of the design matrix is ", rankX))
      # print(paste0("The rank of the augmented matrix is ", rankAug))
      refname <- set[refix]
      compname <- set[i]
      
      ifelse(refix > 0,
             print(paste0("The contrast ", compname," vs ", refname, " is ", ifelse(estimable[i], "", "not "), "uniquely identifiable.")),
             print(paste0("The component ", compname, " is ", ifelse(estimable[i], "", "not "), "uniquely identifiable.")))
      
      
      
    }
    
  }

  return(estimable)
  
}

# set is  a chacter vector of components or treatments - components in multicomponent treatments should be separated by a `+`
cnmaRank <- function(cnet, set, small.values, re = F, verbose = F) {
  
  if(!("netcomb" %in% class(cnet))) {
    
    stop("cnet must be of class netcomb or discomb")
    
  }
  
  ss <- length(set)
  ncomp <- cnet$c
  
  M <- cnet$X.matrix
  
  # Check if contrasts with arbitrary reference identifiable
  canest <- identified(M = M, set = set, verbose = verbose, refix = 1)
  
  if(any(!canest)) {
    
    # check identifiability of components involved in non-identifiable ones
    
    suspect <- set[!canest]
    checkcomps <- unique(unlist(str_split(suspect, "[+]")))
    newcanest <- identified(M, set = checkcomps, verbose = FALSE, refix = 0)
    
    warning(paste("The following component(s) are not identifiable:\n",
                paste(checkcomps[!newcanest], collapse = ", "),
               "\nRerun with `verbose = TRUE` to determine problematic contrasts.\nRanking not done.", collapse = ""))
    
    return(NA)
    
  }
  
  if(re) {
    
    TEtype <- "TE.random"
    seTEtype <- "seTE.random"
    
  } else {
    
    TEtype <- "TE.common"
    seTEtype <- "seTE.common"
    
  }
  
  # First need to get tE and seTE for all desired ones
  TEmat <- matrix(nrow = ss, ncol = ss, dimnames = list(set, set))
  seTEmat <- TEmat
  
  for(i in 1:ss) {
    
    comps <- netcomparison(cnet, treat1 = set[i], treat2 = set)
    
    TEmat[,i] <- comps[[TEtype]]
    seTEmat[,i] <- comps[[seTEtype]]
    
  }
  
  pscores <- cnmaPscores(TEmat,seTEmat, small.values = small.values)
  
  return(pscores)
  
}


cnmaPscores <- function (TE, seTE, small.values = "desirable", trts = NULL)  {
  n <- nrow(TE)
  n.seq <- seq_len(n)
  if (length(trts) != n) {
    if (is.null(colnames(TE))) {
      trts <- paste0("trt", n.seq)
      colnames(TE) <- colnames(seTE) <- trts
    }
    trts <- colnames(TE)
  }
  else colnames(TE) <- trts
  a_mat <- matrix(NA, nrow = n, ncol = n)
  for (i in n.seq) for (j in n.seq) if (i != j) 
    a_mat[i, j] <- TE[i, j]/seTE[i, j]
  pscores <- numeric(n)
  direction <- ifelse(small.values == "undesirable", 1, -1)
  if (n == 1) 
    pscores <- 1
  else {
    for (i in n.seq) pscores[i] <- sum(pnorm(a_mat[i, ] * 
                                               direction), na.rm = TRUE)/(n - 1)
  }
  names(pscores) <- trts
  pscores
}