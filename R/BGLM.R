#' @useDynLib BGLM, .registration = TRUE
NULL

# Functions
#-------------------------------------------------------------------------------

################################# Uni-trait #################################

build_blocks<-function(p,k){
  if (p==k || p<k){
    stop("The number of blocks can not be equal or greater than the covariates size (",p,")")
  }
  sizes<-rep(floor(p/k),k)
  remainder<-p%%k
  if(remainder > 0) sizes[1:remainder]<-sizes[1:remainder] + 1

  blocks<-list(); idx<-1
  for (i in 1:k){
    blocks[[i]]<-idx:(idx + sizes[i] - 1)
    idx<-idx + sizes[i]
  }
  return(blocks)
}

setLT.Fixed<-function(LT,whichNa,y,n,response_type,saveAt,intercept,i)
{
  if (any(is.na(LT$X))) {
    stop("matrix X has NAs")
  }
  if (nrow(as.matrix(LT$X)) != n) {
    stop("Number of rows of matrix X not equal to the number of phenotypes")
  }

  if(intercept[[1]]){
    LT$X = as.matrix(LT$X)
    LT$X_train = LT$X[whichNa,]
  }
  else{
    if(i == intercept[[2]]){
      LT$X = cbind(Intercept = 1,as.matrix(LT$X))
      LT$X_train = LT$X[whichNa,]
    }
    else{ LT$X = as.matrix(LT$X); LT$X_train = LT$X[whichNa,] }
  }

  LT$colNames = colnames(data.frame(LT$X))
  LT$p = ncol(LT$X)

  if(is.null(LT$update)){ LT$update = "scalar" }

  if(LT$update == "blocks"){
    if(is.null(LT$nblocks)){LT$nblocks = 1}
    LT$blocks<-build_blocks(LT$p,LT$nblocks)
  }

  if(response_type%in%c("DLN","gaussian")){
    x2 = apply(LT$X,2L,function(x) sum(x^2))
    if(LT$update == "scalar"){
      LT$x2 = x2
    }
    if(LT$update == "blocks"){
      LT$tX = t(LT$X)
      LT$tXX = crossprod(LT$X)
    }
  }
  if(response_type%in%c("Poisson","PLN")){
    x2 = apply(LT$X,2L,function(x) sum(x^2))
    if(LT$update == "scalar"){
      LT$X2 = apply(LT$X,2L,function(x) x^2)
    }
    if(LT$update == "blocks"){
      LT$tX = t(LT$X)
    }
  }

  if(LT$update == "scalar"){
    LT$priorvarB = rep(1e+10,LT$p)
  }
  if(LT$update == "blocks"){
    LT$InvpriorvarB = diag((1/1e+10),LT$p)
  }
  LT$Bv = rep(0,LT$p)
  LT$post_Bv = rep(0,LT$p)
  LT$post_Bv2 = rep(0,LT$p)

  fname = paste(saveAt,LT$Name,"_parameters.dat",sep = "")
  LT$NamefileOut = fname
  LT$fileOut = file(description = fname,open = "w")
  tmp = LT$colNames
  write(tmp,ncolumns = LT$p,file = LT$fileOut,append = TRUE)

  return(LT)

}

setLT.BRR<-function(LT,whichNa,y,n,R2,thin,nIter,burnIn,nLT,response_type,
                      saveAt,intercept,i)
{

  if (any(is.na(LT$X))) {
    stop("matrix X has NAs")
  }
  if (nrow(as.matrix(LT$X)) != n) {
    stop("Number of rows of matrix X not equal to the number of phenotypes")
  }

  if(intercept[[1]]){
    LT$X = as.matrix(LT$X)
    LT$X_train = LT$X[whichNa,]
  }
  else{
    if(i == intercept[[2]]){
      LT$X = cbind(Intercept = 1,as.matrix(LT$X))
      LT$X_train = LT$X[whichNa,]
    }
    else{ LT$X = as.matrix(LT$X); LT$X_train = LT$X[whichNa,] }
  }

  LT$colNames = colnames(data.frame(LT$X))
  LT$p = ncol(LT$X)

  if(is.null(LT$update)){ LT$update = "scalar" }

  if(LT$update == "blocks"){
    if(is.null(LT$nblocks)){LT$nblocks = 1}
    LT$blocks<-build_blocks(LT$p,LT$nblocks)
  }

  if(response_type%in%c("DLN","gaussian")){
    x2 = apply(LT$X,2L,function(x) sum(x^2))
    if(LT$update == "scalar"){
      LT$x2 = x2
    }
    if(LT$update == "blocks"){
      LT$tX = t(LT$X)
      LT$tXX = crossprod(LT$X)
    }
  }
  if(response_type%in%c("Poisson","PLN")){
    x2 = apply(LT$X,2L,function(x) sum(x^2))
    if(LT$update == "scalar"){
      LT$X2 = apply(LT$X,2L,function(x) x^2)
    }
    if(LT$update == "blocks"){
      LT$tX = t(LT$X)
    }
  }

  if (is.null(LT$R2)) {
    LT$R2 = R2/nLT
  }

  if (is.null(LT$priorvB)){
    LT$priorvB = 5
  }

  MSx<-sum(x2)/n

  if(response_type == "gaussian"){
    if (is.null(LT$priorSB)){
      LT$priorSB = (LT$priorvB+2)*(LT$R2)*var(y,na.rm=TRUE)/(MSx)
    }
  }
  if(response_type == "DLN"){
    if (is.null(LT$priorSB)){
      LT$priorSB = (LT$priorvB+2)*(LT$R2)*var(log(y+1),na.rm=TRUE)/(MSx)
    }
  }
  if(response_type %in% c("Poisson","PLN")){
    if (is.null(LT$priorSB)){
      LT$priorSB = (LT$priorvB+2)*(LT$R2)*mean(log(y+1),na.rm=TRUE)/(MSx/2)
    }
  }

  LT$Bv = rep(0,LT$p)
  LT$varB = 1/rgamma(1,LT$priorvB/2,LT$priorSB/2)
  LT$post_Bv = rep(0,LT$p)
  LT$post_Bv2 = rep(0,LT$p)
  LT$post_varB = 0
  LT$post_varB2 = 0

  fname = paste(saveAt,LT$Name,"_varB.dat",sep = "")
  LT$NamefileOut = fname
  LT$fileOut = file(description = fname,open = "w")

  return(LT)

}

setLT.BayesA<-function(LT,whichNa,y,n,R2,thin,nIter,burnIn,nLT,response_type,
                         saveAt,intercept,i)
{

  if (any(is.na(LT$X))) {
    stop("matrix X has NAs")
  }
  if (nrow(as.matrix(LT$X)) != n) {
    stop("Number of rows of matrix X not equal to the number of phenotypes")
  }

  if(intercept[[1]]){
    LT$X = as.matrix(LT$X)
    LT$X_train = LT$X[whichNa,]
  }
  else{
    if(i == intercept[[2]]){
      LT$X = cbind(Intercept = 1,as.matrix(LT$X))
      LT$X_train = LT$X[whichNa,]
    }
    else{ LT$X = as.matrix(LT$X); LT$X_train = LT$X[whichNa,] }
  }

  LT$colNames = colnames(data.frame(LT$X))
  LT$p = ncol(LT$X)

  if(is.null(LT$update)){ LT$update = "scalar" }

  if(LT$update == "blocks"){
    if(is.null(LT$nblocks)){LT$nblocks = 1}
    LT$blocks<-build_blocks(LT$p,LT$nblocks)
  }

  if(response_type%in%c("DLN","gaussian")){
    x2 = apply(LT$X,2L,function(x) sum(x^2))
    if(LT$update == "scalar"){
      LT$x2 = x2
    }
    if(LT$update == "blocks"){
      LT$tX = t(LT$X)
      LT$tXX = crossprod(LT$X)
    }
  }
  if(response_type%in%c("Poisson","PLN")){
    x2 = apply(LT$X,2L,function(x) sum(x^2))
    if(LT$update == "scalar"){
      LT$X2 = apply(LT$X,2L,function(x) x^2)
    }
    if(LT$update == "blocks"){
      LT$tX = t(LT$X)
    }
  }

  if (is.null(LT$R2)) {
    LT$R2 = R2/nLT
  }

  if (is.null(LT$priorvB)){
    LT$priorvB = 5
  }

  MSx<-sum(x2)/n

  if(response_type == "gaussian"){
    if (is.null(LT$priorSB)){
      LT$priorSB = (LT$priorvB+2)*(LT$R2)*var(y,na.rm=TRUE)/(MSx)
      LT$SB = LT$priorSB
    }
  }
  if(response_type == "DLN"){
    if (is.null(LT$priorSB)){
      LT$priorSB = (LT$priorvB+2)*(LT$R2)*var(log(y+1),na.rm=TRUE)/(MSx)
      LT$SB = LT$priorSB
    }
  }
  if(response_type%in%c("Poisson","PLN")){
    if (is.null(LT$priorSB)){
      LT$priorSB = (LT$priorvB+2)*(LT$R2)*mean(log(y+1),na.rm=TRUE)/(MSx/2)
      LT$SB = LT$priorSB
    }
  }

  if (is.null(LT$shape0)) {
    LT$shape0 = 1.1
  }

  if (is.null(LT$rate0)) {
    LT$rate0 = (LT$shape0 - 1)/LT$SB
  }

  LT$Bv = rep(0,LT$p)
  LT$varB = rep(1/rgamma(1,LT$priorvB/2,LT$priorSB/2),LT$p)
  LT$post_Bv = rep(0,LT$p)
  LT$post_Bv2 = rep(0,LT$p)
  LT$post_varB = rep(0,LT$p)
  LT$post_varB2 = rep(0,LT$p)
  LT$post_SB = 0
  LT$post_SB2 = 0

  fname = paste(saveAt,LT$Name,"_SB.dat",sep = "")
  LT$fileOut = file(description = fname,open = "w")
  LT$NamefileOut = fname

  return(LT)

}

setLT.BL<-function(LT,whichNa,y,n,R2,thin,nIter,burnIn,nLT,response_type,
                     saveAt,intercept,i)
{

  if (any(is.na(LT$X))) {
    stop("matrix X has NAs")
  }
  if (nrow(as.matrix(LT$X)) != n) {
    stop("Number of rows of matrix X not equal to the number of phenotypes")
  }

  if(intercept[[1]]){
    LT$X = as.matrix(LT$X)
    LT$X_train = LT$X[whichNa,]
  }
  else{
    if(i == intercept[[2]]){
      LT$X = cbind(Intercept = 1,as.matrix(LT$X))
      LT$X_train = LT$X[whichNa,]
    }
    else{ LT$X = as.matrix(LT$X); LT$X_train = LT$X[whichNa,] }
  }

  LT$colNames = colnames(data.frame(LT$X))
  LT$p = ncol(LT$X)

  if(is.null(LT$update)){ LT$update = "scalar" }

  if(LT$update == "blocks"){
    if(is.null(LT$nblocks)){LT$nblocks = 1}
    LT$blocks<-build_blocks(LT$p,LT$nblocks)
  }

  if(response_type%in%c("DLN","gaussian")){
    x2 = apply(LT$X,2L,function(x) sum(x^2))
    if(LT$update == "scalar"){
      LT$x2 = x2
    }
    if(LT$update == "blocks"){
      LT$tX = t(LT$X)
      LT$tXX = crossprod(LT$X)
    }
  }

  if(response_type%in%c("Poisson","PLN")){
    x2 = apply(LT$X,2L,function(x) sum(x^2))
    if(LT$update == "scalar"){
      LT$X2 = apply(LT$X,2L,function(x) x^2)
    }
    if(LT$update == "blocks"){
      LT$tX = t(LT$X)
    }
  }

  if (is.null(LT$R2)) {
    LT$R2 = R2/nLT
  }

  MSx<-sum(x2)/n

  if (!is.null(LT$lambda)) {
    if (LT$lambda < 0) {
      stop("lambda should be positive")
    }
  }

  if(response_type %in% c("DLN","PLN","gaussian")){
    if (is.null(LT$lambda)) {
      LT$lambda2 = 2 * (1 - R2)/(LT$R2) * MSx
      LT$lambda = sqrt(LT$lambda2)
    }
    else {
      if (LT$lambda < 0)
        stop("lambda should be positive")
      LT$lambda2 = LT$lambda^2
    }
  }else{
    if (is.null(LT$lambda)) {
      LT$lambda2 = MSx / (LT$R2 * mean(log(y+1),na.rm=TRUE))
      LT$lambda = sqrt(LT$lambda2)
    }
    else {
      if (LT$lambda < 0)
        stop("lambda should be positive")
      LT$lambda2 = LT$lambda^2
    }
  }

  if (is.null(LT$type)) {
    LT$type = "gamma"
  }
  else {
    if (!LT$type %in% c("gamma","FIXED"))
      stop("The prior for lambda^2 should be gamma or a point of mass (i.e.,fixed lambda)")
  }


  if (LT$type == "gamma") {
    if (is.null(LT$shape)) {
      LT$shape = 1.1
    }
    if (is.null(LT$rate)) {
      LT$rate = (LT$shape - 1)/LT$lambda2
    }
  }

  LT$Bv = rep(0,LT$p)
  LT$post_Bv = rep(0,LT$p)
  LT$post_Bv2 = rep(0,LT$p)
  if(response_type%in%c("Poisson","PLN")){tmp = ((mean(log(y+1),na.rm=TRUE) * R2/nLT)/(MSx))}
  if(response_type == "DLN"){tmp = ((var(log(y+1),na.rm=TRUE) * R2/nLT)/(MSx))}
  if(response_type == "gaussian"){tmp = ((var(y,na.rm=TRUE) * R2/nLT)/(MSx))}
  LT$tau2 = rep(tmp,LT$p)
  LT$post_tau2 = 0
  LT$post_lambda = 0

  fname = paste(saveAt,LT$Name,"_lambda.dat",sep = "")
  LT$NamefileOut = fname
  LT$fileOut = file(description = fname,open = "w")
  return(LT)
}

setLT.RKHS=function(LT,whichNa,y,n,R2,nLT,response_type,saveAt,i)
{

  # Checking inputs
  if(is.null(LT$V))
  {
    if(is.null(LT$K)) stop("Kernel for linear term ",i," was not provided,specify it with list(K=?,model='RKHS'),where ? is the kernel matrix")

    if(!is.matrix(LT$K)) stop("Kernel for linear term ",i," should be a matrix,the kernel provided is of class ",class(LT$K))

    LT$K = as.matrix(LT$K)

    if(nrow(LT$K)!=ncol(LT$K)) stop("Kernel for linear term ",i," is not a square matrix")

    tmp = eigen(LT$K,symmetric=TRUE)
    LT$V = tmp$vectors
    LT$d = tmp$values
    rm(tmp)

  }

  # Default value for tolD
  # Only those eigenvectors whose eigenvalues> tolD are kept.
  if (is.null(LT$tolD))
  {
    LT$tolD = 1e-10
  }

  # Removing elements whose eigenvalues < tolD
  tmp = LT$d > LT$tolD
  LT$levelsU = sum(tmp)
  LT$d = LT$d[tmp]

  LT$V_train = LT$V[whichNa,tmp]
  LT$V_test = LT$V[!whichNa,tmp]
  LT$V = LT$V[,tmp]

  if(is.null(LT$update)){ LT$update = "scalar" }

  if (LT$update == "blocks") {

    repeat {
      mensaje<-sprintf("Insert the number of blocks (it should be less than the number of eigenvalues (%d)): ",LT$levelsU)
      nblocks<-readline(prompt = mensaje)
      nblocks_num<-as.numeric(nblocks)

      if (!is.na(nblocks_num) && nblocks_num > 0 && nblocks_num < LT$levelsU) break
    }

    LT$nblocks<-nblocks_num
    LT$blocks<-build_blocks(LT$levelsU,LT$nblocks)
  }

  if(response_type%in%c("DLN","gaussian")){
    if(LT$update == "blocks"){
      LT$tV = t(LT$V)
    }
  }

  if(response_type%in%c("Poisson","PLN")){
    v2 = rep(1,LT$levelsU)
    if(LT$update == "scalar"){
      LT$V2 = apply(LT$V,2L,function(x) x^2)
    }
    if(LT$update == "blocks"){
      LT$tV = t(LT$V)
    }
  }

  # Default degrees of freedom and scale parameter associated with the variance component for marker effect
  if (is.null(LT$priorvU))
  {
    LT$priorvU = 5
  }

  if(is.null(LT$R2))
  {
    LT$R2=R2/nLT
  }

  if (is.null(LT$priorSU))
  {
    if(LT$priorvU<=0) stop("priorvU>0 in RKHS in order to set priorSU");

    if(response_type == "gaussian"){
      LT$priorSU = (LT$priorvU+2)*(LT$R2)*var(y,na.rm=TRUE)/(mean(LT$d))
    }
    if(response_type == "DLN"){
      LT$priorSU = (LT$priorvU+2)*(LT$R2)*var(log(y+1),na.rm=TRUE)/(mean(LT$d))
    }
    if(response_type %in% c("Poisson","PLN")){
      LT$priorSU = (LT$priorvU+2)*(LT$R2)*mean(log(y+1),na.rm=TRUE)/((mean(LT$d))/2)
    }

  }

  LT$u=rep(0,n)

  LT$varU=1/rgamma(1,LT$priorvU/2,LT$priorSU/2)

  LT$uStar=rep(0,LT$levelsU)

  fname=paste(saveAt,LT$Name,"_varU.dat",sep="")
  LT$NamefileOut=fname
  LT$fileOut=file(description=fname,open="w")

  LT$post_varU=0
  LT$post_varU2=0
  LT$post_uStar = rep(0,LT$levelsU)
  LT$post_u = rep(0,nrow(LT$V))
  LT$post_u2 = rep(0,nrow(LT$V))

  return(LT)
}

################################# Multi-trait ##################################

setResCov<-function (n,resCov,ntraits,Sy,R2,saveAt)
{
  message("Initializing resCov")

  if (is.null(resCov$type)) {
    resCov$type<-"UN"
    message("Modelling SigmaE as UNstructured")
  }
  else {
    if (!(resCov$type %in% c("UN","DIAG","FA","REC"))) {
      stop("Error '",resCov$type,"' not implemented (note: evaluation is case sensitive)")
    }
  }
  if (resCov$type == "UN") {
    message("Setting hyperparameters for UNstructured SigmaE")
    if (is.null(resCov$priorv)) {
      resCov$priorv<-ntraits + 1
      message("v was set to ",resCov$priorv)
    }
    if (is.null(resCov$priorS)) {
      resCov$priorS<-(1 - R2) * Sy * (resCov$priorv + ntraits + 1)
      message("S was set to ")
      print(resCov$priorS)
    }
    resCov$SigmaE<-MCMCpack::riwish(v = resCov$priorv,S = resCov$priorS)
    resCov$SigmaEInv<-solve(resCov$SigmaE)
  }

  if (resCov$type == "DIAG") {
    message("Setting hyperparameters for DIAG SigmaE")
    if (is.null(resCov$priorv)) {
      resCov$priorv<-rep(ntraits + 1,ntraits)
      message("v set to ",ntraits + 1," for all the traits")
    }
    if (is.null(resCov$priorS)) {
      resCov$priorS<-(1 - R2) * diag(Sy) * (resCov$priorv + 2)
      message("S was set to ")
      print(resCov$priorS)
    }
    resCov$SigmaE<-diag(1/rgamma(ntraits,resCov$priorv/2,resCov$priorS/2))
    resCov$SigmaEInv<-diag(1/diag(resCov$SigmaE))
  }

  if (resCov$type == "FA") {
    message("Setting hyperparameters for FA SigmaE")
    if (is.null(resCov$M))
      stop("M can not be null")
    if (!is.logical(resCov$M))
      stop("M must be logical matrix (with entries being TRUE/FALSE)")
    if (!is.matrix(resCov$M))
      stop("M must be a matrix")
    if (nrow(resCov$M) != ntraits)
      stop("M must have ",ntraits," rows")
    if (ncol(resCov$M) > ntraits)
      stop("Number of columns of M must be smaller than ",ntraits)
    resCov$nF<-ncol(resCov$M)
    if (is.null(resCov$priorv)) {
      resCov$priorv<-rep(ntraits + 1,ntraits)
      message("v set to ",ntraits + 1," for all the traits")
    }
    if (is.null(resCov$priorS)) {
      resCov$priorS<-(1 - R2) * diag(Sy) * (resCov$priorv + 2)
      message("S was set to ")
      print(resCov$priorS)
    }
    if (is.null(resCov$var)) {
      resCov$var<-100
      message("var was set to 100")
    }

    resCov$SigmaE<-MCMCpack::riwish(v = ntraits + 1,S = diag(resCov$priorS))
    sdU<-sqrt(diag(resCov$SigmaE))
    if (is.null(resCov$W)) {
      FA<-factanal(covmat = resCov$SigmaE,factors = resCov$nF,
                     nstart = 10)
      resCov$W<-matrix(nrow = ntraits,ncol = resCov$nF,
                         0)
      resCov$W[resCov$M]<-(diag(sdU) %*% FA$loadings)[resCov$M]
      resCov$PSI<-(sdU^2) * FA$uniquenesses + 1e-04
    }
    else {
      resCov$PSI<-(sdU^2)/2
    }
    resCov$SigmaE<-tcrossprod(resCov$W) + diag(resCov$PSI)
    resCov$SigmaEInv<-solve(resCov$SigmaE)
    resCov$F<-matrix(nrow = n,ncol = resCov$nF,0)
    resCov$post_W<-matrix(0,nrow = ntraits,ncol = resCov$nF)
    resCov$post_W2<-matrix(0,nrow = ntraits,ncol = resCov$nF)
    resCov$post_PSI<-rep(0,ntraits)
    resCov$post_PSI2<-rep(0,ntraits)

    resCov$fName_W<-paste(saveAt,"W_SigmaE.dat",sep = "")
    resCov$fName_PSI<-paste(saveAt,"PSI_SigmaE.dat",sep = "")
    resCov$f_W<-file(description = resCov$fName_W,open = "w")
    resCov$f_PSI<-file(description = resCov$fName_PSI,open = "w")
  }

  if (resCov$type == "REC") {
    message("Setting hyperparameters for REC SigmaE")
    if (is.null(resCov$M))
      stop("M can not be null")
    if (!is.logical(resCov$M))
      stop("M must be logical matrix (with entries being TRUE/FALSE)")
    if (!is.matrix(resCov$M))
      stop("M must be a matrix")
    if (nrow(resCov$M) != ncol(resCov$M))
      stop("M must be a square matrix")
    if (nrow(resCov$M) != ntraits)
      stop("M must have ",ntraits," rows and columns")
    if (any(diag(resCov$M) == TRUE))
      stop("All diagonal entries of M must be set to FALSE")
    resCov$M[upper.tri(resCov$M)]<-FALSE
    if (is.null(resCov$priorv)) {
      resCov$priorv<-rep(ntraits + 1,ntraits)
      message("v set to ",ntraits + 1," for all the traits")
    }
    if (is.null(resCov$priorS)) {
      resCov$priorS<-(1 - R2) * diag(Sy) * (resCov$priorv +  2)
      message("S was set to ")
      print(resCov$priorS)
    }
    if (is.null(resCov$var)) {
      resCov$var<-100
      message("var was set to 100")
    }
    resCov$W<-matrix(0,nrow = ntraits,ncol = ntraits)
    resCov$PSI<-rep(NA,ntraits)
    for (k in 1:ntraits) {
      resCov$PSI[k]<-1/rgamma(n=1,resCov$priorv[k]/2,resCov$priorS[k]/2)
    }
    resCov$SigmaE<-diag(ntraits) * diag(resCov$PSI) * diag(ntraits)
    resCov$SigmaEInv<-diag(1/diag(resCov$SigmaE))
    resCov$post_W<-matrix(0,nrow = ntraits,ncol = ntraits)
    resCov$post_W2<-matrix(0,nrow = ntraits,ncol = ntraits)
    resCov$post_PSI<-rep(0,ntraits)
    resCov$post_PSI2<-rep(0,ntraits)

    resCov$fName_W<-paste(saveAt,"W_SigmaE.dat",sep = "")
    resCov$fName_PSI<-paste(saveAt,"PSI_SigmaE.dat",sep = "")
    resCov$f_W<-file(description = resCov$fName_W,open = "w")
    resCov$f_PSI<-file(description = resCov$fName_PSI,open = "w")
  }

  resCov$post_SigmaE<-matrix(0,nrow = ntraits,ncol = ntraits)
  resCov$post_SigmaE2<-matrix(0,nrow = ntraits,ncol = ntraits)

  resCov$fName_SigmaE<-paste(saveAt,"SigmaE.dat",sep = "")
  resCov$f_SigmaE<-file(description = resCov$fName_SigmaE,open = "w")
  message("Done")
  return(resCov)
}

setCov.UN<-function (Cov,ntraits,i,mo,saveAt)
{
  message("UNstructured covariance matrix")
  if (is.null(Cov$priorvB)) {
    Cov$priorvB<-ntraits + 1
    message("vB was set to ",Cov$priorvB)
  }
  if (is.null(Cov$priorSB)) {
    Cov$priorSB<-mo * (Cov$priorvB + ntraits + 1)
    message("SB set to ")
    print(Cov$priorSB)
  }
  Cov$SigmaB<-MCMCpack::riwish(v = Cov$priorvB,S = Cov$priorSB)
  Cov$SigmaBInv<-solve(Cov$SigmaB)
  Cov$post_SigmaB<-matrix(0,nrow = ntraits,ncol = ntraits)
  Cov$post_SigmaB2<-matrix(0,nrow = ntraits,ncol = ntraits)

  Cov$fName_SigmaB<-paste(saveAt,"SigmaB_",i,".dat",sep = "")
  Cov$f_SigmaB<-file(description = Cov$fName_SigmaB,open = "w")
  return(Cov)
}

setCov.DIAG<-function (Cov,ntraits,i,mo,saveAt)
{
  message("DIAGonal covariance matrix")
  if (is.null(Cov$priorvB)) {
    Cov$priorvB<-rep(ntraits + 1,ntraits)
    message("vB set to  ",ntraits + 1," for all the traits")
  }
  if (is.null(Cov$priorSB)) {
    Cov$priorSB<-mo * (Cov$priorvB + 2)
    message("SB was set to: ")
    print(Cov$priorSB)
  }
  Cov$SigmaB<-diag(1/rgamma(ntraits,Cov$priorvB/2,Cov$priorSB/2))
  Cov$SigmaBInv<-diag(1/diag(Cov$SigmaB))
  Cov$post_SigmaB<-matrix(0,nrow = ntraits,ncol = ntraits)
  Cov$post_SigmaB2<-matrix(0,nrow = ntraits,ncol = ntraits)

  Cov$fName_SigmaB<-paste(saveAt,"SigmaB_",i,".dat",sep = "")
  Cov$f_SigmaB<-file(description = Cov$fName_SigmaB,open = "w")
  return(Cov)
}

setCov.FA<-function (Cov,ntraits,nD,i,mo,saveAt)
{
  message("FA covariance matrix")
  if (is.null(Cov$priorvB)) {
    Cov$priorvB<-rep(ntraits + 1,ntraits)
    message("vB set to ",ntraits + 1," for all the traits")
  }
  if (is.null(Cov$priorSB)) {
    Cov$priorSB<-mo * (Cov$priorvB + 2)
    message("SB was set to: ")
    print(Cov$priorSB)
  }
  if (is.null(Cov$var)) {
    Cov$var<-100
    message("var was set to 100")
  }
  if (is.null(Cov$varimax)) {
    Cov$varimax = TRUE
    message("Rotation set to varimax")
  }
  if (is.null(Cov$M))
    stop("M can not be null")
  if (!is.logical(Cov$M))
    stop("M must be logical matrix (with entries being TRUE/FALSE)")
  if (!is.matrix(Cov$M))
    stop("M must be a matrix")
  if (nrow(Cov$M) != ntraits)
    stop("M must have ",ntraits," rows")
  if (ncol(Cov$M) > ntraits)
    stop("Number of columns of M must be smaller than ",
         ntraits)
  Cov$nF<-ncol(Cov$M)
  Cov$nD<-nD
  Cov$SigmaB<-MCMCpack::riwish(v = ntraits + 1,S = diag(Cov$priorSB))
  sdU<-sqrt(diag(Cov$SigmaB))
  if (is.null(Cov$W)) {
    FA<-factanal(covmat = Cov$SigmaB,factors = Cov$nF,
                   nstart = 10)
    Cov$W<-matrix(nrow = ntraits,ncol = Cov$nF,0)
    Cov$W[Cov$M]<-(diag(sdU) %*% FA$loadings)[Cov$M]
    Cov$PSI<-(sdU^2) * FA$uniquenesses + 1e-04
  }
  else {
    Cov$PSI<-(sdU^2)/2
  }
  Cov$SigmaB<-tcrossprod(Cov$W) + diag(Cov$PSI)
  Cov$SigmaBInv<-solve(Cov$SigmaB)
  Cov$F<-matrix(nrow = nD,ncol = Cov$nF,0)
  Cov$post_SigmaB<-matrix(0,nrow = ntraits,ncol = ntraits)
  Cov$post_SigmaB2<-matrix(0,nrow = ntraits,ncol = ntraits)
  Cov$post_W<-matrix(0,nrow = ntraits,ncol = Cov$nF)
  Cov$post_W2<-matrix(0,nrow = ntraits,ncol = Cov$nF)
  Cov$post_PSI<-rep(0,ntraits)
  Cov$post_PSI2<-rep(0,ntraits)

  Cov$fName_W<-paste(saveAt,"W_",i,".dat",sep = "")
  Cov$fName_PSI<-paste(saveAt,"PSI_",i,".dat",sep = "")
  Cov$f_W<-file(description = Cov$fName_W,open = "w")
  Cov$f_PSI<-file(description = Cov$fName_PSI,open = "w")
  return(Cov)
}

setCov.REC<-function (Cov,ntraits,i,mo,saveAt)
{
  message("RECursive covariance matrix")
  if (is.null(Cov$priorvB)) {
    Cov$priorvB<-rep(ntraits + 1,ntraits)
    message("vB set to ",ntraits + 1," for all the traits")
  }
  if (is.null(Cov$priorSB)) {
    Cov$priorSB<-mo * (Cov$priorvB + 2)
    message("SB was set to: ")
    print(Cov$priorSB)
  }
  if (is.null(Cov$var)) {
    Cov$var<-100
    message("var was set to 100")
  }
  if (is.null(Cov$M))
    stop("M can not be null")
  if (!is.logical(Cov$M))
    stop("M must be logical matrix (with entries being TRUE/FALSE)")
  if (!is.matrix(Cov$M))
    stop("M must be a matrix")
  if (nrow(Cov$M) != ncol(Cov$M))
    stop("M must be a square matrix")
  if (nrow(Cov$M) != ntraits)
    stop("M must have ",ntraits," rows and columns")
  if (any(diag(Cov$M) == TRUE))
    stop("All diagonal entries of M must be set to FALSE")
  Cov$M[upper.tri(Cov$M)]<-FALSE
  Cov$W<-matrix(0,nrow = ntraits,ncol = ntraits)
  Cov$PSI<-rep(NA,ntraits)
  for (k in 1:ntraits) {
    Cov$PSI[k]<-1/rgamma(n=1,Cov$priorvB[k]/2,Cov$priorSB[k]/2)
  }

  Cov$SigmaB<-MCMCpack::riwish(v = ntraits + 1,S = diag(Cov$priorSB))
  Cov$SigmaBInv<-solve(Cov$SigmaB)
  Cov$post_SigmaB<-matrix(0,nrow = ntraits,ncol = ntraits)
  Cov$post_SigmaB2<-matrix(0,nrow = ntraits,ncol = ntraits)
  Cov$post_W<-matrix(0,nrow = ntraits,ncol = ntraits)
  Cov$post_W2<-matrix(0,nrow = ntraits,ncol = ntraits)
  Cov$post_PSI<-rep(0,ntraits)
  Cov$post_PSI2<-rep(0,ntraits)

  Cov$fName_W<-paste(saveAt,"W_",i,".dat",sep = "")
  Cov$fName_PSI<-paste(saveAt,"PSI_",i,".dat",sep = "")
  Cov$f_W<-file(description = Cov$fName_W,open = "w")
  Cov$f_PSI<-file(description = Cov$fName_PSI,open = "w")
  return(Cov)
}

sample_G0_FA<-function (U,F,M,B,PSI,ntraits,nF,nD,df0 = rep(1,ntraits),
                          S0 = rep(1/100,ntraits),priorVar = 100,varimaxRotate = TRUE)
{
  for (i in 1:nF) {
    tmpY<-U - F[,-i] %*% matrix((B[,-i]),ncol = ntraits)
    rhs<-tmpY %*% matrix(B[,i]/PSI,ncol = 1)
    CInv<-1/(sum((B[,i]^2)/PSI) + 1)
    sol<-CInv * rhs
    SD<-sqrt(CInv)
    F[,i]<-rnorm(n = nD,sd = SD,mean = sol)
  }
  for (i in 1:ntraits) {
    for (j in 1:nF) {
      if (M[i,j]) {
        tmpY<-U[,i] - F[,-j] %*% matrix(B[i,-j],ncol = 1)
        CInv<-1/as.numeric(crossprod(F[,j])/PSI[i] + 1/priorVar)
        rhs<-as.numeric(crossprod(F[,j],tmpY)/PSI[i])
        sol<-CInv * rhs
        SD<-sqrt(CInv)
        B[i,j]<-rnorm(n = 1,mean = sol,sd = SD)
      }
    }
    D<-U[,i] - F %*% B[i,]
    df<-df0[i] + nD
    SS<-S0[i] + crossprod(D)
    PSI[i]<-1/rgamma(n=1,df/2,SS/2)
  }
  if ((nF > 1) && varimaxRotate) {
    B<-varimax(B)$loadings[]
  }
  G<-tcrossprod(B) + diag(PSI)
  out<-list(F = F,PSI = PSI,B = B,G = G)
  return(out)
}

sample_G0_REC<-function (U,M,PSI,ntraits,priorVar = 100,df0 = rep(0,ntraits),
                           S0 = rep(0,ntraits))
{
  B<-matrix(nrow = ntraits,ncol = ntraits,0)
  for (i in 1:ntraits) {
    dimX<-sum(M[i,])
    if (dimX > 0) {
      tmpX<-U[,M[i,]]
      tmpY<-U[,i]
      C<-crossprod(tmpX)/PSI[i] + 1/priorVar
      CInv<-chol2inv(chol(C))
      rhs<-crossprod(tmpX,tmpY)/PSI[i]
      sol<-crossprod(CInv,rhs)
      L<-chol(CInv)
      shock<-crossprod(L,rnorm(dimX))
      tmpB<-as.numeric(sol + shock)
      B[i,M[i,]]<-tmpB
      uStar<-tmpY - matrix(tmpX,ncol = dimX) %*% (tmpB)
      SS<-as.numeric(crossprod(uStar)) + S0[i]
      df<-nrow(U) + df0[i]
      PSI[i]<-1/rgamma(n=1,df/2,SS/2)
    }
    else {
      SS<-as.numeric(crossprod(U[,i])) + S0[i]
      df<-nrow(U) + df0
      PSI[i]<-1/rgamma(n=1,df/2,SS/2)
    }
  }
  tmp<-solve(diag(ntraits) - B)
  G<-tmp %*% diag(PSI) %*% t(tmp)
  out<-list(B = B,PSI = PSI,G = G)
  return(out)
}

setLT.Fixed_mt<-function(LT,n,ntraits,i,saveAt,response_type,NoWhichNa)
{

  message("Setting linear term ",i)
  if (is.null(LT$common)) {
    LT$common<-TRUE
    message("matrix for fixed effects X is the same for all the traits,")
    message("so the same effects are assumed for all the traits")
  }
  else {
    message("matrix for fixed effects X is the same for all the traits,")
    message("so the same effects are assumed for all the traits")
  }
  if (is.null(LT$X))
    stop("X can not be NULL\n")
  if (!is.matrix(LT$X))
    stop("X must be a matrix\n")
  if (any(is.na(LT$X)))
    stop("X has NAs\n")

  LT$X_train = LT$X[NoWhichNa,]
  LT$colNames = colnames(LT$X)

  LT$Cov<-list()
  LT$Cov$SigmaB<-diag(rep(1e+10,ntraits))
  LT$Cov$SigmaBInv<-diag(1/diag(LT$Cov$SigmaB))

  if (qr(LT$X)$rank < ncol(LT$X))
    stop("X is rank deficient")

  if(response_type%in%c("Poisson","PLN")){
    LT$X2 = apply(LT$X,2L,function(x) x^2)
  }
  if(response_type%in%c("DLN","gaussian")){
    LT$x2 = apply(LT$X,2L,function(x) sum(x^2))
  }

  LT$p<-ncol(LT$X)
  LT$Bv<-matrix(0,nrow = LT$p,ncol = ntraits)
  LT$post_Bv<-matrix(0,nrow = LT$p,ncol = ntraits)
  LT$post_Bv2<-matrix(0,nrow = LT$p,ncol = ntraits)

  return(LT)

}

setLT.BRR_mt<-function(LT,n,ntraits,i,Sy,nLT,R2,saveAt,response_type,NoWhichNa)
{

  message("Setting linear term ",i)
  if (is.null(LT$X))
    stop("X can not be NULL\n")
  if (!is.matrix(LT$X))
    stop("X must be a matrix\n")
  if (any(is.na(LT$X)))
    stop("X has NAs\n")

  LT$X_train = LT$X[NoWhichNa,]
  LT$colNames = colnames(LT$X)

  if(response_type%in%c("DLN","gaussian")){
    LT$x2 = apply(LT$X,2L,function(x) sum(x^2))
    MSx<-sum(LT$x2)/n
  }
  if(response_type%in%c("Poisson","PLN")){
    LT$X2 = apply(LT$X,2L,function(x) x^2)
    LT$x2 = apply(LT$X,2L,function(x) sum(x^2))
    MSx<-sum(LT$x2)/(2*n)
  }

  LT$p<-ncol(LT$X)
  LT$Bv<-matrix(0,nrow = LT$p,ncol = ntraits)
  if (is.null(LT$Cov)) {
    LT$Cov<-list()
    LT$Cov$type<-"UN"
  }
  else {
    if (is.null(LT$Cov$type)) {
      LT$Cov$type<-"UN"
    }
    else {
      if (!(LT$Cov$type %in% c("UN","DIAG","FA","REC"))) {
        stop("Error '",LT$Cov$type,"' not implemented (note: evaluation is case sensitive)")
      }
    }
  }
  LT$Cov<-switch(LT$Cov$type,UN = setCov.UN(Cov = LT$Cov,
                   ntraits = ntraits,i = i,mo = (R2/nLT) * Sy/MSx,saveAt = saveAt),
                   DIAG = setCov.DIAG(Cov = LT$Cov,ntraits = ntraits,i = i,
                   mo = (R2/nLT) * diag(Sy)/MSx,saveAt = saveAt),FA = setCov.FA(Cov = LT$Cov,
                   ntraits = ntraits,nD = LT$p,i = i,mo = (R2/nLT) * diag(Sy)/MSx,saveAt = saveAt),
                   REC = setCov.REC(Cov = LT$Cov,ntraits = ntraits,i = i,mo = (R2/nLT) * diag(Sy)/MSx,
                   saveAt = saveAt))

  LT$post_Bv<-matrix(0,nrow = LT$p,ncol = ntraits)
  LT$post_Bv2<-matrix(0,nrow = LT$p,ncol = ntraits)

  return(LT)

}

setLT.RKHS_mt<-function(LT,n,ntraits,i,Sy,nLT,R2,saveAt,response_type,NoWhichNa)
{

  if(is.null(LT$EVD) && is.null(LT$K))
  {
    text<-"Either variance co-variance matrix K or its eigen-value decomposition\n"
    text<-paste(text,"must be provided for linear term ",i,"\n")
    text<-paste(text,"To specify the variance covariance matrix K use:\n")
    text<-paste(text,"list(K=?,model='RKHS'),where ? is the user defined (between subjects) co-variance matrix\n")
    text<-paste(text,"To specify the eigen-value decomposition for K use:\n")
    text<-paste(text,"list(EVD=?,model='RKHS'),where ? is the output from eigen function for a user defined (between subjects) co-variance matrix\n")
    stop(text)
  }

  if((!is.null(LT$K)) && (!is.null(LT$EVD)))
  {
    message("Variance covariance matrix K and its eigen-value decomposition for linear term ",i," was provided")
    message("ONLY EVD will be used")
    LT$K<-NULL
  }

  if((!is.null(LT$K)) && is.null(LT$EVD))
  {
    message("Checking variance co-variance matrix K for linear term ",i)
    if(nrow(LT$K)!=ncol(LT$K)) stop("variance covariance matrix must be square")
    LT$EVD<-eigen(LT$K,symmetric=TRUE)
    message("Ok")
  }

  if(is.null(LT$K) && (!is.null(LT$EVD)))
  {
    message("Checking EVD provided for linear term ",i)
    if(!is.matrix(LT$EVD$vectors)) stop("eigen-vectors must be a matrix\n")
    if(nrow(LT$EVD$vectors)!=ncol(LT$EVD$vectors)) stop("eigen-vectors must be a square matrix\n")
    if(!is.numeric(LT$EVD$values)) stop("eigen-values must be a numeric vector\n")
    message("Ok")
  }

  keep<-LT$EVD$values>1e-10
  LT$EVD$vectors<-LT$EVD$vectors[,keep]
  LT$EVD$values<-LT$EVD$values[keep]

  # X=Gamma*Lambda^{1/2}
  LT$X<-sweep(x=LT$EVD$vectors,MARGIN=2,STATS=sqrt(LT$EVD$values),FUN="*")

  LT<-setLT.BRR_mt(LT=LT,n=n,ntraits=ntraits,i=i,Sy=Sy,nLT=nLT,R2=R2,saveAt=saveAt,
                    response_type=response_type,NoWhichNa=NoWhichNa)

  return(LT)

}

######################## Multi-trait Multi-Environment #########################

sample_G0_FA_MTME<-function (U1,U2,F1,F2,M,B,B2,const1,const2,PSI,
                               ntraits,nF,nD1,nD2,df0 = rep(1,ntraits),
                               S0 = rep(1/100,ntraits),priorVar = 100,
                               varimaxRotate = TRUE)
{
  for (i in 1:nF) {
    tmpY1<-U1 - F1[,-i] %*% matrix((B[,-i]),ncol = ntraits)
    rhs1<-tmpY1 %*% matrix(B[,i]/PSI,ncol = 1)
    CInv<-1/(sum((B[,i]^2)/PSI) + const1)
    sol1<-CInv * rhs1
    SD1<-sqrt(CInv)
    F1[,i]<-rnorm(n = nD1,sd = SD1,mean = sol1)

    tmpY2<-U2 - F2[,-i] %*% matrix((B[,-i]),ncol = ntraits)
    rhs2<-tmpY2 %*% matrix(B[,i]/PSI,ncol = 1)
    CInv2<-1/(sum((B2[,i]^2)/PSI) + const2)
    sol2<-CInv2 * rhs2
    SD2<-sqrt(CInv2)
    F2[,i]<-rnorm(n = nD2,sd = SD2,mean = sol2)
  }
  for (i in 1:ntraits) {
    for (j in 1:nF) {
      if (M[i,j]) {
        tmpY1<-U1[,i] - F1[,-j] %*% matrix(B[i,-j],ncol = 1)
        tmpY2<-U2[,i] - F2[,-j] %*% matrix(B[i,-j],ncol = 1)
        CInv<-1/as.numeric((crossprod(F1[,j])+crossprod(F2[,j]))/PSI[i] + 1/priorVar)
        rhs<-as.numeric((crossprod(F1[,j],tmpY1)+crossprod(F2[,j],tmpY2))/PSI[i])
        sol<-CInv * rhs
        SD<-sqrt(CInv)
        B[i,j]<-rnorm(n = 1,mean = sol,sd = SD)
      }
    }
    D1<-U1[,i] - F1 %*% B[i,]; D2<-U2[,i] - F2 %*% B[i,]
    df<-df0[i] + nD1 + nD2
    SS<-S0[i] + crossprod(D1) + crossprod(D2)
    PSI[i]<-1/rgamma(n=1,df/2,SS/2)
  }
  if ((nF > 1) && varimaxRotate) {
    B<-varimax(B)$loadings[]
  }
  G<-tcrossprod(B) + diag(PSI)
  out<-list(F1 = F1,F2 = F2,PSI = PSI,B = B,G = G)
  return(out)
}

sample_G0_REC_MTME<-function (U1,U2,M,PSI,ntraits,priorVar = 100,
                                df0 = rep(0,ntraits),S0 = rep(0,ntraits))
{
  B<-matrix(nrow = ntraits,ncol = ntraits,0)
  for (i in 1:ntraits) {
    dimX<-sum(M[i,])
    if (dimX > 0) {
      tmpX1<-U1[,M[i,]]; tmpX2<-U2[,M[i,]]
      tmpY1<-U1[,i]; tmpY2<-U2[,i]
      C<-(crossprod(tmpX1)+crossprod(tmpX2))/PSI[i] + 1/priorVar
      CInv<-chol2inv(chol(C))
      rhs<-(crossprod(tmpX1,tmpY1)+crossprod(tmpX2,tmpY2))/PSI[i]
      sol<-crossprod(CInv,rhs)
      L<-chol(CInv)
      shock<-crossprod(L,rnorm(dimX))
      tmpB<-as.numeric(sol + shock)
      B[i,M[i,]]<-tmpB
      uStar1<-tmpY1 - matrix(tmpX1,ncol = dimX) %*% (tmpB)
      uStar2<-tmpY2 - matrix(tmpX2,ncol = dimX) %*% (tmpB)
      SS<-as.numeric(crossprod(uStar1)) + as.numeric(crossprod(uStar2)) + S0[i]
      df<-nrow(U1) + nrow(U2) + df0[i]
      PSI[i]<-1/rgamma(n=1,df/2,SS/2)
    }
    else {
      SS<-as.numeric(crossprod(U1[,i])) + as.numeric(crossprod(U2[,i])) + S0[i]
      df<-nrow(U1) + nrow(U2) + df0
      PSI[i]<-1/rgamma(n=1,df/2,SS/2)
    }
  }
  tmp<-solve(diag(ntraits) - B)
  G<-tmp %*% diag(PSI) %*% t(tmp)
  out<-list(B = B,PSI = PSI,G = G)
  return(out)
}

setLT.MTME<-function(LT,n,ntraits,i,Sy,nLT,R2,saveAt,response_type,NoWhichNa)
{

  if (is.null(LT$X))
    stop("X can not be NULL\n")
  if (!is.matrix(LT$X))
    stop("X must be a matrix\n")
  if (any(is.na(LT$X)))
    stop("X has NAs\n")

  LT$X_train = LT$X[NoWhichNa,]
  LT$colNames = colnames(LT$X)

  if(response_type%in%c("DLN","gaussian")){
    LT$x2 = apply(LT$X,2L,function(x) sum(x^2))
  }
  if(response_type%in%c("Poisson","PLN")){
    LT$X2 = apply(LT$X,2L,function(x) x^2)
    LT$x2 = apply(LT$X,2L,function(x) sum(x^2))
  }

  LT$p<-ncol(LT$X)
  LT$Bv<-matrix(0,nrow = LT$p,ncol = ntraits)
  LT$post_Bv<-matrix(0,nrow = LT$p,ncol = ntraits)
  LT$post_Bv2<-matrix(0,nrow = LT$p,ncol = ntraits)

  return(LT)

}

setCov.MTME<-function(LT,n,nEnvs,i,Sy,nLT,R2,saveAt,response_type,NoWhichNa)
{

  if(response_type%in%c("DLN","gaussian")){
    x2 = apply(LT$X,2L,function(x) sum(x^2))
    MSx<-sum(x2)/n
  }
  if(response_type%in%c("Poisson","PLN")){
    x2 = apply(LT$X,2L,function(x) sum(x^2))
    MSx<-sum(x2)/(2*n)
  }

  if (is.null(LT$Cov)) {
    LT$Cov<-list()
    LT$Cov$type<-"UN"
  }
  else {
    if (is.null(LT$Cov$type)) {
      LT$Cov$type<-"UN"
    }
    else {
      if (!(LT$Cov$type %in% c("UN","DIAG","FA","REC"))) {
        stop("Error '",LT$Cov$type,"' not implemented (note: evaluation is case sensitive)")
      }
    }
  }
  LT$Cov<-switch(LT$Cov$type,UN = setCov.UN(Cov = LT$Cov,
                   ntraits = nEnvs,i = i,mo = (R2/nLT) * Sy/MSx,saveAt = saveAt),
                   DIAG = setCov.DIAG(Cov = LT$Cov,ntraits = nEnvs,i = i,
                   mo = (R2/nLT) * diag(Sy)/MSx,saveAt = saveAt),FA = setCov.FA(Cov = LT$Cov,
                   ntraits = nEnvs,nD = ncol(LT$X),i = i,mo = (R2/nLT) * diag(Sy)/MSx,saveAt = saveAt),
                   REC = setCov.REC(Cov = LT$Cov,ntraits = nEnvs,i = i,mo = (R2/nLT) * diag(Sy)/MSx,
                   saveAt = saveAt))

  tmp = which(names(LT) %in% c("X"))
  LT = LT[-tmp]
  rm(tmp)

  if(LT$Cov$type%in%c("DIAG","UN")){
    LT$Cov<-listr::list_rename(LT$Cov,priorvE = "priorvB",priorSE = "priorSB",SigmaEnv = "SigmaB",
                        SigmaEnvInv = "SigmaBInv",post_SigmaEnv="post_SigmaB",
                        post_SigmaEnv2="post_SigmaB2",fName_SigmaEnv="fName_SigmaB",
                        f_SigmaEnv="f_SigmaB")
  }

  if(LT$Cov$type%in%c("FA","REC")){
    LT$Cov<-listr::list_rename(LT$Cov,priorvE = "priorvB",priorSE = "priorSB",SigmaEnv = "SigmaB",
                        SigmaEnvInv = "SigmaBInv",post_SigmaEnv="post_SigmaB",
                        post_SigmaEnv2="post_SigmaB2")
  }

  return(LT)
}

setCov.MTME_GT<-function(LT,n,ntraits,i,Sy,nLT,R2,saveAt,response_type,NoWhichNa)
{

  message("Setting linear term ",i)
  if(response_type%in%c("DLN","gaussian")){
    x2 = apply(LT$X,2L,function(x) sum(x^2))
    MSx<-sum(x2)/n
  }
  if(response_type%in%c("Poisson","PLN")){
    x2 = apply(LT$X,2L,function(x) sum(x^2))
    MSx<-sum(x2)/(2*n)
  }

  if (is.null(LT$Cov)) {
    LT$Cov<-list()
    LT$Cov$type<-"UN"
  }
  else {
    if (is.null(LT$Cov$type)) {
      LT$Cov$type<-"UN"
    }
    else {
      if (!(LT$Cov$type %in% c("UN","DIAG","FA","REC"))) {
        stop("Error '",LT$Cov$type,"' not implemented (note: evaluation is case sensitive)")
      }
    }
  }
  LT$Cov<-switch(LT$Cov$type,UN = setCov.UN(Cov = LT$Cov,
                   ntraits = ntraits,i = i,mo = (R2/nLT) * Sy/MSx,saveAt = saveAt),
                   DIAG = setCov.DIAG(Cov = LT$Cov,ntraits = ntraits,i = i,
                   mo = (R2/nLT) * diag(Sy)/MSx,saveAt = saveAt),FA = setCov.FA(Cov = LT$Cov,
                   ntraits = ntraits,nD = ncol(LT$X),i = i,mo = (R2/nLT) * diag(Sy)/MSx,saveAt = saveAt),
                   REC = setCov.REC(Cov = LT$Cov,ntraits = ntraits,i = i,mo = (R2/nLT) * diag(Sy)/MSx,
                   saveAt = saveAt))

  tmp = which(names(LT) %in% c("X"))
  LT = LT[-tmp]
  rm(tmp)

  return(LT)
}

setLT.RKHS_mtme<-function(LT,n,ntraits,i,y,Sy,nLT,R2,saveAt,response_type,NoWhichNa)
{

  if(is.null(LT$EVD) && is.null(LT$K))
  {
    text<-"Either variance co-variance matrix K or its eigen-value decomposition\n"
    text<-paste(text,"must be provided for linear term ",i,"\n")
    text<-paste(text,"To specify the variance covariance matrix K use:\n")
    text<-paste(text,"list(K=?,model='RKHS'),where ? is the user defined (between subjects) co-variance matrix\n")
    text<-paste(text,"To specify the eigen-value decomposition for K use:\n")
    text<-paste(text,"list(EVD=?,model='RKHS'),where ? is the output from eigen function for a user defined (between subjects) co-variance matrix\n")
    stop(text)
  }

  if((!is.null(LT$K)) && (!is.null(LT$EVD)))
  {
    message("Variance covariance matrix K and its eigen-value decomposition for linear term ",i," was provided")
    message("ONLY EVD will be used")
    LT$K<-NULL
  }

  if((!is.null(LT$K)) && is.null(LT$EVD))
  {
    message("Checking variance co-variance matrix K for linear term ",i)
    if(nrow(LT$K)!=ncol(LT$K)) stop("variance covariance matrix must be square")
    LT$EVD<-eigen(LT$K,symmetric=TRUE)
    message("Ok")
  }

  if(is.null(LT$K) && (!is.null(LT$EVD)))
  {
    message("Checking EVD provided for linear term ",i)
    if(!is.matrix(LT$EVD$vectors)) stop("eigen-vectors must be a matrix\n")
    if(nrow(LT$EVD$vectors)!=ncol(LT$EVD$vectors)) stop("eigen-vectors must be a square matrix\n")
    if(!is.numeric(LT$EVD$values)) stop("eigen-values must be a numeric vector\n")
    message("Ok")
  }

  # LT$K_diag= diag(LT$nLines)
  keep<-LT$EVD$values>1e-10
  LT$EVD$vectors<-LT$EVD$vectors[,keep]
  LT$EVD$values<-LT$EVD$values[keep]

  if(is.null(LT$GxT) && is.null(LT$GxExT)){
    stop("You need to introduce the incidence matrix for some interaction: Genotype x Traits (GxT)
    or Genotype x Environment x Traits (GxExT)")
  }

  if(!is.null(LT$GxT) && is.null(LT$GxExT)){
    # Number of lines and Environments
    LT$nLines = ncol(LT$GxT$Z); LT$nEnv = n/LT$nLines; LT$nLinesTraits = LT$nLines*ntraits
    LT$nLinesEnvs = LT$nLines*LT$nEnv

    LT$SE1 = diag(0,ntraits)
    LT$U2 = matrix(nrow=LT$nLines*LT$nEnv,ncol=ntraits,0)

    Aux<-sweep(x=LT$EVD$vectors,MARGIN=2,STATS=sqrt(LT$EVD$values),FUN="*")
    LT$X<-LT$GxT$Z%*%Aux

    LT<-setLT.BRR_mt(LT=LT,n=n,ntraits=ntraits,i=i,Sy=Sy,nLT=nLT,R2=R2,saveAt=saveAt,
                     response_type=response_type,NoWhichNa=NoWhichNa)

    if(LT$Cov$type=="FA"){
      LT$F2<-matrix(nrow = LT$nLinesEnvs,ncol = LT$Cov$nF,0)
      LT$W2<-matrix(nrow = ntraits,ncol = ncol(LT$Cov$M),0)
    }
  }

  if(is.null(LT$GxT) && !is.null(LT$GxExT)){
    # Number of lines and Environments
    LT$nLines = 0; LT$nEnv = ncol(LT$GxExT$Z)/LT$GxExT$nLines; LT$nLinesTraits = LT$GxExT$nLines*ntraits
    LT$nLinesEnvs = LT$GxExT$nLines*LT$nEnv

    LT$SE = diag(0,ntraits)
    LT$U1 = matrix(nrow=LT$GxExT$nLines,ncol=ntraits,0)

    LT$X<-sweep(x=LT$EVD$vectors,MARGIN=2,STATS=sqrt(LT$EVD$values),FUN="*")

    LT<-setCov.MTME_GT(LT=LT,n=n,ntraits=ntraits,i=i,Sy=Sy,nLT=nLT,R2=R2,
                       saveAt=saveAt,response_type=response_type,NoWhichNa=NoWhichNa)

    # For environment sample variances
    y_Env = matrix(matrixcalc::vec(t(y)),ncol=LT$nEnv,byrow=F)

    if(response_type == "gaussian"){
      Sy_Env<-var(y_Env,na.rm = T)
      if(any(is.na(Sy_Env))){Sy_Env=diag(1,LT$nEnv)}
    }
    if(response_type == "DLN"){
      Sy_Env<-var(log(y_Env+1),na.rm = T)
      if(any(is.na(Sy_Env))){Sy_Env=diag(1,LT$nEnv)}
    }
    if(response_type%in%c("Poisson","PLN")){
      aux<-colMeans(log(y_Env+1),na.rm=T)
      Sy_Env<-outer(sqrt(aux),sqrt(aux))
      if(any(is.na(Sy_Env))){Sy_Env=diag(1,LT$nEnv)}
      if(!matrixcalc::is.positive.definite(Sy_Env)){
        eigens = eigen(Sy_Env)$values
        MinEigen = min(eigens)
        alpha = -MinEigen + 1e-6
        diag(Sy_Env)<-diag(Sy_Env) + alpha
      }
    }

    # Eigen decomposition of Traits co-variance
    EVD_Traits = eigen(LT$Cov$SigmaB,symmetric=TRUE)

    keep<-EVD_Traits$values>1e-10
    EVD_Traits_vectors<-EVD_Traits$vectors[,keep]
    EVD_Traits_values<-EVD_Traits$values[keep]

    # Interaction ExGxT components
    LT$GxExT$eigenvals_GT<-as.vector(kronecker(LT$EVD$values,EVD_Traits_values,"*"))
    LT$GxExT$eigenvecs_GT<-kronecker(LT$EVD$vectors,EVD_Traits_vectors)

    # X=Gamma*Lambda^{1/2}
    Aux<-LT$GxExT$eigenvecs_GT
    LT$GxExT$X<-Aux

    LT$GxExT<-setCov.MTME(LT=LT$GxExT,n=n,nEnvs=LT$nEnv,i=i,Sy=Sy_Env,nLT=nLT,R2=R2,
                           saveAt=saveAt,response_type=response_type,NoWhichNa=NoWhichNa)

    EVD_Env = eigen(LT$GxExT$Cov$SigmaEnv,symmetric = T)

    keep<-EVD_Env$values>1e-10
    eigenvals_E<-EVD_Env$values[keep]
    eigenvecs_E<-EVD_Env$vectors[,keep]

    LT$GxExT$eigenvecs_EG = kronecker(eigenvecs_E,LT$EVD$vectors)
    LT$GxExT$eigenvals_EG = as.vector(kronecker(eigenvals_E,LT$EVD$values))

    # X=Z*Gamma*Lambda^{1/2}
    LT$GxExT$Aux<-sweep(x=LT$GxExT$eigenvecs_EG,MARGIN=2,STATS=sqrt(LT$GxExT$eigenvals_EG),FUN="*")
    LT$GxExT$X<-LT$GxExT$Z%*%LT$GxExT$Aux

    LT$GxExT<-setLT.MTME(LT=LT$GxExT,n=n,ntraits=ntraits,i=i,Sy=Sy_Env,nLT=nLT,R2=R2,
              saveAt=saveAt,response_type=response_type,NoWhichNa=NoWhichNa)

    if(LT$Cov$type=="FA"){
      LT$Cov$F<-matrix(nrow = LT$nLines,ncol = ncol(LT$Cov$M),0)
      LT$F2<-matrix(nrow = LT$nLinesEnvs,ncol = ncol(LT$Cov$M),0)
      LT$W2<-matrix(nrow = ntraits,ncol = ncol(LT$Cov$M),0)
    }

  }

  if(!is.null(LT$GxT) && !is.null(LT$GxExT)){
    # Number of lines and Environments
    LT$nLines = ncol(LT$GxT$Z); LT$nEnv = n/LT$nLines; LT$nLinesTraits = LT$nLines*ntraits
    LT$nLinesEnvs = LT$nLines*LT$nEnv

    Aux<-sweep(x=LT$EVD$vectors,MARGIN=2,STATS=sqrt(LT$EVD$values),FUN="*")
    LT$X<-LT$GxT$Z%*%Aux

    LT<-setLT.BRR_mt(LT=LT,n=n,ntraits=ntraits,i=i,Sy=Sy,nLT=nLT,R2=R2,saveAt=saveAt,
                     response_type=response_type,NoWhichNa=NoWhichNa)

    # For environment sample variances
    y_Env = matrix(matrixcalc::vec(t(y)),ncol=LT$nEnv,byrow=F)

    if(response_type == "gaussian"){
      Sy_Env<-var(y_Env,na.rm = T)
      if(any(is.na(Sy_Env))){Sy_Env=diag(1,LT$nEnv)}
    }
    if(response_type == "DLN"){
      Sy_Env<-var(log(y_Env+1),na.rm = T)
      if(any(is.na(Sy_Env))){Sy_Env=diag(1,LT$nEnv)}
    }
    if(response_type%in%c("Poisson","PLN")){
      aux<-colMeans(log(y_Env+1),na.rm=T)
      Sy_Env<-outer(sqrt(aux),sqrt(aux))
      if(any(is.na(Sy_Env))){Sy_Env=diag(1,LT$nEnv)}
      if(!matrixcalc::is.positive.definite(Sy_Env)){
        eigens = eigen(Sy_Env)$values
        MinEigen = min(eigens)
        alpha = -MinEigen + 1e-6
        diag(Sy_Env)<-diag(Sy_Env) + alpha
      }
    }

    # Eigen decomposition of Traits co-variance
    EVD_Traits = eigen(LT$Cov$SigmaB,symmetric=TRUE)

    keep<-EVD_Traits$values>1e-10
    EVD_Traits_vectors<-EVD_Traits$vectors[,keep]
    EVD_Traits_values<-EVD_Traits$values[keep]

    # Interaction ExGxT components
    LT$GxExT$eigenvals_GT<-as.vector(kronecker(LT$EVD$values,EVD_Traits_values))
    LT$GxExT$eigenvecs_GT<-kronecker(LT$EVD$vectors,EVD_Traits_vectors)

    # X=Gamma*Lambda^{1/2}
    Aux<-LT$GxExT$eigenvecs_GT
    LT$GxExT$X<-Aux

    LT$GxExT<-setCov.MTME(LT=LT$GxExT,n=n,nEnvs=LT$nEnv,i=i,Sy=Sy_Env,nLT=nLT,R2=R2,
                          saveAt=saveAt,response_type=response_type,NoWhichNa=NoWhichNa)

    EVD_Env = eigen(LT$GxExT$Cov$SigmaEnv,symmetric = T)

    keep<-EVD_Env$values>1e-10
    eigenvals_E<-EVD_Env$values[keep]
    eigenvecs_E<-EVD_Env$vectors[,keep]

    LT$GxExT$eigenvecs_EG = kronecker(eigenvecs_E,LT$EVD$vectors)
    LT$GxExT$eigenvals_EG = as.vector(kronecker(eigenvals_E,LT$EVD$values))

    # X=Z*Gamma*Lambda^{1/2}
    LT$GxExT$Aux<-sweep(x=LT$GxExT$eigenvecs_EG,MARGIN=2,STATS=sqrt(LT$GxExT$eigenvals_EG),FUN="*")
    LT$GxExT$X<-LT$GxExT$Z%*%LT$GxExT$Aux

    LT$GxExT<-setLT.MTME(LT=LT$GxExT,n=n,ntraits=ntraits,i=i,Sy=Sy_Env,nLT=nLT,R2=R2,
                         saveAt=saveAt,response_type=response_type,NoWhichNa=NoWhichNa)

    if(LT$Cov$type=="FA"){
      LT$F2<-matrix(nrow = LT$nLinesEnvs,ncol = LT$Cov$nF,0)
    }

  }

  return(LT)

}

########################### Posterior distributions ############################

# Posterior distribution for latent variable in DLN-UT model
fL<-function(n=NULL,a=NULL,b=NULL,rv=NULL,varE=NULL)
{
  return(truncnorm::rtruncnorm(n,a=a,b=b,mean = rv,sd = sqrt(varE)))
}

# Posterior distribution for latent variable in DLN-MT model
lmulti<-function(n = NULL,ntraits = NULL,rv = NULL,SigmaE = NULL,SigmaEInv = NULL,
                   ay = NULL,by = NULL,type = NULL,Iters_latent = NULL,
                   Burn_latent = NULL,Thin_latent = NULL){
  if((matrixcalc::is.diagonal.matrix(SigmaE))){
    sds<-rep(sqrt(diag(SigmaE)),each = n)
    l<-matrix(truncnorm::rtruncnorm(n = n * ntraits,a = as.vector(ay),b = as.vector(by),
                           mean = as.vector(rv),sd = sds),
                nrow = n,ncol = ntraits,byrow = FALSE)
  }
  else{
    l<-.Call("rtmvnorm_gibbs",rv,ay,by,SigmaE,SigmaEInv,n = n,ntraits = ntraits,
               n_iter = Iters_latent,burn_in = Burn_latent,
               thin = Thin_latent,PACKAGE="BGLM")
  }
  return(l)
}

# Posterior distribution for latent variable in Poisson and Poisson-lognormal UT
# models
fL_P<-function(y_r=NULL,n=NULL,rv=NULL)
{
  return(BayesLogit::rpg(n,y_r,rv))
}

# Posterior distribution for latent variable in Poisson and Poisson-lognormal MT
# models
wmulti<-function(y_r=NULL,n=NULL,ntraits=NULL,rv=NULL)
{
  w<-matrix(BayesLogit::rpg(n * ntraits,as.vector(y_r),as.vector(rv)),
              nrow = n,ncol = ntraits,byrow = FALSE)
  return(w)
}

# Posterior distribution for intercept in DLN-UT model
fB0_DLN<-function(priorvarbeta0=NULL,varE=NULL,rv=NULL,n=NULL)
{
  varbeta0=1/((1/priorvarbeta0)+(n/varE)); Mu_0=varbeta0*((1/varE)*sum(rv))
  return(rnorm(1,Mu_0,sqrt(varbeta0)))
}

# Posterior distribution for intercept in DLN-MT model
fB0_DLN_multi<-function(priorSigma0Inv = NULL,SigmaEInv,rv = NULL,n = NULL) {

  tmp<-(n * SigmaEInv) + priorSigma0Inv
  Sigma0<-solve(tmp)
  sum_term<-SigmaEInv %*% colSums(rv)
  Mu0<-Sigma0 %*% (sum_term)

  return(MASS::mvrnorm(n=1,mu=Mu0,Sigma=Sigma0))
}

# Posterior distribution for intercept in Poisson and Poisson-lognormal UT models
fB0_P<-function(l=NULL,Syr=NULL,priorvarbeta0=NULL,rv=NULL)
{
  varbeta0=1/(sum(l)+(1/priorvarbeta0)); Mu_0=varbeta0*((Syr/2)-(l%*%rv))
  return(rnorm(1,Mu_0,sqrt(varbeta0)))
}

# Posterior distribution for intercept in Poisson and Poisson-lognormal MT models
fB0_P_multi<-function(l = NULL,Syr = NULL,priorSigma0Inv = NULL,rv = NULL){

  tmp<-diag(colSums(l)) + priorSigma0Inv
  Sigma0<-solve(tmp)
  sum_term<-Syr/2 - colSums(l*rv)
  Mu0<-Sigma0 %*% (sum_term)

  return(MASS::mvrnorm(n=1,mu=Mu0,Sigma=Sigma0))
}

# Posterior distribution for intercept in linear-UT model
fB0_G<-function(priorvarbeta0=NULL,varE=NULL,rv=NULL,n=NULL)
{
  varbeta0=1/((1/priorvarbeta0)+(n/varE)); Mu_0=varbeta0*((1/varE)*sum(rv))
  return(rnorm(1,Mu_0,sqrt(varbeta0)))
}

# Posterior distribution for intercept in linear-MT model
fB0_G_multi<-function(priorSigma0Inv = NULL,
                  SigmaEInv,rv = NULL,n = NULL) {

  tmp<-(n * SigmaEInv) + priorSigma0Inv
  Sigma0<-solve(tmp)
  sum_term<-SigmaEInv %*% colSums(rv)
  Mu0<-Sigma0 %*% (sum_term)

  return(MASS::mvrnorm(n=1,mu=Mu0,Sigma=Sigma0))
}

# Posterior Log-likelihood for linear-UT model
fllp_G<-function(rv=NULL,varE=NULL,y=NULL)
{
  return(sum(dnorm(y,y-rv,sqrt(varE),log=TRUE)))
}

# Partial Posterior Log-likelihood for linear-MT model
partial_fllp_G<-function(rv = NULL,SigmaE = NULL)
{
  error = rv
  n<-nrow(error)
  Linv<-solve(chol(SigmaE))
  Pll<--0.5 * n * log(det(SigmaE)) - 0.5 * sum(crossprod(t(error),Linv)^2)

  # Partial Log-likelihood
  return(Pll)
}

# Posterior Log-likelihood for DLN-UT model
fllp_DLN<-function(rv=NULL,varE=NULL,y=NULL)
{
  return(sum(log(plnorm(y+1,rv,sqrt(varE))-plnorm(y,rv,sqrt(varE)))))
}

# Posterior Log-likelihood for DLN-MT model
fllp_DLN_multi<-function(y = NULL,n = NULL,ntraits = NULL,a = NULL,b = NULL,
                     rv = NULL,SigmaE = NULL,U_qmc = NULL)
{
  if ((matrixcalc::is.diagonal.matrix(SigmaE))) {
    sds<-rep(sqrt(diag(SigmaE)),each = n)
    loglik_total<-sum(log(plnorm(as.vector(y)+1,as.vector(rv),sds)-plnorm(as.vector(y),as.vector(rv),sds)))
  }else {
    loglik_total<-.Call("fllp_qmc",a,b,rv,SigmaE,n,ntraits,U_qmc,PACKAGE="BGLM")
  }

  # Log-likelihood
  return(loglik_total)
}

# Posterior Log-likelihood for Poisson-UT model
fllp_P<-function(rv=NULL,y=NULL,r=NULL)
{
  return(sum(dpois(y,exp(rv + log(r)),log=TRUE)))
}

# Posterior Log-likelihood for Poisson-MT model
fllp_P_multi<-function(rv=NULL,y=NULL,r=NULL,ntraits=NULL)
{
  mu = exp(rv + log(r))
  loglik_total = 0
  for (i in 1:ntraits) {
    loglik_total = loglik_total + sum(dpois(y[,i],mu[,i],log=TRUE))
  }

  return(loglik_total)
}

# Posterior Log-likelihood for Poisson-lognormal-MT model
fllp_PLN_multi_R<-function(y=NULL,n=NULL,ntraits=NULL,logy_factorial=NULL,
                     rvPois=NULL,SigmaE=NULL,gh=NULL,Q=NULL,
                     nodes_mult=NULL,gh_weights_mult=NULL)
{
  if (matrixcalc::is.diagonal.matrix(SigmaE)) {
    loglik_total = 0.0
    varE<-diag(SigmaE)
    for (i in 1:ntraits) {
      loglik = .Call("fllp_PLN",y[,i],logy_factorial[,i],rvPois[,i],varE[i],
                     gh$nodes,gh$weights,n,Q,PACKAGE="BGLM")
      loglik_total = loglik_total + loglik
    }
  }else {
    L<-chol(SigmaE); gh_nodes_mult<-sqrt(2) * t(L %*% t(nodes_mult))
    loglik_total = .Call("fllp_PLN_multi",y,logy_factorial,rvPois,
                         gh_nodes_mult,gh_weights_mult,n,Q,ntraits,
                         PACKAGE="BGLM")
  }

  # Log-likelihood
  return(loglik_total)
}

#-------------------------------------------------------------------------------

#' Fit a Bayesian Uni-Trait Multi-Environment Regression Model
#'
#' Fits a Bayesian Uni-trait and Multi-environment model using MCMC sampling.
#' Supports Gaussian, Poisson, DLN and PLN responses.
#'
#' @param ETA List used to specify the linear predictor.
#' By default it is set to NULL, in which case only the intercept is included.
#' @param y Response vector. Can contain NA values.
#' @param response_type Type of response ("DLN","gaussian","Poisson","PLN").
#' For DLN, Poisson or PLN responses, the argument `y` must contain integer values,
#' while for gaussian responses any type of numeric value can be accepted.
#' @param nIter Number of MCMC iterations.
#' @param nBurnin Burn-in iterations.
#' @param nThin Thinning interval.
#' @param priorv Prior degrees of freedom for the error variance.
#' @param priorS Prior scale for the error variance.
#' @param R2 Proportion of variance explained a priori for the error term.
#' @param priorvarbeta0 Choose the variance, σ_0^2, to obtain a flat prior
#' for intercept,β_0∼N(0, σ_0^2). By default (internally) is 1e10.
#' @param rcontrol Specify the k-value used to control the difference between the
#' Conditional Mean and Conditional Variance of the Negative Binomial model.
#' Higher values of k indicate a smaller difference between the mean and variance,
#' approximating the equidispersion of a Poisson model.
#' @param type_prediction Used to specify if you want to predict with the posterior
#' mean or median in DLN model; default is mean.
#' @param intercept A list of length two. The first element is used to specify if
#' you want to include (TRUE) or not (FALSE) the intercept as an isolated term in
#' linear predictor; in case you specify FALSE, the second argument assigns the
#' position, according to the covariates order, where you want to include it.
#' The default values are list(TRUE,0).
#' @param verbose Used to show (TRUE) or not (FALSE) the progress of iterations
#' in the Gibbs Sampler; default is FALSE
#' @param saveAt Used to indicate UTME where to store the samples and to provide
#' a pre-fix to be appended to the names of the file where samples are stored.
#' By default samples are saved in the current working directory and no pre-fix
#' is added to the file names
#' @return A list containing posterior samples and summaries.
#'
#' @export
UTME<-function(
    ETA=NULL,
    y=NULL,
    response_type="DLN",
    nIter=1e3,
    nBurnin=1e2,
    nThin=1,
    R2=0.5,
    priorv=5,
    priorS=NULL,
    priorvarbeta0=NULL,
    rcontrol=1.15,
    type_prediction="mean",
    intercept=list(TRUE,0),
    verbose=F,
    saveAt=""
    ){

  # nIter and nBurnin validation
  if (!is.null(nIter) && !is.null(nBurnin)) {
    if (nIter <= nBurnin) {
      stop("\033[31m\nError: nIter should be greater than nBurnin by at least 1.\033[39m")
    }
  }

  # Response type validation
  if (!(response_type %in% c("DLN","Poisson","PLN","gaussian"))){
    stop("Only DLN,Poisson,PLN,or gaussian responses are allowed
           (note: evaluation is case sensitive)")
  }

  # Type prediction validation
  if (response_type == "DLN" && !(type_prediction %in% c("mean","median"))){
    stop("For the DLN model,only mean and median type predictions are allowed
           (note: evaluation is case sensitive)")
  }

  if (saveAt == "") {
    saveAt = paste(getwd(),"/",sep = "")
  }

  if(verbose){
    pb<-progress::progress_bar$new(
      format = "  Progress [:bar] :percent in :elapsed",
      total = nIter,clear = FALSE,width = 60
    )
  }

  y = as.vector(y)
  n = length(y)
  rowNames = names(y)
  whichNa = which(is.na(y))
  nNa = length(whichNa)
  obs_idx<-rep(TRUE,n)
  if (nNa > 0) obs_idx[whichNa]<-FALSE

  if(intercept[[1]]){ if (is.null(priorvarbeta0)) { priorvarbeta0<-1e+10 } }
  #-----------------------------------------------------------------------------

  # nLT = ifelse(is.null(ETA),0,length(ETA))
  nLT<-ifelse(is.null(ETA),0,length(ETA))
  nLT_total<-ifelse(is.null(ETA),0,sum(sapply(ETA,function(x){
    sublists<-x[sapply(x,is.list)]
    sublists<-sublists[names(sublists)!="Cov"]
    if(length(sublists)==0) 1 else length(sublists)
  })))

  if (!intercept[[1]]) {
    if (intercept[[2]] == 0 || intercept[[2]]>nLT) {
      stop("The intercept in position ",intercept[[2]]," is out of the linear predictor")
    }
  }

  if (nLT > 0) {
    if (is.null(names(ETA))) {
      names(ETA)<-rep("",nLT)
    }
    for (i in 1:nLT) {

      if (names(ETA)[i] == "") {
        ETA[[i]]$Name = paste("ETA_",i,sep = "")
      }
      else {
        ETA[[i]]$Name = paste("ETA_",names(ETA)[i],
                              sep = "")
      }
      # Model validation
      if (!(ETA[[i]]$model %in% c("FIXED","BRR","BayesA","BL","RKHS","RKHS_utme"))) {
        stop("Error in ETA[[",i,"]]"," model ",ETA[[i]]$model,
             " not implemented (note: evaluation is case sensitive)")
      }

      if(ETA[[i]]$model == "RKHS_utme"){
        ETA[[i]]$Cov$type = "UN"
        ntraits = 1; y1 = as.matrix(y)
        if(response_type == "gaussian"){
          Sy<-cov(y1,use = "pairwise.complete.obs")
        }
        if(response_type == "DLN"){
          Sy<-cov(log(y1+1),use = "pairwise.complete.obs")
        }
        if(response_type%in%c("Poisson","PLN")){
          aux<-colMeans(log(y1+1),na.rm=TRUE)
          Sy<-outer(sqrt(aux),sqrt(aux))
          if(!matrixcalc::is.positive.definite(Sy)){
            eigens = eigen(Sy)$values
            MinEigen = min(eigens)
            alpha = -MinEigen + 1e-6
            diag(Sy)<-diag(Sy) + alpha
          }
        }
      }

      # Update validation
      if (!(is.null(ETA[[i]]$update)) && !(ETA[[i]]$update %in% c("scalar","blocks"))){
        stop("Only scalar and blocks updates are allowed (note: evaluation is case sensitive)")
      }

      if(!ETA[[i]]$model%in%c("RKHS","RKHS_utme")){
        if(nNa != 0){
          ETA[[i]]$X_test = as.matrix(ETA[[i]]$X[whichNa,])
          if(!intercept[[1]] && i==intercept[[2]]){
            ETA[[i]]$X_test = cbind(Intercept = 1,ETA[[i]]$X_test)
          }
        }
      }

      ETA[[i]] = switch(ETA[[i]]$model,FIXED = setLT.Fixed(LT = ETA[[i]],whichNa = obs_idx,
                        n = n,y = y,response_type,saveAt = saveAt,intercept = intercept,i = i),
                        BRR = setLT.BRR(LT = ETA[[i]],whichNa = obs_idx,
                        n = n,y = y,R2 = R2,thin = nThin,nIter = nIter,
                        burnIn = nBurnin,nLT_total,response_type,saveAt = saveAt,intercept = intercept,i = i),
                        BayesA = setLT.BayesA(LT = ETA[[i]],whichNa = obs_idx,
                        n = n,y = y,R2 = R2,thin = nThin,nIter = nIter,
                        burnIn = nBurnin,nLT_total,response_type,saveAt = saveAt,intercept = intercept,i = i),
                        BL = setLT.BL(LT = ETA[[i]],whichNa = obs_idx,
                        n = n,y = y,R2 = R2,thin = nThin,nIter = nIter,
                        burnIn = nBurnin,nLT_total,response_type,saveAt = saveAt,intercept = intercept,i = i),
                        RKHS = setLT.RKHS(LT = ETA[[i]],whichNa = obs_idx,
                        n = n,y = y,R2 = R2,nLT_total,response_type,saveAt = saveAt,i = i),
                        RKHS_utme=setLT.RKHS_mtme(LT=ETA[[i]],n=n,ntraits=ntraits,i=i,y=y1,Sy=Sy,
                        nLT=nLT_total,R2=R2,saveAt=saveAt,response_type=response_type,NoWhichNa=obs_idx))
    }
  }

  if(!(response_type == "DLN" && type_prediction == "median")){
    y_pred = rep(0,n)
    y_pred2 = rep(0,n)
  }

  rv_mean = rep(0,n)
  #-----------------------------------------------------------------------------

  # Priors
  #-----------------------------------------------------------------------------

  if(response_type == "gaussian"){
    if (is.null(priorS)){
      priorS = (priorv+2)*(1-R2)*var(y,na.rm=TRUE)
    }
    varE = 1/rgamma(1,priorv/2,priorS/2)
    post_varE = 0
    post_varE2 = 0
    f_varE<-file(description = paste(saveAt,"varE.dat",sep = ""),
                   open = "w")
  }
  if(response_type == "DLN"){
    if (is.null(priorS)){
      priorS = (priorv+2)*(1-R2)*var(log(y+1),na.rm=TRUE)
    }
    varE = 1/rgamma(1,priorv/2,priorS/2)
    post_varE = 0
    post_varE2 = 0
    f_varE<-file(description = paste(saveAt,"varE.dat",sep = ""),
                   open = "w")
  }
  if(response_type == "PLN"){
    if (is.null(priorS)){
      priorS = (priorv+2)*(1-R2)*mean(log(y+1),na.rm=TRUE)/(1/2)
    }
    varE = 1/rgamma(1,priorv/2,priorS/2)
    post_varE = 0
    post_varE2 = 0
    f_varE<-file(description = paste(saveAt,"varE.dat",sep = ""),
                   open = "w")
  }

  #-----------------------------------------------------------------------------

  # Initial values
  #-----------------------------------------------------------------------------
  if (nNa > 0) {
    f_y_posterior<-file(description = paste(saveAt,"y_posterior.dat"),open = "w")
  }

  if(intercept[[1]]){
    beta0 = 0
    post_beta0 = 0
    post_beta02 = 0
    f_beta0<-file(description = paste(saveAt,"beta0.dat",sep = ""),
                    open = "w")
  }
  post_logLik = 0
  post_logLik2 = 0
  f_logLik<-file(description = paste(saveAt,"logLik.dat",sep = ""),
                   open = "w")
  nk = 0
  #-----------------------------------------------------------------------------

  # Algorithm
  #-----------------------------------------------------------------------------
  time<-proc.time()[3]
  yStar = y
  if (nNa > 0) {
    yStar[whichNa] = 0
  }
  if(response_type == "DLN"){
    ay = log(yStar);  by = log(yStar+1)
    rv = rep(0,n)
  }
  if(response_type == "gaussian"){
    rv = yStar - rep(0,n)
  }
  if(response_type == "Poisson"){
    r = (yStar+1) * 10^(rcontrol)
    y_r = yStar+r
    yr = (yStar - r) / 2
    Syr = sum(yStar-r)
    post_r = rep(0,n)
    post_r2 = rep(0,n)
    rv = rep(0,n) - log(r)
  }
  if(response_type == "PLN"){
    gh<-statmod::gauss.quad(20,kind = "hermite")
    Q = length(gh$nodes)
    logy_factorial = lgamma(y[obs_idx]+1)
    U = rep(0,n)
    r = ((1/mean(yStar))+1) * 10^(rcontrol)
    y_r = yStar+r
    yr = (yStar - r) / 2
    Syr = sum(yStar-r)
    post_r = 0
    post_r2 = 0
    rv = rep(0,n) + U - log(r)
  }

  fL_fun<-switch(response_type,
                   DLN = function(...) fL(n=n,a=ay,b=by,rv=rv,varE=varE),
                   Poisson = function(...) fL_P(y_r=y_r,n=n,rv=rv),
                   PLN = function(...) fL_P(y_r=y_r,n=n,rv=rv)
  )

  fB0_fun<-switch(response_type,
                    DLN = function(...) fB0_DLN(priorvarbeta0=priorvarbeta0,varE=varE,rv=rv,n=n),
                    Poisson = function(...) fB0_P(l=l,Syr=Syr,priorvarbeta0=priorvarbeta0,rv=rv),
                    PLN = function(...) fB0_P(l=l,Syr=Syr,priorvarbeta0=priorvarbeta0,rv=rv),
                    gaussian = function(...) fB0_G(priorvarbeta0=priorvarbeta0,varE=varE,rv=rv,n=n)
  )

  fllp_fun<-switch(response_type,
                     DLN = function(...) fllp_DLN(rv=rv[obs_idx],varE=varE,y=y[obs_idx]),
                     PLN = function(...) .Call("fllp_PLN",y[obs_idx],logy_factorial,
                                               rvPois[obs_idx],varE,gh$nodes,gh$weights,
                                               n-nNa,Q,PACKAGE="BGLM"),
                     Poisson = function(...) fllp_P(rv=rv[obs_idx],y=y[obs_idx],r=r[obs_idx]),
                     gaussian = function(...) fllp_G(rv=rv[obs_idx],varE=varE,y=y[obs_idx])
  )

  # GIBBS SAMPLER
  #*****************************************************************************

  message("Generating ",nIter," samples,discarting ",nBurnin,',and using a thinning of ',nThin,
          " for the ",response_type," model")

  for(i in 1:nIter){

    deltaSS = 0
    deltadf = 0

    if(!is.null(fL_fun)) l<-fL_fun()

    if(response_type=="DLN"){rv = l - rv}

    if(intercept[[1]]){
      ifelse(response_type%in%c("DLN","gaussian"),rv<-rv+beta0,rv<-rv-beta0)

      beta0<-fB0_fun()

      ifelse(response_type%in%c("DLN","gaussian"),rv<-rv-beta0,rv<-rv+beta0)
    }

    if (nLT > 0) {
      for (j in 1:nLT) {

        if (ETA[[j]]$model == "FIXED") {
          if(response_type%in%c("DLN","gaussian")){
            # Sampling from full conditional of Bj's
            if(ETA[[j]]$update == "scalar"){
              B = .Call("fBj_DLN_G",ETA[[j]]$Bv,ETA[[j]]$priorvarB,
                        varE,rv,ETA[[j]]$X,ETA[[j]]$x2,n,ETA[[j]]$p,
                        PACKAGE="BGLM")
            }
            if(ETA[[j]]$update == "blocks"){
              B = .Call("fB_DLN_G",ETA[[j]]$Bv,ETA[[j]]$InvpriorvarB,rv,
                        ETA[[j]]$X,ETA[[j]]$tXX,ETA[[j]]$tX,varE,
                        ETA[[j]]$blocks,ETA[[j]]$nblocks,n,
                        PACKAGE="BGLM")
            }
          }

          if(response_type %in% c("Poisson","PLN")){
            # Sampling from full conditional of Bj's
            if(ETA[[j]]$update == "scalar"){
              B = .Call("fBj_P",l,yr,ETA[[j]]$Bv,ETA[[j]]$priorvarB,
                        rv,ETA[[j]]$X,ETA[[j]]$X2,n,ETA[[j]]$p,PACKAGE="BGLM")
            }
            if(ETA[[j]]$update == "blocks"){
              B = .Call("fB_P",l,yr,ETA[[j]]$Bv,ETA[[j]]$InvpriorvarB,
                        rv,ETA[[j]]$X,ETA[[j]]$tX,ETA[[j]]$blocks,
                        ETA[[j]]$nblocks,n,PACKAGE="BGLM")
            }
          }

          ETA[[j]]$Bv = B[[1]]
          rv = B[[2]]
        }

        if (ETA[[j]]$model == "BRR") {
          if(response_type%in%c("DLN","gaussian")){
            # Sampling from full conditional of Bj's
            if(ETA[[j]]$update == "scalar"){
              B = .Call("fBj_DLN_G",ETA[[j]]$Bv,rep(ETA[[j]]$varB,ETA[[j]]$p),
                        varE,rv,ETA[[j]]$X,ETA[[j]]$x2,n,ETA[[j]]$p,
                        PACKAGE="BGLM")
            }
            if(ETA[[j]]$update == "blocks"){
              B = .Call("fB_DLN_G",ETA[[j]]$Bv,diag(1/ETA[[j]]$varB,ETA[[j]]$p),rv,
                        ETA[[j]]$X,ETA[[j]]$tXX,ETA[[j]]$tX,varE,
                        ETA[[j]]$blocks,ETA[[j]]$nblocks,n,PACKAGE="BGLM")
            }
          }

          if(response_type %in% c("Poisson","PLN")){
            # Sampling from full conditional of Bj's
            if(ETA[[j]]$update == "scalar"){
              B = .Call("fBj_P",l,yr,ETA[[j]]$Bv,rep(ETA[[j]]$varB,ETA[[j]]$p),
                        rv,ETA[[j]]$X,ETA[[j]]$X2,n,ETA[[j]]$p,PACKAGE="BGLM")
            }
            if(ETA[[j]]$update == "blocks"){
              B = .Call("fB_P",l,yr,ETA[[j]]$Bv,diag(1/ETA[[j]]$varB,ETA[[j]]$p),
                        rv,ETA[[j]]$X,ETA[[j]]$tX,ETA[[j]]$blocks,
                        ETA[[j]]$nblocks,n,PACKAGE="BGLM")
            }
          }

          ETA[[j]]$Bv = B[[1]]

          # Sampling from full conditional of varB
          VA = ETA[[j]]$priorvB+ETA[[j]]$p
          SA = ETA[[j]]$priorSB+sum((ETA[[j]]$Bv)^2)
          ETA[[j]]$varB = 1/rgamma(1,VA/2,SA/2)

          rv = B[[2]]
        }

        if (ETA[[j]]$model == "BayesA") {
          if(response_type%in%c("DLN","gaussian")){
            # Sampling from full conditional of Bj's
            if(ETA[[j]]$update == "scalar"){
              B = .Call("fBj_DLN_G",ETA[[j]]$Bv,ETA[[j]]$varB,varE,rv,
                        ETA[[j]]$X,ETA[[j]]$x2,n,ETA[[j]]$p,PACKAGE="BGLM")
            }
            if(ETA[[j]]$update == "blocks"){
              B = .Call("fB_DLN_G",ETA[[j]]$Bv,diag(1/ETA[[j]]$varB),rv,
                        ETA[[j]]$X,ETA[[j]]$tXX,ETA[[j]]$tX,varE,
                        ETA[[j]]$blocks,ETA[[j]]$nblocks,n,PACKAGE="BGLM")
            }
          }

          if(response_type %in% c("Poisson","PLN")){
            # Sampling from full conditional of Bj's
            if(ETA[[j]]$update == "scalar"){
              B = .Call("fBj_P",l,yr,ETA[[j]]$Bv,ETA[[j]]$varB,rv,ETA[[j]]$X,
                        ETA[[j]]$X2,n,ETA[[j]]$p,PACKAGE="BGLM")
            }
            if(ETA[[j]]$update == "blocks"){
              B = .Call("fB_P",l,yr,ETA[[j]]$Bv,diag(1/ETA[[j]]$varB),
                        rv,ETA[[j]]$X,ETA[[j]]$tX,ETA[[j]]$blocks,
                        ETA[[j]]$nblocks,n,PACKAGE="BGLM")
            }
          }

          ETA[[j]]$Bv = B[[1]]

          # Sampling from full conditional of varB
          VA = ETA[[j]]$priorvB + 1
          SA = ETA[[j]]$SB + (ETA[[j]]$Bv)^2
          ETA[[j]]$varB = 1/rgamma(ETA[[j]]$p,VA/2,SA/2)

          rv = B[[2]]
          tmpShape = ETA[[j]]$shape0 + ((ETA[[j]]$p * ETA[[j]]$priorvB)/2)
          tmpRate = ETA[[j]]$rate0 + (sum(1/ETA[[j]]$varB)/2)
          ETA[[j]]$SB = rgamma(n = 1,shape = tmpShape,rate = tmpRate)
        }

        if (ETA[[j]]$model == "BL") {

          if(response_type=="Poisson"){
            varBj=ETA[[j]]$tau2
          }else{varBj=varE*ETA[[j]]$tau2}

          if(response_type%in%c("DLN","gaussian")){
            # Sampling from full conditional of Bj's
            if(ETA[[j]]$update == "scalar"){
              B = .Call("fBj_DLN_G",ETA[[j]]$Bv,varBj,varE,rv,ETA[[j]]$X,
                        ETA[[j]]$x2,n,ETA[[j]]$p,PACKAGE="BGLM")
            }
            if(ETA[[j]]$update == "blocks"){
              B = .Call("fB_DLN_G",ETA[[j]]$Bv,diag(1/varBj),rv,
                        ETA[[j]]$X,ETA[[j]]$tXX,ETA[[j]]$tX,varE,
                        ETA[[j]]$blocks,ETA[[j]]$nblocks,n,PACKAGE="BGLM")
            }
          }

          if(response_type == "Poisson"){
            # Sampling from full conditional of Bj's
            if(ETA[[j]]$update == "scalar"){
              B = .Call("fBj_P",l,yr,ETA[[j]]$Bv,varBj,rv,ETA[[j]]$X,
                        ETA[[j]]$X2,n,ETA[[j]]$p,PACKAGE="BGLM")
            }
            if(ETA[[j]]$update == "blocks"){
              B = .Call("fB_P",l,yr,ETA[[j]]$Bv,diag(1/varBj),rv,ETA[[j]]$X,
                        ETA[[j]]$tX,ETA[[j]]$blocks,ETA[[j]]$nblocks,n,
                        PACKAGE="BGLM")
            }
          }

          if(response_type == "PLN"){
            # Sampling from full conditional of Bj's
            if(ETA[[j]]$update == "scalar"){
              B = .Call("fBj_P",l,yr,ETA[[j]]$Bv,varBj,rv,ETA[[j]]$X,
                        ETA[[j]]$X2,n,ETA[[j]]$p,PACKAGE="BGLM")
            }
            if(ETA[[j]]$update == "blocks"){
              B = .Call("fB_P",l,yr,ETA[[j]]$Bv,diag(1/varBj),rv,ETA[[j]]$X,
                        ETA[[j]]$tX,ETA[[j]]$blocks,ETA[[j]]$nblocks,n,
                        PACKAGE="BGLM")
            }
          }

          ETA[[j]]$Bv = B[[1]]
          rv = B[[2]]
          if(response_type=="Poisson"){
            nu = sqrt(ETA[[j]]$lambda^2 / ETA[[j]]$Bv^2)
          }
          else{nu = sqrt(varE * ETA[[j]]$lambda^2 / ETA[[j]]$Bv^2)}
          tmp = NULL
          try(tmp<-statmod::rinvgauss(n = ETA[[j]]$p,mean = nu,shape = ETA[[j]]$lambda2))
          if (!is.null(tmp) && !any(tmp < 0)) {
            if (!any(is.na(sqrt(tmp)))) {
              ETA[[j]]$tau2 = 1/tmp
            }
            else {
              warning(paste("tau2 was not updated in iteration",
                            i,"due to numeric problems with beta\n",
                            sep = " "),immediate. = TRUE)
            }
          }
          else {
            warning(paste("tau2 was not updated  in iteration",
                          i,"due to numeric problems with beta\n",
                          sep = " "),immediate. = TRUE)
          }
          if (ETA[[j]]$type == "gamma") {
            rate = sum(ETA[[j]]$tau2)/2 + ETA[[j]]$rate
            shape = ETA[[j]]$p + ETA[[j]]$shape
            ETA[[j]]$lambda2 = rgamma(rate = rate,shape = shape,
                                      n = 1)
            if (!is.na(ETA[[j]]$lambda2)) {
              ETA[[j]]$lambda = sqrt(ETA[[j]]$lambda2)
            }
            else {
              warning(paste("lambda was not updated in iteration",
                            i,"due to numeric problems with beta\n",
                            sep = " "),immediate. = TRUE)
            }
          }
          deltaSS = deltaSS + sum((ETA[[j]]$Bv/sqrt(ETA[[j]]$tau2))^2)
          deltadf = deltadf + ETA[[j]]$p
        }

        if (ETA[[j]]$model == "RKHS") {
          varU = ETA[[j]]$varU * ETA[[j]]$d

          if(response_type%in%c("DLN","gaussian")){
            # Sampling from full conditional of Ui's
            if(ETA[[j]]$update == "scalar"){
              res = .Call("fUi_DLN_G",ETA[[j]]$uStar,varU,varE,rv,ETA[[j]]$V,
                          n,ETA[[j]]$levelsU,PACKAGE="BGLM")
            }
            if(ETA[[j]]$update == "blocks"){
              res = .Call("fU_DLN_G",ETA[[j]]$uStar,diag(1/varU),rv,ETA[[j]]$V,
                          ETA[[j]]$tV,varE,ETA[[j]]$blocks,ETA[[j]]$nblocks,n,
                          PACKAGE="BGLM")
            }
          }

          if(response_type %in% c("Poisson","PLN")){
            # Sampling from full conditional of Ui's
            if(ETA[[j]]$update == "scalar"){
              res = .Call("fBj_P",l,yr,ETA[[j]]$uStar,varU,rv,ETA[[j]]$V,
                          ETA[[j]]$V2,n,ETA[[j]]$levelsU,PACKAGE="BGLM")
            }
            if(ETA[[j]]$update == "blocks"){
              res = .Call("fB_P",l,yr,ETA[[j]]$uStar,diag(1/varU),rv,
                          ETA[[j]]$V,ETA[[j]]$tV,ETA[[j]]$blocks,
                          ETA[[j]]$nblocks,n,PACKAGE="BGLM")
            }
          }

          ETA[[j]]$uStar = res[[1]]
          ETA[[j]]$u = as.vector(ETA[[j]]$V %*% res[[1]])

          # Sampling from full conditional of varU
          tmp = ETA[[j]]$uStar/sqrt(ETA[[j]]$d)
          VU = ETA[[j]]$priorvU+ETA[[j]]$levelsU
          SU = ETA[[j]]$priorSU+as.numeric(crossprod(tmp))
          ETA[[j]]$varU = 1/rgamma(1,VU/2,SU/2)

          rv = res[[2]]
        }

        if (ETA[[j]]$model=="RKHS_utme") {
          if(!is.null(ETA[[j]]$GxT)){

            # Sampling from posterior of Traits parameters
            if(response_type%in%c("DLN","gaussian")){
              # Sampling from full conditional of Bj's
              B1 = .Call("fBj_DLN_G_mtme",ETA[[j]]$Bv,ETA[[j]]$Cov$SigmaBInv,
                         as.matrix(1/varE),as.matrix(rv),ETA[[j]]$X,ETA[[j]]$x2,
                         n,ETA[[j]]$p,ntraits,PACKAGE="BGLM")
            }

            if(response_type%in%c("Poisson","PLN")){
              # Sampling from full conditional of Bj's
              B1 = .Call("fBj_P_mtme",as.matrix(l),as.matrix(yr),ETA[[j]]$Bv,
                         ETA[[j]]$Cov$SigmaBInv,as.matrix(rv),ETA[[j]]$X,
                         ETA[[j]]$X2,ETA[[j]]$p,ntraits,n,PACKAGE="BGLM")
            }

            ETA[[j]]$Bv = B1[[1]]
            rv = B1[[2]]

          }

          if(!is.null(ETA[[j]]$GxT) && !is.null(ETA[[j]]$GxExT)){
            if(ETA[[j]]$Cov$type %in% c("UN","DIAG")) {
              SE = crossprod(ETA[[j]]$Bv); SE1 = crossprod(ETA[[j]]$GxExT$Bv)
            }

            if(ETA[[j]]$Cov$type %in% c("FA","REC")) {
              U1 = ETA[[j]]$Bv; U2 = ETA[[j]]$GxExT$Bv
            }

            if(ETA[[j]]$Cov$type %in% c("FA")) {
              F1 = ETA[[j]]$Cov$F; F2 = ETA[[j]]$F2; const1=1
              W1 = ETA[[j]]$Cov$W; W2 = W1; const2=1
            }

          }
          if(!is.null(ETA[[j]]$GxT) && is.null(ETA[[j]]$GxExT)){
            if(ETA[[j]]$Cov$type %in% c("UN","DIAG")) {
              SE = crossprod(ETA[[j]]$Bv); SE1 = ETA[[j]]$SE1
            }

            if(ETA[[j]]$Cov$type %in% c("FA","REC")) {
              U1 = ETA[[j]]$Bv; U2 = ETA[[j]]$U2
            }

            if(ETA[[j]]$Cov$type %in% c("FA")) {
              F1 = ETA[[j]]$Cov$F; F2 = ETA[[j]]$F2; const1=1
              W1 = ETA[[j]]$Cov$W; W2 = ETA[[j]]$W2; const2=Inf
            }

          }
          if(is.null(ETA[[j]]$GxT) && !is.null(ETA[[j]]$GxExT)){
            if(ETA[[j]]$Cov$type %in% c("UN","DIAG")) {
              SE = ETA[[j]]$SE; SE1 = crossprod(ETA[[j]]$GxExT$Bv)
            }

            if(ETA[[j]]$Cov$type %in% c("FA","REC")) {
              U1 = ETA[[j]]$U1; U2 = ETA[[j]]$GxExT$Bv
            }

            if(ETA[[j]]$Cov$type %in% c("FA")) {
              F1 = ETA[[j]]$Cov$F; F2 = ETA[[j]]$F2; const1=Inf
              W2 = ETA[[j]]$W2; W1 = W2; const2=1
            }

          }

          Term1 = (ETA[[j]]$Cov$priorvB+ETA[[j]]$nLines+ETA[[j]]$nLinesEnvs)/2
          Term2 = (SE + SE1 + ETA[[j]]$Cov$priorSB)/2
          ETA[[j]]$Cov$SigmaB<-as.matrix(1/rgamma(1,Term1,Term2))

          ETA[[j]]$Cov$SigmaBInv<-1/(ETA[[j]]$Cov$SigmaB)

          # Sampling from posterior of Environment parameters

          if(!is.null(ETA[[j]]$GxExT)){

            Q_old=ETA[[j]]$GxExT$eigenvecs_EG
            Lambda_old = ETA[[j]]$GxExT$eigenvals_EG

            if(response_type%in%c("DLN","gaussian")){
              # Sampling from full conditional of Bj's
              B2 = .Call("fBj_DLN_G_mtme",ETA[[j]]$GxExT$Bv,ETA[[j]]$Cov$SigmaBInv,
                         as.matrix(1/varE),as.matrix(rv),ETA[[j]]$GxExT$X,
                         ETA[[j]]$GxExT$x2,n,ETA[[j]]$GxExT$p,ntraits,
                         PACKAGE="BGLM")
            }

            if(response_type%in%c("Poisson","PLN")){
              # Sampling from full conditional of Bj's
              B2 = .Call("fBj_P_mtme",as.matrix(l),as.matrix(yr),ETA[[j]]$GxExT$Bv,
                         ETA[[j]]$Cov$SigmaBInv,as.matrix(rv),ETA[[j]]$GxExT$X,
                         ETA[[j]]$GxExT$X2,ETA[[j]]$GxExT$p,ntraits,n,
                         PACKAGE="BGLM")
            }

            ETA[[j]]$GxExT$Bv = B2[[1]]
            rv = B2[[2]]

            MEnv = matrix(matrixcalc::vec(t(ETA[[j]]$GxExT$Bv)),ncol = ETA[[j]]$nEnv,byrow = F)
            SE2 = crossprod(MEnv)

            if (ETA[[j]]$GxExT$Cov$type == "UN") {
              ETA[[j]]$GxExT$Cov$SigmaEnv<-MCMCpack::riwish(v = ETA[[j]]$GxExT$Cov$priorvE+ETA[[j]]$nLinesTraits,
                                                    S = SE2 + ETA[[j]]$GxExT$Cov$priorSE)
            }
            if (ETA[[j]]$GxExT$Cov$type == "DIAG") {
              VVB = ETA[[j]]$GxExT$Cov$priorvE+ETA[[j]]$nLinesTraits
              SSB = ETA[[j]]$GxExT$Cov$priorSE+diag(SE2)
              ETA[[j]]$GxExT$Cov$SigmaEnv<-diag(1/rgamma(ETA[[j]]$nEnv,VVB/2,SSB/2))
            }
            if (ETA[[j]]$GxExT$Cov$type == "FA") {
              tmp<-sample_G0_FA(U = MEnv,F = ETA[[j]]$GxExT$Cov$F,
                                  M = ETA[[j]]$GxExT$Cov$M,B = ETA[[j]]$GxExT$Cov$W,PSI = ETA[[j]]$GxExT$Cov$PSI,
                                  ntraits = ETA[[j]]$nEnv,nF = ETA[[j]]$GxExT$Cov$nF,nD = ETA[[j]]$GxExT$Cov$nD,
                                  df0 = ETA[[j]]$GxExT$Cov$priorvE,S0 = ETA[[j]]$GxExT$Cov$priorSE,
                                  priorVar = ETA[[j]]$GxExT$Cov$var,varimaxRotate = ETA[[j]]$GxExT$Cov$varimax)
              ETA[[j]]$GxExT$Cov$F<-tmp$F
              ETA[[j]]$GxExT$Cov$PSI<-tmp$PSI
              ETA[[j]]$GxExT$Cov$W<-tmp$B
              ETA[[j]]$GxExT$Cov$SigmaEnv<-tmp$G
              rm(tmp)
            }
            if (ETA[[j]]$GxExT$Cov$type == "REC") {
              tmp<-sample_G0_REC(U = MEnv,M = ETA[[j]]$GxExT$Cov$M,
                                   PSI = ETA[[j]]$GxExT$Cov$PSI,ntraits = ETA[[j]]$nEnv,
                                   priorVar = ETA[[j]]$GxExT$Cov$var,df0 = ETA[[j]]$GxExT$Cov$priorvE,
                                   S0 = ETA[[j]]$GxExT$Cov$priorSE)
              ETA[[j]]$GxExT$Cov$SigmaEnv<-tmp$G
              ETA[[j]]$GxExT$Cov$W<-tmp$B
              ETA[[j]]$GxExT$Cov$PSI<-tmp$PSI
              rm(tmp)
            }

            EVD_Env = eigen(ETA[[j]]$GxExT$Cov$SigmaEnv,symmetric = T)

            keep<-EVD_Env$values>1e-10
            eigenvals_E<-EVD_Env$values[keep]
            eigenvecs_E<-EVD_Env$vectors[,keep]

            ETA[[j]]$GxExT$eigenvecs_EG = .Call("kron",eigenvecs_E,ETA[[j]]$EVD$vectors,
                                                PACKAGE="BGLM")
            ETA[[j]]$GxExT$eigenvals_EG = as.vector(.Call("kron_vec",eigenvals_E,
                                                          ETA[[j]]$EVD$values,
                                                          PACKAGE="BGLM"))

            Q_new=ETA[[j]]$GxExT$eigenvecs_EG
            Lambda_new = ETA[[j]]$GxExT$eigenvals_EG

            Bv_old<-ETA[[j]]$GxExT$Bv
            Bv_scaled<-sweep(Bv_old,1,sqrt(Lambda_old),FUN="*")
            Bv_proj<-.Call("prod",t(Q_new),.Call("prod",Q_old,Bv_scaled,
                                                 PACKAGE="BGLM"),PACKAGE="BGLM")
            Bv_new<-sweep(Bv_proj,1,1/sqrt(Lambda_new),FUN="*")

            ETA[[j]]$GxExT$Bv<-Bv_new

            rv<-.Call("fXb_multi_prod",ETA[[j]]$GxExT$X,Bv_old,rv,1.0,PACKAGE="BGLM")

            ETA[[j]]$GxExT$Aux<-sweep(x=Q_new,MARGIN=2L,STATS=sqrt(Lambda_new),FUN="*")

            ETA[[j]]$GxExT$X<-.Call("prod",ETA[[j]]$GxExT$Z,ETA[[j]]$GxExT$Aux,
                                    PACKAGE="BGLM")

            ETA[[j]]$GxExT$x2 = as.vector(colSums(ETA[[j]]$GxExT$X^2))
            ETA[[j]]$GxExT$p<-ncol(ETA[[j]]$GxExT$X)

            if(response_type%in%c("Poisson","PLN")){ETA[[j]]$GxExT$X2<-ETA[[j]]$GxExT$X^2}

            rv<-.Call("fXb_multi_prod",ETA[[j]]$GxExT$X,Bv_new,rv,-1.0,
                      PACKAGE="BGLM")
          }
        }

      }
    }

    # Predictions,residual variances and missing values
    if(response_type == "DLN"){
      SAE = priorS + sum(rv^2) + deltaSS
      VAE = priorv + n + deltadf
      varE = 1 / rgamma(1,VAE / 2,SAE / 2)
      rv = l - rv
      if(type_prediction == "mean"){
        y_tr = floor(exp(rv+rnorm(n=n,sd=sqrt(varE))))
      }
      if(nNa>0){
        yStar[whichNa] = floor(exp(rv[whichNa]+rnorm(n=nNa,sd=sqrt(varE))))
        ay[whichNa] = log(yStar[whichNa]);  by[whichNa] = log(yStar[whichNa]+1)

        #Save the posterior samples
        write(yStar[whichNa],ncolumns = nNa,file = f_y_posterior,append = TRUE)
      }
    }

    if(response_type == "Poisson"){
      rvPois = rv + log(r)
      y_tr = exp(rvPois)
      if(nNa>0){
        yStar[whichNa] = rpois(n=nNa,lambda=y_tr[whichNa])

        #Save the posterior samples
        write(yStar[whichNa],ncolumns = nNa,file = f_y_posterior,append = TRUE)
      }
    }

    if(response_type == "PLN"){
      Us = .Call("fUi",l,yr,U,varE,rv,n,PACKAGE="BGLM")
      U = Us[[1]]
      rv = Us[[2]]
      rvPois = rv + log(r) - U
      SAE = priorS + sum(U^2) + deltaSS
      VAE = priorv + n + deltadf
      varE = 1 / rgamma(1,VAE / 2,SAE / 2)
      y_tr = exp(rvPois + varE/2)
      if(nNa>0){
        yStar[whichNa] = rpois(n=nNa,lambda=exp(rvPois[whichNa]+U[whichNa]))

        #Save the posterior samples
        write(yStar[whichNa],ncolumns = nNa,file = f_y_posterior,append = TRUE)
      }
    }

    if(response_type == "gaussian"){
      y_tr = yStar - rv
      SAE = priorS + sum(rv^2) + deltaSS
      VAE = priorv + n + deltadf
      varE = 1 / rgamma(1,VAE / 2,SAE / 2)
      if(nNa>0){
        yStar[whichNa] = y_tr[whichNa] + rnorm(n=nNa,sd=sqrt(varE))
        rv[whichNa] = yStar[whichNa] - y_tr[whichNa]

        #Save the posterior samples
        write(yStar[whichNa],ncolumns = nNa,file = f_y_posterior,append = TRUE)
      }
    }

    # loglik
    loglik<-fllp_fun()

    if(response_type == "Poisson"){
      rv = rv + log(r)
      Lambda = exp(rv)
      r = (Lambda+1) * 10^(rcontrol)
      y_r = yStar+r
      yr = (yStar - r) / 2
      Syr = sum(yStar-r)
      rv = rv - log(r)
    }

    if(response_type == "PLN"){
      rv = rv + log(r)
      r = ((1/(exp(varE)-1))+1) * 10^(rcontrol)
      y_r = yStar+r
      yr = (yStar - r) / 2
      Syr = sum(yStar-r)
      rv = rv - log(r)
    }

    if (i > nBurnin)
    {
      if(i%%nThin==0){
        nk = nk + 1

        if(intercept[[1]]){
          write(beta0,ncolumns = length(beta0),file = f_beta0,append = TRUE,
                sep = " ")
        }

        write(loglik,ncolumns = length(loglik),file = f_logLik,append = TRUE,
              sep = " ")

        if(intercept[[1]]){
          # Mean of Beta0
          post_beta0 = (beta0 + (nk - 1) * post_beta0) / nk
          post_beta02 = (beta0^2 + (nk - 1) * post_beta02) / nk
        }

        if (nLT > 0) {
          for (j in 1:nLT) {
            if (ETA[[j]]$model == "FIXED") {
              # Mean of Beta's
              ETA[[j]]$post_Bv = (ETA[[j]]$Bv + (nk - 1) * ETA[[j]]$post_Bv) / nk
              ETA[[j]]$post_Bv2 = (ETA[[j]]$Bv^2 + (nk - 1) * ETA[[j]]$post_Bv2) / nk

              write(ETA[[j]]$Bv,ncolumns = ETA[[j]]$p,file = ETA[[j]]$fileOut,append = TRUE)
            }

            if (ETA[[j]]$model == "BRR") {
              # Mean of Beta's
              ETA[[j]]$post_Bv = (ETA[[j]]$Bv + (nk - 1) * ETA[[j]]$post_Bv) / nk
              ETA[[j]]$post_Bv2 = (ETA[[j]]$Bv^2 + (nk - 1) * ETA[[j]]$post_Bv2) / nk

              # Mean of varB
              ETA[[j]]$post_varB = (ETA[[j]]$varB + (nk - 1) * ETA[[j]]$post_varB) / nk
              ETA[[j]]$post_varB2 = (ETA[[j]]$varB^2 + (nk - 1) * ETA[[j]]$post_varB2) / nk

              write(ETA[[j]]$varB,file = ETA[[j]]$fileOut,append = TRUE)
            }

            if (ETA[[j]]$model == "BayesA") {
              # Mean of Beta's
              ETA[[j]]$post_Bv = (ETA[[j]]$Bv + (nk - 1) * ETA[[j]]$post_Bv) / nk
              ETA[[j]]$post_Bv2 = (ETA[[j]]$Bv^2 + (nk - 1) * ETA[[j]]$post_Bv2) / nk

              # Mean of varB
              ETA[[j]]$post_varB = (ETA[[j]]$varB + (nk - 1) * ETA[[j]]$post_varB) / nk
              ETA[[j]]$post_varB2 = (ETA[[j]]$varB^2 + (nk - 1) * ETA[[j]]$post_varB2) / nk

              # Mean of SB
              ETA[[j]]$post_SB = (ETA[[j]]$SB + (nk - 1) * ETA[[j]]$post_SB) / nk
              ETA[[j]]$post_SB2 = (ETA[[j]]$SB^2 + (nk - 1) * ETA[[j]]$post_SB2) / nk

              write(ETA[[j]]$SB,file = ETA[[j]]$fileOut,append = TRUE)
            }

            if (ETA[[j]]$model == "BL") {
              # Mean of Beta's
              ETA[[j]]$post_Bv = (ETA[[j]]$Bv + (nk - 1) * ETA[[j]]$post_Bv) / nk
              ETA[[j]]$post_Bv2 = (ETA[[j]]$Bv^2 + (nk - 1) * ETA[[j]]$post_Bv2) / nk

              # Mean of tau2
              ETA[[j]]$post_tau2 = (ETA[[j]]$tau2 + (nk - 1) * ETA[[j]]$post_tau2) / nk

              # Mean of lambda
              ETA[[j]]$post_lambda = (ETA[[j]]$lambda + (nk - 1) * ETA[[j]]$post_lambda) / nk

              write(ETA[[j]]$lambda,file = ETA[[j]]$fileOut,append = TRUE)
            }

            if (ETA[[j]]$model == "RKHS") {

              # Mean of U's
              ETA[[j]]$post_u = (ETA[[j]]$u + (nk - 1) * ETA[[j]]$post_u) / nk
              ETA[[j]]$post_u2 = (ETA[[j]]$u^2 + (nk - 1) * ETA[[j]]$post_u2) / nk

              # Mean of varU
              ETA[[j]]$post_varU = (ETA[[j]]$varU + (nk - 1) * ETA[[j]]$post_varU) / nk
              ETA[[j]]$post_varU2 = (ETA[[j]]$varU^2 + (nk - 1) * ETA[[j]]$post_varU2) / nk

              # Mean of Ustar's
              ETA[[j]]$post_uStar = (ETA[[j]]$uStar + (nk - 1) * ETA[[j]]$post_uStar) / nk

              write(ETA[[j]]$varU,file = ETA[[j]]$fileOut,append = TRUE)
            }

            if (ETA[[j]]$model=="RKHS_utme") {
              #Mean of Beta's
              ETA[[j]]$post_Bv = (ETA[[j]]$Bv + (nk - 1) * ETA[[j]]$post_Bv) / nk
              ETA[[j]]$post_Bv2 = (ETA[[j]]$Bv^2 + (nk - 1) * ETA[[j]]$post_Bv2) / nk

              # Mean of varB
              ETA[[j]]$Cov$post_SigmaB = (ETA[[j]]$Cov$SigmaB + (nk - 1) * ETA[[j]]$Cov$post_SigmaB) / nk
              ETA[[j]]$Cov$post_SigmaB2 = (ETA[[j]]$Cov$SigmaB^2 + (nk - 1) * ETA[[j]]$Cov$post_SigmaB2) / nk

              if (ETA[[j]]$Cov$type %in% c("UN","DIAG")) {
                tmp<-matrixcalc::vech(ETA[[j]]$Cov$SigmaB)
                write(tmp,ncolumns = length(tmp),file = ETA[[j]]$Cov$f_SigmaB,
                      append = TRUE,sep = " ")
                rm(tmp)
              }

              if (ETA[[j]]$Cov$type %in% c("FA","REC")) {
                ETA[[j]]$Cov$post_W<-(ETA[[j]]$Cov$W + (nk - 1) * ETA[[j]]$Cov$post_W) / nk
                ETA[[j]]$Cov$post_W2<-(ETA[[j]]$Cov$W^2 + (nk - 1) * ETA[[j]]$Cov$post_W2) / nk
                ETA[[j]]$Cov$post_PSI<-(ETA[[j]]$Cov$PSI + (nk - 1) * ETA[[j]]$Cov$post_PSI) / nk
                ETA[[j]]$Cov$post_PSI2<-(ETA[[j]]$Cov$PSI^2 + (nk - 1) * ETA[[j]]$Cov$post_PSI2) / nk

                if (sum(ETA[[j]]$Cov$M) > 0) {
                  tmp<-ETA[[j]]$Cov$W
                  write(tmp,ncolumns = length(tmp),file = ETA[[j]]$Cov$f_W,
                        append = TRUE,sep = " ")
                  rm(tmp)
                }
                write(ETA[[j]]$Cov$PSI,ncolumns = length(ETA[[j]]$Cov$PSI),
                      file = ETA[[j]]$Cov$f_PSI,append = TRUE,
                      sep = " ")
              }

              if(!is.null(ETA[[j]]$GxExT)){

                ETA[[j]]$GxExT$post_Bv = (ETA[[j]]$GxExT$Bv + (nk - 1) * ETA[[j]]$GxExT$post_Bv) / nk
                ETA[[j]]$GxExT$post_Bv2 = (ETA[[j]]$GxExT$Bv^2 + (nk - 1) * ETA[[j]]$GxExT$post_Bv2) / nk

                ETA[[j]]$GxExT$Cov$post_SigmaEnv = (ETA[[j]]$GxExT$Cov$SigmaEnv + (nk - 1) * ETA[[j]]$GxExT$Cov$post_SigmaEnv) / nk
                ETA[[j]]$GxExT$Cov$post_SigmaEnv2 = (ETA[[j]]$GxExT$Cov$SigmaEnv^2 + (nk - 1) * ETA[[j]]$GxExT$Cov$post_SigmaEnv2) / nk

                if (ETA[[j]]$GxExT$Cov$type %in% c("UN","DIAG")) {
                  tmp<-matrixcalc::vech(ETA[[j]]$GxExT$Cov$SigmaEnv)
                  write(tmp,ncolumns = length(tmp),file = ETA[[j]]$GxExT$Cov$f_SigmaEnv,
                        append = TRUE,sep = " ")
                  rm(tmp)
                }
                if (ETA[[j]]$GxExT$Cov$type %in% c("FA","REC")) {
                  ETA[[j]]$GxExT$Cov$post_W<-(ETA[[j]]$GxExT$Cov$W + (nk - 1) * ETA[[j]]$GxExT$Cov$post_W) / nk
                  ETA[[j]]$GxExT$Cov$post_W2<-(ETA[[j]]$GxExT$Cov$W^2 + (nk - 1) * ETA[[j]]$GxExT$Cov$post_W2) / nk
                  ETA[[j]]$GxExT$Cov$post_PSI<-(ETA[[j]]$GxExT$Cov$PSI + (nk - 1) * ETA[[j]]$GxExT$Cov$post_PSI) / nk
                  ETA[[j]]$GxExT$Cov$post_PSI2<-(ETA[[j]]$GxExT$Cov$PSI^2 + (nk - 1) * ETA[[j]]$GxExT$Cov$post_PSI2) / nk

                  if (sum(ETA[[j]]$GxExT$Cov$M) > 0) {
                    tmp<-ETA[[j]]$GxExT$Cov$W
                    write(tmp,ncolumns = length(tmp),file = ETA[[j]]$GxExT$Cov$f_W,
                          append = TRUE,sep = " ")
                    rm(tmp)
                  }
                  write(ETA[[j]]$GxExT$Cov$PSI,ncolumns = length(ETA[[j]]$GxExT$Cov$PSI),
                        file = ETA[[j]]$GxExT$Cov$f_PSI,append = TRUE,
                        sep = " ")
                }

              }
            }

          }
        }

        if(response_type %in% c("DLN","PLN","gaussian")){
          # Mean of varE
          post_varE = (varE + (nk - 1) * post_varE) / nk
          post_varE2 = (varE^2 + (nk - 1) * post_varE2) / nk
          write(varE,ncolumns = length(varE),file = f_varE,append = TRUE,
                sep = " ")
        }

        # Mean of loglik
        post_logLik = (loglik + (nk - 1) * post_logLik) / nk
        post_logLik2 = (loglik^2 + (nk - 1) * post_logLik2) / nk

        if(response_type %in% c("Poisson","PLN")){
          post_r = (r + (nk - 1) * post_r) / nk
          post_r2 = (r^2 + (nk - 1) * post_r2) / nk
        }

        # Mean predictions
        if(response_type == "DLN"){
          if(type_prediction == "mean"){
            y_pred = ( y_tr + (nk - 1) * y_pred ) / nk
            y_pred2 = ( y_tr^2 + (nk - 1) * y_pred2 ) / nk
          }
          rv_mean = ( rv + (nk - 1) * rv_mean ) / nk
        }

        if(response_type == "Poisson"){
          y_pred = ( y_tr + (nk - 1) * y_pred ) / nk
          y_pred2 = ( y_tr^2 + (nk - 1) * y_pred2 ) / nk
          rv_mean = ( rvPois + (nk - 1) * rv_mean ) / nk
        }

        if(response_type == "PLN"){
          y_pred = ( y_tr + (nk - 1) * y_pred ) / nk
          y_pred2 = ( y_tr^2 + (nk - 1) * y_pred2 ) / nk
          rv_mean = ( rvPois + (nk - 1) * rv_mean ) / nk
        }

        if(response_type == "gaussian"){
          y_pred = ( y_tr + (nk - 1) * y_pred ) / nk
          y_pred2 = ( y_tr^2 + (nk - 1) * y_pred2 ) / nk
          rv_mean = ( rv + (nk - 1) * rv_mean ) / nk
        }

      }
    }

    if (verbose) {
      pb$tick()
    }

  }
  tmp<-proc.time()[3]
  message("Time = ",round((tmp-time)/60,3)," minutes")
  cat('\n')
  #-----------------------------------------------------------------------------
  if (nNa > 0) {
    close(f_y_posterior); f_y_posterior<-NULL
  }

  if(intercept[[1]]){
    close(f_beta0); f_beta0<-NULL
  }
  close(f_logLik); f_logLik<-NULL
  if(response_type %in% c("DLN","PLN","gaussian")){
    close(f_varE); f_varE<-NULL
  }

  # Bayesian estimations
  #---------------------------------------------------------------------------
  # Predictions
  fit = list()

  if (nLT > 0) {
    for (i in 1:nLT) {
      if (!is.null(ETA[[i]]$fileOut)) {
        flush(ETA[[i]]$fileOut)
        close(ETA[[i]]$fileOut)
        ETA[[i]]$fileOut = NULL
      }
    }
  }

  if(response_type == "DLN"){
    fit$logLikAtPostMean = fllp_DLN(rv=rv_mean[obs_idx],varE=post_varE,y=y[obs_idx])
    fit$pD = -2 * (post_logLik - fit$logLikAtPostMean)
    fit$DIC = fit$pD - 2 * post_logLik
    if(type_prediction == "median"){
      if(nNa == 0){
        fit$y_train = as.vector(ceiling(exp(rv_mean)-1))
        names(fit$y_train) = rowNames
      }else{
        fit$y_train = as.vector(ceiling(exp(rv_mean[obs_idx])-1))
        names(fit$y_train) = rowNames[obs_idx]
        fit$y_test = as.vector(ceiling(exp(rv_mean[!obs_idx]) - 1))
        names(fit$y_test) = rowNames[!obs_idx]
      }
    }else{
      fit$y_train = as.vector(y_pred[obs_idx])
      fit$SD.y_train = as.vector(sqrt(y_pred2[obs_idx] - y_pred[obs_idx]^2))
      if(nNa == 0){
        names(fit$y_train) = rowNames
        names(fit$SD.y_train) = rowNames
      }else{
        names(fit$y_train) = rowNames[obs_idx]
        names(fit$SD.y_train) = rowNames[obs_idx]
        fit$y_test = as.vector(y_pred[!obs_idx])
        fit$SD.y_test = as.vector(sqrt(y_pred2[!obs_idx] - y_pred[!obs_idx]^2))
        names(fit$y_test) = rowNames[!obs_idx]
        names(fit$SD.y_test) = rowNames[!obs_idx]
      }
    }
    SD.varE = sqrt(post_varE2 - post_varE^2)
    fit$varE = post_varE
    fit$SD.varE = SD.varE
  }
  if(response_type == "Poisson"){
    fit$logLikAtPostMean = fllp_P(rv=rv_mean[obs_idx]-log(post_r[obs_idx]),
                                  y=y[obs_idx],r=post_r[obs_idx])
    fit$pD = -2 * (post_logLik - fit$logLikAtPostMean)
    fit$DIC = fit$pD - 2 * post_logLik
    fit$y_train = as.vector(y_pred[obs_idx])
    fit$SD.y_train = as.vector(sqrt(y_pred2[obs_idx] - y_pred[obs_idx]^2))
    if(nNa == 0){
      names(fit$y_train) = rowNames
      names(fit$SD.y_train) = rowNames
    }else{
      names(fit$y_train) = rowNames[obs_idx]
      names(fit$SD.y_train) = rowNames[obs_idx]
      fit$y_test = as.vector(y_pred[!obs_idx])
      fit$SD.y_test = as.vector(sqrt(y_pred2[!obs_idx] - y_pred[!obs_idx]^2))
      names(fit$y_test) = rowNames[!obs_idx]
      names(fit$SD.y_test) = rowNames[!obs_idx]
    }
    fit$r = post_r
    fit$SD.r = sqrt(post_r2 - post_r^2)
  }
  if(response_type == "PLN"){
    fit$logLikAtPostMean = .Call("fllp_PLN",y[obs_idx],logy_factorial,
                                 rv_mean[obs_idx],post_varE,gh$nodes,gh$weights,
                                 n-nNa,Q,PACKAGE="BGLM")
    fit$pD = -2 * (post_logLik - fit$logLikAtPostMean)
    fit$DIC = fit$pD - 2 * post_logLik
    fit$y_train = as.vector(y_pred[obs_idx])
    fit$SD.y_train = as.vector(sqrt(y_pred2[obs_idx] - y_pred[obs_idx]^2))
    if(nNa == 0){
      names(fit$y_train) = rowNames
      names(fit$SD.y_train) = rowNames
    }else{
      names(fit$y_train) = rowNames[obs_idx]
      names(fit$SD.y_train) = rowNames[obs_idx]
      fit$y_test = as.vector(y_pred[!obs_idx])
      fit$SD.y_test = as.vector(sqrt(y_pred2[!obs_idx] - y_pred[!obs_idx]^2))
      names(fit$y_test) = rowNames[!obs_idx]
      names(fit$SD.y_test) = rowNames[!obs_idx]
    }
    SD.varE = sqrt(post_varE2 - post_varE^2)
    fit$varE = post_varE
    fit$SD.varE = SD.varE
    fit$r = post_r
    fit$SD.r = sqrt(post_r2 - post_r^2)
  }
  if(response_type == "gaussian"){
    fit$logLikAtPostMean = fllp_G(rv=rv_mean[obs_idx],varE=post_varE,y=y[obs_idx])
    fit$pD = -2 * (post_logLik - fit$logLikAtPostMean)
    fit$DIC = fit$pD - 2 * post_logLik
    fit$y_train = as.vector(y_pred[obs_idx])
    fit$SD.y_train = as.vector(sqrt(y_pred2[obs_idx] - y_pred[obs_idx]^2))
    if(nNa == 0){
      names(fit$y_train) = rowNames
      names(fit$SD.y_train) = rowNames
    }else{
      names(fit$y_train) = rowNames[obs_idx]
      names(fit$SD.y_train) = rowNames[obs_idx]
      fit$y_test = as.vector(y_pred[!obs_idx])
      fit$SD.y_test = as.vector(sqrt(y_pred2[!obs_idx] - y_pred[!obs_idx]^2))
      names(fit$y_test) = rowNames[!obs_idx]
      names(fit$SD.y_test) = rowNames[!obs_idx]
    }
    SD.varE = sqrt(post_varE2 - post_varE^2)
    fit$varE = post_varE
    fit$SD.varE = SD.varE
  }

  if(intercept[[1]]){
    fit$beta0 = post_beta0
    SD.beta0 = sqrt(post_beta02 - post_beta0^2)
    fit$SD.beta0 = SD.beta0
  }

  fit$post_logLik = post_logLik
  SD.logLik = sqrt(post_logLik2 - post_logLik^2)
  fit$SD.postlogLik = SD.logLik

  if (nLT > 0) {
    for (i in 1:nLT) {
      if (ETA[[i]]$model %in% c("FIXED","BRR","BayesA","BL")) {
        ETA[[i]]$Bv = as.vector(ETA[[i]]$post_Bv)
        ETA[[i]]$SD.Bv = as.vector(sqrt(ETA[[i]]$post_Bv2 - ETA[[i]]$post_Bv^2))
        names(ETA[[i]]$Bv) = ETA[[i]]$colNames
        names(ETA[[i]]$SD.Bv) = ETA[[i]]$colNames
        if(response_type %in% c("DLN","gaussian")){
          tmp = which(names(ETA[[i]]) %in% c("post_Bv","post_Bv2","X","x2","tX",
                                             "tXX","nblocks"))
          ETA[[i]] = ETA[[i]][-tmp]
        }
        if(response_type %in% c("Poisson","PLN")){
          tmp = which(names(ETA[[i]]) %in% c("post_Bv","post_Bv2","X","X2","tX",
                                             "tXX","nblocks"))
          ETA[[i]] = ETA[[i]][-tmp]
        }
      }
      if (ETA[[i]]$model %in% c("BRR","BayesA")) {
        ETA[[i]]$SD.varB = sqrt(ETA[[i]]$post_varB2 - (ETA[[i]]$post_varB^2))
        tmp = which(names(ETA[[i]]) %in% c("post_varB","post_varB2"))
        ETA[[i]] = ETA[[i]][-tmp]
      }
      if (ETA[[i]]$model == c("BayesA")) {
        ETA[[i]]$SB = ETA[[i]]$post_SB
        ETA[[i]]$SD.SB = sqrt(ETA[[i]]$post_SB2 - (ETA[[i]]$post_SB^2))
        tmp = which(names(ETA[[i]]) %in% c("post_SB","post_SB2"))
        ETA[[i]] = ETA[[i]][-tmp]
      }
      if (ETA[[i]]$model == "BL") {
        ETA[[i]]$tau2 = ETA[[i]]$post_tau2
        ETA[[i]]$lambda = ETA[[i]]$post_lambda
        tmp = which(names(ETA[[i]]) %in% c("post_tau2","post_lambda","lambda2"))
        ETA[[i]] = ETA[[i]][-tmp]
      }
      if (ETA[[i]]$model == "RKHS") {
        ETA[[i]]$SD.u = sqrt(ETA[[i]]$post_u2 - ETA[[i]]$post_u^2)
        ETA[[i]]$u = ETA[[i]]$post_u
        ETA[[i]]$uStar = ETA[[i]]$post_uStar
        ETA[[i]]$varU = ETA[[i]]$post_varU
        ETA[[i]]$SD.varU = sqrt(ETA[[i]]$post_varU2 - ETA[[i]]$post_varU^2)
        tmp = which(names(ETA[[i]]) %in% c("post_varU","post_varU2","post_uStar","post_u","post_u2",
                                           "V","V2","v2","tV","tVV","nblocks"))
        ETA[[i]] = ETA[[i]][-tmp]
      }
      if (ETA[[i]]$model=="RKHS_utme") {
        if(!is.null(ETA[[i]]$GxT)){
          ETA[[i]]$u = ETA[[i]]$X%*%ETA[[i]]$post_Bv
          ETA[[i]]$uStar = ETA[[i]]$post_Bv
          ETA[[i]]$SD.uStar = sqrt(ETA[[i]]$post_Bv2 - ETA[[i]]$post_Bv^2)
        }
        if(!(response_type%in%c("Poisson","PLN"))){
          tmp = which(names(ETA[[i]]) %in% c("Bv","post_Bv","post_Bv2","X","x2",
                                             "X_train","p","eigenvals_GT",
                                             "eigenvecs_GT","eigenvals_EG",
                                             "eigenvecs_EG","Aux","nEnv",
                                             "nLines","nLinesTraits","nLinesEnvs",
                                             "SE","U1","SE1","U2"))
          ETA[[i]] = ETA[[i]][-tmp]
        }
        else{
          tmp = which(names(ETA[[i]]) %in% c("Bv","post_Bv","post_Bv2","X","X2",
                                             "X_train","p","eigenvals_GT",
                                             "eigenvecs_GT","eigenvals_EG",
                                             "eigenvecs_EG","Aux","nEnv",
                                             "nLines","nLinesTraits","nLinesEnvs",
                                             "SE","U1","SE1","U2"))
          ETA[[i]] = ETA[[i]][-tmp]
        }
        if(!is.null(ETA[[i]]$GxExT)){
          ETA[[i]]$GxExT$u = ETA[[i]]$GxExT$X%*%ETA[[i]]$GxExT$post_Bv
          ETA[[i]]$GxExT$uStar = ETA[[i]]$GxExT$post_Bv
          ETA[[i]]$GxExT$SD.uStar = sqrt(ETA[[i]]$GxExT$post_Bv2 - ETA[[i]]$GxExT$post_Bv^2)
          if(!(response_type%in%c("Poisson","PLN"))){
            tmp = which(names(ETA[[i]]$GxExT) %in% c("Bv","post_Bv","post_Bv2","X","x2",
                                                     "X_train","p","eigenvals_GT",
                                                     "eigenvecs_GT","eigenvals_EG",
                                                     "eigenvecs_EG","Aux","nEnv",
                                                     "nLines","nLinesTraits","nLinesEnvs"))
            ETA[[i]]$GxExT = ETA[[i]]$GxExT[-tmp]
          }
          else{
            tmp = which(names(ETA[[i]]$GxExT) %in% c("Bv","post_Bv","post_Bv2","X","X2",
                                                     "X_train","p","eigenvals_GT",
                                                     "eigenvecs_GT","eigenvals_EG",
                                                     "eigenvecs_EG","Aux","nEnv",
                                                     "nLines","nLinesTraits","nLinesEnvs"))
            ETA[[i]]$GxExT = ETA[[i]]$GxExT[-tmp]
          }
        }
      }
      if (ETA[[i]]$model =="RKHS_utme") {
        ETA[[i]]$Cov$SigmaB = ETA[[i]]$Cov$post_SigmaB
        ETA[[i]]$Cov$SD.SigmaB = sqrt(ETA[[i]]$Cov$post_SigmaB2 - (ETA[[i]]$Cov$post_SigmaB^2))

        if (ETA[[i]]$Cov$type %in% c("UN","DIAG")) {
          close(ETA[[i]]$Cov$f_SigmaB)
          ETA[[i]]$Cov$f_SigmaB<-NULL
        }
        if (ETA[[i]]$Cov$type %in% c("FA","REC")) {
          ETA[[i]]$Cov$W<-ETA[[i]]$Cov$post_W
          ETA[[i]]$Cov$SD.W<-sqrt(ETA[[i]]$Cov$post_W2 - ETA[[i]]$Cov$post_W^2)
          ETA[[i]]$Cov$PSI<-ETA[[i]]$Cov$post_PSI
          ETA[[i]]$Cov$SD.PSI<-sqrt(ETA[[i]]$Cov$post_PSI2 - ETA[[i]]$Cov$post_PSI^2)

          close(ETA[[i]]$Cov$f_W)
          ETA[[i]]$Cov$f_W<-NULL
          close(ETA[[i]]$Cov$f_PSI)
          ETA[[i]]$Cov$f_PSI<-NULL
        }

        tmp<-which(names(ETA[[i]]$Cov) %in% c("post_SigmaB","post_SigmaB2","SigmaBInv",
                                                "post_W","post_W2","post_PSI","post_PSI2","F"))
        ETA[[i]]$Cov<-ETA[[i]]$Cov[-tmp]
        rm(tmp)

        if(!is.null(ETA[[j]]$GxExT)){

          ETA[[i]]$GxExT$Cov$SigmaEnv = ETA[[i]]$GxExT$Cov$post_SigmaEnv
          ETA[[i]]$GxExT$Cov$SD.SigmaEnv = sqrt(ETA[[i]]$GxExT$Cov$post_SigmaEnv2 - (ETA[[i]]$GxExT$Cov$post_SigmaEnv^2))

          if (ETA[[i]]$GxExT$Cov$type %in% c("UN","DIAG")) {
            close(ETA[[i]]$GxExT$Cov$f_SigmaEnv)
            ETA[[i]]$GxExT$Cov$f_SigmaEnv<-NULL
          }
          if (ETA[[i]]$GxExT$Cov$type %in% c("FA","REC")) {
            ETA[[i]]$GxExT$Cov$W<-ETA[[i]]$GxExT$Cov$post_W
            ETA[[i]]$GxExT$Cov$SD.W<-sqrt(ETA[[i]]$GxExT$Cov$post_W2 - ETA[[i]]$GxExT$Cov$post_W^2)
            ETA[[i]]$GxExT$Cov$PSI<-ETA[[i]]$GxExT$Cov$post_PSI
            ETA[[i]]$GxExT$Cov$SD.PSI<-sqrt(ETA[[i]]$GxExT$Cov$post_PSI2 - ETA[[i]]$GxExT$Cov$post_PSI^2)

            close(ETA[[i]]$GxExT$Cov$f_W)
            ETA[[i]]$GxExT$Cov$f_W<-NULL
            close(ETA[[i]]$GxExT$Cov$f_PSI)
            ETA[[i]]$GxExT$Cov$f_PSI<-NULL
          }

          tmp<-which(names(ETA[[i]]$GxExT$Cov) %in% c("post_SigmaEnv","post_SigmaEnv2","SigmaEnvInv",
                                                        "post_W","post_W2","post_PSI","post_PSI2","F"))
          ETA[[i]]$GxExT$Cov<-ETA[[i]]$GxExT$Cov[-tmp]
          rm(tmp)
        }
      }
    }
  }

  out<-list(y=y,whichNa = whichNa,fit = fit,nIter = nIter,nBurnin = nBurnin,
              nThin = nThin,priorvarbeta0 = priorvarbeta0,R2 = R2,priorv = priorv,
              priorS = priorS,response_type=response_type)

  out$ETA<-ETA
  class(out) = "BGLM"
  return(out)

}

#' Fit a Bayesian Multi-Trait Multi-Environment Regression Model
#'
#' Fits a Bayesian Multi-trait and Multi-environment model using MCMC sampling.
#' Supports Gaussian, Poisson, DLN and PLN responses.
#'
#' @param ETA List used to specify the linear predictor.
#' By default it is set to NULL, in which case only the intercept is included.
#' @param y Response vector. Can contain NA values.
#' @param response_type Type of response ("DLN","gaussian","Poisson","PLN").
#' For DLN, Poisson or PLN responses, the argument `y` must contain integer values,
#' while for gaussian responses any type of numeric value can be accepted.
#' @param nIter Number of MCMC iterations.
#' @param nBurnin Burn-in iterations.
#' @param nThin Thinning interval.
#' @param R2 Proportion of variance explained a priori for the error term.
#' @param priorv Prior degrees of freedom for the error co-variance matrix.
#' @param priorS Prior scale matrix for the error co-variance matrix.
#' @param type Co-variance structure (Unstructured, Diagonal, Factor Analysis
#' or Recursive)
#' @param priorSigma0 Choose the variance, σ_0^2, to obtain a flat prior
#' for the vector of intercepts, β_0∼N(0, Iσ_0^2).By default (internally) is 1e10.
#' @param rcontrol Specify the k-value used to control the difference between the
#' Conditional Mean and Conditional Variance of the Negative Binomial model.
#' Higher values of k indicate a smaller difference between the mean and variance,
#' approximating the equidispersion of a Poisson model
#' @param type_prediction Used to specify if you want to predict with the posterior
#' mean or median in DLN model; default is mean.
#' @param Iters_latent Number of MCMC iterations to obtain Multivariate Truncated
#' Normal samples.
#' @param Burn_latent Burn-in iterations for the Multivariate Truncated Normal
#' samples.
#' @param Thin_latent Thinning interval for the Multivariate Truncated Normal
#' samples.
#' @param n_QMCpoints Number of Sobol points for the Quasi-Monte Carlo method
#' used for the approximation of the log-likelihood of the Multivariate DLN model.
#' @param intercept_mt Used to include (TRUE) or not (FALSE) the intercept in
#' the linear predictor; default is TRUE.
#' @param verbose Used to show (TRUE) or not (FALSE) the progress of iterations
#' in the Gibbs Sampler; default is FALSE
#' @param saveAt Used to indicate UTME where to store the samples and to provide
#' a pre-fix to be appended to the names of the file where samples are stored.
#' By default samples are saved in the current working directory and no pre-fix
#' is added to the file names
#' @return A list containing posterior samples and summaries.
#'
#' @export
MTME<-function(
    ETA=NULL,
    y=NULL,
    response_type="DLN",
    nIter=1e3,
    nBurnin=1e2,
    nThin=1,
    R2=0.5,
    resCov=list(priorv=NULL,priorS=NULL,type="UN"),
    priorSigma0=NULL,
    rcontrol=1.15,
    type_prediction="mean",
    Iters_latent=1,
    Burn_latent=0,
    Thin_latent=1,
    n_QMCpoints=10,
    intercept_mt=T,
    verbose=F,
    saveAt=""){

  # nIter and nBurnin validation
  if (!is.null(nIter) && !is.null(nBurnin)) {
    if (nIter <= nBurnin) {
      stop("\033[31m\nError: nIter should be greater than nBurnin by at least 1.\033[39m")
    }
  }

  # Response type validation
  if (!(response_type %in% c("DLN","Poisson","PLN","gaussian"))){
    stop("Only DLN,Poisson,PLN,or gaussian responses are allowed
           (note: evaluation is case sensitive)")
  }

  # Type prediction validation
  if (response_type == "DLN" && !(type_prediction %in% c("mean","median"))){
    stop("For the DLN model,only mean and median type predictions are allowed
           (note: evaluation is case sensitive)")
  }

  if (saveAt == "") {
    saveAt = paste(getwd(),"/",sep = "")
  }

  if(verbose){
    pb<-progress::progress_bar$new(
      format = "  Progress [:bar] :percent in :elapsed",
      total = nIter,clear = FALSE,width = 60
    )
  }

  if (response_type == "DLN"){
    if (Iters_latent <= Burn_latent) {
      stop("\033[31m\nError: To use the Gibbs Sampler for latent variables,
           Iters_latent should be greater than Burn_latent by at least 1.\033[39m")
    }
  }

  if (!is.matrix(y)){ stop("y must be a matrix\n") }
  ntraits<-ncol(y)

  if (ntraits < 2) { stop("y must have at least 2 columns\n") }

  n = nrow(y)
  rowNames = rownames(y); colNames = colnames(y)
  NoWhichNa = complete.cases(y)
  nNa = length(which(!NoWhichNa))
  nObs = n - nNa

  if(!(response_type == "DLN" && type_prediction == "median")){
    y_pred = matrix(0,nrow = n,ncol = ntraits)
    y_pred2 = matrix(0,nrow = n,ncol = ntraits)
  }

  if(response_type == "gaussian"){
    Sy<-cov(y,use = "pairwise.complete.obs")
  }
  if(response_type == "DLN"){
    Sy<-cov(log(y+1),use = "pairwise.complete.obs")
  }
  if(response_type%in%c("Poisson","PLN")){
    aux<-colMeans(log(y+1),na.rm=TRUE)
    Sy<-outer(sqrt(aux),sqrt(aux))
    if(!matrixcalc::is.positive.definite(Sy)){
      eigens = eigen(Sy)$values
      MinEigen = min(eigens)
      alpha = -MinEigen + 1e-6
      diag(Sy)<-diag(Sy) + alpha
    }
  }

  if (intercept_mt && is.null(priorSigma0)) {
    priorSigma0<-diag(ntraits)*1e+10; priorSigma0Inv<-diag(1/diag(priorSigma0))
  }

  #-----------------------------------------------------------------------------

  nLT<-ifelse(is.null(ETA),0,length(ETA))
  nLT_total<-sum(sapply(ETA,function(x){
    sublists<-x[sapply(x,is.list)]
    sublists<-sublists[names(sublists)!="Cov"]
    if(length(sublists)==0) 1 else length(sublists)
  }))

  if(nLT > 0){
    if (is.null(names(ETA))) {
      names(ETA)<-rep("",nLT)
    }

    for (i in 1:nLT) {

      if (names(ETA)[i] == "") {
        ETA[[i]]$Name = paste("ETA_",i,sep = "")
      }
      else {
        ETA[[i]]$Name = paste("ETA_",names(ETA)[i],
                              sep = "")
      }
      # Model validation
      if (!(ETA[[i]]$model %in% c("FIXED","BRR","RKHS","RKHS_mtme"))) {
        stop("Error in ETA[[",i,"]]"," model ",ETA[[i]]$model,
             " not implemented (note: evaluation is case sensitive)")
      }

      if(!(ETA[[i]]$model%in%c("RKHS","RKHS_mtme"))){
        if(nNa != 0){
          ETA[[i]]$X_test = as.matrix(ETA[[i]]$X[!NoWhichNa,])
        }
      }

      ETA[[i]] = switch(ETA[[i]]$model,FIXED=setLT.Fixed_mt(LT=ETA[[i]],n=n,
                                                            ntraits=ntraits,i=i,saveAt=saveAt,response_type=response_type,NoWhichNa=NoWhichNa),
                        BRR=setLT.BRR_mt(LT=ETA[[i]],n=n,ntraits=ntraits,i=i,
                                         Sy=Sy,nLT=nLT_total,R2=R2,saveAt=saveAt,response_type=response_type,NoWhichNa=NoWhichNa),
                        RKHS=setLT.RKHS_mt(LT=ETA[[i]],n=n,ntraits=ntraits,i=i,Sy=Sy,
                                           nLT=nLT_total,R2=R2,saveAt=saveAt,response_type=response_type,NoWhichNa=NoWhichNa),
                        RKHS_mtme=setLT.RKHS_mtme(LT=ETA[[i]],n=n,ntraits=ntraits,i=i,y=y,Sy=Sy,
                                                  nLT=nLT_total,R2=R2,saveAt=saveAt,response_type=response_type,NoWhichNa=NoWhichNa))
    }
  }

  rv_mean = matrix(0,n,ncol = ntraits)
  #-----------------------------------------------------------------------------

  # Priors
  #-----------------------------------------------------------------------------

  if(!(response_type == "Poisson")){
    resCov<-setResCov(n=n,resCov=resCov,ntraits=ntraits,Sy=Sy,R2=R2,saveAt=saveAt)
    if(response_type == "PLN"){
      resCov<-setResCov(n=n,resCov=resCov,ntraits=ntraits,Sy=2*Sy,R2=R2,saveAt=saveAt)
    }
  }

  #-----------------------------------------------------------------------------

  # Initial values
  #-----------------------------------------------------------------------------
  if(intercept_mt){
    beta0 = rep(0,ntraits)
    post_beta0 = rep(0,ntraits)
    post_beta02 = rep(0,ntraits)
    f_beta0<-file(description = paste(saveAt,"beta0.dat",sep = ""),
                    open = "w")
  }
  post_logLik = 0
  post_logLik2 = 0
  f_logLik<-file(description = paste(saveAt,"logLik.dat",sep = ""),
                   open = "w")

  nk = 0
  #-----------------------------------------------------------------------------

  # Algorithm
  #-----------------------------------------------------------------------------
  time<-proc.time()[3]

  yStar<-matrix(NA,nrow=n,ncol=ntraits)
  yStar<-y[]
  if (nNa > 0) {
    for(k in 1:ntraits)
    {
      yStar[!NoWhichNa,k] = 0
    }
  }

  if(response_type == "DLN"){
    ay = log(yStar);  by = log(yStar+1)
    if (!(matrixcalc::is.diagonal.matrix(resCov$SigmaE))){
      U_qmc = randtoolbox::sobol(n = n_QMCpoints,dim = ntraits,seed=123)
    }
    rv = matrix(0,nrow = n,ncol = ntraits)
  }
  if(response_type == "gaussian"){
    cte = - (nObs * ntraits)/2 * log(2 * pi)
    rv = yStar - matrix(0,nrow = n,ncol = ntraits)
  }
  if(response_type == "Poisson"){
    r = (yStar+1) * 10^(rcontrol)
    y_r = yStar+r
    yr = (yStar - r) / 2
    Syr = colSums(yStar-r)
    post_r = matrix(0,nrow = n,ncol = ntraits)
    post_r2 = matrix(0,nrow = n,ncol = ntraits)
    rv = matrix(0,nrow = n,ncol = ntraits) - log(r)
  }
  if(response_type == "PLN"){
    gh<-statmod::gauss.quad(20,kind = "hermite")
    Q = length(gh$nodes)
    logy_factorial = lgamma(y[NoWhichNa,]+1)
    if(!(resCov$type == "DIAG")){
      logy_factorial<-rowSums(logy_factorial)

      Q<-n_QMCpoints
      # Uniform QMC nodes in [0,1]^ntraits
      sobol_nodes<-randtoolbox::sobol(n = Q,dim = ntraits,seed = 123)

      # Transform to normal standard nodes using the inverse (N(0,1))
      nodes_mult<-qnorm(sobol_nodes)/sqrt(2)

      # Equal weights QMC
      gh_weights_mult<-rep(1 / Q,Q)
    }
    U = matrix(0,nrow = n,ncol = ntraits);
    r = ((1/colMeans(yStar))+1) * 10^(rcontrol)
    y_r = sweep(x=yStar,MARGIN=2,STATS=r,FUN="+")
    yminusr = sweep(x=yStar,MARGIN=2,STATS=r,FUN="-")
    yr = (yminusr) / 2
    Syr = colSums(yminusr)
    post_r = rep(0,ntraits)
    post_r2 = rep(0,ntraits)
    rv = sweep(x=U,MARGIN=2,STATS=log(r),FUN="-")
  }

  fL_fun<-switch(response_type,
                   DLN = function(...) lmulti(n = n,ntraits = ntraits,rv = rv,SigmaE = resCov$SigmaE,
                                              SigmaEInv = resCov$SigmaEInv,ay = ay,
                                              by = by,type = resCov$type,Iters_latent = Iters_latent,
                                              Burn_latent = Burn_latent,Thin_latent = Thin_latent),
                   Poisson = function(...) wmulti(y_r = y_r,n = n,ntraits = ntraits,rv = rv),
                   PLN = function(...) wmulti(y_r = y_r,n = n,ntraits = ntraits,rv = rv)
  )

  fB0_fun<-switch(response_type,
                    DLN = function(...) fB0_DLN_multi(priorSigma0Inv=priorSigma0Inv,
                                                      SigmaEInv=resCov$SigmaEInv,rv=rv,n=n),
                    Poisson = function(...) fB0_P_multi(l=l,Syr=Syr,priorSigma0Inv=priorSigma0Inv,rv=rv),
                    PLN = function(...) fB0_P_multi(l=l,Syr=Syr,priorSigma0Inv=priorSigma0Inv,rv=rv),
                    gaussian = function(...) fB0_G_multi(priorSigma0Inv=priorSigma0Inv,
                                                         SigmaEInv=resCov$SigmaEInv,rv=rv,n=n)
  )

  # GIBBS SAMPLER
  #*****************************************************************************

  message("Generating ",nIter," samples,discarting ",nBurnin,',and using a thinning of ',nThin,
          " for the ",response_type," model")


  for(i in 1:nIter){

    if(!is.null(fL_fun)) l<-fL_fun()

    if(response_type=="DLN"){rv = l - rv}

    if(intercept_mt){
      ifelse(response_type%in%c("DLN","gaussian"),rv<-sweep(x=rv,MARGIN=2,STATS=beta0,FUN="+"),
             rv<-sweep(x=rv,MARGIN=2,STATS=beta0,FUN="-"))

      beta0<-fB0_fun()

      ifelse(response_type%in%c("DLN","gaussian"),rv<-sweep(x=rv,MARGIN=2,STATS=beta0,FUN="-"),
             rv<-sweep(x=rv,MARGIN=2,STATS=beta0,FUN="+"))
    }

    if(nLT > 0){
      for (j in 1:nLT) {

        if (ETA[[j]]$model == "FIXED") {
          if(response_type%in%c("DLN","gaussian")){
            # Sampling from full conditional of Bj's
            B = .Call("fBj_DLN_G_mtme",ETA[[j]]$Bv,ETA[[j]]$Cov$SigmaBInv,
                      resCov$SigmaEInv,rv,ETA[[j]]$X,ETA[[j]]$x2,
                      n,ETA[[j]]$p,ntraits,PACKAGE="BGLM")
          }
          if(response_type%in%c("Poisson","PLN")){
            # Sampling from full conditional of Bj's
            B = .Call("fBj_P_mtme",l,yr,ETA[[j]]$Bv,ETA[[j]]$Cov$SigmaBInv,rv,
                      ETA[[j]]$X,ETA[[j]]$X2,ETA[[j]]$p,ntraits,n,PACKAGE="BGLM")
          }
          ETA[[j]]$Bv = B[[1]]
          rv = B[[2]]
        }

        if (ETA[[j]]$model%in%c("BRR","RKHS")) {
          if(response_type%in%c("DLN","gaussian")){
            # Sampling from full conditional of Bj's
            B = .Call("fBj_DLN_G_mtme",ETA[[j]]$Bv,ETA[[j]]$Cov$SigmaBInv,
                      resCov$SigmaEInv,rv,ETA[[j]]$X,ETA[[j]]$x2,n,ETA[[j]]$p,
                      ntraits,PACKAGE="BGLM")
          }
          if(response_type%in%c("Poisson","PLN")){
            # Sampling from full conditional of Bj's
            B = .Call("fBj_P_mtme",l,yr,ETA[[j]]$Bv,ETA[[j]]$Cov$SigmaBInv,rv,
                      ETA[[j]]$X,ETA[[j]]$X2,ETA[[j]]$p,ntraits,n,PACKAGE="BGLM")
          }

          ETA[[j]]$Bv = B[[1]]
          S4<-crossprod(ETA[[j]]$Bv)

          if (ETA[[j]]$Cov$type == "UN") {
            ETA[[j]]$Cov$SigmaB<-MCMCpack::riwish(v = ETA[[j]]$Cov$priorvB + ETA[[j]]$p,
                                          S = S4 + ETA[[j]]$Cov$priorSB)
          }

          if (ETA[[j]]$Cov$type == "DIAG") {
            VVB = ETA[[j]]$Cov$priorvB+ETA[[j]]$p
            SSB = ETA[[j]]$Cov$priorSB+diag(S4)
            ETA[[j]]$Cov$SigmaB<-diag(1/rgamma(ntraits,VVB/2,SSB/2))
          }

          if (ETA[[j]]$Cov$type == "FA") {
            tmp<-sample_G0_FA(U = ETA[[j]]$Bv,F = ETA[[j]]$Cov$F,
                                M = ETA[[j]]$Cov$M,B = ETA[[j]]$Cov$W,PSI = ETA[[j]]$Cov$PSI,
                                ntraits = ntraits,nF = ETA[[j]]$Cov$nF,nD = ETA[[j]]$Cov$nD,
                                df0 = ETA[[j]]$Cov$priorvB,S0 = ETA[[j]]$Cov$priorSB,
                                priorVar = ETA[[j]]$Cov$var,varimaxRotate = ETA[[j]]$Cov$varimax)
            ETA[[j]]$Cov$F<-tmp$F
            ETA[[j]]$Cov$PSI<-tmp$PSI
            ETA[[j]]$Cov$W<-tmp$B
            ETA[[j]]$Cov$SigmaB<-tmp$G
            rm(tmp)
          }
          if (ETA[[j]]$Cov$type == "REC") {
            tmp<-sample_G0_REC(U = ETA[[j]]$Bv,M = ETA[[j]]$Cov$M,
                                 PSI = ETA[[j]]$Cov$PSI,ntraits = ntraits,
                                 priorVar = ETA[[j]]$Cov$var,df0 = ETA[[j]]$Cov$priorvB,
                                 S0 = ETA[[j]]$Cov$priorSB)
            ETA[[j]]$Cov$SigmaB<-tmp$G
            ETA[[j]]$Cov$W<-tmp$B
            ETA[[j]]$Cov$PSI<-tmp$PSI
            rm(tmp)
          }
          if(ETA[[j]]$Cov$type=="DIAG"){
            ETA[[j]]$Cov$SigmaBInv<-diag(1/diag(ETA[[j]]$Cov$SigmaB))
          }
          else{
            ETA[[j]]$Cov$SigmaBInv<-solve(ETA[[j]]$Cov$SigmaB)
          }
          rv = B[[2]]
        }

        if (ETA[[j]]$model=="RKHS_mtme") {

          if(!is.null(ETA[[j]]$GxT)){

            #Sampling from posterior of Traits parameters
            if(response_type%in%c("DLN","gaussian")){
              # Sampling from full conditional of Bj's
              B1 = .Call("fBj_DLN_G_mtme",ETA[[j]]$Bv,ETA[[j]]$Cov$SigmaBInv,
                         resCov$SigmaEInv,rv,ETA[[j]]$X,ETA[[j]]$x2,n,
                         ETA[[j]]$p,ntraits,PACKAGE="BGLM")
            }
            if(response_type%in%c("Poisson","PLN")){
              # Sampling from full conditional of Bj's
              B1 = .Call("fBj_P_mtme",l,yr,ETA[[j]]$Bv,ETA[[j]]$Cov$SigmaBInv,
                         rv,ETA[[j]]$X,ETA[[j]]$X2,ETA[[j]]$p,ntraits,n,
                         PACKAGE="BGLM")
            }

            ETA[[j]]$Bv = B1[[1]]
            rv = B1[[2]]

          }

          if(!is.null(ETA[[j]]$GxT) && !is.null(ETA[[j]]$GxExT)){
            if(ETA[[j]]$Cov$type %in% c("UN","DIAG")) {
              SE = crossprod(ETA[[j]]$Bv); SE1 = crossprod(ETA[[j]]$GxExT$Bv)
            }

            if(ETA[[j]]$Cov$type %in% c("FA","REC")) {
              U1 = ETA[[j]]$Bv; U2 = ETA[[j]]$GxExT$Bv
            }

            if(ETA[[j]]$Cov$type %in% c("FA")) {
              F1 = ETA[[j]]$Cov$F; F2 = ETA[[j]]$F2; const1=1
              W1 = ETA[[j]]$Cov$W; W2 = W1; const2=1
            }

          }
          if(!is.null(ETA[[j]]$GxT) && is.null(ETA[[j]]$GxExT)){
            if(ETA[[j]]$Cov$type %in% c("UN","DIAG")) {
              SE = crossprod(ETA[[j]]$Bv); SE1 = ETA[[j]]$SE1
            }

            if(ETA[[j]]$Cov$type %in% c("FA","REC")) {
              U1 = ETA[[j]]$Bv; U2 = ETA[[j]]$U2
            }

            if(ETA[[j]]$Cov$type %in% c("FA")) {
              F1 = ETA[[j]]$Cov$F; F2 = ETA[[j]]$F2; const1=1
              W1 = ETA[[j]]$Cov$W; W2 = ETA[[j]]$W2; const2=Inf
            }

          }
          if(is.null(ETA[[j]]$GxT) && !is.null(ETA[[j]]$GxExT)){
            if(ETA[[j]]$Cov$type %in% c("UN","DIAG")) {
              SE = ETA[[j]]$SE; SE1 = crossprod(ETA[[j]]$GxExT$Bv)
            }

            if(ETA[[j]]$Cov$type %in% c("FA","REC")) {
              U1 = ETA[[j]]$U1; U2 = ETA[[j]]$GxExT$Bv
            }

            if(ETA[[j]]$Cov$type %in% c("FA")) {
              F1 = ETA[[j]]$Cov$F; F2 = ETA[[j]]$F2; const1=Inf
              W2 = ETA[[j]]$W2; W1 = W2; const2=1
            }

          }

          if (ETA[[j]]$Cov$type == "UN") {
            ETA[[j]]$Cov$SigmaB<-MCMCpack::riwish(v = ETA[[j]]$Cov$priorvB+ETA[[j]]$nLines+ETA[[j]]$nLinesEnvs,
                                          S = SE + SE1 + ETA[[j]]$Cov$priorSB)
          }
          if (ETA[[j]]$Cov$type == "DIAG") {
            VVB = ETA[[j]]$Cov$priorvB+ETA[[j]]$nLines+ETA[[j]]$nLinesEnvs
            SSB = ETA[[j]]$Cov$priorSB+diag(SE)+diag(SE1)
            ETA[[j]]$Cov$SigmaB<-diag(1/rgamma(ntraits,VVB/2,SSB/2))
          }
          if (ETA[[j]]$Cov$type == "FA") {
            tmp<-sample_G0_FA_MTME(U1 = U1,U2 = U2,F1 = F1,F2 = F2,
                                     M = ETA[[j]]$Cov$M,B = W1,B2 = W2,const1=const1,
                                     const2=const2,PSI = ETA[[j]]$Cov$PSI,
                                     ntraits = ntraits,nF = ETA[[j]]$Cov$nF,nD1 = ETA[[j]]$nLines,
                                     nD2 = ETA[[j]]$nLinesEnvs,
                                     df0 = ETA[[j]]$Cov$priorvB,S0 = ETA[[j]]$Cov$priorSB,
                                     priorVar = ETA[[j]]$Cov$var,varimaxRotate = ETA[[j]]$Cov$varimax)
            ETA[[j]]$Cov$F<-tmp$F1
            ETA[[j]]$F2<-tmp$F2
            ETA[[j]]$Cov$PSI<-tmp$PSI
            ETA[[j]]$Cov$W<-tmp$B
            ETA[[j]]$Cov$SigmaB<-tmp$G
            rm(tmp)
          }
          if (ETA[[j]]$Cov$type == "REC") {
            tmp<-sample_G0_REC_MTME(U1 = U1,U2 = U2,M = ETA[[j]]$Cov$M,
                                      PSI = ETA[[j]]$Cov$PSI,ntraits = ntraits,
                                      priorVar = ETA[[j]]$Cov$var,df0 = ETA[[j]]$Cov$priorvB,
                                      S0 = ETA[[j]]$Cov$priorSB)
            ETA[[j]]$Cov$SigmaB<-tmp$G
            ETA[[j]]$Cov$W<-tmp$B
            ETA[[j]]$Cov$PSI<-tmp$PSI
            rm(tmp)
          }
          if(ETA[[j]]$Cov$type=="DIAG"){
            ETA[[j]]$Cov$SigmaBInv<-diag(1/diag(ETA[[j]]$Cov$SigmaB))
          }
          else{
            ETA[[j]]$Cov$SigmaBInv<-solve(ETA[[j]]$Cov$SigmaB)
          }

          #Sampling from posterior of Environment parameters

          if(!is.null(ETA[[j]]$GxExT)){

            Q_old=ETA[[j]]$GxExT$eigenvecs_EG
            Lambda_old = ETA[[j]]$GxExT$eigenvals_EG

            if(response_type%in%c("DLN","gaussian")){
              # Sampling from full conditional of Bj's
              B2 = .Call("fBj_DLN_G_mtme",ETA[[j]]$GxExT$Bv,ETA[[j]]$Cov$SigmaBInv,
                         resCov$SigmaEInv,rv,ETA[[j]]$GxExT$X,ETA[[j]]$GxExT$x2,
                         n,ETA[[j]]$GxExT$p,ntraits,PACKAGE="BGLM")
            }
            if(response_type%in%c("Poisson","PLN")){
              # Sampling from full conditional of Bj's
              B2 = .Call("fBj_P_mtme",l,yr,ETA[[j]]$GxExT$Bv,ETA[[j]]$Cov$SigmaBInv,
                         rv,ETA[[j]]$GxExT$X,ETA[[j]]$GxExT$X2,ETA[[j]]$GxExT$p,
                         ntraits,n,PACKAGE="BGLM")
            }

            ETA[[j]]$GxExT$Bv = B2[[1]]
            rv = B2[[2]]

            MEnv = matrix(matrixcalc::vec(t(ETA[[j]]$GxExT$Bv)),ncol = ETA[[j]]$nEnv,byrow = F)
            SE2 = crossprod(MEnv)

            if (ETA[[j]]$GxExT$Cov$type == "UN") {
              ETA[[j]]$GxExT$Cov$SigmaEnv<-MCMCpack::riwish(v = ETA[[j]]$GxExT$Cov$priorvE+ETA[[j]]$nLinesTraits,
                                                    S = SE2 + ETA[[j]]$GxExT$Cov$priorSE)
            }
            if (ETA[[j]]$GxExT$Cov$type == "DIAG") {
              VVB = ETA[[j]]$GxExT$Cov$priorvE+ETA[[j]]$nLinesTraits
              SSB = ETA[[j]]$GxExT$Cov$priorSE+diag(SE2)
              ETA[[j]]$GxExT$Cov$SigmaEnv<-diag(1/rgamma(ETA[[j]]$nEnv,VVB/2,SSB/2))
            }
            if (ETA[[j]]$GxExT$Cov$type == "FA") {
              tmp<-sample_G0_FA(U = MEnv,F = ETA[[j]]$GxExT$Cov$F,
                                  M = ETA[[j]]$GxExT$Cov$M,B = ETA[[j]]$GxExT$Cov$W,PSI = ETA[[j]]$GxExT$Cov$PSI,
                                  ntraits = ETA[[j]]$nEnv,nF = ETA[[j]]$GxExT$Cov$nF,nD = ETA[[j]]$GxExT$Cov$nD,
                                  df0 = ETA[[j]]$GxExT$Cov$priorvE,S0 = ETA[[j]]$GxExT$Cov$priorSE,
                                  priorVar = ETA[[j]]$GxExT$Cov$var,varimaxRotate = ETA[[j]]$GxExT$Cov$varimax)
              ETA[[j]]$GxExT$Cov$F<-tmp$F
              ETA[[j]]$GxExT$Cov$PSI<-tmp$PSI
              ETA[[j]]$GxExT$Cov$W<-tmp$B
              ETA[[j]]$GxExT$Cov$SigmaEnv<-tmp$G
              rm(tmp)
            }
            if (ETA[[j]]$GxExT$Cov$type == "REC") {
              tmp<-sample_G0_REC(U = MEnv,M = ETA[[j]]$GxExT$Cov$M,
                                   PSI = ETA[[j]]$GxExT$Cov$PSI,ntraits = ETA[[j]]$nEnv,
                                   priorVar = ETA[[j]]$GxExT$Cov$var,df0 = ETA[[j]]$GxExT$Cov$priorvE,
                                   S0 = ETA[[j]]$GxExT$Cov$priorSE)
              ETA[[j]]$GxExT$Cov$SigmaEnv<-tmp$G
              ETA[[j]]$GxExT$Cov$W<-tmp$B
              ETA[[j]]$GxExT$Cov$PSI<-tmp$PSI
              rm(tmp)
            }

            EVD_Env = eigen(ETA[[j]]$GxExT$Cov$SigmaEnv,symmetric = T)

            keep<-EVD_Env$values>1e-10
            eigenvals_E<-EVD_Env$values[keep]
            eigenvecs_E<-EVD_Env$vectors[,keep]

            ETA[[j]]$GxExT$eigenvecs_EG = .Call("kron",eigenvecs_E,ETA[[j]]$EVD$vectors,
                                                PACKAGE="BGLM")
            ETA[[j]]$GxExT$eigenvals_EG = as.vector(.Call("kron_vec",eigenvals_E,
                                                          ETA[[j]]$EVD$values,
                                                          PACKAGE="BGLM"))

            Q_new=ETA[[j]]$GxExT$eigenvecs_EG
            Lambda_new = ETA[[j]]$GxExT$eigenvals_EG

            Bv_old<-ETA[[j]]$GxExT$Bv
            Bv_scaled<-sweep(Bv_old,1,sqrt(Lambda_old),FUN="*")
            Bv_proj<-.Call("prod",t(Q_new),.Call("prod",Q_old,Bv_scaled,PACKAGE="BGLM"),
                           PACKAGE="BGLM")
            Bv_new<-sweep(Bv_proj,1,1/sqrt(Lambda_new),FUN="*")

            # Update Bv
            ETA[[j]]$GxExT$Bv<-Bv_new

            #rv<-rv + ETA[[j]]$GxExT$X %*% Bv_old
            rv<-.Call("fXb_multi_prod",ETA[[j]]$GxExT$X,Bv_old,rv,1.0,
                      PACKAGE="BGLM")

            ETA[[j]]$GxExT$Aux<-sweep(x=Q_new,MARGIN=2L,STATS=sqrt(Lambda_new),
                                      FUN="*")
            ETA[[j]]$GxExT$X<-.Call("prod",ETA[[j]]$GxExT$Z,ETA[[j]]$GxExT$Aux,
                                    PACKAGE="BGLM")

            ETA[[j]]$GxExT$x2 = as.vector(colSums(ETA[[j]]$GxExT$X^2))
            ETA[[j]]$GxExT$p<-ncol(ETA[[j]]$GxExT$X)

            if(response_type%in%c("Poisson","PLN")){ETA[[j]]$GxExT$X2<-ETA[[j]]$GxExT$X^2}

            #rv<-rv - ETA[[j]]$GxExT$X %*% Bv_new
            rv<-.Call("fXb_multi_prod",ETA[[j]]$GxExT$X,Bv_new,rv,-1.0,
                      PACKAGE="BGLM")
          }
        }

      }
    }

    if (response_type == "PLN") {
      Us = .Call("fUi_mtme",l,yr,U,resCov$SigmaEInv,rv,ntraits,n,PACKAGE="BGLM")
      U = Us[[1]]
      rv = Us[[2]]
      rvPois = sweep(x=rv,MARGIN=2,STATS=log(r),FUN="+") - U
    }

    if(response_type%in%c("DLN","gaussian")){
      e = rv
    }
    if(response_type == "PLN"){
      e = U
    }

    if(response_type %in% c("DLN","PLN","gaussian")){
      CP<-crossprod(e)
      # Sampling from full conditional of SigmaE
      if (resCov$type == "UN") {
        resCov$SigmaE<-MCMCpack::riwish(v = resCov$priorv + n,S = CP + resCov$priorS)
      }
      if (resCov$type == "DIAG") {
        VV = resCov$priorv+n
        SS = resCov$priorS+diag(CP)
        resCov$SigmaE<-diag(1/rgamma(ntraits,VV/2,SS/2))
      }
      if (resCov$type == "FA") {
        tmp<-sample_G0_FA(U = e,F = resCov$F,M = resCov$M,
                            B = resCov$W,PSI = resCov$PSI,ntraits = ntraits,
                            nF = resCov$nF,nD = nrow(e),df0 = resCov$priorv,
                            S0 = resCov$priorS,priorVar = resCov$var,varimaxRotate = resCov$varimax)
        resCov$F<-tmp$F
        resCov$PSI<-tmp$PSI
        resCov$W<-tmp$B
        resCov$SigmaE<-tmp$G
      }
      if (resCov$type == "REC") {
        tmp<-sample_G0_REC(U = e,M = resCov$M,PSI = resCov$PSI,
                             ntraits = ntraits,priorVar = resCov$var,df0 = resCov$priorv,
                             S0 = resCov$priorS)
        resCov$SigmaE<-tmp$G
        resCov$W<-tmp$B
        resCov$PSI<-tmp$PSI
        rm(tmp)
      }
      if(resCov$type == "DIAG"){
        resCov$SigmaEInv<-diag(1/diag(resCov$SigmaE))
      }
      else{
        resCov$SigmaEInv<-solve(resCov$SigmaE)
      }
    }

    # Prediction,missing values and logLik
    if(response_type == "DLN"){
      rv = l - rv
      if(type_prediction == "mean"){
        y_tr = floor(exp(rv+MASS::mvrnorm(n=n,mu=rep(0,ntraits),Sigma=resCov$SigmaE)))
      }
      if(nNa>0){
        lp = rv[!NoWhichNa,]+MASS::mvrnorm(n=nNa,mu=rep(0,ntraits),Sigma=resCov$SigmaE)
        yStar[!NoWhichNa,] = floor(exp(lp))
        ay[!NoWhichNa,]=log(yStar[!NoWhichNa,]); by[!NoWhichNa,]=log(yStar[!NoWhichNa,]+1)
      }
      loglik = fllp_DLN_multi(y = y[NoWhichNa,],n = nObs,ntraits = ntraits,a = ay[NoWhichNa,],
                              b = by[NoWhichNa,],rv = rv[NoWhichNa,],SigmaE = resCov$SigmaE,
                              U_qmc)
    }

    if(response_type == "Poisson"){
      rvPois = rv + log(r)
      y_tr = exp(rvPois)
      if(nNa>0){
        rp = y_tr[!NoWhichNa,]
        yStar[!NoWhichNa,] = matrix(rpois(n=nNa*ntraits,lambda=as.vector(rp)),nrow=nNa,ncol=ntraits)
      }
      loglik = fllp_P_multi(rv=rv[NoWhichNa,],y=y[NoWhichNa,],r=r[NoWhichNa,],ntraits = ntraits)
    }

    if(response_type == "PLN"){
      y_tr = exp(sweep(x=rvPois,MARGIN=2,STATS=diag(resCov$SigmaE)/2,FUN="+"))
      if(nNa>0){
        rp = exp(rvPois[!NoWhichNa,]+U[!NoWhichNa,])
        yStar[!NoWhichNa,] = matrix(rpois(n=nNa*ntraits,lambda=as.vector(rp)),nrow=nNa,ncol=ntraits)
      }
      loglik = fllp_PLN_multi_R(y=y[NoWhichNa,],n=nObs,ntraits=ntraits,logy_factorial=logy_factorial,
                              rvPois=rvPois[NoWhichNa,],SigmaE=resCov$SigmaE,gh=gh,Q=Q,
                              nodes_mult=nodes_mult,gh_weights_mult=gh_weights_mult)
    }

    if(response_type == "gaussian"){
      y_tr = yStar - rv
      if(nNa>0){
        yStar[!NoWhichNa,] = y_tr[!NoWhichNa,]+MASS::mvrnorm(n=nNa,mu=rep(0,ntraits),Sigma=resCov$SigmaE)
        rv[!NoWhichNa,] = yStar[!NoWhichNa,] - y_tr[!NoWhichNa,]
      }
      loglik = cte + partial_fllp_G(rv = rv[NoWhichNa,],SigmaE = resCov$SigmaE)
    }

    if(response_type == "Poisson"){
      rv = rv + log(r)
      Lambda = exp(rv)
      r = (Lambda+1) * 10^(rcontrol)
      y_r = yStar+r
      yr = (yStar - r) / 2
      Syr = colSums(yStar-r)
      rv = rv - log(r)
    }

    if(response_type == "PLN"){
      rv = sweep(x=rv,MARGIN=2,STATS=log(r),FUN="+")
      Lambda = exp(rv)
      r = ((1/(exp(diag(resCov$SigmaE))-1))+1) * 10^(rcontrol)
      y_r = sweep(x=yStar,MARGIN=2,STATS=r,FUN="+")
      yminusr = sweep(x=yStar,MARGIN=2,STATS=r,FUN="-")
      yr = (yminusr) / 2
      Syr = colSums(yminusr)
      rv = sweep(x=rv,MARGIN=2,STATS=log(r),FUN="-")
    }

    if (i > nBurnin)
    {
      if(i%%nThin==0){
        nk = nk + 1

        if(intercept_mt){
          write(beta0,ncolumns = length(beta0),file = f_beta0,append = TRUE,
                sep = " ")

          # Mean of Beta0
          post_beta0 = (beta0 + (nk - 1) * post_beta0) / nk
          post_beta02 = (beta0^2 + (nk - 1) * post_beta02) / nk
        }

        if(nLT > 0){
          for (j in 1:nLT) {
            if(ETA[[j]]$model=="FIXED"){
              ETA[[j]]$post_Bv = (ETA[[j]]$Bv + (nk - 1) * ETA[[j]]$post_Bv) / nk
              ETA[[j]]$post_Bv2 = (ETA[[j]]$Bv^2 + (nk - 1) * ETA[[j]]$post_Bv2) / nk
            }
            if (ETA[[j]]$model%in%c("BRR","RKHS")) {
              # Mean of Beta's
              ETA[[j]]$post_Bv = (ETA[[j]]$Bv + (nk - 1) * ETA[[j]]$post_Bv) / nk
              ETA[[j]]$post_Bv2 = (ETA[[j]]$Bv^2 + (nk - 1) * ETA[[j]]$post_Bv2) / nk

              # Mean of varB
              ETA[[j]]$Cov$post_SigmaB = (ETA[[j]]$Cov$SigmaB + (nk - 1) * ETA[[j]]$Cov$post_SigmaB) / nk
              ETA[[j]]$Cov$post_SigmaB2 = (ETA[[j]]$Cov$SigmaB^2 + (nk - 1) * ETA[[j]]$Cov$post_SigmaB2) / nk

              if (ETA[[j]]$Cov$type %in% c("UN","DIAG")) {
                tmp<-matrixcalc::vech(ETA[[j]]$Cov$SigmaB)
                write(tmp,ncolumns = length(tmp),file = ETA[[j]]$Cov$f_SigmaB,
                      append = TRUE,sep = " ")
                rm(tmp)
              }

              if (ETA[[j]]$Cov$type %in% c("FA","REC")) {
                ETA[[j]]$Cov$post_W<-(ETA[[j]]$Cov$W + (nk - 1) * ETA[[j]]$Cov$post_W) / nk
                ETA[[j]]$Cov$post_W2<-(ETA[[j]]$Cov$W^2 + (nk - 1) * ETA[[j]]$Cov$post_W2) / nk
                ETA[[j]]$Cov$post_PSI<-(ETA[[j]]$Cov$PSI + (nk - 1) * ETA[[j]]$Cov$post_PSI) / nk
                ETA[[j]]$Cov$post_PSI2<-(ETA[[j]]$Cov$PSI^2 + (nk - 1) * ETA[[j]]$Cov$post_PSI2) / nk

                if (sum(ETA[[j]]$Cov$M) > 0) {
                  tmp<-ETA[[j]]$Cov$W
                  write(tmp,ncolumns = length(tmp),file = ETA[[j]]$Cov$f_W,
                        append = TRUE,sep = " ")
                  rm(tmp)
                }
                write(ETA[[j]]$Cov$PSI,ncolumns = length(ETA[[j]]$Cov$PSI),
                      file = ETA[[j]]$Cov$f_PSI,append = TRUE,
                      sep = " ")
              }
            }
            if (ETA[[j]]$model=="RKHS_mtme") {
              #Mean of Beta's
              ETA[[j]]$post_Bv = (ETA[[j]]$Bv + (nk - 1) * ETA[[j]]$post_Bv) / nk
              ETA[[j]]$post_Bv2 = (ETA[[j]]$Bv^2 + (nk - 1) * ETA[[j]]$post_Bv2) / nk

              # Mean of varB
              ETA[[j]]$Cov$post_SigmaB = (ETA[[j]]$Cov$SigmaB + (nk - 1) * ETA[[j]]$Cov$post_SigmaB) / nk
              ETA[[j]]$Cov$post_SigmaB2 = (ETA[[j]]$Cov$SigmaB^2 + (nk - 1) * ETA[[j]]$Cov$post_SigmaB2) / nk

              if (ETA[[j]]$Cov$type %in% c("UN","DIAG")) {
                tmp<-matrixcalc::vech(ETA[[j]]$Cov$SigmaB)
                write(tmp,ncolumns = length(tmp),file = ETA[[j]]$Cov$f_SigmaB,
                      append = TRUE,sep = " ")
                rm(tmp)
              }

              if (ETA[[j]]$Cov$type %in% c("FA","REC")) {
                ETA[[j]]$Cov$post_W<-(ETA[[j]]$Cov$W + (nk - 1) * ETA[[j]]$Cov$post_W) / nk
                ETA[[j]]$Cov$post_W2<-(ETA[[j]]$Cov$W^2 + (nk - 1) * ETA[[j]]$Cov$post_W2) / nk
                ETA[[j]]$Cov$post_PSI<-(ETA[[j]]$Cov$PSI + (nk - 1) * ETA[[j]]$Cov$post_PSI) / nk
                ETA[[j]]$Cov$post_PSI2<-(ETA[[j]]$Cov$PSI^2 + (nk - 1) * ETA[[j]]$Cov$post_PSI2) / nk

                if (sum(ETA[[j]]$Cov$M) > 0) {
                  tmp<-ETA[[j]]$Cov$W
                  write(tmp,ncolumns = length(tmp),file = ETA[[j]]$Cov$f_W,
                        append = TRUE,sep = " ")
                  rm(tmp)
                }
                write(ETA[[j]]$Cov$PSI,ncolumns = length(ETA[[j]]$Cov$PSI),
                      file = ETA[[j]]$Cov$f_PSI,append = TRUE,
                      sep = " ")
              }

              if(!is.null(ETA[[j]]$GxExT)){

                ETA[[j]]$GxExT$post_Bv = (ETA[[j]]$GxExT$Bv + (nk - 1) * ETA[[j]]$GxExT$post_Bv) / nk
                ETA[[j]]$GxExT$post_Bv2 = (ETA[[j]]$GxExT$Bv^2 + (nk - 1) * ETA[[j]]$GxExT$post_Bv2) / nk


                ETA[[j]]$GxExT$Cov$post_SigmaEnv = (ETA[[j]]$GxExT$Cov$SigmaEnv + (nk - 1) * ETA[[j]]$GxExT$Cov$post_SigmaEnv) / nk
                ETA[[j]]$GxExT$Cov$post_SigmaEnv2 = (ETA[[j]]$GxExT$Cov$SigmaEnv^2 + (nk - 1) * ETA[[j]]$GxExT$Cov$post_SigmaEnv2) / nk

                if (ETA[[j]]$GxExT$Cov$type %in% c("UN","DIAG")) {
                  tmp<-matrixcalc::vech(ETA[[j]]$GxExT$Cov$SigmaEnv)
                  write(tmp,ncolumns = length(tmp),file = ETA[[j]]$GxExT$Cov$f_SigmaEnv,
                        append = TRUE,sep = " ")
                  rm(tmp)
                }
                if (ETA[[j]]$GxExT$Cov$type %in% c("FA","REC")) {
                  ETA[[j]]$GxExT$Cov$post_W<-(ETA[[j]]$GxExT$Cov$W + (nk - 1) * ETA[[j]]$GxExT$Cov$post_W) / nk
                  ETA[[j]]$GxExT$Cov$post_W2<-(ETA[[j]]$GxExT$Cov$W^2 + (nk - 1) * ETA[[j]]$GxExT$Cov$post_W2) / nk
                  ETA[[j]]$GxExT$Cov$post_PSI<-(ETA[[j]]$GxExT$Cov$PSI + (nk - 1) * ETA[[j]]$GxExT$Cov$post_PSI) / nk
                  ETA[[j]]$GxExT$Cov$post_PSI2<-(ETA[[j]]$GxExT$Cov$PSI^2 + (nk - 1) * ETA[[j]]$GxExT$Cov$post_PSI2) / nk

                  if (sum(ETA[[j]]$GxExT$Cov$M) > 0) {
                    tmp<-ETA[[j]]$GxExT$Cov$W
                    write(tmp,ncolumns = length(tmp),file = ETA[[j]]$GxExT$Cov$f_W,
                          append = TRUE,sep = " ")
                    rm(tmp)
                  }
                  write(ETA[[j]]$GxExT$Cov$PSI,ncolumns = length(ETA[[j]]$GxExT$Cov$PSI),
                        file = ETA[[j]]$GxExT$Cov$f_PSI,append = TRUE,
                        sep = " ")
                }

              }
            }

          }
        }

        if(response_type %in% c("DLN","PLN","gaussian")){
          # Mean of SigmaE
          resCov$post_SigmaE = (resCov$SigmaE + (nk - 1) * resCov$post_SigmaE) / nk
          resCov$post_SigmaE2 = (resCov$SigmaE^2 + (nk - 1) * resCov$post_SigmaE2) / nk

          tmp<-matrixcalc::vech(resCov$SigmaE)
          write(tmp,ncolumns = length(tmp),file = resCov$f_SigmaE,
                append = TRUE,sep = " ")

          rm(tmp)

          if (resCov$type %in% c("REC","FA")) {
            resCov$post_W<-(resCov$W + (nk - 1) * resCov$post_W) / nk
            resCov$post_W2<-(resCov$W2 + (nk - 1) * resCov$post_W2) / nk
            resCov$post_PSI<-(resCov$PSI + (nk - 1) * resCov$post_PSI) / nk
            resCov$post_PSI2<-(resCov$PSI2 + (nk - 1) * resCov$post_PSI2) / nk

            if (sum(resCov$M) > 0) {
              tmp<-resCov$W
              write(tmp,ncolumns = length(tmp),file = resCov$f_W,
                    append = TRUE,sep = " ")
              rm(tmp)
            }
            write(resCov$PSI,ncolumns = length(resCov$PSI),
                  file = resCov$f_PSI,append = TRUE,sep = " ")
          }
        }

        write(loglik,ncolumns = length(loglik),file = f_logLik,append = TRUE,
              sep = " ")

        # Mean of loglik
        post_logLik = (loglik + (nk - 1) * post_logLik) / nk
        post_logLik2 = (loglik^2 + (nk - 1) * post_logLik2) / nk

        if(response_type%in%c("Poisson","PLN")){
          post_r = (r + (nk - 1) * post_r) / nk
          post_r2 = (r^2 + (nk - 1) * post_r2) / nk
        }

        #mean predictions
        if(response_type == "DLN"){
          if(type_prediction == "mean"){
            y_pred = ( y_tr + (nk - 1) * y_pred ) / nk
            y_pred2 = ( y_tr^2 + (nk - 1) * y_pred2 ) / nk
          }
          rv_mean = ( rv + (nk - 1) * rv_mean ) / nk
        }

        if(response_type == "Poisson"){
          y_pred = ( y_tr + (nk - 1) * y_pred ) / nk
          y_pred2 = ( y_tr^2 + (nk - 1) * y_pred2 ) / nk
          rv_mean = ( rvPois + (nk - 1) * rv_mean ) / nk
        }

        if(response_type == "PLN"){
          y_pred = ( y_tr + (nk - 1) * y_pred ) / nk
          y_pred2 = ( y_tr^2 + (nk - 1) * y_pred2 ) / nk
          rv_mean = ( rvPois + (nk - 1) * rv_mean ) / nk
        }

        if(response_type == "gaussian"){
          y_pred = ( y_tr + (nk - 1) * y_pred ) / nk
          y_pred2 = ( y_tr^2 + (nk - 1) * y_pred2 ) / nk
          rv_mean = ( rv + (nk - 1) * rv_mean ) / nk
        }

      }
    }

    if (verbose) {
      pb$tick()
    }
  }
  tmp2<-proc.time()[3]
  message("Time = ",round((tmp2-time)/60,3)," minutes")
  cat('\n')
  #-----------------------------------------------------------------------------

  if(intercept_mt){
    close(f_beta0)
    f_beta0<-NULL
  }

  close(f_logLik)
  f_logLik<-NULL

  # Bayesian estimations
  #---------------------------------------------------------------------------
  # Predictions
  fit = list()

  if(response_type == "DLN"){
    fit$logLikAtPostMean = fllp_DLN_multi(y=y[NoWhichNa,],n=nObs,ntraits=ntraits,a=ay[NoWhichNa,],
                                          b=by[NoWhichNa,],rv=rv_mean[NoWhichNa,],
                                          SigmaE=resCov$post_SigmaE,U_qmc)
    fit$pD = -2 * (post_logLik - fit$logLikAtPostMean)
    fit$DIC = fit$pD - 2 * post_logLik
    if(type_prediction == "median"){
      if (nNa==0){
        fit$y_train = ceiling(exp(rv_mean)-1)
        rownames(fit$y_train) = rowNames; colnames(fit$y_train) = colNames
      }else{
        fit$y_train = ceiling(exp(rv_mean[NoWhichNa,])-1)
        rownames(fit$y_train) = rowNames[NoWhichNa]; colnames(fit$y_train) = colNames
        fit$y_test = ceiling(exp(rv_mean[!NoWhichNa,]) - 1)
        rownames(fit$y_test) = rowNames[!NoWhichNa]; colnames(fit$y_test) = colNames
      }
    }
    else{
      fit$y_train = y_pred[NoWhichNa,]
      fit$SD.y_train = sqrt(y_pred2[NoWhichNa,] - y_pred[NoWhichNa,]^2)
      if(nNa==0){
        rownames(fit$y_train) = rowNames; colnames(fit$y_train) = colNames
        rownames(fit$SD.y_train) = rowNames; colnames(fit$SD.y_train) = colNames
      }
      else{
        rownames(fit$y_train) = rowNames[NoWhichNa]; colnames(fit$y_train) = colNames
        rownames(fit$SD.y_train) = rowNames[NoWhichNa]; colnames(fit$SD.y_train) = colNames
        fit$y_test = y_pred[!NoWhichNa,]
        fit$SD.y_test = sqrt(y_pred2[!NoWhichNa,] - y_pred[!NoWhichNa,]^2)
        rownames(fit$y_test) = rowNames[!NoWhichNa]; colnames(fit$y_test) = colNames
        rownames(fit$SD.y_test) = rowNames[!NoWhichNa]; colnames(fit$SD.y_test) = colNames
      }
    }
  }

  if(response_type == "Poisson"){
    fit$logLikAtPostMean = fllp_P_multi(rv=rv_mean[NoWhichNa,]-log(post_r[NoWhichNa,]),
                                        y=y[NoWhichNa,],r=post_r[NoWhichNa,],ntraits=ntraits)
    fit$pD = -2 * (post_logLik - fit$logLikAtPostMean)
    fit$DIC = fit$pD - 2 * post_logLik
    fit$y_train = y_pred[NoWhichNa,]
    fit$SD.y_train = sqrt(y_pred2[NoWhichNa,] - y_pred[NoWhichNa,]^2)
    if(nNa==0){
      rownames(fit$y_train) = rowNames; colnames(fit$y_train) = colNames
      rownames(fit$SD.y_train) = rowNames; colnames(fit$SD.y_train) = colNames
    }
    else{
      rownames(fit$y_train) = rowNames[NoWhichNa]; colnames(fit$y_train) = colNames
      rownames(fit$SD.y_train) = rowNames[NoWhichNa]; colnames(fit$SD.y_train) = colNames
      fit$y_test = y_pred[!NoWhichNa,]
      fit$SD.y_test = sqrt(y_pred2[!NoWhichNa,] - y_pred[!NoWhichNa,]^2)
      rownames(fit$y_test) = rowNames[!NoWhichNa]; colnames(fit$y_test) = colNames
      rownames(fit$SD.y_test) = rowNames[!NoWhichNa]; colnames(fit$SD.y_test) = colNames
    }
    fit$r = post_r
    fit$SD.r = sqrt(post_r2 - post_r^2)
  }

  if(response_type == "PLN"){
    fit$logLikAtPostMean = fllp_PLN_multi_R(y=y[NoWhichNa,],n=nObs,ntraits=ntraits,logy_factorial=logy_factorial,
                                          rvPois=rv_mean[NoWhichNa,],SigmaE=resCov$post_SigmaE,gh=gh,Q=Q,
                                          nodes_mult=nodes_mult,gh_weights_mult=gh_weights_mult)
    fit$pD = -2 * (post_logLik - fit$logLikAtPostMean)
    fit$DIC = fit$pD - 2 * post_logLik
    fit$y_train = y_pred[NoWhichNa,]
    fit$SD.y_train = sqrt(y_pred2[NoWhichNa,] - y_pred[NoWhichNa,]^2)
    if(nNa==0){
      rownames(fit$y_train) = rowNames; colnames(fit$y_train) = colNames
      rownames(fit$SD.y_train) = rowNames; colnames(fit$SD.y_train) = colNames
    }
    else{
      rownames(fit$y_train) = rowNames[NoWhichNa]; colnames(fit$y_train) = colNames
      rownames(fit$SD.y_train) = rowNames[NoWhichNa]; colnames(fit$SD.y_train) = colNames
      fit$y_test = y_pred[!NoWhichNa,]
      fit$SD.y_test = sqrt(y_pred2[!NoWhichNa,] - y_pred[!NoWhichNa,]^2)
      rownames(fit$y_test) = rowNames[!NoWhichNa]; colnames(fit$y_test) = colNames
      rownames(fit$SD.y_test) = rowNames[!NoWhichNa]; colnames(fit$SD.y_test) = colNames
    }
    fit$r = post_r
    fit$SD.r = sqrt(post_r2 - post_r^2)
  }

  if(response_type == "gaussian"){
    fit$logLikAtPostMean = cte + partial_fllp_G(rv = rv_mean[NoWhichNa,],
                                                SigmaE = resCov$post_SigmaE)
    fit$pD = -2 * (post_logLik - fit$logLikAtPostMean)
    fit$DIC = fit$pD - 2 * post_logLik
    fit$y_train = y_pred[NoWhichNa,]
    fit$SD.y_train = sqrt(y_pred2[NoWhichNa,] - y_pred[NoWhichNa,]^2)
    if(nNa==0){
      rownames(fit$y_train) = rowNames; colnames(fit$y_train) = colNames
      rownames(fit$SD.y_train) = rowNames; colnames(fit$SD.y_train) = colNames
    }
    else{
      rownames(fit$y_train) = rowNames[NoWhichNa]; colnames(fit$y_train) = colNames
      rownames(fit$SD.y_train) = rowNames[NoWhichNa]; colnames(fit$SD.y_train) = colNames
      fit$y_test = y_pred[!NoWhichNa,]
      fit$SD.y_test = sqrt(y_pred2[!NoWhichNa,] - y_pred[!NoWhichNa,]^2)
      rownames(fit$y_test) = rowNames[!NoWhichNa]; colnames(fit$y_test) = colNames
      rownames(fit$SD.y_test) = rowNames[!NoWhichNa]; colnames(fit$SD.y_test) = colNames
    }
  }

  # Standard deviations for intercept and logLik
  if(intercept_mt){
    fit$beta0 = post_beta0
    fit$SD.beta0 = sqrt(post_beta02 - post_beta0^2)
  }
  fit$post_logLik = post_logLik
  fit$SD.postlogLik = sqrt(post_logLik2 - post_logLik^2)

  if(response_type %in% c("DLN","PLN","gaussian")){
    resCov$SigmaE<-resCov$post_SigmaE
    resCov$SD.SigmaE<-sqrt(resCov$post_SigmaE2 - resCov$post_SigmaE^2)
    close(resCov$f_SigmaE)
    resCov$f_SigmaE<-NULL
    if (resCov$type %in% c("FA","REC")) {
      resCov$W<-resCov$post_W
      resCov$SD.W<-sqrt(resCov$post_W2 - resCov$post_W^2)
      resCov$PSI<-resCov$post_PSI
      resCov$SD.PSI<-sqrt(resCov$post_PSI2 - resCov$post_PSI^2)
      close(resCov$f_W)
      resCov$f_W<-NULL
      close(resCov$f_PSI)
      resCov$f_PSI<-NULL
    }
    tmp<-which(names(resCov) %in% c("post_SigmaE","post_SigmaE2","SigmaEInv",
                                      "post_W","post_W2","post_PSI","post_PSI2","F"))
    resCov<-resCov[-tmp]
    rm(tmp)
  }

  if(nLT > 0){
    for (i in 1:nLT) {
      if (ETA[[i]]$model %in% c("FIXED","BRR")) {
        ETA[[i]]$Bv = ETA[[i]]$post_Bv
        ETA[[i]]$SD.Bv = sqrt(ETA[[i]]$post_Bv2 - ETA[[i]]$post_Bv^2)
        if(!(response_type%in%c("Poisson","PLN"))){
          tmp = which(names(ETA[[i]]) %in% c("post_Bv","post_Bv2","X","x2"))
          ETA[[i]] = ETA[[i]][-tmp]
        }
        else{
          tmp = which(names(ETA[[i]]) %in% c("post_Bv","post_Bv2","X","X2"))
          ETA[[i]] = ETA[[i]][-tmp]
        }
      }

      if (ETA[[i]]$model=="RKHS") {
        ETA[[i]]$u = ETA[[i]]$X%*%ETA[[i]]$post_Bv
        ETA[[i]]$uStar = ETA[[i]]$post_Bv
        ETA[[i]]$SD.uStar = sqrt(ETA[[i]]$post_Bv2 - ETA[[i]]$post_Bv^2)
        if(!(response_type%in%c("Poisson","PLN"))){
          tmp = which(names(ETA[[i]]) %in% c("Bv","post_Bv","post_Bv2","X","x2",
                                             "X_train","p"))
          ETA[[i]] = ETA[[i]][-tmp]
        }
        else{
          tmp = which(names(ETA[[i]]) %in% c("Bv","post_Bv","post_Bv2","X","X2",
                                             "X_train","p"))
          ETA[[i]] = ETA[[i]][-tmp]
        }
      }
      if (ETA[[i]]$model=="RKHS_mtme") {
        if(!is.null(ETA[[i]]$GxT)){
          ETA[[i]]$u = ETA[[i]]$X%*%ETA[[i]]$post_Bv
          ETA[[i]]$uStar = ETA[[i]]$post_Bv
          ETA[[i]]$SD.uStar = sqrt(ETA[[i]]$post_Bv2 - ETA[[i]]$post_Bv^2)
        }
        if(!(response_type%in%c("Poisson","PLN"))){
          tmp = which(names(ETA[[i]]) %in% c("Bv","post_Bv","post_Bv2","X","x2",
                                             "X_train","p","eigenvals_GT",
                                             "eigenvecs_GT","eigenvals_EG",
                                             "eigenvecs_EG","Aux","nEnv",
                                             "nLines","nLinesTraits","nLinesEnvs",
                                             "SE","U1","SE1","U2"))
          ETA[[i]] = ETA[[i]][-tmp]
        }
        else{
          tmp = which(names(ETA[[i]]) %in% c("Bv","post_Bv","post_Bv2","X","X2",
                                             "X_train","p","eigenvals_GT",
                                             "eigenvecs_GT","eigenvals_EG",
                                             "eigenvecs_EG","Aux","nEnv",
                                             "nLines","nLinesTraits","nLinesEnvs",
                                             "SE","U1","SE1","U2"))
          ETA[[i]] = ETA[[i]][-tmp]
        }
        if(!is.null(ETA[[i]]$GxExT)){
          ETA[[i]]$GxExT$u = ETA[[i]]$GxExT$X%*%ETA[[i]]$GxExT$post_Bv
          ETA[[i]]$GxExT$uStar = ETA[[i]]$GxExT$post_Bv
          ETA[[i]]$GxExT$SD.uStar = sqrt(ETA[[i]]$GxExT$post_Bv2 - ETA[[i]]$GxExT$post_Bv^2)
          if(!(response_type%in%c("Poisson","PLN"))){
            tmp = which(names(ETA[[i]]$GxExT) %in% c("Bv","post_Bv","post_Bv2","X","x2",
                                                     "X_train","p","eigenvals_GT",
                                                     "eigenvecs_GT","eigenvals_EG",
                                                     "eigenvecs_EG","Aux","nEnv",
                                                     "nLines","nLinesTraits","nLinesEnvs"))
            ETA[[i]]$GxExT = ETA[[i]]$GxExT[-tmp]
          }
          else{
            tmp = which(names(ETA[[i]]$GxExT) %in% c("Bv","post_Bv","post_Bv2","X","X2",
                                                     "X_train","p","eigenvals_GT",
                                                     "eigenvecs_GT","eigenvals_EG",
                                                     "eigenvecs_EG","Aux","nEnv",
                                                     "nLines","nLinesTraits","nLinesEnvs"))
            ETA[[i]]$GxExT = ETA[[i]]$GxExT[-tmp]
          }
        }
      }
      if (ETA[[i]]$model %in% c("BRR","RKHS")) {
        ETA[[i]]$Cov$SigmaB = ETA[[i]]$Cov$post_SigmaB
        ETA[[i]]$Cov$SD.SigmaB = sqrt(ETA[[i]]$Cov$post_SigmaB2 - (ETA[[i]]$Cov$post_SigmaB^2))

        if (ETA[[i]]$Cov$type %in% c("UN","DIAG")) {
          close(ETA[[i]]$Cov$f_SigmaB)
          ETA[[i]]$Cov$f_SigmaB<-NULL
        }

        if (ETA[[i]]$Cov$type %in% c("FA","REC")) {
          ETA[[i]]$Cov$W<-ETA[[i]]$Cov$post_W
          ETA[[i]]$Cov$SD.W<-sqrt(ETA[[i]]$Cov$post_W2 - ETA[[i]]$Cov$post_W^2)
          ETA[[i]]$Cov$PSI<-ETA[[i]]$Cov$post_PSI
          ETA[[i]]$Cov$SD.PSI<-sqrt(ETA[[i]]$Cov$post_PSI2 - ETA[[i]]$Cov$post_PSI^2)

          close(ETA[[i]]$Cov$f_W)
          ETA[[i]]$Cov$f_W<-NULL
          close(ETA[[i]]$Cov$f_PSI)
          ETA[[i]]$Cov$f_PSI<-NULL
        }
        tmp<-which(names(ETA[[i]]$Cov) %in% c("post_SigmaB","post_SigmaB2","SigmaBInv",
                                                "post_W","post_W2","post_PSI","post_PSI2","F"))
        ETA[[i]]$Cov<-ETA[[i]]$Cov[-tmp]
        rm(tmp)
      }
      if (ETA[[i]]$model =="RKHS_mtme") {
        ETA[[i]]$Cov$SigmaB = ETA[[i]]$Cov$post_SigmaB
        ETA[[i]]$Cov$SD.SigmaB = sqrt(ETA[[i]]$Cov$post_SigmaB2 - (ETA[[i]]$Cov$post_SigmaB^2))

        if (ETA[[i]]$Cov$type %in% c("UN","DIAG")) {
          close(ETA[[i]]$Cov$f_SigmaB)
          ETA[[i]]$Cov$f_SigmaB<-NULL
        }
        if (ETA[[i]]$Cov$type %in% c("FA","REC")) {
          ETA[[i]]$Cov$W<-ETA[[i]]$Cov$post_W
          ETA[[i]]$Cov$SD.W<-sqrt(ETA[[i]]$Cov$post_W2 - ETA[[i]]$Cov$post_W^2)
          ETA[[i]]$Cov$PSI<-ETA[[i]]$Cov$post_PSI
          ETA[[i]]$Cov$SD.PSI<-sqrt(ETA[[i]]$Cov$post_PSI2 - ETA[[i]]$Cov$post_PSI^2)

          close(ETA[[i]]$Cov$f_W)
          ETA[[i]]$Cov$f_W<-NULL
          close(ETA[[i]]$Cov$f_PSI)
          ETA[[i]]$Cov$f_PSI<-NULL
        }

        tmp<-which(names(ETA[[i]]$Cov) %in% c("post_SigmaB","post_SigmaB2","SigmaBInv",
                                                "post_W","post_W2","post_PSI","post_PSI2","F"))
        ETA[[i]]$Cov<-ETA[[i]]$Cov[-tmp]
        rm(tmp)

        if(!is.null(ETA[[j]]$GxExT)){

          ETA[[i]]$GxExT$Cov$SigmaEnv = ETA[[i]]$GxExT$Cov$post_SigmaEnv
          ETA[[i]]$GxExT$Cov$SD.SigmaEnv = sqrt(ETA[[i]]$GxExT$Cov$post_SigmaEnv2 - (ETA[[i]]$GxExT$Cov$post_SigmaEnv^2))

          if (ETA[[i]]$GxExT$Cov$type %in% c("UN","DIAG")) {
            close(ETA[[i]]$GxExT$Cov$f_SigmaEnv)
            ETA[[i]]$GxExT$Cov$f_SigmaEnv<-NULL
          }
          if (ETA[[i]]$GxExT$Cov$type %in% c("FA","REC")) {
            ETA[[i]]$GxExT$Cov$W<-ETA[[i]]$GxExT$Cov$post_W
            ETA[[i]]$GxExT$Cov$SD.W<-sqrt(ETA[[i]]$GxExT$Cov$post_W2 - ETA[[i]]$GxExT$Cov$post_W^2)
            ETA[[i]]$GxExT$Cov$PSI<-ETA[[i]]$GxExT$Cov$post_PSI
            ETA[[i]]$GxExT$Cov$SD.PSI<-sqrt(ETA[[i]]$GxExT$Cov$post_PSI2 - ETA[[i]]$GxExT$Cov$post_PSI^2)

            close(ETA[[i]]$GxExT$Cov$f_W)
            ETA[[i]]$GxExT$Cov$f_W<-NULL
            close(ETA[[i]]$GxExT$Cov$f_PSI)
            ETA[[i]]$GxExT$Cov$f_PSI<-NULL
          }

          tmp<-which(names(ETA[[i]]$GxExT$Cov) %in% c("post_SigmaEnv","post_SigmaEnv2","SigmaEnvInv",
                                                        "post_W","post_W2","post_PSI","post_PSI2","F"))
          ETA[[i]]$GxExT$Cov<-ETA[[i]]$GxExT$Cov[-tmp]
          rm(tmp)
        }
      }
    }
  }

  out<-list(y=y,NoWhichNa=NoWhichNa,fit = fit,nIter = nIter,nBurnin = nBurnin,nThin = nThin,
              priorSigma0 = priorSigma0,R2 = R2,resCov = resCov,
              response_type=response_type,ntraits = ntraits)



  out$ETA<-ETA
  class(out) = "BGLM"
  return(out)

}
