# PCovRMIX #
# Bjarne Daems #

# 1. pre-processing function ####
pre_proc <- function(predictors){
  # inputs #
  # predictors: only pre-process the predictor variables. Set aside the outcome variable Y and don't add that 
  # to the pre-processing function
  
  #Split the dataset into quantitative and qualitative variables
  split <- PCAmixdata::splitmix(predictors)
  
  X1 <- split$X.quanti
  X2 <- split$X.quali
  
  #Scaling and centering the quantitative data (X1)
  center_X1 <- scale(X1, center = T, scale = F)
  
  #Creating and centering the indicator matrix (G) from the qualitative data (X2)
  dummy <- mltools::one_hot(data.table::as.data.table(X2))
  dummymatrix <- as.matrix(dummy)
  
  G <- scale(dummymatrix, center = T, scale = F)
  
  #This version of SD is used over the standard on in the scale function
  #This version uses length(x) instead of length(x)-1 as in Chavent et al. 2017
  sd_own <- function(x){
    sqrt(sum(x^2)/length(x))
  }
  
  sd_X1 <- apply(center_X1, 2, sd_own)
  scale_X1 <- center_X1 %*% diag(1/sd_X1)
  
  #Combine the processed qualitative and quantitative data back into one matrix  
  Z <- cbind(scale_X1,G)
  
  #The proportion of the columns
  prop <- nrow(dummy)/colSums(dummy)
  
  #Creating the matrices of weights of the columns (M) and rows (N)
  M <- diag(c(rep(1, ncol(X1)), nrow(dummy) / colSums(dummy)))
  N <-diag(rep(1/nrow(Z), nrow(Z)))
  
  returnobject <- list(Z=Z, M=M, N=N, prop=prop)
  
  return(returnobject)
}

# 2. pcovr function I've sent before ####
# this function gives the same solution as "pcovr_est" function in PCovR package
pcovr_svd <- function (X, Y, R, alpha) {
  
  N <- nrow(X)
  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  S <- t(X) %*% X
  
  # computing projection matrix #
  if (det(S) < 1e-12) {
    l <- eigen(S)$values
    w <- which(l > 1e-08 * max(l))
    Sh <- eigen(S)$vectors[, w] %*% diag(l[w]^-1) %*% t(eigen(S)$vectors[, 
                                                                         w])
    Hx <- X %*% Sh %*% t(X)
  } else {
    Hx <- X %*% solve(S) %*% t(X)
  }
  
  # block1: predictors 
  # scaling according to alpha
  block1 <- sqrt(alpha / sum(X^2)) * X
  
  # block2: response variable
  # scaling according to alpha
  block2 <- Hx %*% Y
  block2 <- sqrt((1 - alpha)/ sum(Y^2)) * block2
  
  Xstar <- cbind(block1, block2)
  
  # P <- V %*% eigen(t(V) %*% V)$vectors %*% diag(eigen(t(V) %*% V)$values^(-1/2))
  
  svdd <- svd(Xstar)
  
  Te <- svdd$u[,1:R]
  
  Px <- t(Te) %*% X
  Py <- t(Te) %*% Y
  
  Rx2 <- sum((Te %*% Px)^2) * sum(X^2)^-1
  Ry2 <- sum((Te %*% Py)^2) * sum(Y^2)^-1
  
  results <- list(Rx2 = Rx2, Ry2 = Ry2, Te = Te, 
                  Px = Px, Py = Py)
  
  return(results)
}

# 3. pcovrmix function ####
pcovrmix <- function(Z, M, N, Y, prop, R, alpha, qualis){
  
  # input #
  # Z: concatenated quantitative and qualitative predictor variables, pre-processed
  # M: M
  # N: N
  # Y: response variable
  # alpha: alpha weighting parameter
  # qualis: vector specifying the number of levels in each qualitative variable
  
  Z <- as.matrix(Z)
  Y <- as.matrix(Y)
  
  # M matrix re-defined to take account of the continuous outcome
  # M <- diag(c(rep(1, ncol(Y)), diag(M)))
  
  Msqrt <- sqrt(M)
  Nsqrt <- sqrt(N)
  
  # Zpca for the GSVD step
  Zpca <- Nsqrt %*% Z %*% Msqrt
  # continuous columns at this point of unit-norm
  
  # computing the projection matrix Hx #
  # (refer to the report that I sent you to know what the projection matrix is)
  S <- t(Zpca) %*% Zpca
  
  if (det(S)<1e-12){ 
    l <- eigen(S)$values
    w <- which(l>1e-8*max(l))
    Sh <- eigen(S)$vectors[,w] %*% diag(l[w]^-1) %*% t(eigen(S)$vectors[,w])
    Hx <- Zpca %*% Sh %*% t(Zpca)
  } else {
    Hx <- Zpca %*% solve(S) %*% t(Zpca)
  }
  
  # block1 is the block of predictors #
  # you can see that it is being scaled according to alpha 
  block1 <- sqrt(alpha / sum(Zpca^2))*Zpca
  
  alpha_scaling <- sqrt(alpha / sum(Zpca^2))
  
  # block2 is the block of the response variables 
  block2 <- Hx %*% Y
  
  block2 <- sqrt((1-alpha) / sum(Y^2)) * block2
  
  # concatenation of block1 and block2
  Zstar <- cbind(block1,block2)
  
  
  # svd step 
  svdd <- svd(Zstar)
  
  if (R == 1){
    Utilde <- matrix(svdd$u[,1:R], ncol = R)
    Ltilde <- matrix(svdd$d[1:R], ncol = R, nrow = R)
    Vtilde <- matrix(svdd$v[,1:R], ncol = R)
  }
  else{
    Utilde <- svdd$u[,1:R]
    Ltilde <- diag(svdd$d[1:R])
    Vtilde <- svdd$v[,1:R]
  }
  # Te = same as PCovR covariate
  Te <- svdd$u[,1:R]
  #Te <- scale(Te, F, T)
  ## Px and Py are computed according to the paper of Vervloet et al., (2015)
  # Px = loadings in the sense of PCovR. It recovers predictors Z. 
  # Py = regression coefficients
  Px <- t(Te) %*% Zpca
  
  Px <- sqrt(nrow(Z)) * Px %*% solve(Msqrt) # rescaling such that the solution is the same as PCovR's
  
  Py <- t(Te) %*% Y
  
  # W = weights matrix. Z (not Ztilde) * W = Te. 
  # this can be used for applying the model on a test dataset
  W <- MASS::ginv(Z) %*% Te
  
  # Fc = same as factor coordinates for observations in PCAmix
  # or, N^(-1/2) * Utilde * Ltilde (equation 18 in Chavent, 2017)
  # multiplying with alpha_scaling to restore scaling 
  Fc <- solve(Nsqrt) %*% Te %*% (Ltilde / alpha_scaling)
  
  singular_values <- svdd$d
  # total variance in the predictors, scaled according to PCAmix
  # this is only to calculate the percentage of variance explained by each component
  
  # Astar 
  # elements that correspond to quantitative variables are the same as PCAmix$quanti$coord
  Astar <- M %*% solve(Msqrt) %*% Vtilde[1:ncol(Z), ] %*% (Ltilde / alpha_scaling)
  
  # Astar part for the response variable
  Astar_y <- matrix(Vtilde[(ncol(Z)+1):nrow(Vtilde), ], ncol = R) %*% (Ltilde / alpha_scaling)
  
  # Cmat = C matrix (squared loadings) in Chavent (2017)
  # elements that correspond to qualitative variables are the same as PCAmix$quali$contrib
  Astarsq <- Astar^2
  
  
  if (!is.null(qualis)){
    
    Mdiag <- diag(M)
    
    quantis <- length(which(Mdiag == 1))
    
    qualis2 <- cumsum(qualis) + quantis
    
    quali_index <- list()
    
    quali_index[[1]] <- (quantis + 1):qualis2[1]
    
    for (i in 2:length(qualis)){
      quali_index[[i]] <- (qualis2[i-1]+1):qualis2[i]
    }
    
    prop2 <- c(rep(1, quantis), prop)
    
    # Cmat = C matrix in Chavent's paper; provides squared correlation, 
    # which represents the contribution of each variable to each component
    Cmat <- matrix(NA, nrow = length(qualis), ncol = ncol(Astar))
    
    for (i in 1:length(qualis)){
      
      Cmat[i,] <- t(t(Astarsq[quali_index[[i]],]) %*% (1/prop2)[quali_index[[i]]])
      
    }
    
    Cmat <- rbind(matrix(Astarsq[1:quantis,], ncol = R), Cmat)
    
  } else {
    
    Cmat <- Astarsq
    
  }
  
  # concatenating the Astar^2 part of response variable into Cmat,
  # so the contribution of response variable to each component can also be studied
  Cmat <- rbind(Cmat, Astar_y^2)
  
  # loading1 <- Astarsq[3,1] * (1/prop)[1] + Astarsq[4,1] * (1/prop)[2]
  # 
  # loading2 <- Astarsq[5,1] * (1/prop)[3] + Astarsq[6,1] * (1/prop)[4]
  
  # output #
  # Te: covariate scores (or component scores) Same as Te in PCovR
  # Px, Py: loadings and regression coefficients. Same as in PCovR
  # W: weights. Z * W = Te. Same as in PCovR
  # Fc: factor coordinates of observations. Same as in PCAmix
  # Astar: Astar matrix. Same as in PCAmix
  # Cmat: squared loadings matrix. Same as in PCAmix
  # singular_values: raw singular values from SVD
  # alpha_scaling: scaling done to the Z matrix due to the alpha processing for PCovR
  
  results <- list(Te = Te, Fc = Fc, Px=Px, Py=Py, W = W, Astar = Astar, Cmat = Cmat, 
                  singular_values = singular_values,
                  alpha_scaling = alpha_scaling)
  
  return(results)
  
}

# 4. cross-validation for components R ####
# initiated: 22-april-2021 #

# 1. explanation 

# for cross validation, first you can define the range of alpha (and / or with the number of components) 
# that you want to investigate:
#alpha_range <- c(0.3, 0.4, 0.5, 0.6)

R_range <- c(1,2,3,4,5) # (if you also want to CV for R)

# ranges <- expand.grid(alpha_range, R_range)

# colnames(ranges) <- c("alpha", "R")

# defining the matrix that would record the results of CV
cv_results_R <- matrix(NA, ncol = 3, nrow = length(R_range)) 

for (i in 1:length(R_range)){
  
  # this is where the method is being applied, using the CV function (defined below)
  # you can see that the alpha being used as input is different for each iteration i
  scdcv1 <- scdlogr_cv(predictors = dat[,-1], Y = Y, R = R_range[i], 
                       alpha = alpha_ml,
                       nrFolds = 10, seed = i)
  
  # record the CV results on the cv_results matrix
  cv_results_R[i,] <- c(scdcv1$cve, scdcv1$se, R_range[i])
  
  print(i)
}

cv_min_index <- which.min(cv_results_R[,1])

serule <- cv_results_R[cv_min_index,1] + cv_results_R[cv_min_index,2]

scd_possible <- cv_results_R[cv_results_R[,1] <= serule,]

if (is.vector(scd_possible)){
  scd_chosen <- scd_possible
} else {
  scd_chosen <- scd_possible[1,]
}
scd_chosen

# you can then look into the alpha value that results in the smallest CV error,
# or the biggest alpha value within the 1 standard error region from the minimal CV error

# CV function #
scdlogr_cv <- function(predictors, Y, R, alpha, qualis, nrFolds, seed){
  
  # set seed, because we would randomly divide the dataset into folds
  set.seed(seed)
  
  sd_own <- function(x){
    sqrt(sum(x^2)/length(x))
  }
  
  #Pre-process the data so PCovRMIX can use it
  split <- PCAmixdata::splitmix(predictors)
  
  X1 <- split$X.quanti
  X2 <- split$X.quali
  
  #Create the dummy (G) matrix from X2 and scale data
  
  dummy <- mltools::one_hot(data.table::as.data.table(X2))
  dummymatrix <- as.matrix(dummy)
  
  center_X1 <- scale(X1, center = T, scale = F)
  sd_X1 <- apply(center_X1, 2, sd_own)
  scale_X1 <- center_X1 %*% diag(1/sd_X1)
  
  # randomly shuffle the data
  sampled <- sample(length(Y))
  
  quanti <- scale_X1[sampled,]
  
  quali <- dummymatrix[sampled,]
  
  Y <- Y[sampled,]
  
  # Create equally sized folds
  folds <- cut(seq(1,nrow(predictors)),breaks=nrFolds,labels=FALSE)
  
  cve_k <- data.frame(error_y = NA)
  
  R2_k <- data.frame()
  
  for(k in 1:nrFolds){
    
    # defining test and train data  
    test_index <- which(folds==k, arr.ind=TRUE)
    
    X_train.quanti <- quanti[-test_index,]
    X_train.quali <- quali[-test_index,]
    Z_train <- cbind(X_train.quanti, X_train.quali)
    
    X_test.quanti <- quanti[test_index,]
    X_test.quali <- quali[test_index,]
    Z_test <- cbind(X_test.quanti, X_test.quali)
    
    Y_train <- Y[-test_index]
    Y_test <- Y[test_index]
    
    Y_train <- scale(Y_train, T, T)
    Y_test <- scale(Y_test, T, T)
    
    prop <- nrow(X_train.quali)/colSums(X_train.quali)
    M <- diag(c(rep(1, ncol(X1)), nrow(dummy) / colSums(dummy)))
    N <- diag(rep(1/nrow(Z_train), nrow(Z_train)))
    
    # model fitting #
    scd1 <- pcovrmix(Z = Z_train, M = M, N = N, Y = Y_train, prop = prop, R = R, alpha = alpha, qualis = qualis)
    
    # out of sample prediction #
    
    # (here, the criterion being used is the sum of squared errors)
    # (and the inverse-logit function is used because this is a logistic regression model)
    # (you can simply compute the sum of squared errors between the predicted Y and the observed Y)
    
    # (you would need to use the formula: Z_newdata * W * Py)
    pred_train <- Z_train %*% scd1$W %*% scd1$Py 
    
    R2 <-  1- sum((Y_train - pred_train)^2) * sum((Y_train)^2)^-1
    
    pred <- Z_test %*% scd1$W %*% scd1$Py 
    
    SSE <-  sum((Y_test - pred)^2) * sum((Y_test)^2)^-1
    
    cve_k[k,1] <- SSE
  }
  
  cve <- colMeans(cve_k)
  # cve = mean cv error
  
  se <- apply(cve_k,2,sd) / sqrt(nrFolds)
  # se = standard deviation of the mean cv errors per k / sqrt(nrFold)
  
  return(list(cve = cve, se = se, cve_k = cve_k))
}

# 5. cross-validation for alpha ####
dat[,-1]
library(PCAmixdata)
library(RegularizedSCA)
library(dplyr)
data("gironde")
dat <- gironde$housing

dim(dat)

head(dat)

Y <- dat[,1] # set aside the "density" variable as the response variable

Y <- scale(Y, T, T)

# pre-processing the predictors, which are columns 2-4 of dat
processed <- pre_proc(dat[,-1])

head(processed$Z)

x_pc <- prcomp(processed$Z)

#proportion unexplained variance for X
var_Ex <- (sum(x_pc$sdev^2) - sum(x_pc$sdev[1:3]^2)) /  sum(x_pc$sdev^2)

# proportion of unexplained variance by the dominant components = 22.4%

# for estimation of ey, this is the proportion of unexplained variance when y is regressed on X.

reg <- MASS::ginv(t(processed$Z) %*% processed$Z) %*% t(processed$Z) %*% Y

var_ey <- sum((processed$Z %*% reg - Y)^2) / sum(Y^2) # proportion of unexplained variance (1 - R^2)

# combining them together:
a_ml <- sum(processed$Z^2) / (sum(processed$Z^2) + sum(Y^2) * var_Ex / var_ey)  # 1 - alpha = 0.866
# therefore a lot of weight on X
a_ml

#Perform PCovRMIX with the maximum likelihood alpha
mix <- pcovrmix(Z = processed$Z, M = processed$M, N = processed$N, Y = Y, 
                prop = processed$prop, R = 6, alpha = a_ml, qualis = c(2,2))
t <- as.data.frame(mix$Te)
t %>% summarise_if(is.numeric, var)
#scree plot
plot(mix$singular_values^2, xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     type = "b")

#Based on scree plot, pick a number of components R (elbow method or Kaiser's criterion)

#Use R in cross-validation to determine the optimal alpha

alpha_range = seq(0,1,0.05)

cv_results_alpha <- matrix(NA, ncol = 3, nrow = length(alpha_range)) 

for (i in 1:length(alpha_range)){
  
  # this is where the method is being applied, using the CV function (defined below)
  # you can see that the alpha being used as input is different for each iteration i
  scdcv1 <- scdlogr_cv(predictors = dat[,-1], Y = Y, R = 3, 
                       alpha = alpha_range[i],
                       nrFolds = 10, seed = i, qualis = c(2,2))
  
  # record the CV results on the cv_results matrix
  cv_results_alpha[i,] <- c(scdcv1$cve, scdcv1$se, alpha_range[i])
  
  print(i)
}

cv_min_index <- which.min(cv_results_alpha[,1])

serule <- cv_results_alpha[cv_min_index,1] + cv_results_alpha[cv_min_index,2]

scd_possible <- cv_results_alpha[cv_results_alpha[,1] <= serule,]

if (is.vector(scd_possible)){
  scd_chosen <- scd_possible
} else {
  scd_chosen <- scd_possible[1,]
}
cv_results_alpha
scd_chosen <- matrix(scd_chosen, nrow = 1, ncol = 3)
colnames(scd_chosen) <- c("CVE", "SE", "alpha")
scd_chosen
cv_results_alpha[cv_min_index,3]

# you can then look into the alpha value that results in the smallest CV error,
# or the biggest alpha value within the 1 standard error region from the minimal CV error

# CV function #
scdlogr_cv <- function(predictors, Y, R, alpha, qualis, nrFolds, seed){
  
  # set seed, because we would randomly divide the dataset into folds
  set.seed(seed)
  
  sd_own <- function(x){
    sqrt(sum(x^2)/length(x))
  }
  
  #Pre-process the data so PCovRMIX can use it
  split <- PCAmixdata::splitmix(predictors)
  
  X1 <- split$X.quanti
  X2 <- split$X.quali
  
  #Create the dummy (G) matrix from X2 and scale data
  
  dummy <- mltools::one_hot(data.table::as.data.table(X2))
  dummymatrix <- as.matrix(dummy)
  
  center_X1 <- scale(X1, center = T, scale = F)
  sd_X1 <- apply(center_X1, 2, sd_own)
  scale_X1 <- center_X1 %*% diag(1/sd_X1)
  
  # randomly shuffle the data
  sampled <- sample(length(Y))
  
  quanti <- scale_X1[sampled,]
  
  quali <- dummymatrix[sampled,]
  
  Y <- Y[sampled,]
  
  # Create equally sized folds
  folds <- cut(seq(1,nrow(predictors)),breaks=nrFolds,labels=FALSE)
  
  cve_k <- data.frame(error_y = NA)
  
  for(k in 1:nrFolds){
    
    # defining test and train data  
    test_index <- which(folds==k, arr.ind=TRUE)
    
    X_train.quanti <- quanti[-test_index,]
    X_train.quali <- quali[-test_index,]
    Z_train <- cbind(X_train.quanti, X_train.quali)
    
    X_test.quanti <- quanti[test_index,]
    X_test.quali <- quali[test_index,]
    Z_test <- cbind(X_test.quanti, X_test.quali)
    
    Y_train <- Y[-test_index]
    Y_test <- Y[test_index]
    
    Y_train <- scale(Y_train, T, T)
    Y_test <- scale(Y_test, T, T)
    
    prop <- nrow(X_train.quali)/colSums(X_train.quali)
    M <- diag(c(rep(1, ncol(X1)), nrow(dummy) / colSums(dummy)))
    N <- diag(rep(1/nrow(Z_train), nrow(Z_train)))
    
    # model fitting #
    scd1 <- pcovrmix(Z = Z_train, M = M, N = N, Y = Y_train, prop = prop, R = R, alpha = alpha, qualis = qualis)
    
    # out of sample prediction #
    
    # (here, the criterion being used is the sum of squared errors)
    # (and the inverse-logit function is used because this is a logistic regression model)
    # (you can simply compute the sum of squared errors between the predicted Y and the observed Y)
    
    pred <- Z_test %*% scd1$W %*% scd1$Py 
    
    SSE <-  sum((Y_test - pred)^2) * sum((Y_test)^2)^-1
    
    cve_k[k,1] <- SSE
  }
  
  cve <- colMeans(cve_k)
  # cve = mean cv error
  
  se <- apply(cve_k,2,sd) / sqrt(nrFolds)
  # se = standard deviation of the mean cv errors per k / sqrt(nrFold)
  
  return(list(cve = cve, se = se, cve_k = cve_k))
}




# 6. cross-validation final model ####

cv_results <- matrix(NA, ncol = 4, nrow = 1) 

# this is where the method is being applied, using the CV function (defined below)
scdcv1 <- scdlogr_cv(predictors = dat[,-1], Y = Y, R = cv_R, 
                     alpha = cv_alpha,
                     nrFolds = 10, seed = i, qualis = c(2,2))

# record the CV results on the cv_results matrix
cv_results <- c(scdcv1$cve, scdcv1$se, cv_R, cv_alpha)
cv_results  

# 7. testing ####

# 4.1. when alpha = 1, solutions from PCovRmix is the same as solutions from PCAmix
library(PCAmixdata)
library(RegularizedSCA)

data("gironde")
dat <- gironde$housing

dim(dat)

head(dat)

Y <- dat[,1] # set aside the "density" variable as the response variable

Y <- scale(Y, T, T)

# pre-processing the predictors, which are columns 2-4 of dat
processed <- pre_proc(dat[,-1])

head(processed$Z)

# fitting the pcovrmix model
pcovrmixed1 <- pcovrmix(Z = processed$Z, M = processed$M, N = processed$N, Y = Y, 
                  prop = processed$prop, R = 3, alpha = 0.999921, qualis = c(2,2))
head(pcovrmixed1$Px)

pcovrmixed2 <- pcovrmix(Z = processed$Z, M = processed$M, N = processed$N, Y = Y, 
                        prop = processed$prop, R = 3, alpha = 1, qualis = c(2,2))

head(pcovrmixed2$Px)
# fitting the PCAmix model
X.quanti <- splitmix(dat)$X.quanti
X.quali <- splitmix(dat)$X.quali

pcamixed1 <-PCAmix(X.quanti = X.quanti[,-1], X.quali = X.quali,ndim=3)


# factor coordinates:
head(pcovrmixed1$Fc)
head(pcamixed1$ind$coord) # same values

TuckerCoef(pcovrmixed1$Fc, pcamixed1$ind$coord) # perfect correlation

# eigenvalues (explained variance per component):
pcamixed1$eig

pcovrmixed1$singular_values[1:3]^2 / sum(pcovrmixed1$singular_values^2)
# proportion of explained variance per component (or factor coordinate) is the same

round((pcovrmixed1$singular_values / pcovrmixed1$alpha_scaling)^2, 3)
# "eigenvalue" from PCovRmix can be found by scaling like this

# Astar matrix (very similar to loadings):
pcovrmixed1$Astar[1:2,] # first two rows give Astar matrix of PCAmix
pcamixed1$quanti$coord

# squared loadings (contribution of the qualitative variables to the components):
pcovrmixed1$Cmat[3:4,] # last two rows of Cmat
pcamixed1$quali$contrib


# squared loadings quantitative variables this time
pcovrmixed1$Cmat[1:2,]
pcamixed1$quanti$contrib

# studying the Astar matrix #

# factor coordinates scaled
sd_own <- function(x){
  sqrt(sum(x^2)/length(x))
}

sd_Fc <- apply(pcovrmixed1$Fc, 2, sd_own)
scale_Fc <- pcovrmixed1$Fc %*% diag(1/sd_Fc)

# dummy matrix of qualitative variables
dummy <- mltools::one_hot(data.table::as.data.table(dat[,c(3,5)]))
dummymatrix <- as.matrix(dummy)

colMeans(scale_Fc[dummymatrix[,1] == 1,])
colMeans(scale_Fc[dummymatrix[,2] == 1,])
pcovrmixed1$Astar

# studying the Cmat matrix #
pcovrmixed1$Cmat %*% diag(1/colSums(pcovrmixed1$Cmat))

# 4.2. when only quantitative predictors are used, the PCovRmix solution is the same as PCovR
library(PCovR)
data(alexithymia)

X <- as.matrix(alexithymia$X)

Y <- as.matrix(alexithymia$Y)

X <- scale(X, T, T)

Y <- scale(Y,T,T)

head(X)
head(Y)

pcovr1 <- pcovr_svd(X = X, Y = Y, R = 4, alpha = 0.15)
pcovr2 <- pcovr_est(X = X, Y = Y, r = 4, a = 0.15, cross = F)

pcovrmixed2 <- pcovrmix(Z = X, M = diag(20), N = diag(rep(1/nrow(X), nrow(X))), Y = Y, R = 4, 
                       alpha = 0.15, qualis = NULL, prop = NULL)


# covariate scores (component scores)
head(pcovr1$Te)
head(pcovrmixed2$Te)

# Px matrix (loadings)
head(t(pcovr1$Px))
head(t(pcovrmixed2$Px))

TuckerCoef(t(pcovr1$Px), t(pcovrmixed2$Px))

# W weights
TuckerCoef(pcovrmixed2$W, pcovr2$W)

# Py regression coefficients
pcovr1$Py
pcovrmixed2$Py



# 4.3. compare the results between pcovrmix and pcovr, when mixed data is used.
data("gironde")
dat <- gironde$housing

dim(dat)

head(dat)

Y <- dat[,1] # set aside the "density" variable as the response variable

Y <- scale(Y, T, T)

# pre-processing the predictors, which are columns 2-4 of dat
processed <- pre_proc(dat[,-1])

head(processed$Z)

# fitting the pcovrmix model
pcovrmixed3 <- pcovrmix(Z = processed$Z, M = processed$M, N = processed$N, Y = Y, 
                        prop = processed$prop, R = 3, alpha = 0.9, qualis = c(2,2))

# for pcovr model, we simply provide the dummy matrix for the qualitative variables
X.quali <- splitmix(dat)$X.quali

dummy <- mltools::one_hot(data.table::as.data.table(X.quali))
dummymatrix <- as.matrix(dummy)

X.quanti <- processed$Z[,1:2]

X_for_pcovr <- cbind(X.quanti, dummymatrix)

pcovr3 <- pcovr_svd(X = X_for_pcovr, Y = Y, R = 3, alpha = 0.9)

# R squared calculation #
sum((Y - pcovrmixed3$Te %*% pcovrmixed3$Py)^2) / sum(Y^2)
sum((Y - pcovr3$Te %*% pcovr3$Py)^2) / sum(Y^2)






# 8. Family data set ####
#formatting the data set
setwd("D:/Uni/Thesis/R-code")
file <- "500family_mixed.txt"
dat <- read.delim(file, header = TRUE, sep = " ", dec = ".")
datstr <- str(dat)

#Make sure that only continuous (numerical) and categorical (Factor) variables are present
dat[sapply(dat, is.character)] <- lapply(dat[sapply(dat, is.character)], as.factor)
dat[sapply(dat, is.integer)] <- lapply(dat[sapply(dat, is.integer)], as.numeric)
dat <- dat[,-1]
str(mid)

#Determining the outcome variable
Y_fam <- dat[,"Academic.performance"]
Y_fam <- scale(Y_fam, center = T, scale = T)

#Excluding the outcome variable from the predictors
numb <- grep("Academic.performance", colnames(dat))
Z <- dat[,-c(1,numb)]
colnames(Z)
Zproc<- pre_proc(Z)

#Model selection:
#1. cross-validation R
#maximum likelihood alpha
x_pc <- prcomp(Zproc$Z)

#scree plot
plot(x_pc$sdev^2, xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     type = "b")

var_Ex <- (sum(x_pc$sdev^2) - sum(x_pc$sdev[1:3]^2)) /  sum(x_pc$sdev^2)
# for estimation of ey, this is the proportion of unexplained variance when y is regressed on X.
reg <- MASS::ginv(t(Zproc$Z) %*% Zproc$Z) %*% t(Zproc$Z) %*% Y_fam

var_ey <- sum((Zproc$Z %*% reg - Y_fam)^2) / sum(Y_fam^2) 
# combining them together:
a_ml_fam <- sum(Zproc$Z^2) / (sum(Zproc$Z^2) + sum(Y_fam^2) * var_Ex / var_ey)  
a_ml_fam

#Determining the optimal number of components (R) to use 
#Note don't use more components than variables present in the R range

R_range <- c(1,2,3,4,5,6,7,8,9,10)

cv_results_R_fam <- matrix(NA, ncol = 3, nrow = length(R_range)) 

for (i in 1:length(R_range)){
  
  # this is where the method is being applied, using the CV function (defined below)
  # you can see that the alpha being used as input is different for each iteration i
  family_cv_R <- scdlogr_cv(predictors = Z, Y = Y_fam, R = R_range[i], 
                       alpha = a_ml_fam,
                       nrFolds = 10, qualis = c(8, 3, 3, 7, 4, 2, 5, 5, 4, 5, 4, 5, 5, 4, 5, 6),
                       seed = i)
  
  # record the CV results on the cv_results_R matrix
  cv_results_R_fam[i,] <- c(family_cv_R$cve, family_cv_R$se, R_range[i])
  
  print(i)
}

#Determine the best model according to the logic of Hastie et al. (2001)
cv_min_index_fam <- which.min(cv_results_R_fam[,1])

serule_fam <- cv_results_R_fam[cv_min_index_fam,1] + cv_results_R_fam[cv_min_index_fam,2]

scd_possible_fam <- cv_results_R_fam[cv_results_R_fam[,1] <= serule_fam,]

if (is.vector(scd_possible_fam)){
  scd_chosen_fam <- scd_possible_fam
} else {
  scd_chosen_fam <- scd_possible_fam[1,]
}
cv_results_R_fam
R_chosen_fam <-matrix(scd_chosen_fam, nrow = 1, ncol = 3)
R_chosen_fam <- cbind(R_chosen_fam, a_ml_fam)
colnames(R_chosen_fam) <- c("CVE", "SE", "R", "alpha")
R_chosen_fam
cv_R_fam <- strtoi(R_chosen_fam[,3])
cv_R_fam


#Model selection:
#2. cross validation of alpha

x_pc <- prcomp(Zproc$Z)
#First determine the alpha using maximum likelihood
var_Ex <- (sum(x_pc$sdev^2) - sum(x_pc$sdev[1:3]^2)) /  sum(x_pc$sdev^2)
# for estimation of ey, this is the proportion of unexplained variance when y is regressed on X.
reg <- MASS::ginv(t(Zproc$Z) %*% Zproc$Z) %*% t(Zproc$Z) %*% Y_fam
var_ey <- sum((Zproc$Z %*% reg - Y_fam)^2) / sum(Y_fam^2) 
# combining them together:
a_ml_fam <- sum(Zproc$Z^2) / (sum(Zproc$Z^2) + sum(Y_fam^2) * var_Ex / var_ey)  
a_ml_fam

#Run PCovRMIX
mix <- pcovrmix(Zproc$Z, Zproc$M, Zproc$N, Y_fam, Zproc$prop, 6, a_ml_fam, c(8, 3, 3, 7, 4, 2, 5, 5, 4, 5, 4, 5, 5, 4, 5, 6))
install.packages("rpart.plot")
library(rpart.plot)
#Plot the eigenvalues on a screeplot to determine R (Using the elbow method)
png(filename = "PCovRMIX_screeplot_fam", width = 750, height = 600)
rpart.plot(mix$singular_values^2, xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     cex.lab = 1.6,
     cex.axis = 1.5,
     type = "b")
dev.off()
rpart.plot()
#Based on the screeplot, 5 components are chosen as R
R_scree_fam <- 5

#Cross-validation for alpha

alpha_range = seq(0,1,0.01)

cv_results_alpha_fam <- matrix(NA, ncol = 3, nrow = length(alpha_range)) 

for (i in 1:length(alpha_range)){
  
  # this is where the method is being applied, using the CV function (defined below)
  # you can see that the alpha being used as input is different for each iteration i
  family_cv_a <- scdlogr_cv(predictors = Z, Y = Y_fam, R = R_scree_fam, 
                       alpha = alpha_range[i],
                       nrFolds = 10, seed = i, qualis = c(8, 3, 3, 7, 4, 2, 5, 5, 4, 5, 4, 5, 5, 4, 5, 6))
  
  # record the CV results on the cv_results matrix
  cv_results_alpha_fam[i,] <- c(family_cv_a$cve, family_cv_a$se, alpha_range[i])
  
  print(i)
}

cv_min_index_fam <- which.min(cv_results_alpha_fam[,1])

serule_fam <- cv_results_alpha_fam[cv_min_index_fam,1] + cv_results_alpha_fam[cv_min_index_fam,2]

scd_possible_fam <- cv_results_alpha_fam[cv_results_alpha_fam[,1] <= serule_fam,]

if (is.vector(scd_possible_fam)){
  alpha_chosen_fam <- scd_possible_fam
} else {
  alpha_chosen_fam <- scd_possible_fam[1,]
}
cv_results_alpha_fam
alpha_chosen_fam <- matrix(alpha_chosen_fam, nrow = 1, ncol = 3)
alpha_chosen_fam_CVE <- matrix(alpha_chosen_fam[,c(1,2)], nrow = 1, ncol = 2)
alpha_chosen_fam_a <- alpha_chosen_fam[,3]
alpha_chosen_fam <- cbind(alpha_chosen_fam_CVE, R_scree_fam, alpha_chosen_fam_a)
alpha_chosen_fam
colnames(alpha_chosen_fam) <- c("CVE", "SE", "R", "alpha")
alpha_chosen_fam
cv_alpha_fam <- cv_results_alpha_fam[cv_min_index_fam,3]
cv_alpha_fam

#Model parameters for:
#Cross-validation R:
  a_ml_fam
  cv_R_fam
  
#Cross-validation for alpha:
  R_scree_fam
  cv_alpha_fam

#Create a table with both cross-validation results to compare the models
CV_table_fam <- rbind(R_chosen_fam, alpha_chosen_fam)
rownames(CV_table_fam) <- c("CV_R", "CV_alpha")
CV_table_fam <- round(CV_table_fam, 5)
CV_table_fam

#Export the CV_table
CV_table_family500 <- as.data.frame(round(CV_table, 6))
write.table(x = CV_table_family500, file = "CV.table family500.txt ", sep = ",", dec = ".")

#The model obtained from alpha cross-validation gives the best out-of-sample prediction power with
#the lowest CVE: CV_alpha

#Run PCovRMIX with the optimal parameters
mix_fam <- pcovrmix(Zproc$Z, Zproc$M, Zproc$N, Y_fam, Zproc$prop, R = R_scree_fam, alpha = cv_alpha_fam, c(8, 3, 3, 7, 4, 2, 5, 5, 4, 5, 4, 5, 5, 4, 5, 6))

#Determine the accuracy of in-sample prediction by calculating R2
R2y_family <- sum((Y_fam - (Zproc$Z %*% mix_fam$W %*% mix_fam$Py))^2) * sum((Y_fam)^2)^-1
R2x_family <- sum((Zproc$Z - (Zproc$Z %*% mix_fam$W %*% mix_fam$Px))^2) * sum((Zproc$Z)^2)^-1

#Showing the Py matrix
colnames(dat[numb])
colnames(Y) <- colnames(dat[numb])
colnames(mix_fam$Py) <- colnames(Y)
round(mix_fam$Py, 3)

#The variance of each component
varcomp <- colSums(mix_fam$Cmat)

#Showing the C matrix with each element being the proportion of variance contributed
#and accounted for variance
C <- t(mix_fam$Cmat)/rowSums(t(mix_fam$Cmat))
CC <- round((C/rowSums(C))*100, 2)
CC <- t(CC)
rownames(CC) <- c(colnames(Z), colnames(Y))
round(CC, 3)



#Exporting Py, CC and the variance of each component
Py_fam <- as.data.frame(round(mix_fam$Py, 3))
C_fam <- as.data.frame(round(CC, 3))
write.table(Py_fam, file = "Py.family500.txt", sep = ",", dec = ".")
write.table(C_fam, file = "C.family500.txt", sep = ",", dec = ".")


str(dat)





# 9. MIDUS data set ####
#Loading in the file
setwd("D:/Uni/Thesis/R-code")
file <- "MIDUS2_mixed.txt"
#Read the MIDUS2 file
mid <- read.delim(file, header = TRUE, sep = " ", dec = ".")
#Turn all character variables into factors
mid[sapply(mid, is.character)] <- lapply(mid[sapply(mid, is.character)], as.factor)
mid[sapply(mid, is.integer)] <- lapply(mid[sapply(mid, is.integer)], as.numeric)
str(mid)
#Prepare the data for pcovrmix
Y_mid <- mid[,2]
Y_mid <- scale(Y_mid, T, T)
predictors <- mid[,-2]

Zproc_mid <- pre_proc(predictors)

#Model selection:
#1. cross validation R
#maximum likelihood alpha
x_pc <- prcomp(Zproc_mid$Z)

var_Ex <- (sum(x_pc$sdev^2) - sum(x_pc$sdev[1:3]^2)) /  sum(x_pc$sdev^2)
# for estimation of ey, this is the proportion of unexplained variance when y is regressed on X.
reg <- MASS::ginv(t(Zproc_mid$Z) %*% Zproc_mid$Z) %*% t(Zproc_mid$Z) %*% Y_mid
var_ey <- sum((Zproc_mid$Z %*% reg - Y_mid)^2) / sum(Y_mid^2) 
# combining them together:
a_ml_mid <- sum(Zproc_mid$Z^2) / (sum(Zproc_mid$Z^2) + sum(Y_mid^2) * var_Ex / var_ey)  
a_ml_mid

#Determining the optimal number of components (R) to use 
#Note don't use more components than variables present in the R range

R_range <- c(1,2,3,4,5,6,7,8,9,10)

cv_results_R_mid <- matrix(NA, ncol = 3, nrow = length(R_range)) 

for (i in 1:length(R_range)){
  
  # this is where the method is being applied, using the CV function (defined below)
  # you can see that the alpha being used as input is different for each iteration i
  mid_cv <- scdlogr_cv(predictors = predictors, Y = Y_mid, R = R_range[i], 
                       alpha = a_ml_mid,
                       nrFolds = 10, qualis = c(12,2,2,2,2,2),
                       seed = i)
  
  # record the CV results on the cv_results_R matrix
  cv_results_R_mid[i,] <- c(mid_cv$cve, mid_cv$se, R_range[i])
  
  print(i)
}

#Determine the best model according to the logic of Hastie et al. (2001)
cv_min_index_mid <- which.min(cv_results_R_mid[,1])

serule_mid <- cv_results_R_mid[cv_min_index_mid,1] + cv_results_R_mid[cv_min_index_mid,2]

scd_possible_mid <- cv_results_R_mid[cv_results_R_mid[,1] <= serule_mid,]

if (is.vector(scd_possible_mid)){
  scd_chosen_mid <- scd_possible_mid
} else {
  scd_chosen_mid <- scd_possible_mid[1,]
}
scd_chosen_mid
cv_results_R_mid
R_chosen_mid <-matrix(scd_chosen_mid, nrow = 1, ncol = 3)
R_chosen_mid <- cbind(R_chosen_mid, a_ml_mid)
colnames(R_chosen_mid) <- c("CVE", "SE", "R", "alpha")
R_chosen_mid
cv_R_mid <- strtoi(R_chosen_mid[,3])
cv_R_mid


#Model selection:
#2. cross validation alpha

x_pc <- prcomp(Zproc_mid$Z)
#First determine the alpha using maximum likelihood
var_Ex <- (sum(x_pc$sdev^2) - sum(x_pc$sdev[1:3]^2)) /  sum(x_pc$sdev^2)
# for estimation of ey, this is the proportion of unexplained variance when y is regressed on X.
reg <- MASS::ginv(t(Zproc_mid$Z) %*% Zproc_mid$Z) %*% t(Zproc_mid$Z) %*% Y_mid
var_ey <- sum((Zproc_mid$Z %*% reg - Y_mid)^2) / sum(Y_mid^2) 
# combining them together:
a_ml_mid <- sum(Zproc_mid$Z^2) / (sum(Zproc_mid$Z^2) + sum(Y_mid^2) * var_Ex / var_ey)  
a_ml_mid

#Run PCovRMIX
mix_mid_a <- pcovrmix(Zproc_mid$Z, Zproc_mid$M, Zproc_mid$N, Y_mid, Zproc_mid$prop, 6, a_ml_mid, c(12,2,2,2,2,2))

#Plot the eigenvalues on a screeplot to determine R (Using the elbow method)
png(filename = "PCovRMIX_screeplot_mid")
plot(mix_mid_a$singular_values^2, xlab = "Principal Component",
     ylab = "Proportion of Variance Explained",
     type = "b")
dev.off()

#Based on the screeplot, 2 components are chosen as R
R_scree_mid <- 2

#Cross-validation for alpha

alpha_range = seq(0,1,0.01)

cv_results_alpha_mid <- matrix(NA, ncol = 3, nrow = length(alpha_range)) 

for (i in 1:length(alpha_range)){
  
  # this is where the method is being applied, using the CV function (defined below)
  # you can see that the alpha being used as input is different for each iteration i
  mid_cv_a <- scdlogr_cv(predictors = predictors, Y = Y_mid, R = R_scree_mid, 
                       alpha = alpha_range[i],
                       nrFolds = 10, seed = i, qualis = c(12,2,2,2,2,2))
  
  # record the CV results on the cv_results matrix
  cv_results_alpha_mid[i,] <- c(mid_cv_a$cve, mid_cv_a$se, alpha_range[i])
  
  print(i)
}

cv_min_index_a_mid <- which.min(cv_results_alpha_mid[,1])

serule_a_mid <- cv_results_alpha_mid[cv_min_index_a_mid,1] + cv_results_alpha_mid[cv_min_index_a_mid,2]

scd_possible_a_mid <- cv_results_alpha_mid[cv_results_alpha_mid[,1] <= serule_a_mid,]

if (is.vector(scd_possible_a_mid)){
  alpha_chosen_mid <- scd_possible_a_mid
} else {
  alpha_chosen_mid <- scd_possible_a_mid[1,]
}
cv_results_alpha_mid
alpha_chosen_mid <- matrix(alpha_chosen_mid, nrow = 1, ncol = 3)
alpha_chosen_mid_CVE <- matrix(alpha_chosen_mid[,c(1,2)], nrow = 1, ncol = 2)
alpha_chosen_mid_a <- alpha_chosen_mid[,3]


alpha_chosen_mid <- cbind(alpha_chosen_mid_CVE, R_scree_mid, alpha_chosen_mid_a)
colnames(alpha_chosen_mid) <- c("CVE", "SE", "R", "alpha")
alpha_chosen_mid
cv_alpha_mid <- cv_results_alpha_mid[cv_min_index_a_mid,3]
cv_alpha_mid

#Model parameters for:
#Cross-validation R:
  a_ml_mid
  cv_R_mid

#Cross-validation for alpha:
  R_scree_mid
  cv_alpha_mid

#Create a table with both cross-validation results to compare the models
  CV_table_mid <- rbind(R_chosen_mid, alpha_chosen_mid)
  rownames(CV_table_mid) <- c("CV R", "CV alpha")
  CV_table_mid

#Export the table
  CV_table__midus2 <- as.data.frame(round(CV_table_mid, 6))
  write.table(x = CV_table__midus2, file = "CV.table midus2.txt ", sep = ",", dec = ".")
  
#The model obtained from alpha cross-validation gives the best out-of-sample prediction power with
#the lowest CVE. However, with the SE included, the CVEs overlap, which makes it difficult to determine
#which model predicts better. Following Hastie et al. (2001) one should go with the more parsimonious model.
#Which is the CV R model.

#Perform pcovrmix
mix_mid <- pcovrmix(Z = Zproc_mid$Z, M = Zproc_mid$M, N = Zproc_mid$N, Y = Y_mid, prop = Zproc_mid$prop, 
                    R = cv_R_mid, alpha = a_ml_mid, qualis = c(12, 2, 2, 2, 2, 2))

#Determine the in-sample accuracy
R2y_mid <- sum((Y_mid - (Zproc_mid$Z %*% mix_mid$W %*% mix_mid$Py))^2) * sum((Y_mid)^2)^-1
R2x_mid <- sum((Zproc_mid$Z - (Zproc_mid$Z %*% mix_mid$W %*% mix_mid$Px))^2) * sum((Zproc_mid$Z)^2)^-1

#Showing the Py matrix
colnames(Y_mid) <- colnames(mid[2])
colnames(mix_mid$Py) <- colnames(Y_mid)
round(mix_mid$Py, 3)

#The variance of each component
varcomp_mid <- colSums(mix_mid$Cmat)

#Showing the C matrix with each element being the proportion of variance contributed
#and accounted for variance
C_mid <- t(mix_mid$Cmat)/rowSums(t(mix_mid$Cmat))
CC_mid <- round((C_mid/rowSums(C_mid))*100, 2)
CC_mid <- t(CC_mid)
rownames(CC_mid) <- c(colnames(predictors), colnames(Y_mid))
round(CC_mid, 3)

#Export the relevant tables
Py_mid <- as.data.frame(round(mix_mid$Py, 3))
C_mid <- as.data.frame(round(CC_mid, 3))
write.table(Py_mid, file = "Py.MIDUS2.txt", sep = ",", dec = ".")
write.table(C_mid, file = "C.MIDUS2.txt", sep = ",", dec = ".")


# 10. cross-validation PCovR ####

#CV for pcovr will use the same parameter values for R and alpha as PCovRMIX

#pcovrcv <- pcovr_cv(X = X, Y = Y, R = R, 
                    #alpha = alpha,
                    #nrFolds = 10, seed = R)

# record the CV results on the cv_results matrix
#pcovr_cv_results <- c(pcovrcv$cve, pcovrcv$se, R, alpha)

# 2. CV function
pcovr_cv <- function(X, Y, R, alpha, nrFolds, seed){
  
  # set seed, because we would randomly divide the dataset into folds
  set.seed(seed)
  
  # randomly shuffle the data
  sampled <- sample(nrow(X))
  
  X <- X[sampled,]
  
  Y <- Y[sampled,]
  
  # Create equally sized folds
  folds <- cut(seq(1,nrow(X)),breaks=nrFolds,labels=FALSE)
  
  cve_k <- data.frame(error_y = NA)
  
  for(k in 1:nrFolds){
    
    # defining test and train data  
    test_index <- which(folds==k, arr.ind=TRUE)
    
    X_train <- X[-test_index,]
    X_test <- X[test_index, ]
    
    Y_train <- Y[-test_index]
    Y_test <- Y[test_index]
    
    # model fitting #
    org <- pcovr_est(X = X_train, Y = Y_train, 
                     r = R, a = alpha, cross = F)
    # out of sample prediction #
    
    # (here, the criterion being used is the sum of squared errors)
    # (you can simply compute the sum of squared errors between the predicted Y and the observed Y)
    
    pred <- X_test %*% org$W %*% org$Py
    
    SSE <- sum((Y_test - pred)^2) * sum((Y_test)^2)^-1
    
    cve_k[k,1] <- SSE
    
  }
  
  cve <- colMeans(cve_k)
  # cve = mean cv error
  
  se <- apply(cve_k,2,sd) / sqrt(nrFolds)
  # se = standard deviation of the mean cv errors per k / sqrt(nrFold)
  
  return(list(cve = cve, se = se, cve_k = cve_k))
}

pcovr_cv_results
summary(P)

# 11. Comparison 1: PCAmix ####
#Here PCAMIX is performed and the results regressed on the criterion
#The in-sample prediction accuracy will be measured and compared with PCovRMIx results

#compare the Family_500 dataset
split <- splitmix(dat)
X1 <- split$X.quanti
X2 <- split$X.quali
res.pcamix <- PCAmix(X1, X2, rename.level = T)
res.pcamix$eig
plot(res.pcamix$eig[,1], xlab = "Principle Components",
      ylab = "Proportion of variance explained",
     type = "b")

#Based on the PCovRMIX model, the first 5 components are chosen.

PCAMIX_C <- t(res.pcamix$sqload[,c(1:5)])/rowSums(t(res.pcamix$sqload[,c(1:5)]))
PCAMIX_CC <- round((PCAMIX_C/rowSums(PCAMIX_C))*100, 2)
PCAMIX_CC <-t(PCAMIX_CC)
round(PCAMIX_CC, 3)
round(CC, 3)

#Variance of each component
#PCAmix
colSums(res.pcamix$sqload)
#PCovRMIX
colSums(mix_fam$Cmat)

#Scale the factor coordinates of the rows back to component scores according to Chavent et al. (2017)
PCAmix_svalues <- (sqrt(res.pcamix$eig[c(1:5),1])^-1)
PCAmix_T <- res.pcamix$ind$coord[,c(1:5)] %*% diag(PCAmix_svalues)
PCAmix_T
PCAmix_svalues
#Regress these components on the criterion "Academic.performance"
PCAmix_df <- as.data.frame(PCAmix_T)
PCAmix_reg_y_fam <- lm(dat[,24] ~ PCAmix_T, data = dat)
PCAmix_reg_x_fam <- lm(Zproc$Z ~ PCAmix_T, data = dat)

#Comparison of the R2y of the PCovRMIX model and the PCAMIx model
R2y_family
R2y_PCAmix_fam <- summary(PCAmix_reg_y_fam)$adj.r.squared
R2y_PCAmix_fam
#Compared of the R2x of the PCovRMIX model and the PCAMIX model
R2x_family
summary(PCAmix_reg_x_fam)

#Convert all relevant measures in a table
PCAmix_R2y_fam_table <- rbind(R2y_family, R2y_PCAmix_fam)
rownames(PCAmix_R2y_fam_table) <- c("PCovRMIX", "PCAmix")
colnames(PCAmix_R2y_fam_table) <- "R2y"
PCAmix_R2y_fam_table

#Export all relevant data
PCAmix_R2y_comp_fam_table <- as.data.frame(round(PCAmix_R2y_fam_table, 3))
write.table(PCAmix_R2y_comp_fam_table, file = "PCAmix.R2y.Family500.txt", sep = ",", dec = ".")





#Compare the MIDUS_2 dataset
split <- splitmix(mid)
X1 <- split$X.quanti
X2 <- split$X.quali
res.pcamix <- PCAmix(X1, X2, rename.level = T)
res.pcamix$eig
plot(res.pcamix$eig[,1], xlab = "Principle Components",
     ylab = "Proportion of variance explained",
     type = "b")

#Based on the PCovRMIX model, the first component is chosen.
PCAMIX_C <- t(res.pcamix$sqload[,1])/rowSums(t(res.pcamix$sqload[,1]))
PCAMIX_CC <- round((PCAMIX_C/rowSums(PCAMIX_C))*100, 2)
PCAMIX_CC <-t(PCAMIX_CC)
round(PCAMIX_CC, 3)

#Variance of each component
sum(res.pcamix$sqload[,1])
colSums(mix_mid$Cmat)
res.pcamix$eig[,1]
res.pcamix$ind$coord[,1]
#Scale the factor coordinates of the rows back to component scores
PCAmix_svalues <- (sqrt(res.pcamix$eig[1,1])^-1)
PCAmix_T <- res.pcamix$ind$coord[,1] * PCAmix_svalues
PCAmix_T
#Regress these components on the criterion "Academic.performance"
PCAmix_df <- as.data.frame(PCAmix_T)
PCAmix_reg_y_mid <- lm(mid[,2] ~ PCAmix_T, data = mid)
PCAmix_reg_x_mid <- lm(Zproc_mid$Z ~ PCAmix_T, data = mid)

#Comparison of the R2y of the PCovRMIX model and the PCAMIx model
R2y_mid
R2y_PCAmix_mid <- summary(PCAmix_reg_y_mid)$adj.r.squared
R2y_PCAmix_mid

#Compared of the R2x of the PCovRMIX model and the PCAMIX model
R2x_mid
summary(PCAmix_reg_x)

#Convert all relevant measures in a table
PCAmix_R2y_mid_table <- rbind(R2y_mid, R2y_PCAmix_mid)
rownames(PCAmix_R2y_mid_table) <- c("PCovRMIX", "PCAmix")
colnames(PCAmix_R2y_mid_table) <- "R2Y"
PCAmix_R2y_mid_table

#Export all relevant data
PCAmix_R2y_comp_mid_table <- as.data.frame(round(PCAmix_R2y_mid_table, 3))
write.table(PCAmix_R2y_comp_mid_table, file = "PCAmix.R2y.MIDUS2.txt", sep = ",", dec = ".")

# 12. Comparison 2: PcovR ####

install.packages("PCovR")
library(PCovR)


#Family_500 dataset
P_fam <- pcovr_est(X = Zproc$Z, Y = Y_fam, r = 5, a = cv_alpha_fam, cross = T, fold = 10)

#The in-sample prediction accuracy compared
#PCovR
P_fam$Ry2
#PCovRMIX
R2y_family

#out-sample prediction error compared
#PCovR
R_fam = 5
alpha = cv_alpha_fam

pcovrcv_fam <- pcovr_cv(X = Zproc$Z, Y = Y_fam, R = R_fam, 
                    alpha = alpha,
                    nrFolds = 10, seed = R_fam)

# record the CV results on the cv_results matrix
pcovr_cv_results_fam <- c(pcovrcv_fam$cve, pcovrcv_fam$se, R_fam, alpha)
pcovr_cv_results_fam[1]

#PCovRMIX
alpha_chosen_fam[1,1]

#Regression coefficients compared
P_fam$Py
mix_fam$Py

#Store all relevant values in tables
#in-sample prediction 
PCovR_Ry2_fam_table <- rbind(R2y_family, P_fam$Ry2)
rownames(PCovR_Ry2_fam_table) <- c("PCovRMIX", "PCovR")
colnames(PCovR_Ry2_fam_table) <- "Qy2"
PCovR_Ry2_fam_table

#out-sample prediction
PCovR_CVE_fam_table <- rbind(alpha_chosen_fam[1,1], pcovr_cv_results_fam[1])
rownames(PCovR_CVE_fam_table) <- c("PCovRMIX", "PCovR")
PCovR_CVE_fam_table

#Export the tables
PCovR_Ry2_comp_fam_table <- as.data.frame(round(PCovR_Ry2_fam_table, 3))
PCovR_CVE_comp_fam_table <- as.data.frame(round(PCovR_CVE_fam_table, 3))
write.table(PCovR_Ry2_comp_fam_table, file = "PCovR.R2y.Family500.txt", sep = ",", dec = ".")
write.table(PCovR_CVE_comp_fam_table, file = "PCovR.CVE.Family500.txt", sep = ",", dec = ".")



#MIDUS_2 dataset
P_mid <- pcovr_est(X = Zproc_mid$Z, Y = Y_mid, r = 1, a = a_ml_mid, cross = T, fold = 10)
summary(P_mid)

#in-sample prediction accuracy:
#PCovR
P_mid$Ry2

#PCovRMIX
R2y_mid

#out-sample prediction accuracy
#PCovR
R_mid <- 1
pcovrcv_mid <- pcovr_cv(X = Zproc$Z, Y = Y_mid, R = R_mid, 
                    alpha = a_ml_mid,
                    nrFolds = 10, seed = R_mid)

pcovr_cv_mid_results <- c(pcovrcv_mid$cve, pcovrcv_mid$se, R_mid, a_ml_mid)
pcovr_cv_mid_results[1]

#PCovRMIX
R_chosen_mid[1,1]

#regression coefficients compared
P_mid$Py
mix_mid$Py

#Store all relevant values in tables
#in-sample prediction 
PCovR_Ry2_mid_table <- rbind(R2y_mid, P_mid$Ry2)
rownames(PCovR_Ry2_mid_table) <- c("PCovRMIX", "PCovR")
colnames(PCovR_Ry2_mid_table) <- "Qy2"
PCovR_Ry2_mid_table

#out-sample prediction
PCovR_CVE_mid_table <- rbind(R_chosen_mid[1,1], pcovr_cv_mid_results[1])
rownames(PCovR_CVE_mid_table) <- c("PCovRMIX", "PCovR")
PCovR_CVE_mid_table

#Export the tables
PCovR_Ry2_comp_mid_table <- as.data.frame(round(PCovR_Ry2_mid_table, 3))
PCovR_CVE_comp_mid_table <- as.data.frame(round(PCovR_CVE_mid_table, 3))
write.table(PCovR_Ry2_comp_mid_table, file = "PCovR.R2y.MIDUS2.txt", sep = ",", dec = ".")
write.table(PCovR_CVE_comp_mid_table, file = "PCovR.CVE.MIDUS2.txt", sep = ",", dec = ".")
