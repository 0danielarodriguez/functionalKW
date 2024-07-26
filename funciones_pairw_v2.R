#' Computes the p-value of a kruskal-wallis test based on random projections 
#' @param X A set of functional data evaluated on a t 
#' @param g A categorical variable indicating the group of each observation
#' @param t The t in which the functions are evaluated 
#' @examples
#' x <- rnorm(8, 0, 1)
#' y <- rt(6, 1)
#' z <- rnorm(10, 4, 1)
#' xx_coef <- c(x, y, z)
#' g <- c(rep(1, 8), rep(2, 6), rep(3, 10))
#' pairwise.wilcox.test(xx_coef, g)
#' lt <- 250
#' t <- seq(0, 1, length = lt)
#' ene1 <-  ene2 <- ene3 <- 8
#' library(fda.usc)
#' X1 <- rproc2fdata(ene1, t, mu = rep(0, lt), sigma = 1)$data
#' matplot(t(X1), type = "l", col = 3)
#' X2 <- rproc2fdata(ene2, t, mu = rep(2, lt), sigma = "OU")$data
#' matplot(t(X2), type = "l", col = 4, add = TRUE)
#' X3 <- rproc2fdata(ene3, t, mu = rep(0, lt), sigma = "vexponential")$data
#' matplot(t(X3), type = "l", col = 5, add = TRUE)
#' X <- rbind(X1, X2, X3)
#' g <- as.factor(c(rep(1, ene1), rep(2, ene2), rep(3, ene3)))
#' func.kw.test(X, g, t, nn = 30)
#' @return The p-value
#' @export
################ pairwise comparison da muy distinto de la version antigua en el ejemplo ################
func.kw.test <- function(X, g, t, nn = 30){
tablakruskal <- list()
pvalorespairw <- list()
pvalorkruskal <- 0
lt <- length(t)
i <- 1
for( i in 1:nn){
  elembr <-rproc2fdata(1, t, sigma = "brownian")
  delta <- 1/(lt-1)
  xx_coef <- X %*% matrix(elembr$data) * delta
  #matplot(t(matrix(as.numeric(as.matrix(tempmed3)), nrow = 71, ncol = 12)), type = "l")
  #boxplot(tempproy ~ region, data = datos)
  tablakruskal[[i]] <- kruskal.test(xx_coef, g)
  pvalorkruskal[i] <- tablakruskal[[i]]$p.value
  pvalorespairw[[i]] <- pairwise.wilcox.test(xx_coef, g, p.adjust.method= "BY" )$p.value
  nrowpairw <- nrow(pvalorespairw[[i]])
}
adjpvaluesasnumeric <- p.adjust(unlist(pvalorespairw), method = "BY")
K <- length(unique(g))
adjpvaluesaslist <- list()
for(i in 1:nn){
  desde <- (K-1) * (K-1) * (i-1) + 1
  hasta <- (K-1) * (K-1) *i
  adjpvaluesaslist[[i]] <- matrix(adjpvaluesasnumeric[desde:hasta], nrow = K-1)
}
####### chequeo con los pvalores sin ajustar si esta bien la forma de armar y desarmar las matrices##############
#pvaluesasnumeric <- unlist(pvalorespairw) ## este paso es necesario para usar p.adjust
#pvaluesaslist <- list()
#for(i in 1:nn){
#  desde <- (K-1) * (K-1) * (i-1) + 1
#  hasta <- (K-1) * (K-1) * i
#  pvaluesaslist[[i]] <- matrix(pvaluesasnumeric[desde:hasta], nrow = K-1)
#}
#####
### busco el minimo p.adjust para cada coordenada de la matriz de pvalores ajustados
elementos<-function(x){x[lower.tri(x,diag=TRUE)]}
pvalorespwdeck <- lapply(pvalorespairw, elementos)
pvalorespwdeckadj <- lapply(pvalorespwdeck, p.adjust, method = "BY")
minps <- Reduce(pmin, pvalorespwdeckadj)
x<-pvalorespairw[[1]]
rownames(x) <- rownames(pvalorespairw[[1]])
colnames(x) <- colnames(pvalorespairw[[1]])
lowerTriangle(x, diag=TRUE, byrow=FALSE) <- minps
#K/(1:K) * sort(pvalorkruskal) - sort(p.adjust(pvalorkruskal, method = "BY"))
#min(K/(1:K) * sort(pvalorkruskal)) - min(sort(p.adjust(pvalorkruskal, method = "BY")))
ans <- list(p.value = min(p.adjust(pvalorkruskal, method = "BY")),
            p.value.pairw = x )
ans
}


#' Computes adjusted for multiple comparisons p-values of wilcoxon tests based on random projections 
#' @param X A set of functional data evaluated on a t 
#' @param g A categorical variable indicating the group of each observation
#' @param t The t in which the functions are evaluated 
#' @examples
#' lt <- 250
#' tes <- seq(0, 1, length = lt)
#' ene1 <-  ene2 <- ene3 <- 8
#' X1 <- rproc2fdata(ene1, t = tes, mu = rep(0, lt), sigma = 1)$data
#' matplot(t(X1), type = "l", col = 3)
#' X2 <- rproc2fdata(ene2, t = tes, mu = rep(2, lt), sigma = "OU")$data
#' matplot(t(X2), type = "l", col = 4, add = TRUE)
#' X3 <- rproc2fdata(ene3, t = tes, mu = rep(0, lt), sigma = "vexponential")$data
#' matplot(t(X3), type = "l", col = 5, add = TRUE)
#' X <- rbind(X1, X2, X3)
#' g <- as.factor(c(rep(1, ene1), rep(2, ene2), rep(3, ene3)))
#' func.pairwise.wilcox.test(X, g, tes, nn = 30)
#' @return Adjusted p-values for pairwise comparisons
#' @export

func.pairwise.wilcox.test <- function(X, g, t, method = "BY", nn = 30){
  grupos <- unique(g)
  K <- length(grupos)
  pvalores <- matrix(NA, nrow = K, ncol = K)
  i <- 1; j <- 2
  for(i in 2:K) {
    for (j in 1:(i-1)) {
      Xij <- X[g == grupos[i] | g == grupos[j],]
      gij <- g[g == grupos[i] | g == grupos[j]]
      pvalores[i, j] <- func.kw.test(Xij, gij, t, nn)$p.value
    }
  }
  ans <- matrix(p.adjust(pvalores, method = method), nrow = K, ncol = K)
  rownames(ans) <- grupos
  colnames(ans) <- grupos
  ans
}


#' Computes the p-value of an anova test based on random projections 
#' @param X A set of functional data evaluated on a t 
#' @param g A categorical variable indicating the group of each observation
#' @param t The t in which the functions are evaluated 
#' @examples
#' @return The p-value
#' @export

func.aov.test <- function(X, g, t, nn = 30){
  tabla <- list()
  pvalorespairw <- list()
  pvalor <- 0
  lt <- length(t)
  i <- 1
  for( i in 1:nn){
    elembr <-rproc2fdata(1, t, sigma = "brownian")
    delta <- 1/(lt-1)
    xx_coef <- X %*% matrix(elembr$data) * delta
    dataaux <- data.frame(xx_coef = xx_coef, g = as.factor(g))
    tabla[[i]] <- aov(xx_coef ~ g, data = dataaux)
    pvalor[i] <- summary(tabla[[i]])[[1]][1,5]
    pvalorespairw[[i]] <- pairwise.t.test(xx_coef, g, p.adjust.method= "BY" )$p.value
    nrowpairw <- nrow(pvalorespairw[[i]])
  }
  adjpvaluesasnumeric <- p.adjust(unlist(pvalorespairw), method = "BY")
  K <- length(unique(g))
  adjpvaluesaslist <- list()
  for(i in 1:nn){
    desde <- (K-1) * (K-1) * (i-1) + 1
    hasta <- (K-1) * (K-1) *i
    adjpvaluesaslist[[i]] <- matrix(adjpvaluesasnumeric[desde:hasta], nrow = K-1)
  }
  ####### chequeo con los pvalores sin ajustar si esta bien la forma de armar y desarmar las matrices##############
  #pvaluesasnumeric <- unlist(pvalorespairw) ## este paso es necesario para usar p.adjust
  #pvaluesaslist <- list()
  #for(i in 1:nn){
  #  desde <- (K-1) * (K-1) * (i-1) + 1
  #  hasta <- (K-1) * (K-1) * i
  #  pvaluesaslist[[i]] <- matrix(pvaluesasnumeric[desde:hasta], nrow = K-1)
  #}
  #####
  ### busco el minimo p.adjust para cada coordenada de la matriz de pvalores ajustados
  elementos<-function(x){x[lower.tri(x,diag=TRUE)]}
  pvalorespwdeck <- lapply(pvalorespairw, elementos)
  pvalorespwdeckadj <- lapply(pvalorespwdeck, p.adjust, method = "BY")
  minps <- Reduce(pmin, pvalorespwdeckadj)
  x<-pvalorespairw[[1]]
  rownames(x) <- rownames(pvalorespairw[[1]])
  colnames(x) <- colnames(pvalorespairw[[1]])
  lowerTriangle(x, diag=TRUE, byrow=FALSE) <- minps
  #K/(1:K) * sort(pvalorkruskal) - sort(p.adjust(pvalorkruskal, method = "BY"))
  #min(K/(1:K) * sort(pvalorkruskal)) - min(sort(p.adjust(pvalorkruskal, method = "BY")))
  ans <- list(p.value = min(p.adjust(pvalor, method = "BY")),
              p.value.pairw = x )
  ans
}

#' Computes adjusted for multiple comparisons p-values based on random projections 
#' @param X A set of functional data evaluated on a t 
#' @param g A categorical variable indicating the group of each observation
#' @param t The t in which the functions are evaluated 
#' @examples
#'    t <- seq(0, 1, length = 100) 
#'    n <- 20
#'    m <- 20
#'    X1 <- rproc2fdata(n, t, sigma = "brownian")$data
#'    X2 <- rproc2fdata(m, t, rep(0.3, length(t)), sigma = "brownian")$data
#'    X3 <- rproc2fdata(m, t, rep(0.3, length(t)), sigma = "brownian")$data
#'    X <- rbind(X1, X2, X3)
#'    g <- c(rep("a", n), rep("b", m), rep("c", m))
#' @return 
#' @export

func.pairwise.test <- function(X, g, t, method = "BY", nn = 30){
  grupos <- unique(g)
  K <- length(grupos)
  pvalores <- matrix(NA, nrow = K, ncol = K)
  i <- 1; j <- 2
    for(i in 2:K) {
      for (j in 1:(i-1)) {
        Xij <- X[g == grupos[i] | g == grupos[j],]
        gij <- g[g == grupos[i] | g == grupos[j]]
        pvalores[i, j] <- func.aov.test(Xij, gij, t, nn)$p.value
      }
    }
  ans <- matrix(p.adjust(pvalores, method = method), nrow = K, ncol =K)
  rownames(ans) <- grupos
  colnames(ans) <- grupos
  ans
}




#' Computes the p-value of a Wallace test based on random projections 
#' @param X A set of functional data evaluated on a t 
#' @param g A categorical variable indicating the group of each observation
#' @param t The t in which the functions are evaluated 
#' @examples
#' @return The p-value
#' @export

func.wall.test <- function(X, g, t, nn = 30){
  tabla <- list()
  pvalor <- 0
  lt <- length(t)
  i <- 1
  K <- length(unique(g))
  N <- nrow(X)
  for( i in 1:nn){
    elembr <-rproc2fdata(1, t, sigma = "brownian")
    delta <- 1/(lt-1)
    xx_coef <- X %*% matrix(elembr$data) * delta
    #matplot(t(matrix(as.numeric(as.matrix(tempmed3)), nrow = 71, ncol = 12)), type = "l")
    #boxplot(tempproy ~ region, data = datos)
    H <- kruskal.test(xx_coef, g)$statistic
    FR <- ((N - K) * H) / ((K - 1) * (N - 1 - H)) 
    pvalor[i] <- 1 - pchisq((K-1) * FR, K-1)
}
  ans <- list(p.value = min(p.adjust(pvalor, method = "BY")))
  ans
}

#' Computes the p-value of a Wallace test based on random projections for a large number of groups
#' @param X A set of functional data evaluated on a t 
#' @param g A categorical variable indicating the group of each observation
#' @param t The t in which the functions are evaluated 
#' @examples
#' @return The p-value
#' @export

func.wall.test.large.K <- function(X, g, t, nn = 30){
  tabla <- list()
  pvalor <- 0
  lt <- length(t)
  i <- 1
  K <- length(unique(g))
  N <- nrow(X)
  n <- N/K
  for( i in 1:nn){
    elembr <-rproc2fdata(1, t, sigma = "brownian")
    delta <- 1/(lt-1)
    xx_coef <- X %*% matrix(elembr$data) * delta
    #matplot(t(matrix(as.numeric(as.matrix(tempmed3)), nrow = 71, ncol = 12)), type = "l")
    #boxplot(tempproy ~ region, data = datos)
    H <- kruskal.test(xx_coef, g)$statistic
    FR <-
      sqrt(K/(2 * n/(n-1))) * 
      (((N - K) * H) / ((K - 1) * (N - 1 - H)) - 1)
    pvalor[i] <- 2 * (1 - pnorm(abs(FR)))
  }
  ans <- list(p.value = min(p.adjust(pvalor, method = "BY")))
  ans
}

#' Computes the p-value of an anova test based on random projections for a large number of groups
#' @param X A set of functional data evaluated on a t 
#' @param g A categorical variable indicating the group of each observation
#' @param t The t in which the functions are evaluated 
#' @examples
#' @return The p-value
#' @export

func.aov.test.large.K <- function(X, g, t, nn = 30){
  tabla <- list()
  K <- length(unique(g))
  pvalor <- 0
  i <- 1
  lt <- length(t)
  N <- nrow(X)
  n <- N/K
  for( i in 1:nn){
    elembr <-rproc2fdata(1, t, sigma = "brownian")
    delta <- 1/(lt-1)
    xx_coef <- X %*% matrix(elembr$data) * delta
    #matplot(t(matrix(as.numeric(as.matrix(tempmed3)), nrow = 71, ncol = 12)), type = "l")
    #boxplot(tempproy ~ region, data = datos)
    dataaux <- data.frame(xx_coef = xx_coef, g = as.factor(g))
    aovfit <- aov(xx_coef ~ g, data = dataaux)
    Z <-
      sqrt(K) * (summary(aovfit)[[1]][1, 4] - 1) / sqrt(2 * n / (n - 1))
    pvalor[i] <- 2 * (1 - pnorm(abs(Z))) 
  }
  ans <- list(p.value = min(p.adjust(pvalor, method = "BY")))
  ans
}

#' Computes the p-value of an kruskal wallis test for functional data based on random projections for a large number of groups
#' @param X A set of functional data evaluated on a t 
#' @param g A categorical variable indicating the group of each observation
#' @param t The t in which the functions are evaluated 
#' @examples
#' @return The p-value
#' @export

func.kw.test.large.K <- function(X, g, t, nn = 30){
  tabla <- list()
  pvalor <- 0
  lt <- length(t)
  i <- 1
  K <- length(unique(g))
  N <- nrow(X)
  n <- N/K
  for( i in 1:nn){
    elembr <-rproc2fdata(1, t, sigma = "brownian")
    delta <- 1/(lt-1)
    xx_coef <- X %*% matrix(elembr$data) * delta
    #matplot(t(matrix(as.numeric(as.matrix(tempmed3)), nrow = 71, ncol = 12)), type = "l")
    #boxplot(tempproy ~ region, data = datos)
    H <- kruskal.test(xx_coef, g)$statistic
    FKW <- (H/(K-1) - 1) / sqrt(2 * (n-1)/(K*n))
    pvalor[i] <- 2 * (1 - pnorm(abs(FKW))) 
  }
  ans <- list(p.value = min(p.adjust(pvalor, method = "BY")))
  ans
}



########### con la correccion de harrar y gupta #################################

#' Computes the p-value of a Wallace test based on random projections for a large number of groups as in Harrar and Gupta (2007)
#' @param X A set of functional data evaluated on a t 
#' @param g A categorical variable indicating the group of each observation
#' @param t The t in which the functions are evaluated 
#' @example 
#' @return The p-value both for balanced and unbalanced case
#' @export

func.wallHG.test.large.K <- function(X, g, t, nn = 30){
  tabla <- list()
  pvalormismon <- pvalordistinton <-  0
  lt <- length(t)
  i <- 1
  K <- length(unique(g))
  N <- nrow(X)
  n <- N/K
  g <- as.factor(g)
  for( i in 1:nn){
    elembr <-rproc2fdata(1, t, sigma = "brownian")
    delta <- 1/(lt-1)
    xx_coef <- X %*% matrix(elembr$data) * delta
    #matplot(t(matrix(as.numeric(as.matrix(tempmed3)), nrow = 71, ncol = 12)), type = "l")
    #boxplot(tempproy ~ region, data = datos)
    r <- rank(xx_coef)
    r <- (r-mean(r))/sd(r)
    FRmismon <- Estharrar(r, g)
    FRdistinton <- EstharrarNoBaln(r, g) 
    pvalormismon[i] <- 2 * (1 - pnorm(abs(FRmismon)))
    pvalordistinton[i] <- 2 * (1 - pnorm(abs(FRdistinton)))
  }
  ans <- list(p.value.balanced = min(p.adjust(pvalormismon, method = "BY")),
              p.value.unbalanced = min(p.adjust(pvalordistinton, method = "BY")))
  ans
}

#' Computes the p-value of an anova test based on random projections for a large number of groups  as in Harrar and Gupta (2007)
#' @param X A set of functional data evaluated on a t 
#' @param g A categorical variable indicating the group of each observation
#' @param t The t in which the functions are evaluated 
#' @examples
#' @return The p-value both for balanced and unbalanced case
#' @export

func.aovHG.test.large.K <- function(X, g, t, nn = 30){
  tabla <- list()
  pvalormismon <- pvalordistinton <-  0
  lt <- length(t)
  i <- 1
  K <- length(unique(g))
  N <- nrow(X)
  n <- N/K
  g <- as.factor(g)
  for( i in 1:nn){
    elembr <-rproc2fdata(1, t, sigma = "brownian")
    delta <- 1/(lt-1)
    xx_coef <- X %*% matrix(elembr$data) * delta
    #matplot(t(matrix(as.numeric(as.matrix(tempmed3)), nrow = 71, ncol = 12)), type = "l")
    #boxplot(tempproy ~ region, data = datos)
    FRmismon <- Estharrar(xx_coef, g)
    FRdistinton <- EstharrarNoBaln(xx_coef, g) 
    pvalormismon[i] <- 2 * (1 - pnorm(abs(FRmismon)))
    pvalordistinton[i] <- 2 * (1 - pnorm(abs(FRdistinton)))
  }
  ans <- list(p.value.balanced = min(p.adjust(pvalormismon, method = "BY")),
              p.value.unbalanced = min(p.adjust(pvalordistinton, method = "BY")))
  ans
}
