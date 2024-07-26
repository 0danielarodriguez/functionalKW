rm(list=ls())
library(fdANOVA)
library(fda.usc)
#setwd("~/Dropbox/testnopfunc/GITHUB")
source("functionsKW.R")

# # # # # # # # # # # # # # # # # # # # # # # # # # # 
# funciones de la base de fourier
phi <- function(k, t){
  sqrt(2) * sin((k-0.5) * pi * t) / ((k-0.5) * pi)
}


# N es el tamaño hasta donde se trunca la serie
N<-1000
# los tt son fijos equiespaciados
lt <- 250
tes <- seq(0, 1, length = lt)
kk<-1:N
matrixphi<-outer(kk,tes,phi)  

GeneroZ<-function(ZZ,matrixphi,mu)
{apply(matrixphi*ZZ,2,sum)+mu}


#NORMAL
ModNormal<-function(nsamp,matrixphi,mu)
{
  N<-dim(matrixphi)[1]
  Z <- rnorm(N)
  X<-GeneroZ(Z,matrixphi,mu)
  for(i in 2:nsamp){
    Z <- rnorm(N)
    X <- rbind(X, GeneroZ(Z,matrixphi,mu) )
  }
  X
}


#T

ModT<-function(nsamp,df,matrixphi,mu)
{
  N<-dim(matrixphi)[1]
  Z <- rt(N,df)
  X<-GeneroZ(Z,matrixphi,mu)
  for(i in 2:nsamp){
    Z <- rt(N,df)
    X <- rbind(X, GeneroZ(Z,matrixphi,mu)) 
  }
  X
}

# Eliptica  t1^2 *normal 
ModEliptica<-function(nsamp,matrixphi,mu)
{
  N<-dim(matrixphi)[1] 
  Z <- rnorm(N)*rt(1,1)^2
  X<-GeneroZ(Z,matrixphi,mu)
  for(i in 2:nsamp){
    Z <- rnorm(N)*rt(1,1)^2
    X <- rbind(X, GeneroZ(Z,matrixphi,mu) )
  }
  X
}



nrep <- 1000
ene1 <- 15
ene2 <- 15
ene3 <- 15

#grafico de la funcion Delta3
plot(tes*(1-tes))

ces <- seq(0.1, 0.6, 0.1) #excluyo la nula 
delta0<-rep(0,lt)
delta0<-t(as.matrix(delta0))
#delta1 <- outer(ces, rep(1,lt))
#delta2 <- outer(ces, tes)
delta3 <- outer(ces, sqrt(tes)) #en cada fila hay una alternativa con otra c. 


media0<-delta0
medias<-rbind(delta3) #cada fila tiene una media distinta
cces<-c(ces) #son los valores de c para cada delta, se usa para el archivo
lmedias <- dim(medias)[1] #usamos delta3 con diferentes valores de c.

proyecciones<-c(10,30,100)
nproy <- length(proyecciones)

# esque mas de simulacion
#normal
#t1
#E
#cada media 0, y delta3 cons 6 ces distintos.



#eliptica


for(mm in 1:lmedias)
{
  delta<-3
  cc<-cces[mm]
  media<-medias[mm,]
  for(i in 1:nrep){
    print(i)
    set.seed(i)  
    X1 <- ModEliptica(ene1, matrixphi,media0)
    X2 <- ModEliptica(ene2, matrixphi,media)
    X3 <- ModEliptica(ene3, matrixphi,media)
    X <- rbind(X1, X2, X3)
    g <- as.factor(c(rep(1, ene1), rep(2, ene2), rep(3, ene3)))
    pkw <- c()
    for(pp in 1:nproy)
    {
      proy<-proyecciones[pp]
      p <- func.kw.test(X, g, tes,nn=proy)$p.value
      pkw <- c(pkw,p)
    }
    
    paov <- c()
    for(pp in 1:nproy)
    {
      proy<-proyecciones[pp]
      p <- func.aov.test(X, g, tes,nn=proy)$p.value
      paov <- c(paov,p)
    }
    
    p03 <- func.wall.test(X, g, tes)$p.value
    
    
    p04todos <- fanova.tests(x = t(X), group.label = g, test = "ALL")
    cc10<-cc*10
    pvalores<-c(i,pkw,paov,p03,sapply(1:12, function(i)p04todos[[i]][2]) ,  p04todos[[12]][3])
    write.table(pvalores, file = paste0("salidasR1/E1_eliptica","c_",cc10,"delta",delta,".txt"),  col.names = FALSE, row.names = FALSE, append = TRUE)
  }
}
