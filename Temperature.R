rm(list=ls())
library(readr)
library(dplyr)
library(fda.usc)
clima_arg <- read.csv("clima_arg.csv")
source("funciones.R", echo=TRUE)


#prepare data and Plot
dim(clima_arg)
table(clima_arg$ValorMediode)
unique(clima_arg$Estacion)
tempmed <- clima_arg[clima_arg$ValorMediode == "Temperatura (°C)",]
tempmed2 <- tempmed[1:71,-2]
dim(tempmed)
matplot(t(tempmed2), type = "l")

unique(tempmed2$Estacion)
region <- c(rep("noa", 5), "chaco", "mesopotamia", "noa", rep("chaco",  3),
            rep("mesopotamia", 2), "noa", "cuyo", "noa", rep("chaco", 3),
            "mesopotamia", "cuyo", rep("sierraspampeanas", 2),
            rep("pampahumeda",6), "mesopotamia", rep("cuyo", 6),
            rep("pampahumeda", 2), rep("cuyo", 2),  rep("pampahumeda", 20),
            "patagonia", "pampahumeda", rep("patagonia",9))

length(region)
cbind(tempmed2$Estacion, region)
table(region)
tempmed3 <- tempmed2[,-1]
grilla.tes   <- seq(1,12) 

par(mfrow = c(2, 4))
matplot(t(tempmed2[region=="chaco",-1]), type ="l", col = 2, ylim = c(0,27), main = "Chaco",
        xlab = "month", ylab = "Mean temperature")
matplot(t(tempmed2)[-1, region=="mesopotamia"], type ="l", col = 3, ylim = c(0,27),
        main = "Mesopotamia", xlab = "month", ylab = "Mean temperature")
matplot(t(tempmed2)[-1, region=="noa"], type ="l", col = 4, ylim = c(0,27), main = "Noroeste",
        xlab = "month", ylab = "Mean temperature")
matplot(t(tempmed2)[-1, region=="pampahumeda"], type ="l", col = 5, ylim = c(0,27),
        main = "Pampa", xlab = "month", ylab = "Mean temperature")
matplot(t(tempmed2)[-1, region=="sierraspampeanas"], type ="l", col = 6, ylim = c(0,27),
        main = "Sierras Pampeanas", xlab = "month", ylab = "Mean temperature")
matplot(t(tempmed2)[-1, region=="patagonia"], type ="l", col = 7, ylim = c(0,27),
        main = "Patagonia", xlab = "month", ylab = "Mean temperature")
matplot(t(tempmed2)[-1, region=="cuyo"], type ="l", col = 8, ylim = c(0,27),
        main = "Cuyo", xlab = "month", ylab = "Mean temperature")

# Check Normality
tablaanova <- tablakruskal <- list()
pvalorkruskal <- pvaloranova <- 0
nn<-300
comparaciones <- array(NA, dim = c(nn, 6, 6))
i <- 1
shapiros <- 0
for( i in 1:nn){
  elembr <-rproc2fdata(1, grilla.tes, sigma = "brownian")
  #plot(elembr)
  #lines(elembr$argvals, elembr$data, col = 2)
  delta <- 1/12
  tempmed4 <- matrix(as.numeric(as.matrix(tempmed3)),
                     nrow = 71, ncol = 12)
  xx_coef <- tempmed4 %*% matrix(elembr$data) * delta
  #matplot(t(matrix(as.numeric(as.matrix(tempmed3)), nrow = 71, ncol = 12)), type = "l")
  xx_coef
  datos <- data.frame(region = region, tempproy = xx_coef)
  datos
  residuosanova <- aov(tempproy ~ region, data = datos)$residuals
  shapiros[i] <- shapiro.test(residuosanova)$p.value
  #boxplot(tempproy ~ region, data = datos)
  tablakruskal[[i]] <- kruskal.test(tempproy ~ region, data = datos)
  tablaanova[[i]] <- aov(tempproy ~ region, data = datos)
  comparacionesaux <- 
    pairwise.wilcox.test(xx_coef, region, p.adjust.method = "BY")$p.value
  for(l in 1:6) {
    for (j in 1:6) {
      comparaciones[i,j,l] <- comparacionesaux[j, l]
    }
  }
  pvalorkruskal[i] <- tablakruskal[[i]]$p.value
  pvaloranova[i] <- summary(tablaanova[[i]])[[1]][1,5]
}
boxplot(shapiros)

mean(p.adjust(shapiros, method = "fdr") < 0.2)
mean(p.adjust(shapiros, method = "fdr") < 0.05)

mean(shapiros<0.2)
mean(shapiros<0.05)

#compute the test


source("funciones_pairw_v2.R", echo=TRUE)
library(gdata)
X <- tempmed4
g <- region
dim(X)
grid <- 1:12
set.seed(56793)
reskw <- func.kw.test(X, g, grid, nn = 30)
reskw
round(reskw$p.value.pairw , digits = 3)


pvaloresaov <- func.aov.test(X, g, grid, nn = 30)
pvaloresaov

