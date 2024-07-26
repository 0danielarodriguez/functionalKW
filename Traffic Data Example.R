rm(list=ls())
library(dplyr)
datos <- read.csv2("flujo-vehicular-2018.csv",sep=",",header=TRUE)

dim(datos)
names(datos)
head(datos)

datos2 <- datos %>%  filter(mes!="enero"& mes!="febrero"& mes!="diciembre" &
                              mes!="julio" &
                              estacion == "ILLIA" & tipo_vehiculo == "Liviano 2 ejes"
)
#detach(datos)
attach(datos2)
dim(datos2)

horas <- names(table(datos2$hora_inicio))
estaciones <- names(table(datos2$estacion))
dias<- names(table(datos2$dia))
fechas <- names(table(datos2$fecha))
meses<-names(table(datos2$mes))
periodos<-names(table(datos2$periodo))
tipos<-names(table(datos2$tipo_vehiculo))
formas<-names(table(datos2$forma_pago))
meses<-names(table(datos2$mes))
observaciones<-names(table(datos2$observacion))
fechas<-names(table(datos2$fecha))


lestaciones<-length(estaciones)
ldias<-length(dias)
lhoras<-length(horas)
lfechas<-length(fechas)
lmeses<-length(meses)
lperiodos<-length(periodos)
ltipos<-length(tipos)
lformas<-length(formas)
lmeses<-length(meses)
lobservaciones<-length(observaciones)
dim(datos)
dim(datos2)

################# assembling the functional data matrix ######################
i<-1
datosdia1_completo <- datos2[datos2$fecha == fechas[i],]
datosdia1 <- datosdia1_completo %>% 
  group_by(hora_inicio) %>% 
  summarise(cantidad_pasos = sum(cantidad_pasos))
if(length(unique(datosdia1_completo$hora_inicio)) == 24){
  # datos3 <- data.frame(cantidad_pasos = datosdia1$cantidad_pasos,
  #                      hora_inicio = datosdia1$hora_inicio,
  #                      mes = datosdia1_completo$mes[1],
  #                      dia = datosdia1_completo$dia[1])
  X <- datosdia1$cantidad_pasos
  dia <- datosdia1_completo$dia[1]
  mes <- datosdia1_completo$mes[1]
  fecha2 <- datosdia1_completo$fecha[1]
}else{
  print(paste(fechas[i], "Dia con datos incompletos"))
}


for(i in 2:lfechas){
  datosdiai_completo <- datos2[datos2$fecha == fechas[i],]
  datosdiai <- datosdiai_completo %>% 
    group_by(hora_inicio) %>% 
    summarise(cantidad_pasos = sum(cantidad_pasos))
  if(length(unique(datosdiai_completo$hora_inicio)) == 24){
    if(length(unique(datosdiai_completo$mes))==1 &
       length(unique(datosdiai_completo$dia))==1){
        X <- rbind(X,datosdiai$cantidad_pasos)
      dia <- c(dia, datosdiai_completo$dia[i])
      mes <- c(mes, datosdiai_completo$mes[i])
      fecha2 <- c(fecha2, datosdiai_completo$fecha[i])
    }else{
      print(paste(fechas[i], i, "error en data frame"))
    }
  }else{
    print(paste(fechas[i], i,  "Dia con datos incompletos"))
  }
}

length(fecha2)
dia[dia == dias[6]] <- "Sabado"
dia[dia == dias[5]] <- "Miercoles"

# nsamp day
table(dia)


#plot the data

library(fda)
par(mfrow=c(2,4))
fboxlu <- fbplot(t(X[dia == "Lunes",]), main = "Monday", xlab = "", ylab = "",ylim=c(0,9000))
fboxma <- fbplot(t(X[dia == "Martes",]), main = "Tuesday", xlab = "", ylab = "",ylim=c(0,9000))
fboxmi <- fbplot(t(X[dia == "Miercoles",]), main = "Wednesday", xlab = "", ylab = "",ylim=c(0,9000))
fboxju <- fbplot(t(X[dia == "Jueves",]), main = "Thursday", xlab = "", ylab = "",ylim=c(0,9000))
fboxvi <- fbplot(t(X[dia == "Viernes",]), main = "Friday", xlab = "", ylab = "",ylim=c(0,9000))
fboxsa <- fbplot(t(X[dia == "Sabado",]), main = "Saturday", xlab = "", ylab = "",ylim=c(0,9000))
fboxdom <- fbplot(t(X[dia == "Domingo",]), main = "Sunday", xlab = "", ylab = "",ylim=c(0,9000))



#####################################################
# attach the test code

source("funciones_pairw_v2.R")
library(fda.usc)
library(gdata)
Xtodos <- X
grid <- 0:23
g <- dia
set.seed(56793)
reskw<- func.kw.test(Xtodos, g, grid, nn = 30)
round(reskw$p.value.pairw,4)


#without weekend
X <- X[dia!="Sabado" & dia!="Domingo",]
g <- dia[dia!="Sabado" & dia!="Domingo"]
set.seed(56793)
reskw<- func.kw.test(X, g, grid, nn = 30)
round(reskw$p.value.pairw,4)


#outliers

fechasoutlu <- fecha2[dia == "Lunes"][fboxlu$outpoint]
fechasoutma <- fecha2[dia == "Martes"][fboxma$outpoint]
fechasoutmi <- fecha2[dia == "Miercoles"][fboxmi$outpoint]
fechasoutju <- fecha2[dia == "Jueves"][fboxju$outpoint]
fechasoutvi <- fecha2[dia == "Viernes"][fboxvi$outpoint]
X1 <- Xtodos[dia == "Lunes",][-fboxlu$outpoint,]
X2 <- Xtodos[dia == "Martes",][-fboxma$outpoint,]
X3 <- Xtodos[dia == "Miercoles",][-fboxmi$outpoint,]
X4 <- Xtodos[dia == "Jueves",][-fboxju$outpoint,]
X5 <- Xtodos[dia == "Viernes",][-fboxvi$outpoint,]

dia1 <- dia[dia == "Lunes"][-fboxlu$outpoint]
dia2 <- dia[dia == "Martes"][-fboxma$outpoint]
dia3 <- dia[dia == "Miercoles"][-fboxmi$outpoint]
dia4 <- dia[dia == "Jueves"][-fboxju$outpoint]
dia5 <- dia[dia == "Viernes"][-fboxvi$outpoint]

g_so <- c(dia1, dia2, dia3, dia4, dia5)
X_so <- rbind(X1, X2, X3, X4, X5)


set.seed(56793)
reskw<- func.kw.test(X, g, grid, nn = 30)
round(reskw$p.value.pairw,4)


set.seed(56793)
resaov <- func.aov.test(X, g, grid, nn = 30)
round(resaov$p.value.pairw, 4)

set.seed(56793)
resaovso <- func.aov.test(X_so, g_so, grid, nn = 30)
round(resaovso$p.value.pairw, 4)
