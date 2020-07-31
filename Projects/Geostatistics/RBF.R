# clear the workspace, plots and console
rm(list = ls())
if(!is.null(dev.list())) dev.off()
cat("\014")
#restart R
# .rs.restartR()
#libraries
library(sp)
library(geoR)
library(geospt)
#data
load("/home/rcamilo/Documents/geoestadistica/IDSTA.ov.rda")
load("/home/rcamilo/Documents/geoestadistica/croatia.grilla 25s.2008.rda")

IDSTAr.ov <- IDSTA.ov[!is.na(IDSTA.ov$LST2008_07_03),]

#temperatura 
Tmt<-IDSTAr.ov$LST2008_07_03
min(Tmt)
max(Tmt)
#guardar geodata 
croacia.geoR <- as.geodata(IDSTAr.ov, data.col = 'LST2008_07_03') 
E<-croacia.geoR$coords[,1]
N<-croacia.geoR$coords[,2]
data2<-cbind(N, E, Tmt)
temp2<-as.data.frame(data2)
###FUNCIONES DE BASE RADIAL####
#toda la grilla para predecir sobre esta
grillabaser<-as.data.frame(IDSTA.OV) #se vuelve data frame
puntosxyb<-grillabaser[,52:53] #se seleccionan las columnas donde estan las coordenadas x, y
names(puntosxyb) <- c("N","E") #esta en x, y entonces cambiamos a N y E para ppoder hacer las RBF
plot(puntosxyb)

temp2.spb<-temp2[,1:2]#coordenadas
# temp2.spb<- puntosxyb
temp2.spb$Tmt<-Tmt
length(Tmt)
coordinates(temp2.spb)
head(temp2.spb) #estos son mis datos, N, E Y Tmt(temperatura)

#Hacemos la optimzacion para cada una de RBF

###### Multicuadratica########## RMSPE dio=2.494421
op.m <- graph.rbf(Tmt~1, temp2.spb, eta.opt=TRUE, rho.opt=TRUE, n.neigh=9, func="M", 
                  eta.dmax=2, rho.dmax=2, x0=c(0.1,0.1), iter=500)
###### Multicuadratica inversa##########RMSPE dio=2.610865
op.mi <- graph.rbf(Tmt~1, temp2.spb, eta.opt=TRUE, rho.opt=TRUE, n.neigh=9, func="IM", 
                   eta.dmax=2, rho.dmax=2, x0=c(0.1,0.1), iter=500)
###### Spline con tension ##########RMSPE dio=2.399818
op.st <- graph.rbf(Tmt~1, temp2.spb, eta.opt=TRUE, rho.opt=TRUE, n.neigh=9, func="ST",
                   eta.dmax=2, rho.dmax=2, x0=c(0.1,0.1), iter=500)
###### Spline completamente regularizada ##########RMSPE dio=2.333328
op.crs <- graph.rbf(Tmt~1, temp2.spb, eta.opt=TRUE, rho.opt=TRUE, n.neigh=9, func="CRS",
                    eta.dmax=2, rho.dmax=2, x0=c(0.1,0.1), iter=500)
###### Spline capa delgada ##########RMSPE dio=3.053522
op.tps <- graph.rbf(Tmt~1, temp2.spb, eta.opt=TRUE, rho.opt=TRUE, n.neigh=9, func="TPS",
                    eta.dmax=2, rho.dmax=2, x0=c(0.1,0.1), iter=500)
###### Exponencial ##########RMSPE dio=2.611089
op.expon <- graph.rbf(Tmt~1, temp2.spb, eta.opt=TRUE, rho.opt=TRUE, n.neigh=9, func="EXPON",
                      eta.dmax=2, rho.dmax=2, x0=c(0.1,0.1), iter=500)
###### Gaussiana#########RMSPE dio=2.611089
op.gau <- graph.rbf(Tmt~1, temp2.spb, eta.opt=TRUE, rho.opt=TRUE, n.neigh=9, func="GAU",
                    eta.dmax=2, rho.dmax=2, x0=c(0.1,0.1), iter=500)

#multicuadratica
cv.m <- rbf.tcv(Tmt~1, data=as.data.frame( temp2.spb) , eta=0.00001,rho=0, n.neigh=9, func="M") 
criterio.cv(cv.m)
warnings()
#multicuadratica Inversa
cv.mi <- rbf.tcv(Tmt~1, data=as.data.frame( temp2.spb) , eta=2,rho=0, n.neigh=9, func="IM") 
criterio.cv(cv.mi)
#ST
cv.st <- rbf.tcv(Tmt~1, data=as.data.frame( temp2.spb) , eta=0.0002882479,rho=0.07792173, n.neigh=9, func="ST") 
criterio.cv(cv.st)
#CRS
cv.crs <- rbf.tcv(Tmt~1, data=as.data.frame( temp2.spb) , eta=0.0003768368,rho=0.07983177,n.neigh=9, func="CRS") 
criterio.cv(cv.crs)
#TPS
cv.tps <- rbf.tcv(Tmt~1, data=as.data.frame( temp2.spb) , eta=2,rho=0.02031208, n.neigh=9, func="TPS") 
criterio.cv(cv.tps)
#EXPONENCIAL
cv.expon <- rbf.tcv(Tmt~1, data=as.data.frame( temp2.spb) , eta=0.1,rho=0.12, n.neigh=9, func="EXPON") 
criterio.cv(cv.expon)
#GAUSSIANA
cv.gau <- rbf.tcv(Tmt~1, data=as.data.frame( temp2.spb) , eta=0.1,rho=0.12, n.neigh=9, func="GAU") 
criterio.cv(cv.gau)


#se elije CRS
pred.rbf.st<-rbf(Tmt~1, data=as.data.frame( temp2.spb) , eta=0.0003768368,rho=0.07983177,n.neigh=9, newdata=puntosxyb, func="CRS") 
# muestra el mapa de prediccion con la rbf Spline completamente regularizada
puntosxyb$pred <- pred.rbf.st$var1.pred
grafb<-as.data.frame(puntosxyb)
coordinates(grafb) = ~N+E
gridded(grafb)=T

spplot(grafb ,"pred", main="Temperatura media 3 de julio de 2008 Predicciones \nSpline with Tension", col.regions=bpy.colors(101), cuts=100, cex.main=0.5, regions=T, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)",  key.space=list(space="left", cex=0.6))

plot(croacia.geoR)

min(pred.rbf.st$var1.pred)
max(pred.rbf.st$var1.pred)

    