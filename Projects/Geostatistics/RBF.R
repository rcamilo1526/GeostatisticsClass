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
#pc 0
# load("/home/rcamilo/Documents/geoestadistica/IDSTA.ov.rda")
# load("/home/rcamilo/Documents/geoestadistica/croatia.grilla 25s.2008.rda")
#pc 1
load("D:/Documentos/11 semestre/Geoestadística/ProyectoGeoestadistica/code/IDSTA.ov.rda")
load("D:/Documentos/11 semestre/Geoestadística/ProyectoGeoestadistica/code/croatia.grilla 25s.2008.rda")


IDSTAr.ov <- IDSTA.ov[!is.na(IDSTA.ov$LST2008_07_03),]

#temperatura 
Tmt<-IDSTAr.ov$LST2008_07_03
min(Tmt)
max(Tmt)
na
#guardar geodata 
croacia.geoR <- as.geodata(IDSTAr.ov, data.col = 'LST2008_07_03') 
E<-croacia.geoR$coords[,1]
N<-croacia.geoR$coords[,2]
plot(E,N)


data2<-cbind(E, N, Tmt)
temp2<-as.data.frame(data2)
temp2.sp <- SpatialPoints(as.data.frame(temp2)) 
temp2.spt<-as.data.frame(temp2.sp)
temp2.spb<-temp2.spt[,1:2]#coordenadas
temp2.spb$Tmt<-Tmt

###FUNCIONES DE BASE RADIAL####
#toda la grilla para predecir sobre esta
grillabaser<-as.data.frame(IDSTA.OV) #se vuelve data frame
puntosxyb<-grillabaser[,52:53] #se seleccionan las columnas donde estan las coordenadas x, y
puntosxyb<-grillabaser[,c('x','y')]
coordinates(puntosxyb)=~E+N
x11()
plot(puntosxyb)


#Hacemos la optimzacion para cada una de RBF

###### Multicuadratica########## RMSPE dio=2.498376
op.m <- graph.rbf(Tmt~1, temp2.spb, eta.opt=TRUE, rho.opt=TRUE, n.neigh=15, func="M", 
                  eta.dmax=2, rho.dmax=2, x0=c(0.1,0.1), iter=500)
###### Multicuadratica inversa##########RMSPE dio=2.592904
op.mi <- graph.rbf(Tmt~1, temp2.spb, eta.opt=TRUE, rho.opt=TRUE, n.neigh=15, func="IM", 
                   eta.dmax=2, rho.dmax=2, x0=c(0.1,0.1), iter=500)
###### Spline con tension ##########RMSPE dio=2.38963
op.st <- graph.rbf(Tmt~1, temp2.spb, eta.opt=TRUE, rho.opt=TRUE, n.neigh=15, func="ST",
                   eta.dmax=2, rho.dmax=2, x0=c(0.1,0.1), iter=500)
###### Spline completamente regularizada ##########RMSPE dio=2.330511
op.crs <- graph.rbf(Tmt~1, temp2.spb, eta.opt=TRUE, rho.opt=TRUE, n.neigh=15, func="CRS",
                    eta.dmax=2, rho.dmax=2, x0=c(0.1,0.1), iter=800)
###### Spline capa delgada ##########RMSPE dio=2.824071
op.tps <- graph.rbf(Tmt~1, temp2.spb, eta.opt=TRUE, rho.opt=TRUE, n.neigh=15, func="TPS",
                    eta.dmax=2, rho.dmax=2, x0=c(0.1,0.1), iter=500)
###### Exponencial ##########RMSPE dio=2.593215
op.expon <- graph.rbf(Tmt~1, temp2.spb, eta.opt=TRUE, rho.opt=TRUE, n.neigh=15, func="EXPON",
                      eta.dmax=2, rho.dmax=2, x0=c(0.1,0.1), iter=500)
###### Gaussiana#########RMSPE dio=2.593215
op.gau <- graph.rbf(Tmt~1, temp2.spb, eta.opt=TRUE, rho.opt=TRUE, n.neigh=15, func="GAU",
                    eta.dmax=2, rho.dmax=2, x0=c(0.1,0.1), iter=500)

#multicuadratica
cv.m <- rbf.tcv(Tmt~1, data=as.data.frame( temp2.spb) , eta=0.00001,rho=0, n.neigh=9, func="M") 
criterio.cv(cv.m)
warnings()
#multicuadratica Inversa
cv.mi <- rbf.tcv(Tmt~1, data=as.data.frame( temp2.spb) , eta=2,rho=0, n.neigh=9, func="IM") 
criterio.cv(cv.mi)
#ST
cv.st <- rbf.tcv(Tmt~1, data=as.data.frame( temp2.spb) , eta=0.00027482429,rho=0.1054797, n.neigh=9, func="ST") 
criterio.cv(cv.st)
#CRS
cv.crs <- rbf.tcv(Tmt~1, data=as.data.frame( temp2.spb) , eta=0.000406909,rho=0.1075205,n.neigh=9, func="CRS") 
criterio.cv(cv.crs)
#TPS
cv.tps <- rbf.tcv(Tmt~1, data=as.data.frame( temp2.spb) , eta=2,rho=1.999798, n.neigh=9, func="TPS") 
criterio.cv(cv.tps)
#EXPONENCIAL
cv.expon <- rbf.tcv(Tmt~1, data=as.data.frame( temp2.spb) , eta=0.1,rho=0.1, n.neigh=9, func="EXPON") 
criterio.cv(cv.expon)
#GAUSSIANA
cv.gau <- rbf.tcv(Tmt~1, data=as.data.frame( temp2.spb) , eta=0.1,rho=0.1, n.neigh=9, func="GAU") 
criterio.cv(cv.gau)


#se elije CRS
pred.rbf.st<-rbf(Tmt~1, data=as.data.frame( temp2.spb), eta=0.000406909,rho=0.1075205,n.neigh=10, newdata=puntosxyb, func="CRS") 
# muestra el mapa de prediccion con la rbf Spline completamente regularizada
puntosxyb$pred <- pred.rbf.st$var1.pred
grafb<-as.data.frame(puntosxyb)
coordinates(grafb) = ~E+N
gridded(grafb)=T
x11()
spplot(grafb ,"pred", main="Temperatura media 3 de julio de 2008 Predicciones \nSpline completamente regularizada", col.regions=bpy.colors(101), cuts=100, cex.main=0.5, regions=T, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)",  key.space=list(space="left", cex=0.6))

plot(croacia.geoR)

min(pred.rbf.st$var1.pred)#24.06542
max(pred.rbf.st$var1.pred)#34.51251
mean(pred.rbf.st$var1.pred)#29.39422
min(Tmt)
max(Tmt)
mean(Tmt)
