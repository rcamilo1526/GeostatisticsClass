# clear the workspace, plots and console
rm(list = ls())
if(!is.null(dev.list())) dev.off()
cat("\014")
#restart R
.rs.restartR()
#libraries
library(sp)
library(geoR)
library(geospt)
#data
load("/home/rcamilo/Documents/geoestadistica/IDSTA.ov.rda")
load("/home/rcamilo/Documents/geoestadistica/croatia.grilla 25s.2008.rda")

IDSTAr.ov <- IDSTA.ov[!is.na(IDSTA.ov$LST2008_07_03),]

options(max.print = 40)

Hrdem<-IDSTAr.ov$HRdem
Hrdsea<-IDSTAr.ov$HRdsea
Hrtwi<-IDSTAr.ov$HRtwi
Lat<-IDSTAr.ov$Lat
Lon<-IDSTAr.ov$Lon
#temperatura 
Tmt<-IDSTAr.ov$LST2008_07_03

#guardar geodata 
croacia.geoR <- as.geodata(IDSTAr.ov, data.col = 'LST2008_07_03') 
E<-croacia.geoR$coords[,1]
N<-croacia.geoR$coords[,2]

#no normalidad anamorfosis gaussiana
library(RGeostats)
library(gstat)
data<-cbind(N, E, Tmt)
temp<-as.data.frame(data)
#01. Fit Gaussian Anamorphosis 
mdb = temp[,c("N","E","Tmt")] #must be 3 columns: x,y,z
mdb.rgdb = db.create(mdb,ndim=2,autoname=F)
mdb.herm = anam.fit(mdb.rgdb,name="Tmt",type="gaus")
mdb.hermtrans = anam.z2y(mdb.rgdb,names="Tmt",anam=mdb.herm)
Tmt.trans = mdb.hermtrans@items$Gaussian.Tmt 


data2<-cbind(N, E, Tmt.trans)
temp2<-as.data.frame(data2)
temp2.sp <- SpatialPoints(as.data.frame(temp2)) 
#no normalidad anamorfosis gaussiana
library(RGeostats)
library(gstat)
data<-cbind(N, E, Tmt)
temp<-as.data.frame(data)
#01. Fit Gaussian Anamorphosis 
mdb = temp[,c("N","E","Tmt")] #must be 3 columns: x,y,z
mdb.rgdb = db.create(mdb,ndim=2,autoname=F)
mdb.herm = anam.fit(mdb.rgdb,name="Tmt",type="gaus")
mdb.hermtrans = anam.z2y(mdb.rgdb,names="Tmt",anam=mdb.herm)
Tmt.trans = mdb.hermtrans@items$Gaussian.Tmt 


#ANALISIS DE TENDENCIA

modelo1<-lm(Tmt.trans~Hrdem+Hrdsea+Hrtwi+N+E)
summary(modelo1)
step.model <- stepAIC(modelo1, direction = "both", trace = FALSE)
summary(step.model)
# residuals(step.model)
#Hay tendencia, se toman los residuos
tp.geo <- croacia.geoR
tp.geo$data <- residuals(step.model)
temp.sin<-residuals(step.model)
#Con residuos
head(Tmt.trans)
head(temp.sin)

head(Tmt.trans)
head(temp.sin)

grillabaser<-as.data.frame(IDSTA.OV)

puntosxyr<-grillabaser[,52:53]
names(puntosxyr) <- c("N","E")

puntosxyt<-grillabaser[,52:53]
names(puntosxyt) <- c("N","E")

#se vuelve data frame
puntosxyb<-grillabaser[,52:53] #se seleccionan las columnas donde estan las coordenadas x, y
names(puntosxyb) <- c("N","E") #esta en x, y entonces cambiamos a N y E para ppoder hacer las RBF


temp2.spt<-as.data.frame(temp2.sp)

temp2.spr<-temp2.spt[,1:2]
coordinates(temp2.spr) = ~N+E
temp2.spr$z<-temp.sin

temp2.spo<-temp2.spt[,1:2]
coordinates(temp2.spo) = ~N+E
temp2.spo$z<-Tmt.trans

temp2.spb<-temp2.spt[,1:2]
coordinates(temp2.spb) = ~N+E
temp2.spb$z<-Tmt

head(temp2.spt[,1:2]) 

#02. Variogram fit to gaussian transformed data using gstat 
zn.svgm <- variogram(z~1,temp2.spr)
vgm1 = vgm(0.3,"Exp",19000, 0.3)
plot(zn.svgm,vgm1)
ko.gstat <- gstat(id="Tmt.trans", formula=Tmt.trans~1, model=vgm1, data=temp2.sp,nmin=2, nmax=30)
########################################################
####################krigeado REsiduos###################
########################################################
pron.pts.kur <- krige(z~N+E, temp2.spr, SpatialPoints(puntosxyr), vgm1, nmax=10)

#04. Back transform using RGeostats
LST1.btsr = cbind(coordinates(pron.pts.kur),pron.pts.kur$var1.pred)
LST1.bts.dbr = db.create(LST1.btsr,autoname = F)
tempdbr = anam.y2z(LST1.bts.dbr,names="V3",anam = mdb.herm)

#CROSSS
KUr.esf.cv.z <- krige.cv(z~N+E, temp2.spr, vgm1, nmin=0, nmax=10)    


resultados.cvr.z <- rbind(criterio.cv(KUr.esf.cv.z))
rownames(resultados.cvr.z) <- c("KUr.esf.cv.z")
resultados.cvr.z

#Prediction map
puntosxyr$pred.LST <- tempdbr@items$Raw.V3
puntosxyr$var.pred<-pron.pts.kur$var1.var

grafr<-as.data.frame(puntosxyr)
coordinates(grafr) = ~N+E
gridded(grafr)=T

min(grafr$var.pred)
max(graf$var.pred)


spplot(grafr, "pred.LST", main="Temperatura media 25 de enero de 2008 Predicciones \nKriging simple", col.regions=bpy.colors(101), cuts=100, cex.main=0.5, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=0.6))


spplot(grafr ,"pred.LST", main="Temperatura media 25 de enero de 2008 Predicciones \nKriging Universal orden 1 Residuos", col.regions=bpy.colors(101), cuts=100, cex.main=0.5, regions=T, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)",  key.space=list(space="left", cex=0.6))

spplot(grafr,"var.pred", main="Temperatura media 25 de enero de 2008 Varianzas \nKriging Universal orden 1 Residuos", col.regions=bpy.colors(101), cuts=100, cex.main=0.5, regions=T, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)",  key.space=list(space="left", cex=0.6))

###########################################
##########Kriging TRANSFORMADA#############
###########################################
pron.pts.kuo <- krige(z~N+E, temp2.spo, SpatialPoints(puntosxyt), vgm1, nmax=10)

#04. Back transform using RGeostats
LST1.btso = cbind(coordinates(pron.pts.kuo),pron.pts.kuo$var1.pred)
LST1.bts.dbo = db.create(LST1.btso,autoname = F)
tempdbo = anam.y2z(LST1.bts.dbo,names="V3",anam = mdb.herm)

#CROSSS
KUo.esf.cv.z <- krige.cv(z~N+E, temp2.spo, vgm1, nmin=0, nmax=10)    


resultados.cvo.z <- rbind(criterio.cv(KUo.esf.cv.z))
rownames(resultados.cvo.z) <- c("KUo.esf.cv.z")
resultados.cvo.z

#Prediction map
puntosxyt$pred.LST <- tempdbo@items$Raw.V3
puntosxyt$var.pred<-pron.pts.kuo$var1.var

grafo<-as.data.frame(puntosxyt)
coordinates(grafo) = ~N+E
gridded(grafo)=T

min(grafo$pred.LST)
max(grafo$pred.LST)

min(grafo$var.pred)
max(grafo$var.pred)




spplot(grafo ,"pred.LST", main="Temperatura media 25 de enero de 2008 Predicciones \nKriging Universal Orden 1 ", col.regions=bpy.colors(101), cuts=100, cex.main=0.5, regions=T, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)",  key.space=list(space="left", cex=0.6))

spplot(grafo,"var.pred", main="Temperatura media 25 de enero de 2008 Varianzas \nKriging Universal Orden 1 ", col.regions=bpy.colors(101), cuts=100, cex.main=0.5, regions=T, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)",  key.space=list(space="left", cex=0.6))





