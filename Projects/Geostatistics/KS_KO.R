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
#load("/home/rcamilo/Documents/geoestadistica/IDSTA.ov.rda")
# load("/home/rcamilo/Documents/geoestadistica/croatia.grilla 25s.2008.rda")
#pc 1
load("D:/Documentos/11 semestre/Geoestadística/ProyectoGeoestadistica/code/IDSTA.ov.rda")
load("D:/Documentos/11 semestre/Geoestadística/ProyectoGeoestadistica/code/croatia.grilla 25s.2008.rda")

IDSTA.OV

library(rgdal)
IDSTA.OV@bbox
rasterToPolygons(IDSTA.OV)
IDSTAr.ov <- IDSTA.ov[!is.na(IDSTA.ov$LST2008_07_03),]

#temperatura 
Tmt<-IDSTAr.ov$LST2008_07_03
min(Tmt)
max(Tmt)
#guardar geodata 
croacia.geoR <- as.geodata(IDSTAr.ov, data.col = 'LST2008_07_03') 
E<-croacia.geoR$coords[,1]
N<-croacia.geoR$coords[,2]
plot(E,N)

#no normalidad  entonces anamorfosis gaussiana
library(RGeostats)
library(gstat)
data<-cbind(E, N, Tmt)
temp<-as.data.frame(data)
#01. Fit Gaussian Anamorphosis 
mdb = temp[,c("E","N","Tmt")] #must be 3 columns: x,y,z
mdb.rgdb = db.create(mdb,ndim=2,autoname=F)
mdb.herm = anam.fit(mdb.rgdb,name="Tmt",type="gaus")
mdb.hermtrans = anam.z2y(mdb.rgdb,names="Tmt",anam=mdb.herm)
Tmt.trans = mdb.hermtrans@items$Gaussian.Tmt 



library(gstat)
library(sp)
data2<-cbind(E, N, Tmt.trans)
mu<-mean(Tmt.trans)
temp2s<-as.data.frame(data2)
coordinates(temp2s) = ~E+N
grillabaser<-as.data.frame(IDSTA.OV) #se vuelve data frame
borde<-grillabaser[,52:53]
names(borde) <- c("E","N")
coordinates(borde) = ~E+N
# m = vgm(0.6,"Exp",35000, 0.4)
detach("package:sgeostat")
# 
# zn.svgm <- variogram(Tmt.trans~1,temp2s)
# vgm1 = vgm(0.544,"Exp",19000, 0.01)
# plot(zn.svgm,m)


# m <- vgm1
m = vgm(2.1,"Exp",32000, 0.02)
# Kriging Ordinario
pron.pts.ko<-krige(Tmt.trans~1,temp2s,borde,model=m) 
ko.bts = cbind(coordinates(pron.pts.ko),pron.pts.ko@data)
ko.bts.db = db.create(ko.bts,autoname = F)
ko.tempdb = anam.y2z(ko.bts.db,names="var1.pred",anam = mdb.herm)
ko.tempdbvar = anam.y2z(ko.bts.db,names="var1.var",anam = mdb.herm)
min(ko.tempdb@items$Raw.var1.pred)
max(ko.tempdb@items$Raw.var1.pred)
mean(ko.tempdb@items$Raw.var1.pred)
summary(ko.tempdbvar@items$R)
#Prediction map
grilla_ko<-grillabaser[,52:53]
grilla_ko$pred<-ko.tempdb@items$Raw.var1.pred
grilla_ko$var<-ko.tempdbvar@items$var1.var
names(grilla_ko) <- c("E","N","pre","var")
coordinates(grilla_ko)= ~E+N
gridded(grilla_ko)=T
x11()
spplot(grilla_ko, "pre", main="Temperatura media terrestre \nPredicciones Kriging Ordinario", col.regions=bpy.colors(150), cuts=15, cex.main=0.2, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)",  key.space=list(space="right", cex=0.6))

spplot(grilla_ko, "var", main="Temperatura media terrestre \nVarianzas Kriging Ordinario", col.regions=bpy.colors(150), cuts=15, cex.main=0.2, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)",  key.space=list(space="right", cex=0.6))

typeof(pron.pts.ko)
# Kriging Simple

ms = vgm (2.1,"Exp",44000, 0.01)
pron.pts.ks<-krige(Tmt.trans~1,temp2s,borde,model=ms,beta=mu) 

ks.bts = cbind(coordinates(pron.pts.ks),pron.pts.ks@data)
ks.bts.db = db.create(ks.bts,autoname = F)
ks.tempdb = anam.y2z(ks.bts.db,names="var1.pred",anam = mdb.herm)
ks.tempdbvar = anam.y2z(ks.bts.db,names="var1.var",anam = mdb.herm)
min(ks.tempdb@items$Raw.var1.pred)
max(ks.tempdb@items$Raw.var1.pred)
mean(ks.tempdb@items$Raw.var1.pred)
min(ks.tempdbvar@items$var1.var)
max(ks.tempdbvar@items$var1.var)
mea
min(Tmt)
max(Tmt)
#Prediction map
grilla_ks<-grillabaser[,52:53]
grilla_ks$pred<-ks.tempdb@items$Raw.var1.pred
grilla_ks$var<-ks.tempdbvar@items$var1.var
names(grilla_ks) <- c("E","N","pre","var")
coordinates(grilla_ks)= ~E+N
gridded(grilla_ks)=T
x11()
spplot(grilla_ks, "pre", main="Temperatura media terrestre\nPredicciones Kriging Simple", col.regions=bpy.colors(150), cuts=15, cex.main=0.2, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)",  key.space=list(space="right", cex=0.6))

spplot(grilla_ks, "var", main="Temperatura media terrestre \nVarianza Kriging Simple", col.regions=bpy.colors(150), cuts=15, cex.main=0.2, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)",  key.space=list(space="right", cex=0.6))


KO.esf.cv.z <- krige.cv(Tmt.trans~1,temp2s,borde,model=m)     
KS.esf.cv.z <- krige.cv(Tmt.trans~1,temp2s,borde,model=ms,beta=mu) 

Tmt<
# Validación

resultados.cv.z <- rbind(criterio.cv(KO.esf.cv.z), criterio.cv(KS.esf.cv.z))
rownames(resultados.cv.z) <- c("KO.esf.cv.z", "KS.esf.cv.z")
resultados.cv.z


#Krigin universal 1

x <- krige(Tmt.trans~E+N,temp2s,borde,model=m, block = c(40,40))
spplot(x["var1.pred"], main = "universal kriging predictions")

x <- krige.cv(Tmt.trans~E+N,temp2s,borde,model=m)
criterio.cv(x)


# contorno del mapa de croacia 
df.pts<-spsample(IDSTA.OV, n=9, type="regular") 
loc1 <- df.pts@coords 
colnames(loc1) <- c("x","y")
plot(IDSTA.OV)
x11()
plot(loc1)

library(maptools)
Polygon <- readShapePoly("D:/Documentos/11 semestre/Geoestadística/ProyectoGeoestadistica/original_data/croatia.shp")
shape1<-readOGR("D:/Documentos/11 semestre/Geoestadística/ProyectoGeoestadistica/original_data/croatia.shp")
#### Creación de la grilla con los puntos de predicción y ploteo
network.design()
vgm2 <-  vgm(2.1,"Exp",32000, 0.02)
set.seed(123)
x11()
NDP1 <- network.design(x~1,vgm.model=vgm2, npoints=100,boundary=shape1,nmax=6, beta=mu, type="stratified")
NDP2 <- network.design(x~1, vgm.model=vgm2,npoints=400, boundary=shape1,nmax=5, beta=mu,  type="stratified")
x11()
NDP3 <- network.design(x~1, vgm.model=vgm2,npoints=900, boundary=shape1,nmax=5, beta=mu, type="stratified")
x11()
NDP4 <- network.design(x~1,vgm.model=vgm2,npoints=1600,boundary=shape1,nmax=5, beta=mu, type="stratified")

Networks.P <- rbind(NDP1,NDP2,NDP3,NDP4)
colnames(Networks.P) <- c("ASEPE")
library(xtable)
xtable(Networks.P)

