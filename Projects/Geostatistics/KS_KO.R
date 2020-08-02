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
data<-cbind(N, E, Tmt)
temp<-as.data.frame(data)
#01. Fit Gaussian Anamorphosis 
mdb = temp[,c("N","E","Tmt")] #must be 3 columns: x,y,z
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
m = vgm(0.6,"Exp",35000, 0.4)
detach("package:sgeostat")
 

# Kriging Ordinario

pron.pts.ko<-krige(Tmt.trans~1,temp2s,borde,model=m) 
ko.bts = cbind(coordinates(pron.pts.ko),pron.pts.ko@data)
ko.bts.db = db.create(ko.bts,autoname = F)
ko.tempdb = anam.y2z(ko.bts.db,names="var1.pred",anam = mdb.herm)
ko.tempdbvar = anam.y2z(ko.bts.db,names="var1.var",anam = mdb.herm)
min(ko.tempdb@items$Raw.var1.pred)
max(ko.tempdb@items$Raw.var1.pred)

#Prediction map
grilla_ko<-grillabaser[,52:53]
grilla_ko$pred<-ko.tempdb@items$Raw.var1.pred
grilla_ko$var<-ko.tempdbvar@items$Raw.var1.var
names(grilla_ko) <- c("E","N","pre","var")
coordinates(grilla_ko)= ~E+N

spplot(grilla_ko, "pre", main="Precipitación Media Anual \nPredicciones Kriging Ordinario", col.regions=bpy.colors(100), cuts=10, cex.main=0.2, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)",  key.space=list(space="right", cex=0.6))

spplot(grilla_ko, "var", main="Precipitación Media Anual \nErrores Kriging Ordinario", col.regions=bpy.colors(100), cuts=10, cex.main=0.2, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)",  key.space=list(space="right", cex=0.6))


# Kriging Simple

pron.pts.ks<-krige(Tmt.trans~1,temp2s,borde,model=m,beta=mu) 

ks.bts = cbind(coordinates(pron.pts.ks),pron.pts.ks@data)
ks.bts.db = db.create(ks.bts,autoname = F)
ks.tempdb = anam.y2z(ks.bts.db,names="var1.pred",anam = mdb.herm)
ks.tempdbvar = anam.y2z(ks.bts.db,names="var1.var",anam = mdb.herm)
min(ks.tempdb@items$Raw.var1.pred)
max(ks.tempdb@items$Raw.var1.pred)

#Prediction map
grilla_ks<-grillabaser[,52:53]
grilla_ks$pred<-ks.tempdb@items$Raw.var1.pred
grilla_ks$var<-ks.tempdbvar@items$var1.var
names(grilla_ks) <- c("E","N","pre","var")
coordinates(grilla_ks)= ~E+N

spplot(grilla_ks, "pre", main="Precipitación Media Anual \nPredicciones Kriging Ordinario", col.regions=bpy.colors(100), cuts=10, cex.main=0.2, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)",  key.space=list(space="right", cex=0.6))

spplot(grilla_ks, "var", main="Precipitación Media Anual \nErrores Kriging Ordinario", col.regions=bpy.colors(100), cuts=10, cex.main=0.2, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)",  key.space=list(space="right", cex=0.6))


KO.esf.cv.z <- krige.cv(Tmt.trans~1,temp2s,borde,model=m)     
KS.esf.cv.z <- krige.cv(Tmt.trans~1,temp2s,borde,model=m,beta=mu) 

# Validación

resultados.cv.z <- rbind(criterio.cv(KO.esf.cv.z), criterio.cv(KS.esf.cv.z))
rownames(resultados.cv.z) <- c("KO.esf.cv.z", "KS.esf.cv.z")
resultados.cv.z


#Krigin universal 1

x <- krige(Tmt.trans~E+N,temp2s,borde,model=m, block = c(40,40))
spplot(x["var1.pred"], main = "universal kriging predictions")

x <- krige.cv(Tmt.trans~E+N,temp2s,borde,model=m)
criterio.cv(x)
