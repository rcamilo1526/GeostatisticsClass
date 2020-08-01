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

temp2.spb<-temp2[,1:2]#coordenadas
# temp2.spb<- puntosxyb
temp2.spb$Tmt<-Tmt
length(Tmt)
coordinates(temp2.spb)
head(temp2.spb) #estos son mis datos, N, E Y Tmt(temperatura)

data2<-cbind(E, N, Tmt)
temp2<-as.data.frame(data2)
temp2.sp <- SpatialPoints(as.data.frame(temp2)) 
temp2.spt<-as.data.frame(temp2.sp)
temp2.spb<-temp2.spt[,1:2]#coordenadas
temp2.spb$Tmt<-Tmt
plot(temp2)temp2.sp <- SpatialPoints(as.data.frame(temp2)) 
umtna<-"+NA"
proj4string(temp2.sp) <- CRS(umtna)

#################################
######## KRIGING ORDINARIO#######
#################################

#02. Variogram fit to gaussian transformed data using gstat 
zn.svgm <- variogram(Tmt.trans~1,temp2.sp)
vgm1 = vgm(0.544,"Exp",19000, 0.01)
plot(zn.svgm,vgm1)
ko.gstat <- gstat(id="Tmt.trans", formula=Tmt.trans~1, model=vgm1, data=temp2.sp,nmin=2, nmax=30)

#03. Kriging predictions with grid (gstat)
grillasin<-as.data.frame(IDSTA.OV)
coordinates(grillasin) = ~x+y
gridded(grillasin)=T
croatia.ks.ok = predict(ko.gstat,grillasin) 

#04. Back transform using RGeostats
ko.bts = cbind(coordinates(croatia.ks.ok),croatia.ks.ok@data)
ko.bts.db = db.create(ko.bts,autoname = F)
ko.tempdb = anam.y2z(ko.bts.db,names="Tmt.trans.pred",anam = mdb.herm)
ko.tempdbvar = anam.y2z(ko.bts.db,names="Tmt.trans.var",anam = mdb.herm)
#Prediction map
IDSTA.OV@data$pred.ko.temp <- ko.tempdb@items$Raw.Tmt.trans.pred
spplot(IDSTA.OV,"pred.ko.temp")
#Var map
IDSTA.OV@data$var.ko.temp <- ko.tempdbvar@items$Raw.Tmt.trans.var
spplot(IDSTA.OV,"var.ko.temp")
IDSTA.OV@data$var.ko.temp.sqrt <- sqrt(ko.tempdbvar@items$Raw.Tmt.trans.var)


##############################
######## KRIGING SIMPLE#######
##############################

mu<-mean(temp2.sp$Tmt.trans)
ks.gstat <- gstat(id="Tmt.trans", formula=Tmt.trans~1, model=vgm1, data=temp2.sp,nmin=2, nmax=30,beta=mu)

#03. Kriging predictions with grid (gstat)
grillasin<-as.data.frame(IDSTA.OV)
coordinates(grillasin) = ~x+y
gridded(grillasin)=T
croatia.ks.ok = predict(ks.gstat,grillasin) 

#04. Back transform using RGeostats
gridded(grillasin)=T
ks.bts = cbind(coordinates(croatia.ks.ok),croatia.ks.ok@data)
ks.bts.db = db.create(ks.bts,autoname = F)
ks.tempdb = anam.y2z(ks.bts.db,names="Tmt.trans.pred",anam = mdb.herm)
ks.tempdbvar = anam.y2z(ks.bts.db,names="Tmt.trans.var",anam = mdb.herm)
#Prediction map
IDSTA.OV@data$pred.ks.temp <-ks.tempdb@items$Raw.Tmt.trans.pred
spplot(IDSTA.OV,"pred.ks.temp")
#Var map
IDSTA.OV@data$var.ks.temp.sqrt <- sqrt(ks.tempdbvar@items$Raw.Tmt.trans.var)
spplot(IDSTA.OV,"var.ks.temp.sqrt")


IDSTA.OV@data$pred.ks.temp.res <- abs(Tmt - IDSTA.OV$pred.ko.temp)
spplot(IDSTA.OV, "pred.ks.temp.res", main="RES ABS \nKriging simple", col.regions=bpy.colors(100), cuts=10, cex.main=0.7, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=0.6))
spplot(IDSTA.OV, "var.ks.temp.sqrt", main="SQRT(VAR) \nKriging simple", col.regions=bpy.colors(100), cuts=10, cex.main=0.7, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=0.6))



#SALIDA GRAFICA

l2 = list("sp.points", So, pch = 3, col = "grey")
rw.colors <- colorRampPalette(c("red", "yellow"))
p1 <- 
  spplot(IDSTA.OV, "pred.ks.temp", main="Temperatura media 25 de enero de 2008 Predicciones \nKriging simple", col.regions=bpy.colors(100), cuts=10, cex.main=0.7, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=0.6))
p2 <- 
  spplot(IDSTA.OV, "var.ks.temp.sqrt", main="Temperatura media 25 de enero de 2008 Errores \nKriging simple", col.regions=bpy.colors(100), cuts=10, cex.main=0.7, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=0.6))
p3 <- 
  spplot(IDSTA.OV, "pred.ko.temp", main="Temperatura media 25 de enero de 2008 Predicciones \nKriging ordinario", col.regions=bpy.colors(100), cuts=10, cex.main=0.7, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=0.6))
p4 <- 
  spplot(IDSTA.OV, "var.ko.temp.sqrt", main="Temperatura media 25 de enero de 2008 Errores \nKriging ordinario", col.regions=bpy.colors(100), cuts=10, cex.main=0.7, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=0.6))
print(p1, split = c(1, 1, 2, 2), more = T)
print(p2, split = c(2, 1, 2, 2), more = T)
print(p3, split = c(1, 2, 2, 2), more = T)
print(p4, split = c(2, 2, 2, 2), more = T)

writeAsciiGrid(IDSTA.OV@data$var.ks.temp, "ko.asc")
write.asciigrid(pred.ks.temp, "ko.asc")