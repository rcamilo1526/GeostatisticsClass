# Limpia el espacio de trabajo, las parcelas y la consola
rm(list = ls())
if(!is.null(dev.list())) dev.off()
cat("\014") 

#Cargar librerias
library(leaps)
library(spgwr)
library(adespatial)
library(raster)
library(tmap)
library(lmtest)
library(RColorBrewer)
library(classInt)
library(spdep)
library(sphet)
library(pgirmess)
library(rgdal)
library(spatialreg)
library(ggplot2)
library(sp)
library(dismo)
library(prettymapr)
library(lattice)
library(xtable)
library(car)
library(tidyverse)
library(stargazer)
library(MASS)
#cargar funciones
load("C:/Users/rmartin/Documents/11 semestre/Geoestadística/FuncDatosArea/ImagenDatosDeArea")
set.seed(123)
setwd("C:/Users/rmartin/Documents/11 semestre/Geoestadística/ProyectoDatosArea")

#puntos
balPt<-readOGR('./baltimore','baltim')
#poligonos de voronoi
balPol <- voronoi(balPt)
x11()
plot(balPol, col='#c7d6fc')
plot(balPt,add=T,col="#3f9945",pch = 18)


bal.poly <- balPol
coords <- coordinates(bal.poly)
##Variable a explicar
PRICE <- as.data.frame(bal.poly)$PRICE
x <- c(PRICE)
x1<-(x-mean(x))/sd(x) 

# writeOGR(balPol, dsn = "shapefile", layer = "balPol",driver = "ESRI Shapefile" )

#corregir BMENT
balPol$BMENT<-ifelse(balPol$BMENT==0,0,1)
bal.poly$BMENT<-balPol$BMENT
# x11()
# op=par(mfrow=c(1,2))
# plot(balPol)
# plot(coords)


rn <- sapply(slot(balPol, "polygons"), function(x) slot(x, "ID"))
k1 <- knn2nb(knearneigh(coords))
all.linked <- max(unlist(nbdists(k1, coords)))
bal.nb.0.all <- dnearneigh(coords, 0, all.linked, row.names=rn)
summary(bal.nb.0.all, coords)
plot(balPol, border="gray7") #mirar las conexiones de los vecinos
plot(bal.nb.0.all, coords, add=TRUE, col="darkorchid")
title(main=paste("Distance based neighbours 0-", format(all.linked), " distance units", sep=""))

#Mapa con etiquetas del ID
x11()
par(mai=c(0,0,0,0))
set.seed(1)
plot(balPol, col=sample(rainbow(500))[1:49])
xy <- coords
points(xy, cex=0.3, pch=20, col='white')
text(balPol, 'STATION', cex=0.52, col='black')


#######################################################
# GRAFICO VARIABLE ENDOGENA
#######################################################

x11()
n=length(bal.poly)
l2 = list("SpatialPolygonsRescale", layout.north.arrow(), offset =  c(0,0), scale = 0.5)   # Layout
l3 = list("SpatialPolygonsRescale", layout.scale.bar(), offset = c(1,1),  scale = 0.5, fill=c("black"))
l4 = list("sp.text", c(6.6,13.7), "0")
l5 = list("sp.text", c(7.1,13.7), "500 m")
#display.brewer.all()
cols <- brewer.pal(9, "YlGn")
pal <- colorRampPalette(cols)
x11()
spplot(bal.poly["PRICE"], 
       scales=list(draw=TRUE), 
       col.regions=pal(20),  
       sp.layout=list(l2,l3,l4,l5),
       main = "House sales price Baltimore")

# rainbow(n, start=.7, end=.1), heat.colors(n), terrain.colors(n), topo.colors(n), cm.colors(n)
plotvar <- round(balPol@data$PRICE,1)    
nclr <- 5
plotclr <- brewer.pal(nclr,"YlOrBr")
class <- classIntervals(plotvar, nclr, style="quantile")
colcode <- findColours(class, plotclr)
plot(bal.poly)   
plot(bal.poly, col=colcode, add=T)
title(main="Tasa de Criminalidad", sub="Cuantiles (Igual-Frecuencia) Intervalos de Clase")
# mtext("Cuantiles (Igual-Frecuencia) Intervalos de Clase", side=1)
legend(locator(1), legend=names(attr(colcode, "table")), fill=attr(colcode, "palette"), cex=0.8, bty="n")
northarrow(c(10.1,14.5),0.2,cex=0.6)
scalebar(c(6.1,10.7),length=1,unit="m",division.cex=0.6)

#######################################################
# CALCULO PESOS ESPACIALES
#######################################################
#EFECTO REINA 
bal.nb.q <- poly2nb(balPol) 
bal.lags.q <- nblag(bal.nb.q, 9)

#EFECTO TORRE
bal.nb.r <- poly2nb(balPol, queen = F)
bal.lags.r <- nblag(bal.nb.r, 9) 

#####pesos K-Vecino####
#k.nb <- read.gwt2nb("atl_hom.gwt") 
IDs <- row.names(as(bal.poly, "data.frame"))
bal.kn1<-knn2nb(knearneigh(coords, k=1), row.names=IDs)
bal.kn2<-knn2nb(knearneigh(coords, k=2), row.names=IDs)
bal.kn3<-knn2nb(knearneigh(coords, k=3), row.names=IDs)
bal.kn4<-knn2nb(knearneigh(coords, k=4), row.names=IDs)
bal.kn5<-knn2nb(knearneigh(coords, k=5), row.names=IDs)
bal.kn6<-knn2nb(knearneigh(coords, k=6), row.names=IDs)
bal.kn7<-knn2nb(knearneigh(coords, k=7), row.names=IDs)  
bal.kn8<-knn2nb(knearneigh(coords, k=8), row.names=IDs)
bal.kn9<-knn2nb(knearneigh(coords, k=9), row.names=IDs)

#============================================================
# Criterios basados en gr?ficas
#============================================================
x11()
op=par(mfrow=c(2,2))
trinb=tri2nb(coords)
delaunay <-nb2listw(trinb, style="W")
plot(bal.poly,border="gray")
plot(trinb,coords,add=T,col="blue")
title(main="Triangulación Delaunay")
soinb=graph2nb(soi.graph(trinb,coords))
esf.influencia <-nb2listw(soinb, style="W")
plot(bal.poly,border="gray")
plot(soinb,coords,add=T,col="green")
title(main="Esfera de influencia")
gabrielnb=graph2nb(gabrielneigh(coords),sym=TRUE)
gabriel <-nb2listw(gabrielnb, style="W")
plot(bal.poly,border="gray")
plot(gabrielnb,coords,add=T,col="red")
title(main="Gráfica de Gabriel")
relativenb=graph2nb(relativeneigh(coords),sym=TRUE)
vec.relative <-nb2listw(relativenb, style="W")
plot(bal.poly,border="gray")
plot(relativenb,coords,add=T,col="orange")
title(main="Vecinos relativos")
par(op)

x11()
sp.cr <- sp.correlogram(trinb, balPol$PRICE, order=9, method="corr", style="W", zero.policy=T)
cor <- sp.correlogram(trinb, balPol$PRICE, order=9, method="I", style="W", zero.policy=T)
plot(cor, main="PRICE")

#==========================================================================
############### SELECCI?N DE MATRICES DE VECINDAD POR (PCNM) ##############
#==========================================================================
#principal coordinates of neighbour matrices (PCNM, Borcard and Legendre (2002))
#=============================================================
############### SELECCI?N DE LA MATRIZ DE PESOS ##############
#=============================================================


#Queen
q.lw.1<-nb2listw(bal.nb.q, style="W")
summary(test.W(balPol@data$PRICE,bal.nb.q))
q.lw.2<-nb2listw(bal.lags.q[[2]], style="W")
summary(test.W(balPol@data$PRICE,bal.lags.q[[2]]))
q.lw.3<-nb2listw(bal.lags.q[[3]], style="W")
summary(test.W(balPol@data$PRICE,bal.lags.q[[3]]))
q.lw.4<-nb2listw(bal.lags.q[[4]], style="W")
summary(test.W(balPol@data$PRICE,bal.lags.q[[4]]))
q.lw.5<-nb2listw(bal.lags.q[[5]], style="W")
summary(test.W(balPol@data$PRICE,bal.lags.q[[5]]))
q.lw.6<-nb2listw(bal.lags.q[[6]], style="W")
summary(test.W(balPol@data$PRICE,bal.lags.q[[6]]))
q.lw.7<-nb2listw(bal.lags.q[[7]], style="W")
summary(test.W(balPol@data$PRICE,bal.lags.q[[7]]))
q.lw.8<-nb2listw(bal.lags.q[[8]], style="W")
summary(test.W(balPol@data$PRICE,bal.lags.q[[8]]))
q.lw.9<-nb2listw(bal.lags.q[[9]], style="W")
summary(test.W(balPol@data$PRICE,bal.lags.q[[9]]))

#Rook
r.lw.1<-nb2listw(bal.nb.r, style="W")
summary(test.W(balPol@data$PRICE,bal.nb.r))
r.lw.2<-nb2listw(bal.lags.r[[2]], style="W")
summary(test.W(balPol@data$PRICE,bal.lags.r[[2]]))
r.lw.3<-nb2listw(bal.lags.r[[3]], style="W")
summary(test.W(balPol@data$PRICE,bal.lags.r[[3]]))
r.lw.4<-nb2listw(bal.lags.r[[4]], style="W")
summary(test.W(balPol@data$PRICE,bal.lags.r[[4]]))
r.lw.5<-nb2listw(bal.lags.r[[5]], style="W")
summary(test.W(balPol@data$PRICE,bal.lags.r[[5]]))
r.lw.6<-nb2listw(bal.lags.r[[6]], style="W")
summary(test.W(balPol@data$PRICE,bal.lags.r[[6]]))
r.lw.7<-nb2listw(bal.lags.r[[7]], style="W")
summary(test.W(balPol@data$PRICE,bal.lags.r[[7]]))
r.lw.8<-nb2listw(bal.lags.r[[8]], style="W")
summary(test.W(balPol@data$PRICE,bal.lags.r[[8]]))
r.lw.9<-nb2listw(bal.lags.r[[9]], style="W")
summary(test.W(balPol@data$PRICE,bal.lags.r[[9]]))

#pesos k-vecinos
k.lw.1<-nb2listw(bal.kn1, style="W")
summary(test.W(balPol@data$PRICE,bal.kn1))
k.lw.2<-nb2listw(bal.kn2, style="W")
summary(test.W(balPol@data$PRICE,bal.kn2))
k.lw.3<-nb2listw(bal.kn3, style="W")
summary(test.W(balPol@data$PRICE,bal.kn3))
k.lw.4<-nb2listw(bal.kn4, style="W")
summary(test.W(balPol@data$PRICE,bal.kn4))
k.lw.5<-nb2listw(bal.kn5, style="W")
summary(test.W(balPol@data$PRICE,bal.kn5))
k.lw.6<-nb2listw(bal.kn6, style="W")
summary(test.W(balPol@data$PRICE,bal.kn6))
k.lw.7<-nb2listw(bal.kn7, style="W")
summary(test.W(balPol@data$PRICE,bal.kn7))
k.lw.8<-nb2listw(bal.kn8, style="W")
summary(test.W(balPol@data$PRICE,bal.kn8))
k.lw.9<-nb2listw(bal.kn9, style="W")
summary(test.W(balPol@data$PRICE,bal.kn9))


1-2*100/211
5.2#graficos
delaunay.lw<-nb2listw(trinb, style="W")
summary(test.W(balPol@data$PRICE,trinb))
esfera.inf.lw<-nb2listw(soinb, style="W")
summary(test.W(balPol@data$PRICE,soinb))
gabriel.lw<-nb2listw(gabrielnb, style="W")
summary(test.W(balPol@data$PRICE,gabrielnb))
vec.relative.lw<-nb2listw(relativenb, style="W")
summary(test.W(balPol@data$PRICE,relativenb))

0.512001474/13.13
moran.test(balPol@data$PRICE,delaunay.lw , zero.policy=T, alternative = "two.sided")
moran.test(balPol@data$PRICE, nb2listw(bal.lags.q[[1]], style="W", zero.policy=T), zero.policy=T, alternative = "two.sided")

Pesos.list<-list(reina1=q.lw.1,reina2=q.lw.2,reina3=q.lw.3,reina4=q.lw.4,reina5=q.lw.5,
                 reina6=q.lw.6,reina7=q.lw.7,reina8=q.lw.8,reina9=q.lw.9,
                 torre1=r.lw.1,torre2=r.lw.2,torre3=r.lw.3,torre4=r.lw.4,torre5=r.lw.5,
                 torre6=r.lw.6,torre7=r.lw.7,torre8=r.lw.8,torre9=r.lw.9,
                 kvecinos1=k.lw.1,kvecinos2=k.lw.2,kvecinos3=k.lw.3,kvecinos4=k.lw.4,
                 kvecinos5=k.lw.5,kvecinos6=k.lw.6,kvecinos7=k.lw.7,kvecinos8=k.lw.8,kvecinos9=k.lw.9,
                 gabriel=gabriel.lw,delaunay=delaunay.lw,esfera.inf=esfera.inf.lw,vec.relativos=vec.relative.lw)

class(Pesos.list)

nbw <- length(Pesos.list)
# al 5% de significanc?a
1 - (1 - 0.05)^(nbw)
# la optimizaci?n para la selecci?n se realiza maximizando el R cuadrado ajustado (R2 Adjust) o minimizando 
# la autocorrelaci?n espacial residual.
set.seed(123)
W_sel <- listw.select(balPol@data$PRICE, Pesos.list, MEM.autocor = "all", p.adjust = TRUE, nperm = 50)
W_sel$candidates
W_sel$best.id

# Para copiar y editar en Latex
xtable(W_sel$candidates, digits = 6)
W_sel$best$MEM.select
# El criterio que m?s valor tiene para la selecci?n definitiva de la matriz de contig?idad
# es el R cuadrado ajustado (RAdjust 2), que minimiza la autocorrelaci?n. Por esto la matriz
# elegida ser? K vecinos con 2 vecinos m?s cercanos.
Best.SWM <- Pesos.list[W_sel$best.id]   


PRICE <- as.data.frame(bal.poly)$PRICE
x <- c(PRICE)

# guardar como matriz
W <- as.matrix(as_dgRMatrix_listw(k.lw.2BIN))
W



##matriz seleccionada lunay##

bal.nb.q <- poly2nb(balPol) 
bal.nb=tri2nb(coords)
bal.lw <-nb2listw(bal.nb, style="W")
write.nb.gal(bal.nb, 'shapefile/balPol.gal')

###########################################################
###### AUTOCORRELACION ####################################
###########################################################
var <- as.data.frame(bal.poly)['PRICE']
xv <- as.numeric(as.character(unlist(c(var))))
x1v<-(xv-mean(xv))/sd(xv)  
set.seed(123)
moran.test(x1v, bal.lw, zero.policy=T, alternative="two.sided")

#### I DE MORAN GLOBAL PARA LA VARIABLE DEPENDIENTE

set.seed(123)
moran.test(x1, bal.lw, zero.policy=T, alternative="two.sided")

# Gr?fico de dispersi?n del ?ndice de Moran
png("D:/Documentos/11 semestre/Geoestadística/ProyectoDatosArea/Dispersion/PRICE.png",
    width = 2200, height = 1600,units="px",res=300)

mp<-moran.plot(x, bal.lw, main="Gráfico de Dispersión de Moran",  ylab="W_PRICE", xlab="PRICE")      # Cambiar por "x1", para la estandarizaci?n
dev.off()

# Correlograma Moran a partir de matriz contiguidad espacial
sp.cr <- sp.correlogram(bal.nb, x, order=9, method="C", style="W", zero.policy=T)
x11()
plot(sp.cr,main='Correlograma cr')
png(paste("D:/Documentos/11 semestre/Geoestadística/ProyectoDatosArea/Correlograma/correl_PRICE.png"),
    width = 2200, height = 1600,units="px",res=300)
cor.s <- sp.correlogram(trinb, x, order=9, method="I", style="W", zero.policy=T)
#x11()
plot(cor.s,main='Correlograma PRICE')
dev.off()

vars<-c('NROOM','DWELL','NBATH','PATIO','FIREPL','AC','BMENT','NSTOR','GAR','AGE','CITCOU','LOTSZ','SQFT')

table_bivariate(bivariadoGlobal$BMENT)
hist(PRICE)
#guarda histogramas
a=palette(rainbow(13)) 
col=1
for (i in vars)
{
  i<-'PRICE'
  j<-i
  var <- as.data.frame(bal.poly)[j]
  xv <- as.numeric(as.character(unlist(c(var))))
  png(paste("D:/Documentos/11 semestre/Geoestadística/ProyectoDatosArea/Histogramas/hist",i,".png"),
      width = 2200, height = 1600,units="px",res=300)
  hist(xv,main=paste("Histogram",i),col=a[col])
  dev.off()
  col=col+1
}

x11()

png(paste("D:/Documentos/11 semestre/Geoestadística/ProyectoDatosArea/Histogramas/hist-boxplot","PRICE",".png"),width = 3000, height = 1000,units="px",res=300)
op=par(mfrow=c(1,3))
hist(balPol$PRICE,main=paste("Histograma","PRICE"),col='#2d572c')
boxplot(balPol$PRICE,main=paste("BoxPlot","PRICE"),col='#2d572c')
qqPlot(balPol$PRICE,main=paste("QqPlot","PRICE"),col='#2d572c')
par(op)
dev.off()
hist(1:5, col="cornflowerblue",breaks=1:1)
a=palette(rainbow(13)) 
png(paste("D:/Documentos/11 semestre/Geoestadística/ProyectoDatosArea/Histogramas/hist","N1",".png"),width = 3000, height = 1000,units="px",res=300)
op=par(mfrow=c(1,3))

hist(balPol$NROOM,main=paste("Histograma",'NROOM'),col=a[1],xlab='min: 3    |    median: 5.182    |    mean: 5.182 \n 1stQ: 5    |    3rdQ: 6    |    max: 10')
hist(balPol$NBATH,main=paste("Histograma",'NBATH'),col=a[2],xlab='min: 1    |    median: 1.5    |    mean: 1.565 \n 1stQ: 1    |    3rdQ: 2    |    max: 5')
hist(balPol$NSTOR,main=paste("Histograma",'NSTOR'),col=a[3],xlab='min: 1    |    median: 2    |    mean: 1.904 \n 1stQ: 2    |    3rdQ: 2    |    max: 3')
par(op)
dev.off()
png(paste("D:/Documentos/11 semestre/Geoestadística/ProyectoDatosArea/Histogramas/hist","N2",".png"),width = 3000, height = 1000,units="px",res=300)
op=par(mfrow=c(1,3))
# summary(balPol$SQFT)
hist(balPol$GAR,main=paste("Histograma",'GAR'),col=a[4],xlab='min: 0    |    median: 0    |    mean: 2.344 \n 1stQ: 0    |    3rdQ: 0    |    max: 3')
hist(balPol$AGE,main=paste("Histograma",'AGE'),col=a[5],xlab='min: 0    |    median: 25    |    mean: 30.26 \n 1stQ: 0    |    3rdQ: 20    |    max: 148')
hist(balPol$LOTSZ,main=paste("Histograma",'LOTSZ'),col=a[6],xlab='min: 5.7    |    median: 56.12    |    mean: 70.89 \n 1stQ: 20.69    |    3rdQ: 84    |    max: 400.37')
par(op)
dev.off()
png(paste("D:/Documentos/11 semestre/Geoestadística/ProyectoDatosArea/Histogramas/hist","N3",".png"),width = 3000, height = 1000,units="px",res=300)
op=par(mfrow=c(1,3))
hist(balPol$SQFT,main=paste("Histograma",'SQFT'),col=a[7],xlab='min: 5.76    |    median: 13.44    |    mean: 16.23 \n 1stQ: 11.02    |    3rdQ: 19.6    |    max: 47.61')
hist(balPol$DWELL,main=paste("Histograma",'DWELL'),col=a[8],breaks=2,xlab= '1:111\n 0:98')
hist(balPol$PATIO,main=paste("Histograma",'PATIO'),col=a[1],breaks=2,xlab= '1:29\n 0:180')
par(op)
dev.off()
length(balPol$CITCOU[balPol$CITCOU == 1])
length(balPol$CITCOU[balPol$CITCOU == 0])
png(paste("D:/Documentos/11 semestre/Geoestadística/ProyectoDatosArea/Histogramas/hist","N4",".png"),width = 3000, height = 1000,units="px",res=300)
op=par(mfrow=c(1,4))
hist(balPol$FIREPL,main=paste("Histograma",'FIREPL'),col=a[2],breaks=2,xlab= '1:49\n 0:160')
hist(balPol$AC,main=paste("Histograma",'AC'),col=a[3],breaks=2,xlab= '1:50\n 0:159')
hist(balPol$BMENT,main=paste("Histograma",'BMENT'),col=a[4],breaks=2,xlab= '1:180\n 0:29')
hist(balPol$CITCOU,main=paste("Histograma",'CITCOU'),col=a[5],breaks=2,xlab= '1:126\n 0:83')
par(op)
dev.off()

#guarda coeficientes y graficas variables independietnes
iMoranGlobal<-list()
bivariadoGlobal<-list()
correlbiva<-list()
vars<-c('NROOM','DWELL','NBATH','PATIO','FIREPL','AC','BMENT','NSTOR','GAR','AGE','CITCOU','LOTSZ','SQFT')
X11()
png(paste("D:/Documentos/11 semestre/Geoestadística/ProyectoDatosArea/Correlograma/corr","N2",".png"),width = 3000, height = 3000,units="px",res=300)
op=par(mfrow=c(2,2))
for (i in vars)
{
  j<-i
  var <- as.data.frame(bal.poly)[j]
  xv <- as.numeric(as.character(unlist(c(var))))
  x1v<-(xv-mean(xv))/sd(xv) 

  moran.cluster(x1, bal.lw, zero.policy = T, bal.poly, significant=T)
  set.seed(123)
  iMoranGlobal[[i]]<-moran.test(x1v, bal.lw, zero.policy=T, alternative="two.sided")
  bivariadoGlobal[[i]]<-moranbi1.test(x=x1,y=x1v,bal.lw,zero.policy =T,randomisation =T,alternative="two.sided",adjust.n=TRUE)
  # 
  # png(paste("D:/Documentos/11 semestre/Geoestadística/ProyectoDatosArea/Dispersion/disp",i,".png"),
  #   width = 2200, height = 1600,units="px",res=300)
  moranbi.plot(x1,x1v,quiet =F,zero.policy =T,listw=bal.lw, xlab = "PRICE", ylab = j,main=paste("PRICE_",i))
  # dev.off()
  # png(paste("D:/Documentos/11 semestre/Geoestadística/ProyectoDatosArea/Correlograma/correl",i,".png"),
  #     width = 2200, height = 1600,units="px",res=300)
  correlbiva[[i]] <- spcorrelogram.bi(bal.nb, x1, x1v, order=9, method="I", style="W", zero.policy=T)
  plot(correlbiva[[i]],main=paste("PRICE_",i))
  # dev.off()
}
par(op)
dev.off()
#crear tabla i de moran global
z<-c()
p<-c()
I<-c()
E<-c()
V<-c()
for (i in vars)
{
  j<-i
  var <- as.data.frame(bal.poly)[j]
  xv <- as.numeric(as.character(unlist(c(var))))
  x1v<-(xv-mean(xv))/sd(xv)  
  set.seed(123)
  z<-c(z,moran.test(x1v, bal.lw, zero.policy=T, alternative="two.sided")$statistic)
  p<-c(p,moran.test(x1v, bal.lw, zero.policy=T, alternative="two.sided")$p.value)
  I<-c(I,moran.test(x1v, bal.lw, zero.policy=T, alternative="two.sided")$estimate["Moran I statistic"])
  E<-c(E,moran.test(x1v, bal.lw, zero.policy=T, alternative="two.sided")$estimate["Expectation"])
  V<-c(V,moran.test(x1v, bal.lw, zero.policy=T, alternative="two.sided")$estimate["Variance"])
}
moranglob.df <- data.frame("VAR" = vars, "Z" = z,"P-value"=p,"I de Moran"=I,"E(I)"=E,"V(I)"=V)
xtable(moranglob.df , digits = 6)

#crear tabla i de moran global
z<-c()
p<-c()
I<-c()
E<-c()
V<-c()
for (i in vars)
{
  j<-i
  var <- as.data.frame(bal.poly)[j]
  xv <- as.numeric(as.character(unlist(c(var))))
  x1v<-(xv-mean(xv))/sd(xv)  
  set.seed(123)
  z<-c(z,moranbi1.test(x=x1,y=x1v,bal.lw,zero.policy =T,randomisation =T,alternative="two.sided",adjust.n=TRUE)$statistic["Bivariate Moran Z(I) statistic"])
  p<-c(p,moranbi1.test(x=x1,y=x1v,bal.lw,zero.policy =T,randomisation =T,alternative="two.sided",adjust.n=TRUE)$p.value[1])
  I<-c(I,moranbi1.test(x=x1,y=x1v,bal.lw,zero.policy =T,randomisation =T,alternative="two.sided",adjust.n=TRUE)$estimate['Bivariate Moran I statistic'])
  E<-c(E,moranbi1.test(x=x1,y=x1v,bal.lw,zero.policy =T,randomisation =T,alternative="two.sided",adjust.n=TRUE)$estimate["Expectation"])
  V<-c(V,moranbi1.test(x=x1,y=x1v,bal.lw,zero.policy =T,randomisation =T,alternative="two.sided",adjust.n=TRUE)$estimate["Variance"])
}
morangBIglob.df <- data.frame("VAR" = vars, "Z" = z,"P-value"=p,"I de Moran Bivariado"=I,"E(I)"=E,"V(I)"=V)
xtable(morangBIglob.df , digits = 6)


#tablas del bivariado
table_bivariate <- function(name_0) {
  y0<-c("Z(I)","P-value","I de Moran","E(I)","V(I)")
  z<-name_0$statistic["Bivariate Moran Z(I) statistic"]
  p<-name_0$p.value[1]
  I<-name_0$estimate['Bivariate Moran I statistic']
  E<-name_0$estimate["Expectation"]
  V<-name_0$estimate["Variance"]
  y1<-c(z,p,I,E,V)
  table.df <- data.frame(y0,y1)
  xtable(data.frame(table.df$y0,table.df$y1), digits = 6)
}
moranBIglob.df <- data.frame("VAR" = vars, "Z" = z,"P-value"=p,"I de Moran"=I,"E(I)"=E,"V(I)"=V)
xtable(moranglob.df , digits = 6)


table_bivariate(bivariadoGlobal$DWELL)



#### I DE MORAN GLOBAL PARA LAS VARIABLES INDEPENDIENTES
iMoranGlobal
# I Moran bivariado global 
bivariadoGlobal
# Getis-Ord global G statistic
globalG.test(x, listw=bal.lw)   


# LISA Cluster Map VARIABLE DEPENDIENTE
x11()
moran.cluster(x1, bal.lw, zero.policy = T, balPol, significant=T)

# LISA Cluster Map VARIABLE INDEPENDIENTE
x11()
moran.cluster(x1, bal.lw, zero.policy = T, bal.poly, significant=T)
moran.cluster(x1.OWNH, pesosw, zero.policy = T, col.poly, significant=T)
moran.cluster(x1.POP65, pesosw, zero.policy = T, col.poly, significant=T)
moran.cluster(x1.UNEMP, pesosw, zero.policy = T, col.poly, significant=T)

# I de moran local para la variable independiente 

set.seed(123)
PRICE_LOCAL<-localmoran(x1,bal.lw,zero.policy = T)
local.df<-as.data.frame(PRICE_LOCAL)

sign<-local.df[['Pr(z > 0)']]< 0.05
localsign<-local.df[sign,]

length(local.df[['Pr(z > 0)']][local.df[['Pr(z > 0)']] <0.05])
#G de getis local 

lg<-localG(x1,listw=bal.lw,zero.policy = T)
localsign$getis<-c(lg)[sign]
xtable(localsign$getis)
###########################################################
###### AUTOCORRELACION BIVARIADA ##################
###########################################################
xtable(localsign) c