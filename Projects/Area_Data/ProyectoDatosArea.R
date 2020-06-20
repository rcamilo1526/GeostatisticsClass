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
#cargar funciones
load("C:/Users/rmartin/Documents/11 semestre/Geoestadística/FuncDatosArea/ImagenDatosDeArea")
set.seed(123)
setwd("C:/Users/rmartin/Documents/11 semestre/Geoestadística/ProyectoDatosArea")

#puntos
balPt<-readOGR('./baltimore','baltim')

x11()
plot(balPol, col='#c7d6fc')
plot(balPt,add=T,col="#3f9945",pch = 18)
#poligonos de voronoi
balPol <- voronoi(balPt)

bal.poly <- balPol
coords <- coordinates(bal.poly)
##Variable a explicar
PRICE <- as.data.frame(bal.poly)$PRICE
x <- c(PRICE)
x1<-(x-mean(x))/sd(x) 

writeOGR(bal.poly, dsn = "shapefile", layer = "map",
         driver = "ESRI Shapefile" )

A_uno <- function(x) {
  if (x > 0) {1} else {0}
}
bment2<-c()
for (i in c(bal.poly$BMENT)){
  bment2<-c(bment2,A_uno(i))
}
bal.poly$BMENT<-bment2

x11()
op=par(mfrow=c(1,2))
plot(balPol)
plot(coords)


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

sp.cr <- sp.correlogram(bal.kn9, balPol$PRICE, order=9, method="corr", style="W", zero.policy=T)
cor <- sp.correlogram(bal.kn9, balPol$PRICE, order=9, method="I", style="W", zero.policy=T)
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

#graficos
delaunay.lw<-nb2listw(trinb, style="W")
summary(test.W(balPol@data$PRICE,trinb))
esfera.inf.lw<-nb2listw(soinb, style="W")
summary(test.W(balPol@data$PRICE,soinb))
gabriel.lw<-nb2listw(gabrielnb, style="W")
summary(test.W(balPol@data$PRICE,gabrielnb))
vec.relative.lw<-nb2listw(relativenb, style="W")
summary(test.W(balPol@data$PRICE,relativenb))

moran.test(balPol@data$PRICE, nb2listw(bal.nb.q, style="W", zero.policy=T), zero.policy=T, alternative = "two.sided")
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



#Mejor matriz de pesos reina de orden 1

##matriz seleccionada REINA 1##

bal.nb.q <- poly2nb(balPol) 

bal.lw<-nb2listw(bal.nb.q, style="W")

###########################################################
###### AUTOCORRELACION ####################################
###########################################################
var <- as.data.frame(bal.poly)['PATIO']
xv <- as.numeric(as.character(unlist(c(var))))
x1v<-(xv-mean(xv))/sd(xv)  
set.seed(123)
moran.test(x1v, bal.lw, zero.policy=T, alternative="two.sided")

#### I DE MORAN GLOBAL PARA LA VARIABLE DEPENDIENTE

set.seed(123)
moran.test(x1, bal.lw, zero.policy=T, alternative="two.sided")

# Gr?fico de dispersi?n del ?ndice de Moran
png("D:/Documentos/11 semestre/Geoestadística/ProyectoDatosArea/Dispersion/test.png",
    width = 2200, height = 1600,units="px",res=300)

mp<-moran.plot(xv, bal.lw, main="Gráfico de Dispersión de Moran",  ylab="W_PRICE", xlab="PRICE")      # Cambiar por "x1", para la estandarizaci?n
dev.off()

# Correlograma Moran a partir de matriz contiguidad espacial
sp.cr <- sp.correlogram(bal.nb.q, x, order=9, method="C", style="W", zero.policy=T)
x11()
plot(sp.cr,main='Correlograma cr')
png(paste("D:/Documentos/11 semestre/Geoestadística/ProyectoDatosArea/Correlograma/correl_PRICE.png"),
    width = 2200, height = 1600,units="px",res=300)
cor.s <- sp.correlogram(bal.nb.q, x, order=9, method="I", style="W", zero.policy=T)
#x11()
plot(cor.s,main='Correlograma PRICE')
dev.off()

vars<-c('NROOM','DWELL','NBATH','PATIO','FIREPL','AC','BMENT','NSTOR','GAR','AGE','CITCOU','LOTSZ','SQFT')

table_bivariate(bivariadoGlobal$SQFT)

#guarda histogramas
a=palette(rainbow(13)) 
col=1
for (i in vars)
{
  j<-i
  var <- as.data.frame(bal.poly)[j]
  xv <- as.numeric(as.character(unlist(c(var))))
  png(paste("D:/Documentos/11 semestre/Geoestadística/ProyectoDatosArea/Histogramas/hist",i,".png"),
      width = 2200, height = 1600,units="px",res=300)
  hist(xv,main=paste("Histogram",i),col=a[col])
  dev.off()
  col=col+1
}


#guarda coeficientes y graficas variables independietnes
iMoranGlobal<-list()
bivariadoGlobal<-list()
correlbiva<-list()
for (i in vars)
{
  j<-i
  var <- as.data.frame(bal.poly)[j]
  xv <- as.numeric(as.character(unlist(c(var))))
  x1v<-(xv-mean(xv))/sd(xv)  
  set.seed(123)
  iMoranGlobal[[i]]<-moran.test(x1v, bal.lw, zero.policy=T, alternative="two.sided")
  bivariadoGlobal[[i]]<-moranbi1.test(x=x1,y=x1v,bal.lw,zero.policy =T,randomisation =T,
                                      alternative="two.sided",adjust.n=TRUE)
  
  png(paste("D:/Documentos/11 semestre/Geoestadística/ProyectoDatosArea/Dispersion/disp",i,".png"),
    width = 2200, height = 1600,units="px",res=300)
  moranbi.plot(x1,x1v,quiet =F,zero.policy =T,listw=bal.lw, xlab = "PRICE", ylab = j,main=paste("PRICE_",i))
  dev.off()
  png(paste("D:/Documentos/11 semestre/Geoestadística/ProyectoDatosArea/Correlograma/correl",i,".png"),
      width = 2200, height = 1600,units="px",res=300)
  correlbiva[[i]] <- spcorrelogram.bi(bal.nb.q, x1, x1v, order=9, 
                                method="I", style="W", zero.policy=T)
  plot(correlbiva[[i]],main=paste("PRICE_",i))
  dev.off()
}
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

table_bivariate(bivariadoGlobal$DWELL)



#### I DE MORAN GLOBAL PARA LAS VARIABLES INDEPENDIENTES
iMoranGlobal
# I Moran bivariado global 
bivariadoGlobal

# LISA Cluster Map VARIABLE DEPENDIENTE
x11()
moran.cluster(x1, pesosw, zero.policy = T, col.poly, significant=T)

# LISA Cluster Map VARIABLE INDEPENDIENTE
x11()
moran.cluster(x1, bal.lw, zero.policy = T, bal.poly, significant=T)
moran.cluster(x1.OWNH, pesosw, zero.policy = T, col.poly, significant=T)
moran.cluster(x1.POP65, pesosw, zero.policy = T, col.poly, significant=T)
moran.cluster(x1.UNEMP, pesosw, zero.policy = T, col.poly, significant=T)

# I de moran local para la variable independiente 

set.seed(123)
PRICE_LOCAL<-localmoran(x1,bal.lw,zero.policy = T)
PRICE_LOCAL

#G de getis local 

lg<-localG(x1,listw=pesosw,zero.policy = T)
lg

###########################################################
###### AUTOCORRELACION BIVARIADA ##################
###########################################################




############### MODELOS ###########3
#PRICE	sales price of house in $1,000 (MLS)
#NROOM	number of rooms
#DWELL	1 if detached unit, 0 otherwise
#NBATH	number of bathrooms
#PATIO	1 if patio, 0 otherwise
#FIREPL	1 if fireplace, 0 otherwise
#AC	1 if air conditioning, 0 otherwise
#BMENT	1 if basement, 0 otherwise
#NSTOR	number of stories
#GAR	number of car spaces in garage (0 = no garage)
#AGE	age of dwelling in years
#CITCOU	1 if dwelling is in Baltimore County, 0 otherwise
#LOTSZ	lot size in hundreds of square feet
#SQFT	interior living space in hundreds of square feet
attach(balPol@data)
mod.lin <- lm(PRICE ~ NROOM+DWELL+NBATH+PATIO+FIREPL+AC+BMENT+NSTOR+GAR+AGE+CITCOU+LOTSZ+SQFT, data = bal.poly@data)
summary(mod.lin)
models <- regsubsets(PRICE~., data = balPol@data, nvmax = 13)
summary(models)


lm.morantest(mod.lin, pesosw, alternative = "greater")
LM_W1 <- lm.LMtests(mod.lin, listw = pesosw, test = "all")
t(sapply(LM_W1, function(x) unlist(x[1:3])))


LMC_AC <- localmoran.bi(PRICE, AC, pesosw, zero.policy =T)
LMCH <- localmoran.bi(PRICE, AGE, pesosw, zero.policy =T)


bal.poly$lm_res <- residuals(mod.lin)
# Mapa residuales: Opci?n 1
spplot(bal.poly["lm_res"], col.regions = rev(terrain.colors(20)))

# Mapa residuales: Opci?n 2
pal2 <- colorRampPalette(c("red3", "wheat1", "blue3"))
spplot(bal.poly,"lm_res", col.regions = pal2(20))

# Mapa residuales: Opci?n 3
my.palette <- brewer.pal(n = 9, name = "YlOrRd")
spplot(bal.poly, "lm_res", col.regions = my.palette, cuts = 8, col = "transparent")

# Mapa residuales: Opci?n 4
breaks.ci <- classIntervals(bal.poly$lm_res, n = 9, style = "quantile", intervalClosure = "right")$brks
# se amplian los limites inferior y superior en: .Machine$double.eps*5000000000000,  con el fin que un poligono no quede siempre en blanco.
breaks.ci[1] <- breaks.ci[1] - .Machine$double.eps*5000000000000
breaks.ci[length(breaks.ci)] <- breaks.ci[length(breaks.ci)] + .Machine$double.eps*5000000000000
spplot(bal.poly, "lm_res", col = "transparent", col.regions = my.palette,  at = breaks.ci)


# Validaci?n supuestos
bptest(mod.lin)
resettest(mod.lin)
raintest(mod.lin)
shapiro.test(residuals(mod.lin))
vif(mod.lin)



# Modelo Spatial Lag:                                                                   
col.lag.sm <- lagsarlm(PRICE ~ NROOM+DWELL+NBATH+PATIO+FIREPL+AC+BMENT+NSTOR+GAR+AGE+CITCOU+LOTSZ+SQFT, data=bal.poly@data, listw=pesosw) # CRIME ~1
summary(col.lag.sm, Nagelkerke=T)
predict.sarlm(col.lag.sm)
AIC(col.lag.sm)
deviance.sarlm(col.lag.sm)
residuals.sarlm(col.lag.sm)
coef.sarlm(col.lag.sm)
fitted.sarlm(col.lag.sm)
bptest.sarlm(col.lag.sm)
hetero.plot <- function(model) {
  plot(residuals(model) ~ fitted(model))
  abline(h=0, lty="dotted")
  lines(lowess(fitted(model), residuals(model)), col="red")
}
hetero.plot(col.lag.sm)

#Nagelkerke NJD (1991) A note on a general definition of the coefficient of determination. Biometrika 78: 691?C692. 
NK <- function(obj, y) { 
  n <- length(obj$residuals) 
  nullLL <- logLik(lm(y ~ 1)) 
  c(1 - exp(-(2/n)*(logLik(obj) - nullLL))) 
} 
NK(col.lag.sm,PRICE)

#bc <- boxCox(CRIME ~ INC + HOVAL, data = columbus, lambda = seq(-1, 2, length = 20))
#bc$x[which.max(bc$y)]

# Modelo Spatial Error:
col.error.sm <- errorsarlm(PRICE ~ NROOM+DWELL+NBATH+PATIO+FIREPL+AC+BMENT+NSTOR+GAR+AGE+CITCOU+LOTSZ+SQFT, data=bal.poly@data, listw=pesosw)
summary(col.error.sm, Nagelkerke=T)
NK(col.error.sm,PRICE)
predict.sarlm(col.error.sm)
AIC(col.error.sm)
deviance.sarlm(col.error.sm)
residuals.sarlm(col.error.sm)
coef.sarlm(col.error.sm)
fitted.sarlm(col.error.sm)
bptest.sarlm(col.error.sm)
hetero.plot(col.error.sm)

# Comparaci?n
anova(col.lag.sm,col.error.sm)    
