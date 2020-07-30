library(sp)
load("/home/rcamilo/Documents/geoestadistica/IDSTA.ov.rda")
load("/home/rcamilo/Documents/geoestadistica/croatia.grilla 25s.2008.rda")
plot(IDSTA.OV)
library(geospt)
# load("~/Geoestadistica/Geoestad?stica/NOSOTROS/variables de codigo.RData")

#Eliminar dato vacio
options(max.print = 99999)
IDSTAr.ov <- IDSTA.ov[!is.na(IDSTA.ov$LST2008_02_26),]
# IDSTAr.ov<-IDSTA.ov[-c(54, 106, 117, 130),]

Hrdem<-IDSTAr.ov$HRdem
Hrdsea<-IDSTAr.ov$HRdsea
Hrtwi<-IDSTAr.ov$HRtwi
Lat<-IDSTAr.ov$Lat
Lon<-IDSTAr.ov$Lon
Tmt<-IDSTAr.ov$LST2008_02_26
shapiro.test(Tmt)
library(geoR)
croacia.geoR <- as.geodata(IDSTAr.ov, data.col = 13) 
E<-croacia.geoR$coords[,2]
N<-croacia.geoR$coords[,1]
names(IDSTAr.ov)
      # #Grafico en 2D con estes y nortes

   # Objeto del tipo geodata (coordenadas y datos)
x11()
points.geodata(croacia.geoR, x.leg=400000, y.leg=4830000,  pt.div="quintile")
# Graficar los datos (puntos) y ver simbolos graduados, adem?s la opci?n "add.to.plot" es ?til para adicionar opciones de la instrucci?n "plot", como en este caso la grilla.
points.geodata(croacia.geoR, x.leg=400000, y.leg=4830000, main=c("Gr?fico de Intensidades", "de Precipitaci?n"), col.main=3, pt.div="quintile", add.to.plot = TRUE, panel.first = grid())  
x11()
plot(croacia.geoR)
length(N)
length(Tmt)
# Para hacer gr?ficos en 3 dimensiones                                    
library(scatterplot3d)               
s3d<-
x11()
scatterplot3d(N, E, Tmt, angle=45, col.main=4, xlab="Coordenada X", ylab="Coordenada Y", zlab="Temperatura", pch=19, color="firebrick3")

IDSTAr.ov[,c('HRdem','HRdsea')]
# Resumen de Estad?sticas descriptivaas.

summary(IDSTAr.ov)
summary(croacia.geoR$coords)

#### C?lculo de la desviaci?n est?ndar de la variable y ploteo de gr?ficos exploratorios
sd(Tmt)
x11()
par(mfrow=c(nrow=1, ncol=3))
hist(Tmt, freq=F, breaks=10, xlab="Temperatura (?C)", ylab="Frecuencia", col="orange", main ="Histograma")
curve(dnorm(x, mean=mean(Tmt), sd=sd(Tmt)), add=T, col="red")                  # X no funciona, solo minuscula
qqnorm(Tmt, ylab="Temperatura (?C)", xlab="Cuantiles te?ricos", col="orange")
qqline(Tmt)                                                        #Agregarle al Q-Q Plot la l?nea media
boxplot(Tmt, main="BOX-PLOT", notch=F, horizontal=T, xlab="Temperatura (?C)",col="orange")
points(mean(Tmt), y=1, pch=1, cex=2)
par(mfrow=c(nrow=1, ncol=1))
min(Tmt)

shapiro.test(Tmt)
n.q <- (Tmt-mean(Tmt))/sd(Tmt)            #Funci?n prec
n.o <- order(n.q)                            # Lista con las posiciones de los datos ordenados
n.qo <- n.q[n.o]                             # Vector de cuantiles estandarizados y ordenados. Necesarios para prueba K-S
ks.test(n.qo, pnorm)                         # Le hago prueba K-S para saber si los datos provienen de una normal
ks.test(as.numeric(scale(sort(Tmt))),pnorm)
detach("package:geoR")


##An?lisis de tendencia
library(MASS)
modelo1<-lm(Tmt~Hrdem+Hrdsea+Hrtwi+N+E)
summary(modelo1)
stargazer(modelo1,single.row = T)
library(texreg)
step.model <- stepAIC(modelo1, direction = "both", trace = FALSE)
stargazer(step.model,single.row = T)
summary(step.model)
residuals(modelo1)
#creo otro dataframe con la variable transformada
data2<-cbind(N, E, Tmt)
temp2<-as.data.frame(data2)

library(spatial)
#aqui cargamos de nuevo la base da datos para poder tener una lista
tp.ls<-surf.ls(3,N,E,Tmt)
tp.trsurf<-trmat(tp.ls, -150, 120, 0, 190, 100) #dio el algo en cada pixel
#al profe no le agrada ni 5
# Representaci?n de superficies de tendencia
par(pty="s",mar=c(2,2,2,2)) #para la salida
contour(tp.trsurf) #gr?fico de contornos como topograf?a 
points(Lon,Lat,pch=20) #superponemos los puntos
par(mar=c(0,0,0,0)) #si quiero que el borde sea completamente nulo
image(tp.trsurf)
points(N,E,pch=20)
par(mfrow=c(1,1)) #perspectiva inclinado
persp(tp.trsurf)
persp(tp.trsurf,theta=30,phi=30,col=2,ltheta=-20,shade=0.25,xlab="N",ylab="E")
#estamos viendo que tiene valores altos en una direcci?n y valores bajos en otra 

##anisotropia
library(intamap)
coordinates(data.geo) <- ~E+N
estimateAnisotropy(data.geo,"head") 
#no hay anisotrop?a 

var4 <- variog4(tp.geo, max.dist=dmax)
x11()
plot(var4) #Parece que si hay pero pues no
install.packages('intamap')

datares<-cbind(E, N, tp.geo$data)
data.geo<-as.data.frame(datares)

library(gstat)
library(geoR)
library(sgeostat)
library(geospt)

data2<-cbind(N, E, Tmt)
temp2<-as.data.frame(data2)
data2.geo<-as.data.frame(croacia.geoR)
class(data2.geo) <- c("data.frame")
data.or<-data2.geo[,c(2,1,3)]

names(data.or) <- c("x", "y", "t")
dmax<-max(dist(croacia.geoR$coords))/2 #
point <- point(data.or)
pair <- pair(point,num.lags=30,maxdist=dmax)

#intento # 492
newgeodata<-croacia.geoR

v <- est.variograms(point,pair,"t",trim=0.1) 
v <- variog(newgeodata, coords = newgeodata$coords , data = newgeodata$data,max.dist = dmax,
            vec = "default", breaks = "default",trend = "cte", lambda = 1,
            option = c("bin", "cloud", "smooth"),estimator.type = c("classical", "modulus"))


####Gr?ficos de semivariogramas experimentales
x11()
dev.off()
x11()

layout(Conf2x2)
op=par(mfrow=c(2,2))
plot(v$bins,v$classic,pch=19, col="cyan3", main="Modelo experimental cl?sico", ylab="Semivarianza", xlab="Distancia")
lines(v$bins,v$classic,col=1,lty=1)
plot(v$bins,v$robust,pch=19,col="deeppink3", main="Modelo experimental robusto", ylab="Semivarianza", xlab="Distancia")
lines(v$bins,v$robust,col=1,lty=1)
plot(v$bins,v$med,pch=19,col="darkorchid3", main="Modelo experimental mediana", ylab="Semivarianza", xlab="Distancia")
lines(v$bins,v$med,col=1,lty=1)
plot(v$bins,v$trimmed.mean,pch=19,col="chartreuse3", main="Modelo experimental media recortada", ylab="Semivarianza", xlab="Distancia")
lines(v$bins,v$trimmed.mean,col=1,lty=1)
par(op)
##Se supone que ahora vamos a hacer los te?ricos
vrobus<-v$bins
x11()

eyefit(v, silent = TRUE)
abline(h=var(v),lty=2)

##################################################################
###################    CLASICO   #################################
#################################################################
data.frame(v$bins)
# exp.ml.c<-likfit(geodata = newgeodata , nugget = 0.1,  ini = c(0.8,40000),fix.nugget = T)
# plot(v$bins, v$classic,col=1, main="Modelo experimental cl?sico", ylim=c(0,13),ylab="Semivarianza", xlab="Distancia", pch=19)
# lines(exp.ml.c,max.dist=dmax,lwd=2,col='red')
exp.ml.c<-likfit(geodata = newgeodata , nugget = 0.1 , ini = c(10,40000),fix.nugget = T)
lines(exp.ml.c,max.dist=dmax,lwd=2,col='red')
sph.ml.c<-likfit(geodata = newgeodata, nugget = 0.01, ini = c(8,60000),cov.model="sph",fix.nugget = T)
lines(sph.ml.c,max.dist=dmax,lwd=2,col='blue')
mat.ml.c<-likfit(geodata = newgeodata,nugget = 1, ini = c(20,200000),cov.model="mat",kappa=1.5,fix.nugget = T)
lines(mat.ml.c,max.dist=dmax,lwd=2,col='green')
cir.ml.c<-likfit(geodata = newgeodata, nugget = 0.8,ini = c(4,80000),cov.model="cir",fix.nugget = T)
lines(cir.ml.c,max.dist=dmax,lwd=2,col='yellow')
pow.ml.c<-likfit(geodata = newgeodata,nugget = 5, ini = c(5,75000),cov.model="powered.exponential",kappa=1.75,fix.nugget = T)
lines(pow.ml.c,max.dist=dmax,lwd=2,col='orange')

plot(v$bins, v$classic,col=1, main="Modelo experimental clásico", ylab="Semivarianza", xlab="Distancia", pch=19,ylim=c(0,17))
lines(exp.ml.c,max.dist=dmax,lwd=2,col='red')
lines(sph.ml.c,max.dist=dmax,lwd=2,col='blue')
lines(mat.ml.c,max.dist=dmax,lwd=2,col='green')
lines(cir.ml.c,max.dist=dmax,lwd=2,col='yellow')
lines(pow.ml.c,max.dist=dmax,lwd=2,col='orange')
legend(locator(1), legend=c('Exponencial','Esf?rico','Matern','Cir','Pow'),lwd=c(2,2,2,2),col=c('red','blue','green','yellow','orange'))


#Aic para todos los variogramas anteriores
exp.ml.c$AIC
sph.ml.c$AIC
mat.ml.c$AIC
cir.ml.c$AIC
pow.ml.c$AIC

#A continuacion se ajustan algunos semivariogramas teoricos al semivariograma experimental  
length(v$bins)
dir.hor <- seq(0, 0, length.out=30)
dir.ver <- seq(0, 0, length.out=30)
id <- seq (length.out=30)
id <- rep("var1",30)

##ESFERICO CLASICO
y <- data.frame (v$n, v$bins,v$classic,dir.hor,dir.ver,id)
names(y) <- c("np", "dist", "gamma", "dir.hor","dir.ver","id")
class(y) <- c("variogram","gstatVariogram","data.frame")
croacia.geoR  # Objeto del tipo geodata (coordenadas y datos)
vgm('Sph')
#Sph.ml <-likfit(preci.geoR, ini = c(250, 3), cov.model="sph")  # metodo MV
Sph.wls.0 <- fit.variogram(y, vgm(10, "Sph", 50000), fit.method = 2)  
dist.s <- v$bins
Sph.ML <- variogramLine(vgm(10.08552, "Sph", 11258.21), min=0, dist_vector=dist.s)
Sph.ols <- fit.variogram(y, vgm(250, "Sph", 3),fit.method = 6) # metodo 6 MCO
#######################################################
################### ROBUSTO ##########################
######################################################
"Exp", "Sph", "Gau", "Mat"
exp.ml.r<-likfit(geodata = newgeodata , nugget = 0.1 , ini = c(10,40000),fix.nugget = T)
sph.ml.r<-likfit(geodata = newgeodata, nugget = 0.01, ini = c(0.6,50000),cov.model="sph",fix.nugget = T)
mat.ml.r<-likfit(geodata = newgeodata,nugget = 0.25, ini = c(0.6,35000),cov.model="mat",kappa=1.5,fix.nugget = T)
cir.ml.r<-likfit(geodata = newgeodata, nugget = 0.001,ini = c(0.6,31000),cov.model="cir",fix.nugget = T)
pow.ml.r<-likfit(geodata = newgeodata,nugget = 0.05, ini = c(0.6,20000),cov.model="powered.exponential",kappa=1.75,fix.nugget = T)

#Aic para todos los variogramas anteriores
exp.ml.r$AIC
sph.ml.r$AIC
mat.ml.r$AIC
cir.ml.r$AIC
pow.ml.r$AIC


plot(v$bins, v$robust, ylim=c(0,15),col=1, main="Modelo experimental robusto", ylab="Semivarianza", xlab="Distancia", pch=19)
lines ( Sph.ML, lwd=2,col='red')
lines(exp.ml.r,max.dist=dmax,lwd=2,col='red')
lines(sph.ml.r,max.dist=dmax,lwd=2,col='blue')
lines(mat.ml.r,max.dist=dmax,lwd=2,col='green')
lines(cir.ml.r,max.dist=dmax,lwd=2,col='yellow')
lines(pow.ml.r,max.dist=dmax,lwd=2,col='orange')
legend(locator(1), legend=c('Exponencial','Esf?rico','Matern','Cir','Pow'),lwd=c(2,2,2,2),col=c('red','blue','green','yellow','orange'))

##########################################################################
###################### MEDIANA###################################
##################################################################
exp.ml.m<-likfit(geodata = newgeodata ,nugget = 0.1 , ini = c(10,40000),fix.nugget = T)
sph.ml.m<-likfit(geodata = newgeodata, nugget = 0.01, ini = c(0.6,50000),cov.model="sph",fix.nugget = T)
mat.ml.m<-likfit(geodata = newgeodata,nugget = 0.25, ini = c(0.6,35000),cov.model="mat",kappa=1.5,fix.nugget = T)
cir.ml.m<-likfit(geodata = newgeodata, nugget = 0.001,ini = c(0.6,31000),cov.model="cir",fix.nugget = T)
pow.ml.m<-likfit(geodata = newgeodata,nugget = 0.05, ini = c(0.6,20000),cov.model="powered.exponential",kappa=1.75,fix.nugget = T)

#Aic para todos los variogramas anteriores
exp.ml.m$AIC
sph.ml.m$AIC
mat.ml.m$AIC
cir.ml.m$AIC
pow.ml.m$AIC


plot(v$bins, v$med, ylim=c(0,15),col=1, main="Modelo experimental mediana", ylab="Semivarianza", xlab="Distancia", pch=19)
lines(exp.ml.m,max.dist=dmax,lwd=2,col='red')
lines(sph.ml.m,max.dist=dmax,lwd=2,col='blue')
lines(mat.ml.m,max.dist=dmax,lwd=2,col='green')
lines(cir.ml.m,max.dist=dmax,lwd=2,col='yellow')
lines(pow.ml.m,max.dist=dmax,lwd=2,col='orange')
legend(locator(1), legend=c('Exponencial','Esf?rico','Matern','Cir','Pow'),lwd=c(2,2,2,2),col=c('red','blue','green','yellow','orange'))


##########################################################################
###################### MEDIA RECORTADA ###################################
##################################################################
exp.ml.mr<-likfit(geodata = newgeodata , nugget = 0.1 , ini = c(10,40000),fix.nugget = T)
sph.ml.mr<-likfit(geodata = newgeodata, nugget = 0.01, ini = c(0.6,55000),cov.model="sph",fix.nugget = T)
mat.ml.mr<-likfit(geodata = newgeodata,nugget = 0.25, ini = c(0.6,35000),cov.model="mat",kappa=1.5,fix.nugget = T)
cir.ml.mr<-likfit(geodata = newgeodata, nugget = 0.001,ini = c(0.6,28000),cov.model="cir",fix.nugget = T)
pow.ml.mr<-likfit(geodata = newgeodata,nugget = 0.05, ini = c(0.6,20000),cov.model="powered.exponential",kappa=1.75,fix.nugget = T)

#Aic para todos los variogramas anteriores
exp.ml.mr$AIC
sph.ml.mr$AIC
mat.ml.mr$AIC
cir.ml.mr$AIC
pow.ml.mr$AIC


plot(v$bins, v$med, ylim=c(0,15),col=1, main="Modelo experimental media recortada", ylab="Semivarianza", xlab="Distancia", pch=19)
lines(exp.ml.mr,max.dist=dmax,lwd=2,col='red')
lines(sph.ml.mr,max.dist=dmax,lwd=2,col='blue')
lines(mat.ml.mr,max.dist=dmax,lwd=2,col='green')
lines(cir.ml.mr,max.dist=dmax,lwd=2,col='yellow')
lines(pow.ml.mr,max.dist=dmax,lwd=2,col='orange')
legend(locator(1), legend=c('Exponencial','Esf?rico','Matern','Cir','Pow'),lwd=c(2,2,2,2),col=c('red','blue','green','yellow','orange'))

#se corre el modelo escog?do con ML WLS, RML

variogramaclasico<-variog( newgeodata,max.dist=dmax, option = "cloud")
exp.ml<-likfit(geodata = newgeodata , nugget = 0.122,  ini = c(0.65,18000),fix.nugget = T)
exp.rml<-likfit(geodata = newgeodata , nugget = 0.1,  ini = c(0.65,19000),fix.nugget = T,method='RML')
exp.wls<-variofit(vario = variogramaclasico, nugget = 0.08,  ini = c(0.62,18000),fix.nugget = T,weights="npairs")

plot(v$bins, v$robust, ylim=c(0,0.8), pch=19, xlab="Distancia", ylab="Semivarianza")
lines(exp.ml,max.dist=dmax,lwd=2, col='red')
lines(exp.rml,max.dist=dmax,lwd=2,col='blue')
lines(exp.wls,max.dist=dmax,lwd=2,col='green')
legend(locator(1),legend=c('ML','RML','WLS'),lwd=c(2,2,2),col=c('red','blue','green'))

#############################
#######MAPA BASE#############
##############################
#PUNTOS ENCIMA DEL SHAPE
IDSTA.ll <- spTransform(IDSTA.OV, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
points(IDSTA.ll)


#install devtools
install.packages("devtools")
#install ggmap from dev space
devtools::install_github("dkahle/ggmap")
devtools::install_github("dkahle/ggmap",force = TRUE)

library(ggmap)
# install devtools
install.packages("devtools")
o# install ggmap from dev space
devtools::install_github("dkahle/ggmap", force = T)
library(devtools)
register_google(key = "*", write=T, account_type = "standard" )
has_google_key()
key="*"
google_key()
gc <- geocode("Croatia, Croatia")
myMap <- get_map(location = "Croatia",
                 source = "google",
                 maptype = "terrain", crop = FALSE,
                 zoom = 6)


myMap <- get_map(location = "Boulder, Colorado",
                 source = "google",
                 maptype = "terrain", crop = FALSE,
                 zoom = 6)
ggmap()
# plot map
map <- get_map(location = 'Australia', zoom = 6)
ggmap(myMap)
geocode("waco, texas", urlonly = TRUE)
ggmap_show_api_key()


ggmap(
  ggmap = get_map(
    "Dayton",
    zoom = 13, scale = "auto",
    maptype = "satellite",
    source = "google"),
  extent = "device",
  legend = "topright"
)



temp2.sp <- SpatialPoints(as.data.frame(temp2)) 
umtna<-"+NA"
proj4string(temp2.sp) <- CRS(umtna)

#################################
######## KRIGING ORDINARIO#######
#################################

#02. Variogram fit to gaussian transformed data using gstat 
zn.svgm <- variogram(Tmt~1,temp2.sp)
vgm1 = vgm(0.544,"Exp",19000, 0.01)
plot(zn.svgm,vgm1)
ko.gstat <- gstat(id="Tmt", formula=Tmt~1, model=vgm1, data=temp2.sp,nmin=2, nmax=30)

#03. Kriging predictions with grid (gstat)
grillasin<-as.data.frame(IDSTA.OV)
coordinates(grillasin) = ~x+y
gridded(grillasin)=T
croatia.ks.ok = predict(ko.gstat,grillasin) 

#04. Back transform using RGeostats
ko.bts = cbind(coordinates(croatia.ks.ok),croatia.ks.ok@data)
ko.bts.db = db.create(ko.bts,autoname = F)
ko.tempdb = anam.y2z(ko.bts.db,names="Tmt.pred",anam = mdb.herm)
ko.tempdbvar = anam.y2z(ko.bts.db,names="Tmt.var",anam = mdb.herm)
#Prediction map
IDSTA.OV@data$pred.ko.temp <- ko.tempdb@items$Raw.Tmt.pred
spplot(IDSTA.OV,"pred.ko.temp")
#Var map
IDSTA.OV@data$var.ko.temp <- ko.tempdbvar@items$Raw.Tmt.var
spplot(IDSTA.OV,"var.ko.temp")
IDSTA.OV@data$var.ko.temp.sqrt <- sqrt(ko.tempdbvar@items$Raw.Tmt.var)


##############################
######## KRIGING SIMPLE#######
##############################

mu<-mean(temp2.sp$Tmt)
ks.gstat <- gstat(id="Tmt", formula=Tmt~1, model=vgm1, data=temp2.sp,nmin=2, nmax=30,beta=mu)

#03. Kriging predictions with grid (gstat)
grillasin<-as.data.frame(IDSTA.OV)
coordinates(grillasin) = ~x+y
gridded(grillasin)=T
croatia.ks.ok = predict(ks.gstat,grillasin) 

#04. Back transform using RGeostats
gridded(grillasin)=T
ks.bts = cbind(coordinates(croatia.ks.ok),croatia.ks.ok@data)
ks.bts.db = db.create(ks.bts,autoname = F)
ks.tempdb = anam.y2z(ks.bts.db,names="Tmt.pred",anam = mdb.herm)
ks.tempdbvar = anam.y2z(ks.bts.db,names="Tmt.var",anam = mdb.herm)
#Prediction map
IDSTA.OV@data$pred.ks.temp <-ks.tempdb@items$Raw.Tmt.pred
spplot(IDSTA.OV,"pred.ks.temp")
#Var map
IDSTA.OV@data$var.ks.temp.sqrt <- sqrt(ks.tempdbvar@items$Raw.Tmt.var)
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

write.asciigrid(IDSTA.OV@data$var.ks.temp, "ko.asc")
write.asciigrid(pred.ks.temp, "ko.asc")




#########################################################################################
##############                  M?TODOS DETERMINISTICOS:                  ###############
##############                  FUNCIONES DE BASE RADIAL                  ###############
#########################################################################################

coordNE<-as.data.frame(croacia.geoR$coords[,1:2])
temperat<-IDSTAr.ov$LST2008_01_25
coordinates(temp) <- ~N+E


grillabaser<-as.data.frame(IDSTA.OV)
puntosxy<-grillabaser[,57:58]  #grillabaser[,57:58] 
names(grillasin)<-c("x","y")
coordinates(grillasin)<- ~x+y
z=temp$Tmt


###### Multicuadratica##########
pred.rbf.m <- rbf(Tmt~1, temp , eta=0.00001,rho=0, newdata=puntosxy, n.neigh=9, func="M")
coordinates(pred.rbf.m) = c("x", "y")
gridded(pred.rbf.m) <- TRUE

# muestra el mapa de prediccion con la rbf multicuadratica
IDSTA.OV@data$pred.rbf.m <- pred.rbf.m$var1.pred
spplot(IDSTA.OV, "pred.rbf.m", main="Temperatura media 25 de enero de 2008 Predicciones \nMulticuadratica", col.regions=bpy.colors(100), cuts=10, cex.main=0.7, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=0.6))


###### Multicuadratica inversa##########
pred.rbf.mi3 <- rbf(Tmt~1, temp , eta=2,rho=0, newdata=puntosxy, n.neigh=9, func="IM")
coordinates(pred.rbf.mi) = c("x", "y")
gridded(pred.rbf.mi) <- TRUE

# muestra el mapa de prediccion con la rbf multicuadratica Inversa
IDSTA.OV@data$pred.rbf.mi <- pred.rbf.mi$var1.pred
spplot(IDSTA.OV, "pred.rbf.mi", main="Temperatura media 25 de enero de 2008 Predicciones \nMUlticuadratica Inversa", col.regions=bpy.colors(100), cuts=10, cex.main=0.7, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=0.6))

###### Spline con tension ##########
pred.rbf.st <- rbf(Tmt~1, temp , eta=0.003804628,rho=0.5246531, newdata=puntosxy, n.neigh=9, func="ST")
coordinates(pred.rbf.st) = c("x", "y")
gridded(pred.rbf.st) <- TRUE

# muestra el mapa de prediccion con la rbf Spline con tension
IDSTA.OV@data$pred.rbf.st <- pred.rbf.st$var1.pred
spplot(IDSTA.OV, "pred.rbf.st", main="Temperatura media 25 de enero de 2008 Predicciones \nSpline con tension", col.regions=bpy.colors(100), cuts=10, cex.main=0.7, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=0.6))


###### Spline completamente regularizada ##########
pred.rbf.crs <- rbf(Tmt~1, temp , eta=0.0036514479,rho=0.3888446, newdata=puntosxy, n.neigh=9, func="CRS")
coordinates(pred.rbf.crs) = c("x", "y")
gridded(pred.rbf.crs) <- TRUE

# muestra el mapa de prediccion con la rbf Spline completamente regularizada
IDSTA.OV@data$pred.rbf.crs <- pred.rbf.crs$var1.pred
spplot(IDSTA.OV, "pred.rbf.crs", main="Temperatura media 25 de enero de 2008 Predicciones \nSCRS", col.regions=bpy.colors(100), cuts=10, cex.main=0.7, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=0.6))

###### Spline capa delgada ##########
pred.rbf.tps <- rbf(Tmt~1, temp , eta=2,rho=2, newdata=puntosxy, n.neigh=9, func="TPS")
coordinates(pred.rbf.tps) = c("x", "y")
gridded(pred.rbf.tps) <- TRUE

# muestra el mapa de prediccion con la rbf Spline capa delgada
IDSTA.OV@data$pred.rbf.tps <- pred.rbf.tps$var1.pred
spplot(IDSTA.OV, "pred.rbf.tps", main="Temperatura media 25 de enero de 2008 Predicciones \nSpline capa delgada", col.regions=bpy.colors(100), cuts=10, cex.main=0.7, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=0.6))

###### Exponencial ##########
pred.rbf.expon <- rbf(Tmt~1, temp , eta=0.5,rho=0.5, newdata=puntosxy, n.neigh=9, func="EXPON")
coordinates(pred.rbf.expon) = c("x", "y")
gridded(pred.rbf.expon) <- TRUE

# muestra el mapa de prediccion con la rbf exponencial
IDSTA.OV@data$pred.rbf.expon <- pred.rbf.expon$var1.pred
spplot(IDSTA.OV, "pred.rbf.expon", main="Temperatura media 25 de enero de 2008 Predicciones \nExponencial", col.regions=bpy.colors(100), cuts=10, cex.main=0.7, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=0.6))

###### Gaussiana ##########
pred.rbf.gau <- rbf(Tmt~1, temp , eta=0.5,rho=0.5, newdata=puntosxy, n.neigh=9, func="GAU")
coordinates(pred.rbf.gau) = c("x", "y")
gridded(pred.rbf.gau) <- TRUE

# muestra el mapa de prediccion con la rbf gausiana
IDSTA.OV@data$pred.rbf.gau <- pred.rbf.gau$var1.pred
spplot(IDSTA.OV, "pred.rbf.gau", main="Temperatura media 25 de enero de 2008 Predicciones \nGaussiana", col.regions=bpy.colors(100), cuts=10, cex.main=0.7, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=0.6))

#SALIDA GRAFICA PARTE 1

l2 = list("sp.points", So, pch = 3, col = "grey")
rw.colors <- colorRampPalette(c("red", "yellow"))

p1 <-spplot(IDSTA.OV, "pred.rbf.m", main="Temperatura media 25 de enero de 2008 Predicciones \nMulticuadratica", col.regions=bpy.colors(100), cuts=10, cex.main=0.7, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=0.6))
p2 <- spplot(IDSTA.OV, "pred.rbf.mi", main="Temperatura media 25 de enero de 2008 Predicciones \nMulticuadratica Inversa", col.regions=bpy.colors(100), cuts=10, cex.main=0.7, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=0.6))
p3 <- spplot(IDSTA.OV, "pred.rbf.st", main="Temperatura media 25 de enero de 2008 Predicciones \nSpline con tension", col.regions=bpy.colors(100), cuts=10, cex.main=0.7, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=0.6))
p4 <- 
  spplot(IDSTA.OV, "pred.rbf.crs", main="Temperatura media 25 de enero de 2008 Predicciones \nCRS", col.regions=bpy.colors(100), cuts=10, cex.main=0.7, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=0.6))
print(p1, split = c(1, 1, 2, 2), more = T)
print(p2, split = c(2, 1, 2, 2), more = T)
print(p3, split = c(1, 2, 2, 2), more = T)
print(p4, split = c(2, 2, 2, 2), more = T)

#SALIDA GRAFICA PARTE 2

l2 = list("sp.points", So, pch = 3, col = "grey")
rw.colors <- colorRampPalette(c("red", "yellow"))

p1 <-spplot(IDSTA.OV, "pred.rbf.tps", main="Temperatura media 25 de enero de 2008 Predicciones \nSpline capa delgada", col.regions=bpy.colors(100), cuts=10, cex.main=0.7, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=0.6))
p2 <-
  spplot(IDSTA.OV, "pred.rbf.expon", main="Temperatura media 25 de enero de 2008 Predicciones \nExponencial", col.regions=bpy.colors(100), cuts=10, cex.main=0.7, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=0.6))
p3 <-spplot(IDSTA.OV, "pred.rbf.gau", main="Temperatura media 25 de enero de 2008 Predicciones \nGaussiana", col.regions=bpy.colors(100), cuts=10, cex.main=0.7, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)", key.space=list(space="right", cex=0.6))
print(p1, split = c(1, 1, 2, 2), more = T)
print(p2, split = c(2, 1, 2, 2), more = T)
print(p3, split = c(1.5, 2, 1.5, 2), more = T)

sqrt(sum(na.omit(Tmt-IDSTA.OV$pred.rbf.m))^2)/length(na.omit(Tmt-IDSTA.OV$pred.rbf.m))

#Multicuadratica
RMSPE(y_pred=na.omit(IDSTA.OV$pred.rbf.m) ,y_true=Tmt)
#Multicuadratica inversa 
RMSPE(y_pred=na.omit(IDSTA.OV$pred.rbf.mi) ,y_true=na.omit(Tmt))
#Spline con tension 
RMSPE(y_pred=na.omit(IDSTA.OV$pred.rbf.st) ,y_true=Tmt)
#CRS
RMSPE(y_pred=na.omit(IDSTA.OV$pred.rbf.crs) ,y_true=Tmt)
#TPS 
RMSPE(y_pred=na.omit(IDSTA.OV$pred.rbf.tps) ,y_true=Tmt)
#Exponencial
RMSPE(y_pred=na.omit(IDSTA.OV$pred.rbf.expon) ,y_true=Tmt)
#Gaussiana
RMSPE(y_pred=na.omit(IDSTA.OV$pred.rbf.gau) ,y_true=Tmt)


mean(IDSTA.OV$pred.ks.temp)
min(IDSTA.OV$pred.ks.temp)
max(IDSTA.OV$pred.ks.temp)

mean(IDSTA.OV$pred.ko.temp)
min(IDSTA.OV$pred.ko.temp)
max(IDSTA.OV$pred.ko.temp)

mean(IDSTA.OV$pred.rbf.crs)
min(IDSTA.OV$pred.rbf.crs)
max(IDSTA.OV$pred.rbf.crs)


min(IDSTA.OV$var.ks.temp.sqrt)
max(IDSTA.OV$var.ks.temp.sqrt)


min(IDSTA.OV$var.ko.temp.sqrt)
max(IDSTA.OV$var.ko.temp.sqrt)


###########validaci?n cruzada

KO.exp.cv.z <- krige.cv(Tmt~1, ~ E+N, model=vgm1, data=as.data.frame(temp2.sp), nmin=2, nmax=30)  
#esto es para validaci?n cruzada
KS.exp.cv.z <- krige.cv(Tmt~1, ~ E+N, model=vgm1, data=as.data.frame(temp2.sp), beta=mu, nmin=2, nmax=30)

# Correr primero la funci?n criterio.cv

resultados.cv.z <- rbind(criterio.cv(KO.exp.cv.z), criterio.cv(KS.exp.cv.z)) #criterio lo program? el profe 
#predicci?n, varianza, observado, errores, algo, el registro, coordenadas x y y
rownames(resultados.cv.z) <- c("KO.esf.cv.z", "KS.esf.cv.z")
resultados.cv.z

RMSPE <-  sqrt(sum(IDSTA.OV$pred.rbf.expon^2)/nrow(na.omit(aa)))
RMSPE

RMSPE <-  
  sqrt(sum(na.exclude(Tmt-IDSTA.OV$pred.rbf.expon))^2)/length(na.exclude(Tmt-IDSTA.OV$pred.rbf.expon))
RMSPE

sqrt((sum(temp6$Tmt-temp6$V4)^2)/length(temp6$Tmt-temp6$V4))


RMSPE <-  sqrt(sum(na.exclude(Tmt-IDSTA.OV$pred.rbf.mi))^2)/length(Tmt)
RMSPE



str(na.omit(Tmt-IDSTA.OV$pred.rbf.expon))

aa<-IDSTA.OV$LST2008_01_25
na.exclude(aa)
expon<-IDSTA.OV$pred.rbf.expon
expon[is.na(expon)] <- 0
r<-(expon-aa)
r2<-r/aa
mean(sqrt(r2))


str(na.omit(IDSTA.OV$LST2008_01_25))
str(na.exclude(IDSTA.OV$pred.rbf.expon))
aa[is.na(aa)] <- 0

RMSPE(y_pred=na.exclude(IDSTA.OV$pred.rbf.expon) ,y_true=na.exclude(IDSTA.OV$pred.ks.temp))

RMSPE(y_pred=na.exclude((IDSTA.OV$pred.rbf.crs)) ,y_true=na.exclude(Tmt))


RMSPE(y_pred=temp6$V4 ,y_true=temp6$Tmt)
RMSPE(y_pred=na.omit(IDSTA.OV$pred.rbf.tps) ,y_true=temp6$Tmt)



data5<-cbind(IDSTA.OV$pred.rbf.mi, Tmt,IDSTA.OV$LST2008_01_25,pred.rbf.mi$var1.pred )
temp5<-as.data.frame(data5)
temp6<-na.omit(temp5)

temp6$Tmt
-temp6$V4

