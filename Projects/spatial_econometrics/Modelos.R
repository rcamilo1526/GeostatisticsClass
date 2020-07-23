# Limpia el espacio de trabajo, las parcelas y la consola
rm(list = ls())
if(!is.null(dev.list())) dev.off()
cat("\014") 

#Cargar librerias
# library(leaps)
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
library(prettymapr)
library(lattice)
library(car)
library(tidyverse)
library(MASS)
#Latex
library(xtable)
library(texreg)
library(stargazer)
#cargar funciones
load("C:/Users/rmartin/Documents/11 semestre/Geoestadística/FuncDatosArea/ImagenDatosDeArea")
set.seed(123)
setwd("C:/Users/rmartin/Documents/11 semestre/Geoestadística/ProyectoDatosArea")

#puntos
balPt<-readOGR('./baltimore','baltim_clear')

library(dismo)
balPol <- voronoi(balPt)
detach(package:dismo)
bal.poly <- balPol
plot(balPol)
coords <- coordinates(balPol)

#prepare data
PRICE <- as.data.frame(balPol)$PRICE
x <- c(PRICE)
# N<-filter(balPol@data, PRICE<=84)
hist(PRICE,breaks=20)
x1<-(x-mean(x))/sd(x) 
boxplot(PRICE)
# writeOGR(bal.poly, dsn = "shapefile", layer = "map",
# driver = "ESRI Shapefile" )
#corregir BMENT
balPol$BMENT<-ifelse(balPol$BMENT==0,0,1)
bal.poly$BMENT<-balPol$BMENT
#log doesnt works
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
summary(PRICE)
#weights selected delunay
bal.nb <- poly2nb(balPol)
bal.nb <- tri2nb(coords)
bal.lw <-nb2listw(bal.nb, style="W")
bal.lags <- nblag(bal.nb, 9)
# set.seed()
# N<-filter(balPol@data, PRICE<=85)
# coords <- coordinates(N)
# balPol<-N
# vars<-c('NROOM','DWELL','NBATH','PATIO','FIREPL','AC','BMENT','NSTOR','GAR','AGE','CITCOU','LOTSZ','SQFT')
attach(balPol@data)
lm_fit <- lm(PRICE~ DWELL+NBATH+PATIO+FIREPL+BMENT+AC+GAR+CITCOU+LOTSZ, data =balPol)
bptest(lm_fit)
resettest(lm_fit)
raintest(lm_fit)
shapiro.test(residuals(lm_fit))
library(car)
vif(lm_fit)
length(N$PRICE)
hetero.plot <- function(model) {
  plot(residuals(model) ~ model$yhat)
  abline(h=0, lty="dotted")
  lines(lowess(model$yhat, residuals(model)), col="red")
}
hetero.plot(model.h)
lm_fit$fitted.values
#####ANAFORMOSIS GAUSSIANA
library(RGeostats)
library(gstat)
library(sp)

attach(balPol@data)
cx<-balPol$X
cy<-balPol$Y

#01. Fit Gaussian Anamorphosis variable dependiente
mdb = balPol[,c("X","Y","PRICE")] #must be 3 columns: x,y,z
mdb.rgdb = db.create(mdb,ndim=2,autoname=F)
mdb.herm = anam.fit(mdb.rgdb,name="PRICE",type="gaus")
mdb.hermtrans = anam.z2y(mdb.rgdb,names="PRICE",anam=mdb.herm)
PRICE.trans = mdb.hermtrans@items$Gaussian.PRICE
hist(PRICE.trans)
shapiro.test(PRICE.trans)

lm_fitr <- lm(PRICE.trans ~ DWELL+NBATH.trans+PATIO+FIREPL+AC+BMENT+GAR.trans+CITCOU+LOTSZ.trans, data=balPol)
summary(lm_fitr)
bptest(lm_fit)
resettest(lm_fit)
raintest(lm_fit)
shapiro.test(residuals(lm_fit))
vif(lm_fit)
hist(PRICE.trans,breaks=20)

plot(balPol$PRICE,balPol$res)
balPol$res=residuals(lm_fit)
#01. Fit Gaussian Anamorphosis variables independientes
mdb2 = balPol[,c("X","Y","NBATH")] #must be 3 columns: x,y,z
mdb.rgdb2 = db.create(mdb2,ndim=2,autoname=F)
mdb.herm2 = anam.fit(mdb.rgdb2,name="NBATH",type="gaus")
mdb.hermtrans2 = anam.z2y(mdb.rgdb2,names="NBATH",anam=mdb.herm2)
NBATH.trans = mdb.hermtrans2@items$Gaussian.NBATH
hist(NBATH.trans)
shapiro.test(log(NBATH.trans))
####################01. Fit Gaussian Anamorphosis variables independientes
mdb2 = balPol[,c("X","Y","GAR")] #must be 3 columns: x,y,z
mdb.rgdb2 = db.create(mdb2,ndim=2,autoname=F)
mdb.herm2 = anam.fit(mdb.rgdb2,name="GAR",type="gaus")
mdb.hermtrans2 = anam.z2y(mdb.rgdb2,names="GAR",anam=mdb.herm2)
GAR.trans = mdb.hermtrans2@items$Gaussian.GAR
hist(GAR.trans)
shapiro.test(GAR.trans)
####################01. Fit Gaussian Anamorphosis variables independientes
mdb2 = balPol[,c("X","Y","LOTSZ")] #must be 3 columns: x,y,z
mdb.rgdb2 = db.create(mdb2,ndim=2,autoname=F)
mdb.herm2 = anam.fit(mdb.rgdb2,name="LOTSZ",type="gaus")
mdb.hermtrans2 = anam.z2y(mdb.rgdb2,names="LOTSZ",anam=mdb.herm2)
LOTSZ.trans = mdb.hermtrans2@items$Gaussian.LOTSZ
hist(LOTSZ.trans)
shapiro.test(LOTSZ.trans)
# DWELL+PATIO+FIREPL+
####################01. Fit Gaussian Anamorphosis variables independientes
mdb2 = balPol[,c("X","Y","DWELL")] #must be 3 columns: x,y,z
mdb.rgdb2 = db.create(mdb2,ndim=2,autoname=F)
mdb.herm2 = anam.fit(mdb.rgdb2,name="DWELL",type="gaus")
mdb.hermtrans2 = anam.z2y(mdb.rgdb2,names="DWELL",anam=mdb.herm2)
DWELL.trans = mdb.hermtrans2@items$Gaussian.DWELL
hist(DWELL.trans)
shapiro.test(DWELL.trans)
####################01. Fit Gaussian Anamorphosis variables independientes
mdb2 = balPol[,c("X","Y","PATIO")] #must be 3 columns: x,y,z
mdb.rgdb2 = db.create(mdb2,ndim=2,autoname=F)
mdb.herm2 = anam.fit(mdb.rgdb2,name="PATIO",type="gaus")
mdb.hermtrans2 = anam.z2y(mdb.rgdb2,names="PATIO",anam=mdb.herm2)
PATIO.trans = mdb.hermtrans2@items$Gaussian.PATIO
hist(PATIO.trans)
shapiro.test(PATIO.trans)
####################01. Fit Gaussian Anamorphosis variables independientes
mdb2 = balPol[,c("X","Y","FIREPL")] #must be 3 columns: x,y,z
mdb.rgdb2 = db.create(mdb2,ndim=2,autoname=F)
mdb.herm2 = anam.fit(mdb.rgdb2,name="FIREPL",type="gaus")
mdb.hermtrans2 = anam.z2y(mdb.rgdb2,names="FIREPL",anam=mdb.herm2)
FIREPL.trans = mdb.hermtrans2@items$Gaussian.FIREPL
hist(FIREPL.trans)
shapiro.test(FIREPL.trans)
#AC+BMENT+CITCOU
####################01. Fit Gaussian Anamorphosis variables independientes
mdb2 = balPol[,c("X","Y","AC")] #must be 3 columns: x,y,z
mdb.rgdb2 = db.create(mdb2,ndim=2,autoname=F)
mdb.herm2 = anam.fit(mdb.rgdb2,name="AC",type="gaus")
mdb.hermtrans2 = anam.z2y(mdb.rgdb2,names="AC",anam=mdb.herm2)
AC.trans = mdb.hermtrans2@items$Gaussian.AC
hist(AC.trans)
shapiro.test(AC.trans)
####################01. Fit Gaussian Anamorphosis variables independientes
mdb2 = balPol[,c("X","Y","BMENT")] #must be 3 columns: x,y,z
mdb.rgdb2 = db.create(mdb2,ndim=2,autoname=F)
mdb.herm2 = anam.fit(mdb.rgdb2,name="BMENT",type="gaus")
mdb.hermtrans2 = anam.z2y(mdb.rgdb2,names="BMENT",anam=mdb.herm2)
BMENT.trans = mdb.hermtrans2@items$Gaussian.BMENT
hist(BMENT.trans)
shapiro.test(BMENT.trans)
####################01. Fit Gaussian Anamorphosis variables independientes
mdb2 = balPol[,c("X","Y","CITCOU")] #must be 3 columns: x,y,z
mdb.rgdb2 = db.create(mdb2,ndim=2,autoname=F)
mdb.herm2 = anam.fit(mdb.rgdb2,name="CITCOU",type="gaus")
mdb.hermtrans2 = anam.z2y(mdb.rgdb2,names="CITCOU",anam=mdb.herm2)
CITCOU.trans = mdb.hermtrans2@items$Gaussian.CITCOU
hist(CITCOU.trans)
shapiro.test(CITCOU.trans)


lm_fit_ag <- lm(PRICE.trans ~ DWELL+NBATH+PATIO+FIREPL+AC+BMENT+GAR+CITCOU+LOTSZ,data=balPol@data)
# summary(lm_fit_ag)
shapiro.test(residuals(lm_fit_ag))# no hay normalidad

attach(balPol@data)

####################################################################################
# SPATIAL MODELS
####################################################################################
# Modelo Spatial Lag
col.lag.sm <- lagsarlm(PRICE~ DWELL+NBATH+PATIO+FIREPL+AC+BMENT+GAR+CITCOU+LOTSZ, data = balPol, listw=bal.lw)
predict(col.lag.sm)
hetero.plot(col.lag.sm)
class(col.lag.sm)
summary(col.lag.sm, Nagelkerke=T)
shapiro.test(residuals(col.lag.sm))
bptest.sarlm(col.lag.sm)
library(car)
vif(col.lag.sm)
# Modelo Spatial Error:PRICE.trans ~ DWELL.trans+PATIO.trans+FIREPL.trans+AC.trans+BMENT.trans+GAR.trans+CITCOU.trans+LOTSZ.trans
col.error.sm <- errorsarlm(PRICE ~ DWELL+NBATH+PATIO+FIREPL+AC+BMENT+GAR+CITCOU+LOTSZ, data=balPol@data, listw=bal.lw)
summary(col.error.sm, Nagelkerke=T)
shapiro.test(residuals(col.error.sm))
bptest.sarlm(col.error.sm)

#SLX
mod.SLX <- lmSLX(PRICE ~ DWELL+NBATH+PATIO+FIREPL+AC+BMENT+GAR+CITCOU+LOTSZ, data=balPol@data, listw=bal.lw)
summary(mod.SLX)
shapiro.test(residuals(mod.SLX))
bptest(mod.SLX)
# SARAR
col.sarar <- sacsarlm(PRICE ~ DWELL+NBATH+PATIO+FIREPL+AC+BMENT+GAR+CITCOU+LOTSZ, data=balPol@data, listw=bal.lw)    # CRIME ~1
summary(col.sarar, Nagelkerke=T)
shapiro.test(residuals(col.sarar))
bptest.sarlm(col.sarar)
attributes(columlagsd)
columlagsd[['rho']]
#durbin
columlagsd <- lagsarlm(PRICE ~ DWELL+NBATH+PATIO+FIREPL+AC+BMENT+GAR+CITCOU+LOTSZ, data=balPol@data, listw=bal.lw, type="mixed") 
summary(columlagsd, Nagelkerke=T)
shapiro.test(columlagsd$residuals)
bptest.sarlm(columlagsd)
#GLM
# Se rezagan las variables, para ello se copia la informaci?n en col.poly1
# PRICE ~ DWELL+NBATH+PATIO+FIREPL+AC+BMENT+GAR+CITCOU+LOTSZ
col.poly1 <- as.data.frame(balPol)
col.poly1$WX <- lag.listw(bal.lw, col.poly1$X)
col.poly1$WY <- lag.listw(bal.lw, col.poly1$Y)
col.poly1$WDWELL <- lag.listw(bal.lw, col.poly1$DWELL)
col.poly1$WNBATH <- lag.listw(bal.lw, col.poly1$NBATH)
col.poly1$WPATIO <- lag.listw(bal.lw, col.poly1$PATIO)
col.poly1$WFIREPL<- lag.listw(bal.lw, col.poly1$FIREPL)
col.poly1$WAC<- lag.listw(bal.lw, col.poly1$AC)
col.poly1$WBMENT <- lag.listw(bal.lw, col.poly1$BMENT)
col.poly1$WGAR <- lag.listw(bal.lw, col.poly1$GAR)
col.poly1$WCITCOU <- lag.listw(bal.lw, col.poly1$CITCOU)
col.poly1$WLOTSZ <- lag.listw(bal.lw, col.poly1$LOTSZ)

# X <- model.matrix(CRIME~INC + HOVAL, data=as.data.frame(col.poly))
# WX <- create_WX(X,a.lw)
a.lw1 <- nb2listw(bal.lags[[1]], style="W")
a.lw2 <- nb2listw(bal.lags[[2]], style="W")
a.lw3 <- nb2listw(bal.lags[[3]], style="W")
a.lw4 <- nb2listw(bal.lags[[4]], style="W")
a.lw5 <- nb2listw(bal.lags[[5]], style="W")
a.lw6 <- nb2listw(bal.lags[[6]], style="W",zero.policy =T)
a.lw7 <- nb2listw(bal.lags[[7]], style="W",zero.policy =T)
col.poly1$W1PRICE <- lag.listw(a.lw1,col.poly1$PRICE)
col.poly1$W2PRICE <- lag.listw(a.lw2,col.poly1$PRICE)
col.poly1$W3PRICE <- lag.listw(a.lw3,col.poly1$PRICE)
col.poly1$W4PRICE <- lag.listw(a.lw4,col.poly1$PRICE)
col.poly1$W5PRICE <- lag.listw(a.lw5,col.poly1$PRICE)
col.poly1$W6PRICE <- lag.listw(a.lw6,col.poly1$PRICE, zero.policy =T)

t.glm<-glm(PRICE ~ W1PRICE+WDWELL+WNBATH+WPATIO+WFIREPL+WAC+WBMENT+WGAR+WCITCOU+WLOTSZ+X+Y+WX+WY,data=col.poly1, family=poisson(link="log"), na.action=na.exclude)
summary(t.glm, Nagelkerke=T)
residuals.glm <- balPol$PRICE-fitted.values(t.glm)
shapiro.test(residuals.glm)
bptest(t.glm)

pseudoR2.glm <- cor(exp(predict(t.glm)),balPol$PRICE)^2
AIC(t.glm)
hetero.plot <- function(model) {
  plot(residuals(model) ~ fitted(model))
  abline(h=0, lty="dotted")
  lines(lowess(fitted(model), residuals(model)), col="red")
}
NK <- function(obj, y) { 
  n <- length(obj$residuals) 
  nullLL <- logLik(lm(y ~ 1)) 
  c(1 - exp(-(2/n)*(logLik(obj) - nullLL))) 
} 
#https://www.rdocumentation.org/packages/sphet/versions/1.7/topics/gstslshet

#Linear espacial heterocedastico
model.h.sarar<- gstslshet(formula = PRICE~ DWELL+NBATH+PATIO+FIREPL+BMENT+AC+GAR+CITCOU+LOTSZ,  
                     listw = bal.lw, data = balPol ,sarar=T)
summary(model.h.sarar, Nagelkerke=T)
impacts(model.h.sarar)
summary(model.hac, Nagelkerke=T)
#bptest.sarlm(model.hac)
hetero.plot <- function(model) {
  plot(residuals(model) ~ model$yhat)
  abline(h=0, lty="dotted")
  lines(lowess(model$yhat, residuals(model)), col="red")
}
attributes(model.h)
hetero.plot(model.h)
Pdist <- read.gwt2dist("D:/Documentos/11 semestre/Geoestadística/ProyectoDatosArea/baltimore/baltim_clear.gwt", region.id = NULL, skip = 1)
attributes(Pdist)
##model HAC
res <- stslshac(formula = PRICE~ DWELL+NBATH+PATIO+FIREPL+BMENT+AC+GAR+CITCOU+LOTSZ,  listw = bal.lw, data = balPol, 
                distance = Pdist, type = 'Triangular')
summary(res)
hetero.plot(res)
bptest(res)


hetero.plot(col.lag.sm)
shapiro.test(res$residuals)

################################################3
#########QUITAR ATIPICOS#######################
############################################

# mod <- lm(PRICE ~ 1,data=balPol@data)
# cooksd <- cooks.distance(mod)
# plot(cooksd, pch="*", cex=2)
# abline(h = 4*mean(cooksd, na.rm=T), col="red")
# text(x=1:length(cooksd)+1, y=cooksd,
#      labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),
#                    names(cooksd),""),
#      col="red")

############################
##  An?lisis de Impactos  ##
############################mos

col.new <- balPol
col.new@data <- col.poly1

# Cambiando price of station from 3.5 to 68
col.new@data[col.new@data$STATION == "16","PRICE"] <- 68.0

# Los valores de las predicciones originales
orig.pred <- res$yhat

# Los valores predichos con la nueva tasa de criminalidad en el barrio 10
col_nbq1
new.pred <- stslshac(formula = PRICE~ DWELL+NBATH+PATIO+FIREPL+BMENT+AC+GAR+CITCOU+LOTSZ,  listw = bal.lw, 
                            data = col.new, 
                            distance = Pdist, type = 'Triangular')$yhat
# Las diferencias entre las predicciones
effect.10 <- new.pred - orig.pred                  
el <- data.frame(name = col.new$STATION, dif_pred_PRICE = effect.10)
el$dif_pred_PRICE
col.new$ef10 <- effect.10
# Ordenando los barrios por el valor absoluto del cambio en la predicci?n del CRIME
el <- el[rev(order(abs(el$dif_pred_PRICE))), ]
xtable(el[1:10, ])  #muestra los 10 primeros barrios

# Mapear estos cambios es tambien importante:

breaks <- c(min(col.new$ef10), -0.05, 0.05, max(col.new$ef10))
labels <- c("Efecto negativo (< -.05)", "Sin efecto (-.05 a .05)", 
            "Efecto positivo (> .05)")
# faltaba all.inside =T para evitar un poligono en blanco asociado al valor m?ximo  
np <- findInterval(col.new$ef10, breaks,all.inside =T) 
colors <- c("red", "orange", "green")
x11()
# Dibujando el mapa
plot(col.new, col = colors[np])
mtext("Efectos de un cambio en el punto 16, (fijando PRICE=68)\n sobre los valores predichos en un modelo Durbin Espacial", 
      side = 3, line = 1)
legend("topleft", legend = labels, fill = colors, bty = "n")
# Tambi?n podr?amos mapear la magnitud de los cambios causados a ra?z de la disminuci?n de la criminalidad en el barrio 10.
x11()
pal5 <- brewer.pal(6, "Spectral")
cats5 <- classIntervals(col.new$ef10, n = 5, style = "jenks")
colors5 <- findColours(cats5, pal5)
plot(col.new, col = colors5)
legend(locator(1), legend = round(cats5$brks, 2), fill = pal5, bty = "n")
mtext("Efectos de un cambio en el punto 16, (fijando PRICE=68)\n sobre los valores de las 
      predicciones en un modelo Durbin Espacial", side = 3, line = 1)

# While these maps and the effects for individual counties may be useful they only tell us about a 
# change in one place. The impacts() function gives us something like OLS regression coefficients for
# a spatial lag model. The logic of the impacts() function is similar to the code above, it tells you
# the direct (local), indirect (spill-over), and total effect of a unit change in each of the 
# predictor variables. The changes reported by impacts are the global average impact:


W <- as(as_dgRMatrix_listw(bal.lw), "CsparseMatrix")
trMatc <- trW(W, type="mult")
trMC <- trW(W, type="MC")
res$coefficients["Wy"]
impacts(res, listw=bal.lw)
attributes(res)
# Impacto Direct Wy:
SrW.I <- solve(diag(209)-as.numeric(res$coefficients["Wy"])*as.matrix(W))%*%diag(209)*as.numeric(res$coefficients["Wy"])
sum(diag(SrW.I))/209

# Impacto Indirect INC:
(sum(SrW.I) - sum(diag(SrW.I)))/209



# Impacto Total INC:
rep(1,209)%*%(SrW.I)%*%matrix(rep(1,209),ncol=1)/209              # otra forma mas simple es: sum(SrW.I)/49




summary(impacts(res, tr=trMC,R=200), zstats=TRUE, short=TRUE)
summary(impacts(res, tr=trMatc,R=200), zstats=TRUE, short=TRUE)

########################################################################
##########       Regresi?n Geogr?ficamente Ponderada       #############
########################################################################

# https://gis.stackexchange.com/questions/241127/how-to-plot-output-from-gwr-in-r
# http://geokitchen.blogspot.com.co/2012/09/r-geographically-weighted-regression.html

adapt <- gwr.sel(PRICE~ DWELL+NBATH+PATIO+FIREPL+BMENT+AC+GAR+CITCOU+LOTSZ+ WX+WY, data=col.poly1, coords=cbind(col.poly1$X,col.poly1$Y)) 
gwr_fit <- gwr(PRICE~ DWELL+NBATH+PATIO+FIREPL+BMENT+AC+GAR+CITCOU+LOTSZ+WX+WY, data=col.poly1, coords=cbind(col.poly1$X,col.poly1$Y), adapt = 1, hatmatrix = TRUE)
gwr_fit
gwr1 <- gwr.sel(PRICE~ DWELL+NBATH+PATIO+FIREPL+BMENT+AC+GAR+CITCOU+LOTSZ+WX+WY, data=col.poly1, coords=cbind(col.poly1$X,col.poly1$Y), gweight=gwr.bisquare)
# Creando un objecto con el valor de Quasi-global R2
globalR2 <- (1 - (gwr_fit$results$rss/gwr_fit$gTSS))
summary(gwr_fit)
results<-as.data.frame(gwr_fit$SDF)
head(results)

head(gwr_fit$SDF)
# Create a colour palette 
lm.palette <- colorRampPalette(c("red","orange", "white"), space = "rgb")
colours = c("dark blue", "blue", "red", "dark red")
bal.poly$Intercept_coef <- gwr_fit$SDF@data[,2]
bal.poly$DWELL_coef <- gwr_fit$SDF@data$DWELL
bal.poly$NBATH_coef <- gwr_fit$SDF@data$NBATH
bal.poly$PATIO_coef <- gwr_fit$SDF@data$PATIO
bal.poly$FIREPL_coef <- gwr_fit$SDF@data$FIREPL
bal.poly$AC_coef <- gwr_fit$SDF@data$AC
bal.poly$BMENT_coef <- gwr_fit$SDF@data$BMENT
bal.poly$GAR_coef <- gwr_fit$SDF@data$GAR
bal.poly$CITCOU_coef <- gwr_fit$SDF@data$CITCOU
bal.poly$LOTSZ_coef <- gwr_fit$SDF@data$LOTSZ
bal.poly$localR2 <- gwr_fit$SDF@data$localR2
bal.poly$pred <- gwr_fit$SDF@data$pred
Conf3x1 = matrix(c(1:3), nrow=1, byrow=TRUE)
png(paste("D:/Documentos/11 semestre/Geoestadística/ProyectoDatosArea/GWR/hist","N0",".png"),width = 3000, height = 900,units="px",res=300)
LR<-spplot(bal.poly, "localR2", key.space = "right", col.regions = lm.palette(20), cuts = 7, 
       main = "PRICE_localR2")
INT<-spplot(bal.poly, "Intercept_coef", key.space = "right", col.regions = lm.palette(20), cuts = 7, 
       main = "Intercept_coef")
print(LR, split = c(1,1,2,1), more = TRUE)
print(INT, split = c(2,1,2,1), more = FALSE)
dev.off()
png(paste("D:/Documentos/11 semestre/Geoestadística/ProyectoDatosArea/GWR/hist","N1",".png"),width = 4500, height = 900,units="px",res=300)
dwell<-spplot(bal.poly, "DWELL_coef", col.regions=lm.palette(20), cex=0.6, cuts=7,main = "DWELL Coef") 
nbath<-spplot(bal.poly, "NBATH_coef", col.regions=lm.palette(20), cex=0.6, cuts=7,main = "NBATH Coef") 
patio <- spplot(bal.poly, "PATIO_coef", col.regions=lm.palette(20), cex=0.6, cuts=7,main = "PATIO Coef") 
print(dwell, split = c(1,1,3,1), more = TRUE)
print(nbath, split = c(2,1,3,1), more = TRUE)
print(patio, split = c(3,1,3,1), more = FALSE)
dev.off()
png(paste("D:/Documentos/11 semestre/Geoestadística/ProyectoDatosArea/GWR/hist","N2",".png"),width = 4500, height = 900,units="px",res=300)
f<-spplot(bal.poly, "FIREPL_coef", col.regions=lm.palette(20), cex=0.6, cuts=7,main = "FIREPL Coef")
ac<-spplot(bal.poly, "AC_coef", col.regions=lm.palette(20), cex=0.6, cuts=7,main = "AC Coef") 
bm <- spplot(bal.poly, "BMENT_coef", col.regions=lm.palette(20), cex=0.6, cuts=7,main = "BMENT Coef")  
print(f, split = c(1,1,3,1), more = TRUE)
print(ac, split = c(2,1,3,1), more = TRUE)
print(bm, split = c(3,1,3,1), more = FALSE)
dev.off()
png(paste("D:/Documentos/11 semestre/Geoestadística/ProyectoDatosArea/GWR/hist","N3",".png"),width = 4500, height = 900,units="px",res=300)
gar<-spplot(bal.poly, "GAR_coef", col.regions=lm.palette(20), cex=0.6, cuts=7,main = "GAR Coef") 
cit<-spplot(bal.poly, "CITCOU_coef", col.regions=lm.palette(20), cex=0.6, cuts=7,main = "CITCOU Coef") 
lot<-spplot(bal.poly, "LOTSZ_coef", col.regions=lm.palette(20), cex=0.6, cuts=7,main = "LOTSZ Coef") 
print(gar, split = c(1,1,3,1), more = TRUE)
print(cit, split = c(2,1,3,1), more = TRUE)
print(lot, split = c(3,1,3,1), more = FALSE)
dev.off()
names(gwr_fit$SDF)


# Makes the data into a "mappable" format 
map = bal.poly

# Lists the "mappable" data 
names(map) 

#calculate t-value
t = gwr_fit$SDF$DWELL / gwr_fit$SDF$DWELL_se  
map@data$t = t 
colours2=c("green","red","green") 

#estimated GWR t- values, red indicates a relationship that is not significant
tdwell<-spplot(map, "t", col.regions=lm.palette(20), cex=c(0.6,1,0.6), main = "t - Value DWELL") 
t = gwr_fit$SDF$NBATH / gwr_fit$SDF$NBATH_se
map@data$t = t 
tNBA<-spplot(map, "t", col.regions=lm.palette(20), cex=c(0.6,1,0.6), main = "t - Value NBATH")
t = gwr_fit$SDF$PATIO / gwr_fit$SDF$PATIO_se
map@data$t = t 
tPATIO<-spplot(map, "t", col.regions=lm.palette(20), cex=c(0.6,1,0.6), main = "t - Value PATIO")
png(paste("D:/Documentos/11 semestre/Geoestadística/ProyectoDatosArea/GWR/hist","t1",".png"),width = 4500, height = 900,units="px",res=300)
print(tdwell, split = c(1,1,3,1), more = TRUE)
print(tNBA, split = c(2,1,3,1), more = TRUE)
print(tPATIO, split = c(3,1,3,1), more = FALSE)
dev.off()

t = gwr_fit$SDF$FIREPL / gwr_fit$SDF$FIREPL_se
map@data$t = t 
tFIR<-spplot(map, "t", col.regions=lm.palette(20), cex=c(0.6,1,0.6), main = "t - Value FIREPL")
t = gwr_fit$SDF$AC / gwr_fit$SDF$AC_se
map@data$t = t 
tAC<-spplot(map, "t", col.regions=lm.palette(20), cex=c(0.6,1,0.6), main = "t - Value AC")
t = gwr_fit$SDF$BMENT / gwr_fit$SDF$BMENT_se
map@data$t = t 
tBMENT<-spplot(map, "t", col.regions=lm.palette(20), cex=c(0.6,1,0.6), main = "t - Value BMENT")
png(paste("D:/Documentos/11 semestre/Geoestadística/ProyectoDatosArea/GWR/hist","t2",".png"),width = 4500, height = 900,units="px",res=300)
print(tFIR, split = c(1,1,3,1), more = TRUE)
print(tAC, split = c(2,1,3,1), more = TRUE)
print(tBMENT, split = c(3,1,3,1), more = FALSE)
dev.off()

t = gwr_fit$SDF$GAR / gwr_fit$SDF$GAR_se
map@data$t = t 
tGAR<-spplot(map, "t", col.regions=lm.palette(20), cex=c(0.6,1,0.6), main = "t - Value GAR")
t = gwr_fit$SDF$CITCOU / gwr_fit$SDF$CITCOU_se
map@data$t = t 
tCIT<-spplot(map, "t", col.regions=lm.palette(20), cex=c(0.6,1,0.6), main = "t - Value CITCOU")
t = gwr_fit$SDF$LOTSZ/ gwr_fit$SDF$LOTSZ_se
map@data$t = t 
tLOT<-spplot(map, "t", col.regions=lm.palette(20), cex=c(0.6,1,0.6), main = "t - Value BMENT")
png(paste("D:/Documentos/11 semestre/Geoestadística/ProyectoDatosArea/GWR/hist","t3",".png"),width = 4500, height = 900,units="px",res=300)
print(tGAR, split = c(1,1,3,1), more = TRUE)
print(tCIT, split = c(2,1,3,1), more = TRUE)
print(tLOT, split = c(3,1,3,1), more = FALSE)
dev.off()
## predicted and standard error
map@data$pred <- gwr_fit$SDF@data$pred
map@data$pred.se <- gwr_fit$SDF@data$pred.se
map@data$localR2 <- gwr_fit$SDF@data$localR2
png(paste("D:/Documentos/11 semestre/Geoestadística/ProyectoDatosArea/GWR/hist","N4",".png"),width = 4500, height = 900,units="px",res=300)
PR<-spplot(map, "pred",col.regions=lm.palette(20), main = "Predicted Value")
ST<-spplot(map, "pred.se", col.regions=lm.palette(20), main = "Standard Error")
BUB<-bubble(gwr_fit$SDF, "localR2",fill = F, main = "local R2 Coef")
print(PR, split = c(1,1,3,1), more = TRUE)
print(ST, split = c(2,1,3,1), more = TRUE)
print(BUB, split = c(3,1,3,1), more = FALSE)
dev.off()

#######Global tests of geographical weighted regressions ###########
##Four related test statistics for comparing OLS and GWR models based on papers by Brunsdon, 
#Fotheringham and Charlton (1999) and Leung et al (2000), and a development from the GWR book (2002).
attach(balPol)

##Brunsdon, Fotheringham & Charlton (1999) ANOVA
BFC99.gwr.test(gwr_fit)
#Brunsdon, Fotheringham & Charlton (2002, pp. 91-2)
BFC02.gwr.test(gwr_fit)

# anova
anova(gwr_fit)

col.poly2 <- col.poly1
coordinates(col.poly2) <- c("X", "Y")
xx <- gwr(CRIME~ INC + HOVAL + X + Y + WX + WY, col.poly2, bandwidth = 9.454624, hatmatrix=TRUE)
xx

# local R2
map@data$localR2b <- xx$SDF@data$localR2
spplot(map, "localR2b", key.space="right", col.regions = lm.palette(20), cuts=7, main = "Local R2")

