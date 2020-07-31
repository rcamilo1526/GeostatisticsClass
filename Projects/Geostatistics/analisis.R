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
# load("/home/rcamilo/Documents/geoestadistica/croatia.grilla 25s.2008.rda")
 
IDSTAr.ov <- IDSTA.ov[!is.na(IDSTA.ov$LST2008_07_03),]

options(max.print = 30)

#guardar de la base 
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

#ANALISIS DE NORMALIDAD
    
  #graficos
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
  #pruebas
    shapiro.test(Tmt)
    ks.test(as.numeric(scale(sort(Tmt))),pnorm)

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

#normal con anamorfosis
    #graficos
      sd(Tmt.trans)
      x11()
      par(mfrow=c(nrow=1, ncol=3))
      hist(Tmt.trans, freq=F, breaks=10, xlab="Temperatura (?C)", ylab="Frecuencia", col="orange", main ="Histograma")
      curve(dnorm(x, mean=mean(Tmt.trans), sd=sd(Tmt.trans)), add=T, col="red")                  # X no funciona, solo minuscula
      qqnorm(Tmt.trans, ylab="Temperatura (?C)", xlab="Cuantiles te?ricos", col="orange")
      qqline(Tmt.trans)                                                        #Agregarle al Q-Q Plot la l?nea media
      boxplot(Tmt.trans, main="BOX-PLOT", notch=F, horizontal=T, xlab="Temperatura (?C)",col="orange")
      points(mean(Tmt.trans), y=1, pch=1, cex=2)
      par(mfrow=c(nrow=1, ncol=1))
    #pruebas
      shapiro.test(Tmt.trans)
      ks.test(as.numeric(scale(sort(Tmt.trans))),pnorm)

#ANALISIS DE TENDENCIA
      
  modelo1<-lm(Tmt.trans~Hrdem+Hrdsea+Hrtwi+N+E)
  summary(modelo1)
  step.model <- stepAIC(modelo1, direction = "both", trace = FALSE)
  summary(step.model)
  # residuals(step.model)
  #Hay tendencia, se toman los residuos
    tp.geo <- croacia.geoR
    tp.geo$data <- residuals(step.model)

    #comprobrando el grado de tendencia presente en el modelo 
    modelot2<-lm(Tmt.trans~Hrdem+Hrdsea+Hrtwi+N+E+N^2+E^2+N*E+Hrdem*N+Hrdsea*N+Hrtwi*N+Hrdem*E+Hrdsea*E+Hrtwi*E+Hrdem*E^2+Hrdsea*E^2+Hrtwi*E^2+Hrdem*N^2+Hrdsea*N^2+Hrtwi*N^2)
    summary(modelot2)
    modelot3<-lm(Tmt.trans~Hrdem+Hrdsea+Hrtwi+N+E+N^2+E^2+N^3+E^3+N*E+((N^2)*(E^2)))
    summary(modelot3)
    modelot4<-lm(Tmt.trans~Hrdem+Hrdsea+Hrtwi+N+E+N^2+E^2+N^3+E^3+N*E+((N^2)*(E^2))+((N^3)*(E^3)))
    summary(modelot4)
    modelot4<-lm(Tmt.trans~Hrdem+Hrdsea+Hrtwi+N+E+N^2+E^2+N^3+E^3+N*E+((N^2)*(E^2))+((N^3)*(E^3)))
    summary(modelot4)
    modelot5<-lm(Tmt.trans~Hrdem+Hrdsea+Hrtwi+N+E+N^2+E^2+N^3+E^3+N^4+E^4+N*E+((N^2)*(E^2))+((N^3)*(E^3))+((N^4)*(E^4))+((N^5)*(E^5))+((N^6)*(E^6))+((N^7)*(E^7)))
    summary(modelot5)
#ANALISIS DE ANISOTROPIA
    
  library(intamap)
  datares<-cbind(E, N, tp.geo$data)
  data.geo<-as.data.frame(datares)
  coordinates(data.geo) <- ~E+N
  data.geo$value=data.geo$V3
  estimateAnisotropy(data.geo) 
  #SIN ANISOTROPIA
  dmax<-max(dist(tp.geo$coords))/2
  var.an <- variog4(tp.geo, max.dist=dmax)
  plot(var.an) 

  
#SEMIVARIOGRAMAS
  library(gstat)
  library(geoR)
  library(sgeostat)
  library(geospt)
  
  
  #EXPERIMENTALES
  data2.geo<-as.data.frame(tp.geo)
  class(data2.geo) <- c("data.frame")
  data.or<-data2.geo[,c(2,1,3)]
  
  names(data.or) <- c("x", "y", "t")
  point <- point(data.or)
  pair <- pair(point,num.lags=30,maxdist=dmax)
  

  v <- est.variograms(point,pair,"t",trim=0.1) 
  
  ####Gráficos de semivariogramas experimentales
  x11()
  
  layout(Conf2x2)
  op=par(mfrow=c(2,2))
  plot(v$bins,v$classic,pch=19, col="cyan3", main="Modelo experimental clásico", ylab="Semivarianza", xlab="Distancia")
  lines(v$bins,v$classic,col=1,lty=1)
  plot(v$bins,v$robust,pch=19,col="deeppink3", main="Modelo experimental robusto", ylab="Semivarianza", xlab="Distancia")
  lines(v$bins,v$robust,col=1,lty=1)
  plot(v$bins,v$med,pch=19,col="darkorchid3", main="Modelo experimental mediana", ylab="Semivarianza", xlab="Distancia")
  lines(v$bins,v$med,col=1,lty=1)
  plot(v$bins,v$trimmed.mean,pch=19,col="chartreuse3", main="Modelo experimental media recortada", ylab="Semivarianza", xlab="Distancia")
  lines(v$bins,v$trimmed.mean,col=1,lty=1)
  par(op)
  
  #TEORICO
  data2<-cbind(N, E, Tmt.trans)
  temp2<-as.data.frame(data2)
  newgeodata<-as.geodata(temp2)

  ##################################################################
  ###################    CLASICO   #################################
  #################################################################
  
  exp.ml.c<-likfit(geodata = newgeodata , nugget = 0.01,  ini = c(0.6,25000),fix.nugget = T)
  sph.ml.c<-likfit(geodata = newgeodata, nugget = 0.01, ini = c(0.6,38000),cov.model="sph",fix.nugget = T)
  mat.ml.c<-likfit(geodata = newgeodata,nugget = 0.25, ini = c(0.6,30000),cov.model="mat",kappa=1.5,fix.nugget = T)
  cir.ml.c<-likfit(geodata = newgeodata, nugget = 0.001,ini = c(0.6,31000),cov.model="cir",fix.nugget = T)
  pow.ml.c<-likfit(geodata = newgeodata,nugget = 0.158, ini = c(0.7,25000),cov.model="powered.exponential",kappa=1.75,fix.nugget = T)
  
  #Aic para todos los variogramas anteriores
  sv.clas<-c(
  exp.ml.c$AIC,
  sph.ml.c$AIC,
  mat.ml.c$AIC,
  cir.ml.c$AIC,
  pow.ml.c$AIC
  )
  
  
  plot(v$bins, v$classic,col=1, main="Modelo experimental clásico", ylab="Semivarianza", xlab="Distancia", pch=19)
  lines(exp.ml.c,max.dist=dmax,lwd=2,col='red')
  lines(sph.ml.c,max.dist=dmax,lwd=2,col='blue')
  lines(mat.ml.c,max.dist=dmax,lwd=2,col='green')
  lines(cir.ml.c,max.dist=dmax,lwd=2,col='yellow')
  lines(pow.ml.c,max.dist=dmax,lwd=2,col='orange')
  legend(locator(1), legend=c('Exponencial','Esférico','Matern','Cir','Pow'),lwd=c(2,2,2,2),col=c('red','blue','green','yellow','orange'))
  
  
  #######################################################
  ################### ROBUSTO ##########################
  ######################################################
  
  exp.ml.r<-likfit(geodata = newgeodata , nugget = 0.01,  ini = c(0.8,26000),fix.nugget = T)
  sph.ml.r<-likfit(geodata = newgeodata, nugget = 0.01, ini = c(0.6,35000),cov.model="sph",fix.nugget = T)
  mat.ml.r<-likfit(geodata = newgeodata,nugget = 0.28, ini = c(0.6,35000),cov.model="mat",kappa=1.5,fix.nugget = T)
  cir.ml.r<-likfit(geodata = newgeodata, nugget = 0.001,ini = c(0.6,30900),cov.model="cir",fix.nugget = T)
  pow.ml.r<-likfit(geodata = newgeodata,nugget = 0.15, ini = c(0.7,25000),cov.model="powered.exponential",kappa=1.75,fix.nugget = T)
  
  #Aic para todos los variogramas anteriores
  sv.rob<-c(exp.ml.r$AIC,
  sph.ml.r$AIC,
  mat.ml.r$AIC,
  cir.ml.r$AIC,
  pow.ml.r$AIC)
  
  
  plot(v$bins, v$robust,col=1, main="Modelo experimental robusto", ylab="Semivarianza", xlab="Distancia", pch=19)
  lines(exp.ml.r,max.dist=dmax,lwd=2,col='red')
  lines(sph.ml.r,max.dist=dmax,lwd=2,col='blue')
  lines(mat.ml.r,max.dist=dmax,lwd=2,col='green')
  lines(cir.ml.r,max.dist=dmax,lwd=2,col='yellow')
  lines(pow.ml.r,max.dist=dmax,lwd=2,col='orange')
  legend(locator(1), legend=c('Exponencial','Esférico','Matern','Cir','Pow'),lwd=c(2,2,2,2),col=c('red','blue','green','yellow','orange'))
  
  ##########################################################################
  ###################### MEDIANA###################################
  ##################################################################
  exp.ml.m<-likfit(geodata = newgeodata ,  nugget = 0.01,  ini = c(0.8,24000),fix.nugget = T)
  sph.ml.m<-likfit(geodata = newgeodata, nugget = 0.01, ini = c(0.6,35000),cov.model="sph",fix.nugget = T)
  mat.ml.m<-likfit(geodata = newgeodata,nugget = 0.26, ini = c(0.6,35000),cov.model="mat",kappa=1.5,fix.nugget = T)
  cir.ml.m<-likfit(geodata = newgeodata, nugget = 0.001,ini = c(0.6,31100),cov.model="cir",fix.nugget = T)
  pow.ml.m<-likfit(geodata = newgeodata,nugget = 0.155, ini = c(0.7,24000),cov.model="powered.exponential",kappa=1.75,fix.nugget = T)
  
  #Aic para todos los variogramas anteriores
  sv.median<-c(
  exp.ml.m$AIC,
  sph.ml.m$AIC,
  mat.ml.m$AIC,
  cir.ml.m$AIC,
  pow.ml.m$AIC
  )
  
  
  plot(v$bins, v$med,col=1, main="Modelo experimental mediana", ylab="Semivarianza", xlab="Distancia", pch=19)
  lines(exp.ml.m,max.dist=dmax,lwd=2,col='red')
  lines(sph.ml.m,max.dist=dmax,lwd=2,col='blue')
  lines(mat.ml.m,max.dist=dmax,lwd=2,col='green')
  lines(cir.ml.m,max.dist=dmax,lwd=2,col='yellow')
  lines(pow.ml.m,max.dist=dmax,lwd=2,col='orange')
  legend(locator(1), legend=c('Exponencial','Esférico','Matern','Cir','Pow'),lwd=c(2,2,2,2),col=c('red','blue','green','yellow','orange'))
  
  
  ##########################################################################
  ###################### MEDIA RECORTADA ###################################
  ##################################################################
  
  exp.ml.mr<-likfit(geodata = newgeodata , nugget = 0.01,  ini = c(0.8,22000),fix.nugget = T)
  sph.ml.mr<-likfit(geodata = newgeodata, nugget = 0.01, ini = c(0.6,39000),cov.model="sph",fix.nugget = T)
  mat.ml.mr<-likfit(geodata = newgeodata,nugget = 0.27, ini = c(0.6,34000),cov.model="mat",kappa=1.5,fix.nugget = T)
  cir.ml.mr<-likfit(geodata = newgeodata, nugget = 0.001,ini = c(0.6,31400),cov.model="cir",fix.nugget = T)
  pow.ml.mr<-likfit(geodata = newgeodata,nugget = 0.153, ini = c(0.7,24000),cov.model="powered.exponential",kappa=1.75,fix.nugget = T)
  
  #Aic para todos los variogramas anteriores
  sv.medrec<-c(
  exp.ml.mr$AIC,
  sph.ml.mr$AIC,
  mat.ml.mr$AIC,
  cir.ml.mr$AIC,
  pow.ml.mr$AIC
  )
  
  plot(v$bins, v$med,col=1, main="Modelo experimental media recortada", ylab="Semivarianza", xlab="Distancia", pch=19)
  lines(exp.ml.mr,max.dist=dmax,lwd=2,col='red')
  lines(sph.ml.mr,max.dist=dmax,lwd=2,col='blue')
  lines(mat.ml.mr,max.dist=dmax,lwd=2,col='green')
  lines(cir.ml.mr,max.dist=dmax,lwd=2,col='yellow')
  lines(pow.ml.mr,max.dist=dmax,lwd=2,col='orange')
  legend(locator(1), legend=c('Exponencial','Esférico','Matern','Cir','Pow'),lwd=c(2,2,2,2),col=c('red','blue','green','yellow','orange'))
  
  
  data.frame(sv.clas,sv.rob,sv.median,sv.medrec)
  #se corre el modelo escog?do con ML WLS, RML
  # clasico nugget = 0.01,  ini = c(0.8,22000)
  variogramaclasico<-variog( newgeodata,max.dist=dmax, option = "cloud")
  exp.ml<-likfit(geodata = newgeodata , nugget = 0.01,  ini = c(0.8,23000),fix.nugget = T)
  exp.rml<-likfit(geodata = newgeodata , nugget = 0.01,  ini = c(0.8,22000),fix.nugget = T,method='RML')
  exp.wls<-variofit(vario = variogramaclasico, nugget = 0.1,  ini = c(1.25,22000),fix.nugget = T,weights="npairs")
  
  plot(v$bins, v$robust, pch=19, xlab="Distancia", ylab="Semivarianza")
  lines(exp.ml,max.dist=dmax,lwd=2, col='red')
  lines(exp.rml,max.dist=dmax,lwd=2,col='blue')
  lines(exp.wls,max.dist=dmax,lwd=2,col='green')
  legend(locator(1),legend=c('ML','RML','WLS'),lwd=c(2,2,2),col=c('red','blue','green'))
  

