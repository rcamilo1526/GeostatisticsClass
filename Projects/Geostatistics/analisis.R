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
load("D:/Documentos/11 semestre/Geoestad�stica/ProyectoGeoestadistica/code/IDSTA.ov.rda")
load("D:/Documentos/11 semestre/Geoestad�stica/ProyectoGeoestadistica/code/croatia.grilla 25s.2008.rda")
 
IDSTAr.ov <- IDSTA.ov[!is.na(IDSTA.ov$LST2008_07_03),]

#guardar de la base 
  Hrdem<-IDSTAr.ov$HRdem
  Hrdsea<-IDSTAr.ov$HRdsea
  Hrtwi<-IDSTAr.ov$HRtwi
  Lat<-IDSTAr.ov$Lat
  Lon<-IDSTAr.ov$Lon
  #temperatura 
    Tmt<-IDSTAr.ov$LST2008_07_03
    # library(stargazer)
    # 
    # stargazer(data.frame(Hrdem=Hrdem,Hrdsea=Hrdsea,Hrtwi=Hrtwi,Tmt=Tmt,N=N,E=E))

#guardar geodata 
  croacia.geoR <- as.geodata(IDSTAr.ov, data.col = 'HRdem') 
  
    E<-croacia.geoR$coords[,1]
    N<-croacia.geoR$coords[,2]
    # Objeto del tipo geodata (coordenadas y datos)
    x11()
    contourLines(N,E,Hrdem)
    x11()
    points.geodata(croacia.geoR, x.leg=400000, y.leg=4830000,  pt.div="quintile")
    # Graficar los datos (puntos) y ver simbolos graduados, adem�s la opci�n "add.to.plot" es �til para adicionar opciones de la instrucci�n "plot", como en este caso la grilla.
    points.geodata(croacia.geoR, x.leg=400000, y.leg=4830000, main=c("Gr�fico de Intensidades", "de Precipitaci�n"), col.main=3, pt.div="quintile", add.to.plot = TRUE, panel.first = grid())  

    # Para hacer gr�ficos en 3 dimensiones                                    
    library(scatterplot3d)               
      x11()
    scatterplot3d(N, E, Tmt, angle=45, col.main=4, xlab="Coordenada X", ylab="Coordenada Y", zlab="Temperatura", pch=19, color="firebrick3")
    
#ANALISIS DE NORMALIDAD
    
  #graficos
    sd(Tmt)
    x11()
    par(mfrow=c(nrow=1, ncol=3))
    hist(Tmt, freq=F, breaks=10, xlab="Temperatura (�C)", ylab="Frecuencia", col="orange", main ="Histograma")
    curve(dnorm(x, mean=mean(Tmt), sd=sd(Tmt)), add=T, col="red")                  # X no funciona, solo minuscula
    qqnorm(Tmt, ylab="Temperatura (�C)", xlab="Cuantiles te�ricos", col="orange")
    qqline(Tmt)                                                        #Agregarle al Q-Q Plot la l?nea media
    boxplot(Tmt, main="BOX-PLOT", notch=F, horizontal=T, xlab="Temperatura (�C)",col="orange")
    points(mean(Tmt), y=1, pch=1, cex=2)
    par(mfrow=c(nrow=1, ncol=1))
  #pruebas
    shapiro.test(Tmt)
    ks.test(as.numeric(scale(sort(Tmt))),pnorm)

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

#normal con anamorfosis
    #graficos
      sd(Tmt.trans)
      x11()
      par(mfrow=c(nrow=1, ncol=3))
      hist(Tmt.trans, freq=F, breaks=10, xlab="Temperatura (�C)", ylab="Frecuencia", col="orange", main ="Histograma")
      curve(dnorm(x, mean=mean(Tmt.trans), sd=sd(Tmt.trans)), add=T, col="red")                  # X no funciona, solo minuscula
      qqnorm(Tmt.trans, ylab="Temperatura (�C)", xlab="Cuantiles te�ricos", col="orange")
      qqline(Tmt.trans)                                                        #Agregarle al Q-Q Plot la l?nea media
      boxplot(Tmt.trans, main="BOX-PLOT", notch=F, horizontal=T, xlab="Temperatura (�C)",col="orange")
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
  
  library(stargazer)
  stargazer(modelo1,step.model,single.row = T)
  shapiro.test(residuals(step.model))
  ks.test(as.numeric(scale(sort(residuals(step.model)))),pnorm)
  #No hay tendencia, se toma la transformada
    tp.geo <- croacia.geoR
    tp.geo$data<-Tmt.trans
    # tp.geo$data <- residuals(step.model)

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
    
  #superficie tendencia
    library(spatial)
    #aqui cargamos de nuevo la base da datos para poder tener una lista
    tp.ls<-surf.ls(3,N,E,Tmt.trans)
    tp.trsurf<-trmat(tp.ls, -150, 120, 0, 190, 100) #dio el algo en cada pixel
    #al profe no le agrada ni 5
    
    # Representaci�n de superficies de tendencia
    x11()
    dev.off()
    par(pty="s",mar=c(2,2,2,2)) #para la salida
    contour(tp.trsurf) #gr�fico de contornos como topograf�a 
    points(Lon,Lat,pch=20) #superponemos los puntos
    par(mar=c(0,0,0,0)) #si quiero que el borde sea completamente nulo
    image(tp.trsurf)
    points(N,E,pch=20)
    par(mfrow=c(1,1)) #perspectiva inclinado
    persp(tp.trsurf)
    persp(tp.trsurf,theta=60,phi=30,col=2,ltheta=-20,shade=0.25,xlab="N",ylab="E")
    #estamos viendo que tiene valores altos en una direcci�n y valores bajos en otra 
    
    
#ANALISIS DE ANISOTROPIA
    
  library(intamap)
  datares<-cbind(E, N, Tmt.trans)
  data.geo<-as.data.frame(datares)
  coordinates(data.geo) <- ~E+N
  data.geo$value=data.geo$Tmt.trans
  estimateAnisotropy(data.geo) 
  #SIN ANISOTROPIA
  dmax<-max(dist(croacia.geoR$coords))/2
  var.an <- variog4(tp.geo, max.dist=dmax)
  x11()
  plot(var.an) 

  
#SEMIVARIOGRAMAS
  library(gstat)
  library(geoR)
  library(sgeostat)
  library(geospt)
  dmax<-max(dist(croacia.geoR$coords))/2
  
  ############Con tendencia#########
  

  #EXPERIMENTALES
  
  tp.geo <- croacia.geoR
  tp.geo$data<-Tmt.trans-mean(Tmt.trans)
  
  data2.geo<-as.data.frame(tp.geo)
  class(data2.geo) <- c("data.frame")
  data.or<-data2.geo[,c(2,1,3)]
  
  names(data.or) <- c("x", "y", "t")
  point <- point(data.or)
  pair <- pair(point,num.lags=30,maxdist=dmax)
  

  v <- est.variograms(point,pair,"t",trim=0.1) 
  
  ####KS
  x11()
  
  layout(Conf2x2)
  op=par(mfrow=c(2,2))
  plot(v$bins,v$classic,pch=19, col="cyan3", main="Modelo experimental cl�sico", ylab="Semivarianza", xlab="Distancia")
  lines(v$bins,v$classic,col=1,lty=1)
  plot(v$bins,v$robust,pch=19,col="deeppink3", main="Modelo experimental robusto", ylab="Semivarianza", xlab="Distancia")
  lines(v$bins,v$robust,col=1,lty=1)
  plot(v$bins,v$med,pch=19,col="darkorchid3", main="Modelo experimental mediana", ylab="Semivarianza", xlab="Distancia")
  lines(v$bins,v$med,col=1,lty=1)
  plot(v$bins,v$trimmed.mean,pch=19,col="chartreuse3", main="Modelo experimental media recortada", ylab="Semivarianza", xlab="Distancia")
  lines(v$bins,v$trimmed.mean,col=1,lty=1)
  par(op)
  
  #TEORICO
  data2<-cbind(N, E, Tmt.trans-mean(Tmt.trans))
  temp2<-as.data.frame(data2)
  newgeodata<-as.geodata(temp2)

  ##################################################################
  ###################    CLASICO   #################################
  #################################################################
  #beta es c1 init[0], sigma es c0 nugget y phi es a init[1]5.5,"Exp",20000, 0nugget = 0.02,  ini = c(1.6,42000)
  exp.ml.c<-likfit(geodata = newgeodata ,nugget = 0.02,  ini = c(1.6,41000),fix.nugget = T)
  sph.ml.c<-likfit(geodata = newgeodata, nugget = 0.02, ini = c(1.2,64900),cov.model="sph",fix.nugget = T)
  mat.ml.c<-likfit(geodata = newgeodata,nugget = 0.3, ini = c(2,67000),cov.model="mat",kappa=1.5,fix.nugget = T)
  cir.ml.c<-likfit(geodata = newgeodata, nugget = 0.002,ini = c(1.2,51500),cov.model="cir",fix.nugget = T)
  pow.ml.c<-likfit(geodata = newgeodata,nugget = 0.1295,ini = c(0.09330,27600),cov.model="powered.exponential",kappa=1.75,fix.nugget = T)
  lines(pow.ml.c,max.dist=dmax,lwd=2,col='orange')
   summary(pow.ml.c)
  #Aic para todos los variogramas anteriores
  sv.clas<-c(
    exp.ml.c$AIC,
    sph.ml.c$AIC,
    mat.ml.c$AIC,
    cir.ml.c$AIC,
    pow.ml.c$AIC
  )
  
  x11()
  plot(v$bins, v$classic,col=1,ylim=c(0,3), main="Modelo experimental cl�sico", ylab="Semivarianza", xlab="Distancia", pch=19)
  lines(exp.ml.c,max.dist=dmax,lwd=2,col='blue')
  lines(sph.ml.c,max.dist=dmax,lwd=2,col='green')
  lines(mat.ml.c,max.dist=dmax,lwd=2,col='yellow')
  lines(cir.ml.c,max.dist=dmax,lwd=2,col='red')
  lines(pow.ml.c,max.dist=dmax,lwd=2,col='orange')
  legend(locator(1), legend=c('Exponencial','Esf�rico','Matern','Cir','Pow'),lwd=c(2,2,2,2),col=c('blue','green','yellow','red','orange'))
  
  
  #######################################################
  ################### ROBUSTO ##########################
  ######################################################
  
  #beta es c1 init[0], sigma es c0 nugget y phi es a init[1]5.5,"Exp",20000, 0
  exp.ml.r<-likfit(geodata = newgeodata ,nugget = 0.02,  ini = c(1.6,46000),fix.nugget = T)
  sph.ml.r<-likfit(geodata = newgeodata, nugget = 0.02, ini = c(1.2,63000),cov.model="sph",fix.nugget = T)
  mat.ml.r<-likfit(geodata = newgeodata,nugget = 0.25, ini = c(5,50000),cov.model="mat",kappa=1.5,fix.nugget = T)
  cir.ml.r<-likfit(geodata = newgeodata, nugget = 0.002,ini = c(1.2,49000),cov.model="cir",fix.nugget = T)
  pow.ml.r<-likfit(geodata = newgeodata,nugget = 0.1295,ini = c(0.09330,27550),cov.model="powered.exponential",kappa=1.75,fix.nugget = T)
  summary(exp.ml.c)
  #Aic para todos los variogramas anteriores
  sv.robust<-c(
    exp.ml.r$AIC,
    sph.ml.r$AIC,
    mat.ml.r$AIC,
    cir.ml.r$AIC,
    pow.ml.r$AIC
  )
  
  # x11()
  plot(v$bins, v$robust,col=1,ylim=c(0,3), main="Modelo experimental robusto", ylab="Semivarianza", xlab="Distancia", pch=19)
  lines(exp.ml.r,max.dist=dmax,lwd=2,col='blue')
  lines(sph.ml.r,max.dist=dmax,lwd=2,col='green')
  lines(mat.ml.r,max.dist=dmax,lwd=2,col='yellow')
  lines(cir.ml.r,max.dist=dmax,lwd=2,col='red')
  lines(pow.ml.r,max.dist=dmax,lwd=2,col='orange')
  legend(locator(1), legend=c('Exponencial','Esf�rico','Matern','Cir','Pow'),lwd=c(2,2,2,2),col=c('blue','green','yellow','red','orange'))
  
  ##########################################################################
  ###################### MEDIANA###################################
  ##################################################################
  #beta es c1 init[0], sigma es c0 nugget y phi es a init[1]5.5,"Exp",20000, 0
  exp.ml.m<-likfit(geodata = newgeodata ,nugget = 0.02,  ini = c(1.6,43000),fix.nugget = T)
  sph.ml.m<-likfit(geodata = newgeodata, nugget = 0.021, ini = c(1.2,65000),cov.model="sph",fix.nugget = T)
  mat.ml.m<-likfit(geodata = newgeodata,nugget = 0.3, ini = c(2,62000),cov.model="mat",kappa=1.5,fix.nugget = T)
  cir.ml.m<-likfit(geodata = newgeodata, nugget = 0.002,ini = c(1.2,47000),cov.model="cir",fix.nugget = T)
  pow.ml.m<-likfit(geodata = newgeodata,nugget = 0.1291,ini = c(0.09310,27570),cov.model="powered.exponential",kappa=1.75,fix.nugget = T)
  summary(exp.ml.c)
  #Aic para todos los variogramas anteriores
  sv.med<-c(
    exp.ml.m$AIC,
    sph.ml.m$AIC,
    mat.ml.m$AIC,
    cir.ml.m$AIC,
    pow.ml.m$AIC
  )
  
  # x11()
  plot(v$bins, v$med,col=1,ylim=c(0,3), main="Modelo experimental mediana", ylab="Semivarianza", xlab="Distancia", pch=19)
  lines(exp.ml.m,max.dist=dmax,lwd=2,col='blue')
  lines(sph.ml.m,max.dist=dmax,lwd=2,col='green')
  lines(mat.ml.m,max.dist=dmax,lwd=2,col='yellow')
  lines(cir.ml.m,max.dist=dmax,lwd=2,col='red')
  lines(pow.ml.m,max.dist=dmax,lwd=2,col='orange')
  legend(locator(1), legend=c('Exponencial','Esf�rico','Matern','Cir','Pow'),lwd=c(2,2,2,2),col=c('blue','green','yellow','red','orange'))
  
  ##########################################################################
  ###################### MEDIA RECORTADA ###################################
  ##################################################################
  
  #beta es c1 init[0], sigma es c0 nugget y phi es a init[1]5.5,"Exp",20000, 0
  exp.ml.t<-likfit(geodata = newgeodata ,nugget = 0.02,  ini = c(1.6,48100),fix.nugget = T)
  sph.ml.t<-likfit(geodata = newgeodata, nugget = 0.02, ini = c(1.2,63900),cov.model="sph",fix.nugget = T)
  mat.ml.t<-likfit(geodata = newgeodata,nugget = 0.31, ini = c(0.197,61500),cov.model="mat",kappa=1.5,fix.nugget = T)
  cir.ml.t<-likfit(geodata = newgeodata, nugget = 0.0021,ini = c(1.2,47000),cov.model="cir",fix.nugget = T)
  pow.ml.t<-likfit(geodata = newgeodata,nugget = 0.12694,ini = c(0.09310,27500),cov.model="powered.exponential",kappa=1.75,fix.nugget = T)
  summary(exp.ml.c)
  #Aic para todos los variogramas anteriores
  sv.trimmed.mean<-c(
    exp.ml.t$AIC,
    sph.ml.t$AIC,
    mat.ml.t$AIC,
    cir.ml.t$AIC,
    pow.ml.t$AIC
  )
  
  # x11()
  plot(v$bins, v$trimmed.mean,col=1,ylim=c(0,3), main="Modelo experimental median recortada", ylab="Semivarianza", xlab="Distancia", pch=19)
  lines(exp.ml.t,max.dist=dmax,lwd=2,col='blue')
  lines(sph.ml.t,max.dist=dmax,lwd=2,col='green')
  lines(mat.ml.t,max.dist=dmax,lwd=2,col='yellow')
  lines(cir.ml.t,max.dist=dmax,lwd=2,col='red')
  lines(pow.ml.t,max.dist=dmax,lwd=2,col='orange')
  legend(locator(1), legend=c('Exponencial','Esf�rico','Matern','Cir','Pow'),lwd=c(2,2,2,2),col=c('blue','green','yellow','red','orange'))
  data.frame(sv.clas,sv.robust,sv.med,sv.trimmed.mean)
  library(xtable)
  xtable(data.frame(sv.clas,sv.robust,sv.med,sv.trimmed.mean))
  #se corre el modelo escog?do con ML WLS, RML
  # clasico nugget = 0.01,  ini = c(0.8,22000)
  variogramaclasico<-variog( newgeodata,max.dist=dmax, option = "cloud")
  exp.ml<-likfit(geodata = newgeodata , nugget = 0.02,  ini = c(2,41000),fix.nugget = T)
  exp.rml<-likfit(geodata = newgeodata , nugget = 0.01,  ini = c(2.1,44000),fix.nugget = T,method='RML')
  exp.wls<-variofit(vario = variogramaclasico, nugget = 0.02,  ini = c(2,30000),fix.nugget = T,weights="npairs",cov.model = 'exp')
  exp.rml$AIC
  exp.ml$AIC
  exp.wls$AIC
  # X11()
  plot(v$bins, v$classic, pch=19, xlab="Distancia", ylab="Semivarianza",main= "Modelo experimental Cl�sico",ylim=c(0,3))
  lines(exp.ml,max.dist=dmax,lwd=2, col='red')
  lines(exp.rml,max.dist=dmax,lwd=2,col='blue')
  lines(exp.wls,max.dist=dmax,lwd=2,col='green')
  legend(locator(1),legend=c('ML','RML','WLS'),lwd=c(2,2,2),col=c('red','blue','green'))
  
  
  ############KO################
  
  tp.geo2<-croacia.geoR
  tp.geo2$data<-Tmt.trans
  
  #EXPERIMENTALES
  data2.geo<-as.data.frame(tp.geo2)
  class(data2.geo) <- c("data.frame")
  data.or<-data2.geo[,c(2,1,3)]
  
  names(data.or) <- c("x", "y", "t")
  point <- point(data.or)
  pair <- pair(point,num.lags=30,maxdist=dmax)
  
  # library(xtable)
  # library(xlsx)
  v <- est.variograms(point,pair,"t",trim=0.1)
  ex<-as.data.frame(v[,c(2,3,7)])
  ex<-ex[order(ex$bins),]
  names(ex)<-c('h','s','p')
  xtable(ex)
  write.xlsx2(ex,'D:/Documentos/11 semestre/Geoestad�stica/Tareas geoestad�stica/Ejercicios inventados/datos_sv.xlsx')
  # ####Gr�ficos de semivariogramas experimentales
  x11()
  
  layout(Conf2x2)
  op=par(mfrow=c(2,2))
  plot(v$bins,v$classic,pch=19, col="cyan3", main="Modelo experimental cl�sico", ylab="Semivarianza", xlab="Distancia")
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
  #beta es c1 init[0], sigma es c0 nugget y phi es a init[1]5.5,"Exp",20000, 0nugget = 0.02,  ini = c(1.6,42000)
  exp.ml.c<-likfit(geodata = newgeodata ,nugget = 0.02,  ini = c(1.6,42000),fix.nugget = T)
  sph.ml.c<-likfit(geodata = newgeodata, nugget = 0.02, ini = c(1.2,65000),cov.model="sph",fix.nugget = T)
  mat.ml.c<-likfit(geodata = newgeodata,nugget = 0.3, ini = c(2,64000),cov.model="mat",kappa=1.5,fix.nugget = T)
  cir.ml.c<-likfit(geodata = newgeodata, nugget = 0.002,ini = c(1.2,52000),cov.model="cir",fix.nugget = T)
  pow.ml.c<-likfit(geodata = newgeodata,nugget = 0.1319,ini = c(0.09310,27500),cov.model="powered.exponential",kappa=1.75,fix.nugget = T)
  summary(exp.ml.c)
  #Aic para todos los variogramas anteriores
  sv.clas<-c(
    exp.ml.c$AIC,
    sph.ml.c$AIC,
    mat.ml.c$AIC,
    cir.ml.c$AIC,
    pow.ml.c$AIC
  )

  x11()
  plot(v$bins, v$classic,col=1,ylim=c(0,3), main="Modelo experimental cl�sico", ylab="Semivarianza", xlab="Distancia", pch=19)
  lines(exp.ml.c,max.dist=dmax,lwd=2,col='blue')
  lines(sph.ml.c,max.dist=dmax,lwd=2,col='green')
  lines(mat.ml.c,max.dist=dmax,lwd=2,col='yellow')
  lines(cir.ml.c,max.dist=dmax,lwd=2,col='red')
  lines(pow.ml.c,max.dist=dmax,lwd=2,col='orange')
  legend(locator(1), legend=c('Exponencial','Esf�rico','Matern','Cir','Pow'),lwd=c(2,2,2,2),col=c('blue','green','yellow','red','orange'))
  
  
  #######################################################
  ################### ROBUSTO ##########################
  ######################################################
  
  #beta es c1 init[0], sigma es c0 nugget y phi es a init[1]5.5,"Exp",20000, 0
  exp.ml.r<-likfit(geodata = newgeodata ,nugget = 0.02,  ini = c(1.6,46000),fix.nugget = T)
  sph.ml.r<-likfit(geodata = newgeodata, nugget = 0.02, ini = c(1.2,62000),cov.model="sph",fix.nugget = T)
  mat.ml.r<-likfit(geodata = newgeodata,nugget = 0.25, ini = c(5,54000),cov.model="mat",kappa=1.5,fix.nugget = T)
  cir.ml.r<-likfit(geodata = newgeodata, nugget = 0.002,ini = c(1.2,48000),cov.model="cir",fix.nugget = T)
  pow.ml.r<-likfit(geodata = newgeodata,nugget = 0.1295,ini = c(0.09330,27500),cov.model="powered.exponential",kappa=1.75,fix.nugget = T)
  summary(exp.ml.c)
  #Aic para todos los variogramas anteriores
  sv.robust<-c(
    exp.ml.r$AIC,
    sph.ml.r$AIC,
    mat.ml.r$AIC,
    cir.ml.r$AIC,
    pow.ml.r$AIC
  )
  
  x11()
  plot(v$bins, v$robust,col=1,ylim=c(0,3), main="Modelo experimental robusto", ylab="Semivarianza", xlab="Distancia", pch=19)
  lines(exp.ml.r,max.dist=dmax,lwd=2,col='blue')
  lines(sph.ml.r,max.dist=dmax,lwd=2,col='green')
  lines(mat.ml.r,max.dist=dmax,lwd=2,col='yellow')
  lines(cir.ml.r,max.dist=dmax,lwd=2,col='red')
  lines(pow.ml.r,max.dist=dmax,lwd=2,col='orange')
  legend(locator(1), legend=c('Exponencial','Esf�rico','Matern','Cir','Pow'),lwd=c(2,2,2,2),col=c('blue','green','yellow','red','orange'))
  
  ##########################################################################
  ###################### MEDIANA###################################
  ##################################################################
  #beta es c1 init[0], sigma es c0 nugget y phi es a init[1]5.5,"Exp",20000, 0
  exp.ml.m<-likfit(geodata = newgeodata ,nugget = 0.02,  ini = c(1.6,44000),fix.nugget = T)
  sph.ml.m<-likfit(geodata = newgeodata, nugget = 0.021, ini = c(1.2,61000),cov.model="sph",fix.nugget = T)
  mat.ml.m<-likfit(geodata = newgeodata,nugget = 0.3, ini = c(2,62000),cov.model="mat",kappa=1.5,fix.nugget = T)
  cir.ml.m<-likfit(geodata = newgeodata, nugget = 0.002,ini = c(1.2,47500),cov.model="cir",fix.nugget = T)
  pow.ml.m<-likfit(geodata = newgeodata,nugget = 0.1291,ini = c(0.09310,27500),cov.model="powered.exponential",kappa=1.75,fix.nugget = T)
  summary(exp.ml.c)
  #Aic para todos los variogramas anteriores
  sv.med<-c(
    exp.ml.m$AIC,
    sph.ml.m$AIC,
    mat.ml.m$AIC,
    cir.ml.m$AIC,
    pow.ml.m$AIC
  )
  
  x11()
  plot(v$bins, v$med,col=1,ylim=c(0,3), main="Modelo experimental mediana", ylab="Semivarianza", xlab="Distancia", pch=19)
  lines(exp.ml.m,max.dist=dmax,lwd=2,col='blue')
  lines(sph.ml.m,max.dist=dmax,lwd=2,col='green')
  lines(mat.ml.m,max.dist=dmax,lwd=2,col='yellow')
  lines(cir.ml.m,max.dist=dmax,lwd=2,col='red')
  lines(pow.ml.m,max.dist=dmax,lwd=2,col='orange')
  legend(locator(1), legend=c('Exponencial','Esf�rico','Matern','Cir','Pow'),lwd=c(2,2,2,2),col=c('blue','green','yellow','red','orange'))
  
  ##########################################################################
  ###################### MEDIA RECORTADA ###################################
  ##################################################################
  
  #beta es c1 init[0], sigma es c0 nugget y phi es a init[1]5.5,"Exp",20000, 0
  exp.ml.t<-likfit(geodata = newgeodata ,nugget = 0.02,  ini = c(1.6,48000),fix.nugget = T)
  sph.ml.t<-likfit(geodata = newgeodata, nugget = 0.02, ini = c(1.2,64000),cov.model="sph",fix.nugget = T)
  mat.ml.t<-likfit(geodata = newgeodata,nugget = 0.31, ini = c(0.197,61000),cov.model="mat",kappa=1.5,fix.nugget = T)
  cir.ml.t<-likfit(geodata = newgeodata, nugget = 0.0021,ini = c(1.2,46900),cov.model="cir",fix.nugget = T)
  pow.ml.t<-likfit(geodata = newgeodata,nugget = 0.1269,ini = c(0.09310,27500),cov.model="powered.exponential",kappa=1.75,fix.nugget = T)
  summary(exp.ml.c)
  #Aic para todos los variogramas anteriores
  sv.trimmed.mean<-c(
    exp.ml.t$AIC,
    sph.ml.t$AIC,
    mat.ml.t$AIC,
    cir.ml.t$AIC,
    pow.ml.t$AIC
  )
  
  x11()
  plot(v$bins, v$trimmed.mean,col=1,ylim=c(0,3), main="Modelo experimental median recortada", ylab="Semivarianza", xlab="Distancia", pch=19)
  lines(exp.ml.t,max.dist=dmax,lwd=2,col='blue')
  lines(sph.ml.t,max.dist=dmax,lwd=2,col='green')
  lines(mat.ml.t,max.dist=dmax,lwd=2,col='yellow')
  lines(cir.ml.t,max.dist=dmax,lwd=2,col='red')
  lines(pow.ml.t,max.dist=dmax,lwd=2,col='orange')
  legend(locator(1), legend=c('Exponencial','Esf�rico','Matern','Cir','Pow'),lwd=c(2,2,2,2),col=c('blue','green','yellow','red','orange'))
  
  library(xtable)
  xtable(data.frame(sv.clas,sv.robust,sv.med,sv.trimmed.mean))
  #se corre el modelo escog?do con ML WLS, RML
  # clasico nugget = 0.01,  ini = c(0.8,22000)
  variogramaclasico<-variog( newgeodata,max.dist=dmax, option = "cloud")
  exp.ml<-likfit(geodata = newgeodata , nugget = 0.02,  ini = c(1.6,45000),fix.nugget = T)
  exp.rml<-likfit(geodata = newgeodata , nugget = 0.02,  ini = c(1.6,42000),fix.nugget = T,method='RML')
  exp.wls<-variofit(vario = variogramaclasico, nugget = 0.02,  ini = c(2.1,32000),fix.nugget = T,weights="npairs",cov.model = 'exp')
  exp.rml$AIC
  exp.ml$AIC
  exp.wls$AIC
  X11()
  plot(v$bins, v$classic, pch=19, xlab="Distancia", ylab="Semivarianza",main= "Modelo experimental Cl�sico",ylim=c(0,3))
  lines(exp.ml,max.dist=dmax,lwd=2, col='red')
  lines(exp.rml,max.dist=dmax,lwd=2,col='blue')
  lines(exp.wls,max.dist=dmax,lwd=2,col='green')
  legend(locator(1),legend=c('ML','RML','WLS'),lwd=c(2,2,2),col=c('red','blue','green'))
  
  
  
  
  
  
  #############################################33
  