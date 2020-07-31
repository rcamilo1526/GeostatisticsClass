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
load("/home/rcamilo/Documents/geoestadistica/IDSTA.ov.rda")
load("/home/rcamilo/Documents/geoestadistica/croatia.grilla 25s.2008.rda")

IDSTAr.ov <- IDSTA.ov[!is.na(IDSTA.ov$LST2008_07_03),]

#temperatura 
Tmt<-IDSTAr.ov$LST2008_07_03
min(Tmt)
max(Tmt)
#guardar geodata 
croacia.geoR <- as.geodata(IDSTAr.ov, data.col = 'LST2008_07_03') 
E<-croacia.geoR$coords[,1]
N<-croacia.geoR$coords[,2]
data2<-cbind(N, E, Tmt)
temp2<-as.data.frame(data2)
#IDW
p.optimo <- function(p, formula, locations, data, newdata, nmax, nmin, maxdist, var.reg){
  idw.pred <- as.data.frame(matrix(NA,nrow= nrow(data), ncol=4))
  colnames(idw.pred) <- c("x","y","var1.pred","var1.var")
  for(i in 1:(nrow(data))){
    idw.pred[i,] <- idw(formula, locations, data[-i,], newdata[i,], nmax, nmin, maxdist, idp=p)
  } 
  RMSPE <-  sqrt(sum((idw.pred$var1.pred-var.reg)^2)/nrow(data))
  RMSPE
}
coordinates(IDSTA.OV)
grillabaser<-as.data.frame(IDSTA.OV) #se vuelve data frame
puntosxyi<-grillabaser[,52:53] #se seleccionan las columnas donde estan las coordenadas x, y
names(puntosxyi) <- c("N","E") #esta en x, y entonces cambiamos a N y E para ppoder hacer las RBF
IDW <-  temp2

###CALCULAR OPTIMO
P2 <- optimize(p.optimo, c(0,10), formula=Tmt~1, locations=~N+E, data=IDW, newdata=IDW, nmax=10, nmin=10, maxdist=Inf, var.reg=Tmt)
cat("Parámetro óptimo IDW: ", "\n", "P       = ", P2$minimum, "\n", "RMSPE   = ", P2$objective, "\n")

P.opt <- optim( par=1 , fn=p.optimo, gr= "Nelder-Mead", method = "Nelder-Mead", formula=Tmt~1, locations=~N+E, data=IDW, newdata=IDW, nmax=10, nmin=10, maxdist=Inf, var.reg=Tmt)
P.opt <- optim( par=1 , fn=p.optimo, method = "L-BFGS-B", formula=Tmt~1, locations=~N+E, data=IDW, newdata=IDW, nmax=10, nmin=10, maxdist=Inf, var.reg=Tmt)

p1 <- p.optimo(p=1, formula=Tmt~1, locations=~N+E, data=IDW, newdata=IDW, nmax=10, nmin=10, maxdist= Inf, var.reg=Tmt)
p2 <- p.optimo(p=2, formula=Tmt~1, locations=~N+E, data=IDW, newdata=IDW, nmax=10, nmin=10, maxdist=Inf, var.reg=Tmt)
p3 <- P2$objective
p4 <- p.optimo(p=3, formula=Tmt~1, locations=~N+E, data=IDW, newdata=IDW, nmax=10, nmin=10, maxdist=Inf, var.reg=Tmt)
p5 <- p.optimo(p=4, formula=Tmt~1, locations=~N+E, data=IDW, newdata=IDW, nmax=10, nmin=10, maxdist=Inf, var.reg=Tmt)
p6 <- p.optimo(p=5, formula=Tmt~1, locations=~N+E, data=IDW, newdata=IDW, nmax=10, nmin=10, maxdist=Inf, var.reg=Tmt)

RMSPE <- c(p1, p2, p3, p4, p5, p6)
plot(c(1,2,P2$minimum,3,4,5),RMSPE, main="Gráfico de optimización del Parámetro (P)\n Distancia Inversa Ponderada", ylab="RMSPE", xlab="p ?ptimo = 1.903656", type="l")
####Interpolacion###########

idw.p <- idw(Tmt~1,~ N+E, IDW, puntosxyi, nmax=10, nmin=10, idp=1.577517)
grafi<-as.data.frame(idw.p)
coordinates(grafi) = ~N+E
gridded(grafi)=T

min(grafi$var1.pred)
max(grafi$var1.pred)


spplot(grafi ,"var1.pred", main="Interpolaciones de Distancia Inversa\n Ponderada de la Precipitación\np=1.5775117", col.regions=bpy.colors(101), cuts=100, cex.main=0.5, regions=T, scales = list(draw =T), xlab="Este (m)", ylab = "Norte (m)",  key.space=list(space="left", cex=0.6))

