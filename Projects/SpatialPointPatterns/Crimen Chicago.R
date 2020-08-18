library(splancs)
library(spatstat)
data(chicago)
x11()
plot(chicago, main="Chicago Crimes", col="grey", cols=c("red", "blue", "black", "blue", "red", "blue", "blue"),
     chars=c(16,2,22,17,24,15,6), leg.side="left", show.window=FALSE)

w <- owin(c(min(chicago$data$x),max(chicago$data$x)),c(min(chicago$data$y),max(chicago$data$y)))
Crime.Chi <- ppp(chicago$data$x,chicago$data$y,marks=chicago$data$marks,window=w,fatal=T)

### a. CRS 

x11()
plot(x = chicago$data$x, y = chicago$data$y, pch=16,cex=0.5, col = "red"
     , ylim = c(min(chicago$data$y)-0.1,max(chicago$data$y)+0.1)
     , xlim = c(min(chicago$data$x)-0.1,max(chicago$data$x)+0.1)
     ,xlab = '',ylab='',axes=FALSE)
par(new=TRUE)
plot(quadratcount(Crime.Chi,nx = 6, ny = 3),add=T)

x11()
plot(x = chicago$data$x, y = chicago$data$y, pch=16,cex=0.5, col = "red"
     , ylim = c(min(chicago$data$y)-0.1,max(chicago$data$y)+0.1)
     , xlim = c(min(chicago$data$x)-0.1,max(chicago$data$x)+0.1)
     ,xlab = '',ylab='',axes=FALSE)
par(new=TRUE)
plot(quadratcount(Crime.Chi,nx = 8, ny = 4),add=T)

x11()
plot(x = chicago$data$x, y = chicago$data$y, pch=16,cex=0.5, col = "red"
     , ylim = c(min(chicago$data$y)-0.1,max(chicago$data$y)+0.1)
     , xlim = c(min(chicago$data$x)-0.1,max(chicago$data$x)+0.1)
     ,xlab = '',ylab='',axes=FALSE)
par(new=TRUE)
plot(quadratcount(Crime.Chi,nx = 10, ny = 5),add=T)

#monte carlo

quadrat.test(Crime.Chi,nx = 6, ny = 3,method = "MonteCarlo", nsim=199)

quadrat.test(Crime.Chi,nx = 8, ny = 4,method = "MonteCarlo", nsim=199)

quadrat.test(Crime.Chi,nx = 10, ny = 5,method = "MonteCarlo", nsim=199)


##chi cuadrado

quadrat.test(Crime.Chi,nx = 6, ny = 3)

quadrat.test(Crime.Chi,nx = 8, ny = 4)

quadrat.test(Crime.Chi,nx = 10, ny = 5)

#llambda

lambda <- function(x, y) 20 * (chicago$data$x + chicago$data$y)


### c. FUNCIONES G Y F

X11()
plot(Gest(Crime.Chi),main='Función Gest')

x11()
par(mfrow=c(1,2))
plot(Fest(Crime.Chi),main='Función Fest')
plot(F,type="l",lwd=2.5,main='Construcción de la función Fest')

#=======================================
# Construcción Función G: evento-evento
#=======================================

Muestra<- as.data.frame(chicago$data)
library(spdep)
coordinates(Muestra) <- ~x+y
coords1 <- coordinates(Muestra)    
IDs <- row.names(as(Muestra, "data.frame"))

M_kn1 <- knn2nb(knearneigh(coords1, k=1), row.names=IDs)

#par(mar = c(0, 0, 0, 0), pty = "s")
x11()
plot(M_kn1, coords1, col=4, lwd=2)
points(coords1, pch = 19, col = "Red")
text(coords1[,1]+2,coords1[,2]+2,IDs,cex=1.2)
length((IDs))
M <- as.matrix(dist(coords1))+diag(NA,116,116)
D <- apply(M,1,min,na.rm=T)

I10  <- mean(ifelse(D<10,1,0))
I20 <- mean(ifelse(D<20,1,0))
I30 <- mean(ifelse(D<30,1,0))
I40 <- mean(ifelse(D<40,1,0))
I50 <- mean(ifelse(D<50,1,0))
I60 <- mean(ifelse(D<60,1,0))
I70 <- mean(ifelse(D<70,1,0))

G <- data.frame(Distance=c(0,10,20,30,40,50,60,70),Gest=c(0,I10,I20,I30,I40,I50,I60,I70))
plot(G,type="l", lwd=2.5,main='Construcción de la función Gest')


#======================================
# Construcción Función F: punto-evento
#======================================


set.seed(125)
ptox <- runif(110,Muestra@bbox[1,1],Muestra@bbox[1,2])
ptoy <- runif(110,Muestra@bbox[2,1],Muestra@bbox[2,2])

plot(M_kn1, coords1, col=2, lwd=2)
points(coords1, pch = 19, col = "Red")
text(coords1[,1]+2,coords1[,2]+2,IDs,cex=1.2,col="Red")
Puntos <- data.frame(x=ptox,y=ptoy)
points(Puntos, pch = 19, col = "Black")
text(Puntos[,1]+2,Puntos[,2]+2,rownames(Puntos),cex=1.2)

library(SpatialTools)
D1 <- dist2(as.matrix(Puntos),coords1)

X <- apply(D1,1,function(x) {which(x == min(x))})
D2 <- apply(D1,1,min)

segments(Puntos[,1],Puntos[,2],coords1[X,1],coords1[X,2],col=4, 
         lwd=2.5)

In20  <- mean(ifelse(D2<20,1,0))
In40 <- mean(ifelse(D2<40,1,0))
In50 <- mean(ifelse(D2<50,1,0))
In60 <- mean(ifelse(D2<60,1,0))
In70 <- mean(ifelse(D2<70,1,0))
In80 <- mean(ifelse(D2<80,1,0))
In90 <- mean(ifelse(D2<90,1,0))
In110 <- mean(ifelse(D2<110,1,0))

F <- data.frame(Distance=c(0,20,40,50,60,70,80,90,110),Fest=c(0,In20,In40,In50,In60,In70,In80,In90,In110))
x11()
plot(F,type="l",lwd=2.5,main='Construcción de la función Fest')
par(mfrow=c(2,2))


# D intensidad
x11()
par(mfrow=c(1,3))
Intensidad <- density(Crime.Chi)
plot(Intensidad,main='Densidad')
plot(Crime.Chi, add = TRUE, cex = 0.5)
persp(Intensidad,theta=30,phi=20,main='Perspectiva')
contour(Intensidad,main='Contorno')


#Ancho de banda

xbound<-c(min(Crime.Chi$x),min(Crime.Chi$x),max(Crime.Chi$x),max(Crime.Chi$x))
ybound<-c(min(Crime.Chi$y),max(Crime.Chi$y),max(Crime.Chi$y),min(Crime.Chi$y))
boundary <- as.points(xbound,ybound)
polymap(boundary, col="gray")
Pts <- as.points(Crime.Chi$x, Crime.Chi$y)

Mse2d<-mse2d(Pts,boundary,nsmse=50, range=6000)
plot(Mse2d$h,Mse2d$mse, type="l")

bs<- bw.smoothppp(Crime.Chi,hmin=20000,hmax=100000,kernel="gaussian")
x11()
plot(bs)
