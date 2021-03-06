# clear the workspace, plots and console
rm(list = ls())
if(!is.null(dev.list())) dev.off()
cat("\014")
        ########################################################################################
        ###########                    Ejercicio Montecarlo Prevosti               #############
        ########################################################################################

# Instrucci�n generaci�n matrices cuando se ingresan valores de una parte tri�ngular
matriz <- function(vect.tri){
n.filas<-function(vect.tri){k<-1
                            while (sum(seq(1:k))<=length(vect.tri))
                            k<-k+1
                            k
                            }
n.filas1<-n.filas(vect.tri)
X <- vect.tri
m.i <- matrix(0, n.filas1, n.filas1)
m.i[lower.tri(m.i)] <- NA
m.i[is.na(m.i)] <- X
ma.i <- t(m.i)+m.i
print(ma.i)
}

# FUNCI�N MANTEL p.prueba
p.prueba <- function(B,matriz.sim, matriz.dist){                      #system.time({
set.seed(127)
nf <- nrow(matriz.sim)
low.tri <- lower.tri(matriz.dist)
prod.datos <- sum(matriz.sim[low.tri]*matriz.dist[low.tri])
conteo <- replicate(B, {
                    filas.perm <- sample(1:nf)
                    ifelse(sum(matriz.sim[low.tri]*matriz.dist[filas.perm,filas.perm][low.tri])>=prod.datos,1,0)
                    }
                    )
pvalor <- (sum(conteo)+1)/(length(conteo)+1)
pvalor
                                                                      #})
}

# Funci�n Mantel p.Mantel (Propuesto)
p.Mantel <- function(B,matriz.sim, matriz.dist, Gr�fico=F){           #system.time({
set.seed(127)
nf <- nrow(matriz.sim)
low.tri <- lower.tri(matriz.dist)
t.tilde <- sum(matriz.sim[low.tri]*matriz.dist[low.tri])
t.estrella<-numeric(B)
for (i in 1:B){
filas.perm <- sample(1:nf)
t.estrella[i]<-sum(matriz.sim[low.tri]*matriz.dist[filas.perm,filas.perm][low.tri])
}
pvalor <- (sum(t.estrella>=t.tilde)+1)/(B+1)

if (Gr�fico) {
       plot(density(t.estrella),type="l", col=2,
       main="Distribuci�n Aproximada del Estad�stico de Mantel",
       xlab="t.estrella",ylab="N�mero de Permutaciones",
       sub=paste("Est. Observado (t.tilde)=",round(t.tilde,3),
        ":", "Cola Superior P-value=",round(pvalor,4),":", B, "Permutaciones"),
       lwd=1.5, col.main="blue", col.lab="blue", cex.lab=1.2, cex.sub=0.7)
       abline(v=t.tilde, lty=2, col=3)
   }
list(T=t.tilde, T.Simetrico=t.tilde/2, Cola.Superior.P=pvalor)
                                                                      #})
}
# Esta funci�n es m�s eficiente en tiempo, y adem�s ofrece la posibilidad de gr�ficar la funci�n de densidad del test y de mostrar algunos resultados asociados al mismo.


         #######################################################################################
         ########################             Ingreso Datos           ##########################
         #######################################################################################

# La matrices de barreras consideran 1 si hay de por medio una barrera (mar o monta�a) entre los dos lugares y cero si no hay barreras.

length(disesp)
disclim <- read.table('distclim.txt', header = T,sep = ",", dec = ".")
disesp <- read.table('distesp.txt', header = T, sep = ",", dec = ".")
TEM<-disclim$TEM
EA<-disesp$EA

tabla.TEM<-matriz(TEM)
tabla.EA<-matriz(EA)

# Test de Mantel
p.prueba(10000, tabla.TEM,tabla.EA)



# evidencia estad�stica para rechazar "Ho: No hay relaci�n entre elementos en las matrices de
# distancias gen�ticas y de barreras". Y en los casos del Estrecho de Gibraltar y de los
# Pirineos no se rechaza Ho.


# Soluci�n a partir de la funci�n para el Test de Mantel "propuesto"
par(mfrow = c(4,3))
X11()
p.Mantel(10000, tabla.TEM,tabla.EA, Gr�fico=T)



                             ##############################################
                             # Trabajo con library ecodist test de Mantel #
                             ##############################################


# Esta libreria tambien permite hacer el test de Mantel, trabaja internamente las dos posibilidades de ">=" y "<=", para relaciones directas e inversas respectivamente, y adem�s ofrece alternativas para la hip�tesis nula (en �ste ejercicio ser�a r=0). Finalmente, la instrucci�n mgram permite obtener un mantelgrama, algo as� como un correlograma con el test de Mantel, y este relaciona las diferentes distancias con los valores asociados de Mantel. La instrucci�n en �ste caso es "mantel" y dentro se referencian las dos matrices a comparar, luego se establece el n�mero de permutaciones, es importante que estas matrices esten como distancias (es decir, en �ste caso se usar�a: as.dist().)
# Los resultados obtenidos mediante esta instrucci�n corroboran los resultados obtenidos anteriormente. Y con los del ejercicio que ya habiamos realizado con los datos de la deriva continental, haciendo la observaci�n que en dicho caso se trabajaba con "<=" (debido a la relaci�n inversa) y con esta funci�n "mantel", de la libreria "ecodist" es m�s simple.
library(ecodist)     
# Disponible para versi�n: R.2.7.0
dist.TEM<-as.dist(tabla.TEM)
dist.EA<-as.dist(tabla.EA)

t.est.T5<-mantel(dist.EA~dist.TEM,nperm=10000)                # Ho: r<=0 (pval1), r>=0 (pval2), r=0 (pval3)

t.mgram.T5<-mgram(dist.EA, dist.TEM, nperm=0)
x11()
plot(t.mgram.T5)



library(ape)
set.seed(127)
x11()
mantel.test(tabla.TEM,tabla.EA, graph = TRUE)


