
# Establecer directorio de trabajo: el mismo donde desempaquetamos los archivos descargados
setwd("C:/Users/rcami/Documents/11 Semestre/Geoestadística/Clase 04/NUTS")

# Limpia el espacio de trabajo, las parcelas y la consola
rm(list = ls())
if(!is.null(dev.list())) dev.off()
cat("\014") 

# Paquetes para trabajar con mapas y visualizar datos en mapas
library(rgdal)
library(sp)

list.files('C:/Users/rcami/Documents/11 Semestre/Geoestadística/Clase 04/NUTS', pattern='\\.shp$')
file.exists('C:/Users/rcami/Documents/11 Semestre/Geoestadística/Clase 04/NUTS/NUTS_RG_01M_2013.shp')

# Mapa de importaci?n 1 - nivel de pobreza, proyecci?n correcta (source: Centro Central de Documentaci?n Geod?sica y Cartogr?fica, 
# http://www.codgik.gov.pl/index.php/darmowe-dane/prg.html)
mapa1 <- readOGR(".", "NUTS_RG_01M_2013")
# Este mapa es preciso y est? bien descrito, pero las coordenadas est?n codificadas de una manera diferente a la que necesitamos.
# Deber?amos recalcularlos en grados de longitud y latitud.
mapa1 <- spTransform(mapa1, "+proj=longlat")
plot(mapa1)


setwd("D:/MGN2018_00_COLOMBIA/WGS84_MGN2019_00_COLOMBIA/ADMINISTRATIVO")
deptos <- readOGR(".", "MGN_DPTO_POLITICO")
plot(deptos)
