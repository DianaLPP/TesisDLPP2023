
rm(list=ls())
setwd("C:\\Users\\Diana\\Dropbox\\Diana\\CodigoR")

load("pois_5difnodos_15sim.RData")
load("poly_5difnodos_15sim.RData")
load("WS_5difnodos_15sim.RData")
load("Cualidades de red, extra.RData")



d100<- c(diam_pois_100_1, diam_pois_100_2, diam_pois_100_3, 
         diam_poly_100_1, diam_poly_100_2, diam_poly_100_3,
         diam_ws_100_1, diam_ws_100_2, diam_ws_100_3)
d500<- c(diam_pois_500_1, diam_pois_500_2, diam_pois_500_3, 
         diam_poly_500_1, diam_poly_500_2, diam_poly_500_3,
         diam_ws_500_1, diam_ws_500_2, diam_ws_500_3)
d1000<- c(diam_pois_1000_1, diam_pois_1000_2, diam_pois_1000_3, 
          diam_poly_1000_1, diam_poly_1000_2, diam_poly_1000_3,
          diam_ws_1000_1, diam_ws_1000_2, diam_ws_1000_3)
d10000<- c(diam_pois_10000_1, diam_pois_10000_2, diam_pois_10000_3, 
           diam_poly_10000_1, diam_poly_10000_2, diam_poly_10000_3,
           diam_ws_10000_1, diam_ws_10000_2, diam_ws_10000_3)
d50000<- c(diam_pois_50000_1, diam_pois_50000_2, diam_pois_50000_3, 
           diam_poly_50000_1, diam_poly_50000_2, diam_poly_50000_3,
           diam_ws_50000_1, diam_ws_50000_2, diam_ws_50000_3)

m100<- c(med_pois_100_1, med_pois_100_2, med_pois_100_3, 
         med_poly_100_1, med_poly_100_2, med_poly_100_3,
         med_ws_100_1, med_ws_100_2, med_ws_100_3)
m500<- c(med_pois_500_1, med_pois_500_2, med_pois_500_3, 
         med_poly_500_1, med_poly_500_2, med_poly_500_3,
         med_ws_500_1, med_ws_500_2, med_ws_500_3)
m1000<- c(med_pois_1000_1, med_pois_1000_2, med_pois_1000_3, 
          med_poly_1000_1, med_poly_1000_2, med_poly_1000_3,
          med_ws_1000_1, med_ws_1000_2, med_ws_1000_3)
m10000<- c(med_pois_10000_1, med_pois_10000_2, med_pois_10000_3, 
           med_poly_10000_1, med_poly_10000_2, med_poly_10000_3,
           med_ws_10000_1, med_ws_10000_2, med_ws_10000_3)
m50000<- c(med_pois_50000_1, med_pois_50000_2, med_pois_50000_3, 
           med_poly_50000_1, med_poly_50000_2, med_poly_50000_3,
           med_ws_50000_1, med_ws_50000_2, med_ws_50000_3)

v100<- c(var_pois_100_1, var_pois_100_2, var_pois_100_3, 
         var_poly_100_1, var_poly_100_2, var_poly_100_3,
         var_ws_100_1, var_ws_100_2, var_ws_100_3)
v500<- c(var_pois_500_1, var_pois_500_2, var_pois_500_3, 
         var_poly_500_1, var_poly_500_2, var_poly_500_3,
         var_ws_500_1, var_ws_500_2, var_ws_500_3)
v1000<- c(var_pois_1000_1, var_pois_1000_2, var_pois_1000_3, 
          var_poly_1000_1, var_poly_1000_2, var_poly_1000_3,
          var_ws_1000_1, var_ws_1000_2, var_ws_1000_3)
v10000<- c(var_pois_10000_1, var_pois_10000_2, var_pois_10000_3, 
           var_poly_10000_1, var_poly_10000_2, var_poly_10000_3,
           var_ws_10000_1, var_ws_10000_2, var_ws_10000_3)
v50000<- c(var_pois_50000_1, var_pois_50000_2, var_pois_50000_3, 
           var_poly_50000_1, var_poly_50000_2, var_poly_50000_3,
           var_ws_50000_1, var_ws_50000_2, var_ws_50000_3)

c100<- c(coef_pois_100_1, coef_pois_100_2, coef_pois_100_3, 
         coef_poly_100_1, coef_poly_100_2, coef_poly_100_3,
         coef_ws_100_1, coef_ws_100_2, coef_ws_100_3)
c500<- c(coef_pois_500_1, coef_pois_500_2, coef_pois_500_3, 
         coef_poly_500_1, coef_poly_500_2, coef_poly_500_3,
         coef_ws_500_1, coef_ws_500_2, coef_ws_500_3)
c1000<- c(coef_pois_1000_1, coef_pois_1000_2, coef_pois_1000_3, 
          coef_poly_1000_1, coef_poly_1000_2, coef_poly_1000_3,
          coef_ws_1000_1, coef_ws_1000_2, coef_ws_1000_3)
c10000<- c(coef_pois_10000_1, coef_pois_10000_2, coef_pois_10000_3, 
           coef_poly_10000_1, coef_poly_10000_2, coef_poly_10000_3,
           coef_ws_10000_1, coef_ws_10000_2, coef_ws_10000_3)
c50000<- c(coef_pois_50000_1, coef_pois_50000_2, coef_pois_50000_3, 
           coef_poly_50000_1, coef_poly_50000_2, coef_poly_50000_3,
           coef_ws_50000_1, coef_ws_50000_2, coef_ws_50000_3)

a100<- c(assor_pois_100_1, assor_pois_100_2, assor_pois_100_3, 
         assor_poly_100_1, assor_poly_100_2, assor_poly_100_3,
         assor_ws_100_1, assor_ws_100_2, assor_ws_100_3)
a500<- c(assor_pois_500_1, assor_pois_500_2, assor_pois_500_3, 
         assor_poly_500_1, assor_poly_500_2, assor_poly_500_3,
         assor_ws_500_1, assor_ws_500_2, assor_ws_500_3)
a1000<- c(assor_pois_1000_1, assor_pois_1000_2, assor_pois_1000_3, 
          assor_poly_1000_1, assor_poly_1000_2, assor_poly_1000_3,
          assor_ws_1000_1, assor_ws_1000_2, assor_ws_1000_3)
a10000<- c(assor_pois_10000_1, assor_pois_10000_2, assor_pois_10000_3, 
           assor_poly_10000_1, assor_poly_10000_2, assor_poly_10000_3,
           assor_ws_10000_1, assor_ws_10000_2, assor_ws_10000_3)
a50000<- c(assor_pois_50000_1, assor_pois_50000_2, assor_pois_50000_3, 
           assor_poly_50000_1, assor_poly_50000_2, assor_poly_50000_3,
           assor_ws_50000_1, assor_ws_50000_2, assor_ws_50000_3)

md100<- c(mdist_pois_100_1, mdist_pois_100_2, mdist_pois_100_3, 
          mdist_poly_100_1, mdist_poly_100_2, mdist_poly_100_3,
          mdist_ws_100_1, mdist_ws_100_2, mdist_ws_100_3)
md500<- c(mdist_pois_500_1, mdist_pois_500_2, mdist_pois_500_3, 
          mdist_poly_500_1, mdist_poly_500_2, mdist_poly_500_3,
          mdist_ws_500_1, mdist_ws_500_2, mdist_ws_500_3)
md1000<- c(mdist_pois_1000_1, mdist_pois_1000_2, mdist_pois_1000_3, 
           mdist_poly_1000_1, mdist_poly_1000_2, mdist_poly_1000_3,
           mdist_ws_1000_1, mdist_ws_1000_2, mdist_ws_1000_3)
md10000<- c(mdist_pois_10000_1, mdist_pois_10000_2, mdist_pois_10000_3, 
            mdist_poly_10000_1, mdist_poly_10000_2, mdist_poly_10000_3,
            mdist_ws_10000_1, mdist_ws_10000_2, mdist_ws_10000_3)
md50000<- c(mdist_pois_50000_1, mdist_pois_50000_2, mdist_pois_50000_3, 
            mdist_poly_50000_1, mdist_poly_50000_2, mdist_poly_50000_3,
            mdist_ws_50000_1, mdist_ws_50000_2, mdist_ws_50000_3)

e100<- c(eigen_pois_100_1, eigen_pois_100_2, eigen_pois_100_3, 
         eigen_poly_100_1, eigen_poly_100_2, eigen_poly_100_3,
         eigen_ws_100_1, eigen_ws_100_2, eigen_ws_100_3)
e500<- c(eigen_pois_500_1, eigen_pois_500_2, eigen_pois_500_3, 
         eigen_poly_500_1, eigen_poly_500_2, eigen_poly_500_3,
         eigen_ws_500_1, eigen_ws_500_2, eigen_ws_500_3)
e1000<- c(eigen_pois_1000_1, eigen_pois_1000_2, eigen_pois_1000_3, 
          eigen_poly_1000_1, eigen_poly_1000_2, eigen_poly_1000_3,
          eigen_ws_1000_1, eigen_ws_1000_2, eigen_ws_1000_3)
e10000<- c(eigen_pois_10000_1, eigen_pois_10000_2, eigen_pois_10000_3, 
           eigen_poly_10000_1, eigen_poly_10000_2, eigen_poly_10000_3,
           eigen_ws_10000_1, eigen_ws_10000_2, eigen_ws_10000_3)
e50000<- c(eigen_pois_50000_1, eigen_pois_50000_2, eigen_pois_50000_3, 
           eigen_poly_50000_1, eigen_poly_50000_2, eigen_poly_50000_3,
           eigen_ws_50000_1, eigen_ws_50000_2, eigen_ws_50000_3)





colores<-c(rgb(0,0,1,0.5),rgb(1,0,0,0.5),rgb(1,0,1,0.5),
           rgb(0.18, 0.55, 0.18,0.5),rgb(178/255, 142/255, 0,0.5)) #definimos colores con tranparencia


#DIAMETRO
par(mar = c(5, 5, 4, 6))
plot(1:9,ylim=c(0,45),xlab = "", ylab = "", main="Diámetro",xaxt="n", t="n")  #no graficamos ahora
axis(1,at=c(2,5,8),labels = c("Poisson", "Polylogaritmica", "Watts Strogatz"),tck=FALSE)
polygon(c(0.5,3.5,3.5,0.5),c(-1,0,47,47),col=rgb(.9,0.9,0.9,.5),border=NA)
polygon(c(6.5,9.9,9.5,6.5),c(-1,0,47,47),col=rgb(.8,0.8,0.8,.5),border=NA)
lines(d50000,col=colores[1],pch=19,cex=1.2,t="o") 
points(d10000,col=colores[2],pch=19,cex=1.1,t="o")
points(d1000,col=colores[3],pch=19,cex=1,t="o")
points(d500,col=colores[4],pch=19,cex=.9,t="o")
points(d100,col=colores[5],pch=19,cex=.8,t="o")
legend("topright",inset = c(-0.26, 0),pch=19,xpd = TRUE, c("50000","10000","1000","500","100"),
       col= colores,lwd=2,cex=.8,bty="n",title = "Núm. Nodos")




#MEDIA
par(mar = c(5, 5, 4, 6))
plot(1:9,ylim=c(2,4.5),xlab = "", ylab = "", main="Media de los grados",xaxt="n", t="n")  #no graficamos ahora
axis(1,at=c(2,5,8),labels = c("Poisson", "Polylogaritmica", "Watts Strogatz"),tck=FALSE)
polygon(c(0.5,3.5,3.5,0.5),c(0,0,5,5),col=rgb(.9,0.9,0.9,.5),border=NA)
polygon(c(6.5,9.9,9.5,6.5),c(0,0,5,5),col=rgb(.8,0.8,0.8,.5),border=NA)
lines(m50000,col=colores[1],pch=19,cex=1.2,t="o") 
points(m10000,col=colores[2],pch=19,cex=1.1,t="o")
points(m1000,col=colores[3],pch=19,cex=1,t="o")
points(m500,col=colores[4],pch=19,cex=.9,t="o")
points(m100,col=colores[5],pch=19,cex=.8,t="o")
legend("topright",inset = c(-0.26, 0),pch=19,xpd = TRUE, c("50000","10000","1000","500","100"),
       col= colores,lwd=2,cex=.8,bty="n",title = "Núm. Nodos")




#VARIANZA
par(mar = c(5, 5, 4, 6))
plot(1:9,ylim=c(2,5),xlab = "", ylab = "", main="Varianza de los grados",xaxt="n", t="n")  #no graficamos ahora
axis(1,at=c(2,5,8),labels = c("Poisson", "Polylogaritmica", "Watts Strogatz"),tck=FALSE)
polygon(c(0.5,3.5,3.5,0.5),c(0,0,45,45),col=rgb(.9,0.9,0.9,.5),border=NA)
polygon(c(6.5,9.9,9.5,6.5),c(0,0,45,45),col=rgb(.8,0.8,0.8,.5),border=NA)
lines(v50000,col=colores[1],pch=19,cex=1.2,t="o")  
points(v10000,col=colores[2],pch=19,cex=1.1,t="o")
points(v1000,col=colores[3],pch=19,cex=1,t="o")
points(v500,col=colores[4],pch=19,cex=.9,t="o")
points(v100,col=colores[5],pch=19,cex=.8,t="o")
legend("topright",inset = c(-0.26, 0),pch=19,xpd = TRUE, c("50000","10000","1000","500","100"),
       col= colores,lwd=2,cex=.8,bty="n",title = "Núm. Nodos")




#COEFICIENTE DE AGRUPAMIENTO
par(mar = c(5, 5, 4, 6))
plot(1:9,ylim=c(0,0.07),xlab = "", ylab = "", main="Coeficiente de agrupamiento",xaxt="n", t="n")  #no graficamos ahora
axis(1,at=c(2,5,8),labels = c("Poisson", "Polylogaritmica", "Watts Strogatz"),tck=FALSE)
polygon(c(0.5,3.5,3.5,0.5),c(-1,-1,4.5,4.5),col=rgb(.9,0.9,0.9,.5),border=NA)
polygon(c(6.5,9.9,9.5,6.5),c(-1,-1,4.5,4.5),col=rgb(.8,0.8,0.8,.5),border=NA)
lines(c50000,col=colores[1],pch=19,cex=1.2,t="o")  
points(c10000,col=colores[2],pch=19,cex=1.1,t="o")
points(c1000,col=colores[3],pch=19,cex=1,t="o")
points(c500,col=colores[4],pch=19,cex=.9,t="o")
points(c100,col=colores[5],pch=19,cex=.8,t="o")
legend("topright",inset = c(-0.26, 0),pch=19,xpd = TRUE, c("50000","10000","1000","500","100"),
       col= colores,lwd=2,cex=.8,bty="n",title = "Núm. Nodos")





#ASORTATIVIDAD
par(mar = c(5, 5, 4, 6))
plot(1:9,ylim=c(-0.35,0.15),xlab = "", ylab = "", main="Asortatividad",xaxt="n", t="n")  #no graficamos ahora
axis(1,at=c(2,5,8),labels = c("Poisson", "Polylogaritmica", "Watts Strogatz"),tck=FALSE)
polygon(c(0.5,3.5,3.5,0.5),c(-1,-1,4,4),col=rgb(.9,0.9,0.9,.5),border=NA)
polygon(c(6.5,9.9,9.5,6.5),c(-1,-1,4,4),col=rgb(.8,0.8,0.8,.5),border=NA)
lines(a50000,col=colores[1],pch=19,cex=1.2,t="o")  
points(a10000,col=colores[2],pch=19,cex=1.1,t="o")
points(a1000,col=colores[3],pch=19,cex=1,t="o")
points(a500,col=colores[4],pch=19,cex=.9,t="o")
points(a100,col=colores[5],pch=19,cex=.8,t="o")
legend("topright",inset = c(-0.26, 0),pch=19,xpd = TRUE, c("50000","10000","1000","500","100"),
       col= colores,lwd=2,cex=.8,bty="n",title = "Núm. Nodos")




#Distancia media
par(mar = c(5, 5, 4, 6))
plot(1:9,ylim=c(0,9),xlab = "", ylab = "", main="Distancia Media",xaxt="n", t="n")  #no graficamos ahora
axis(1,at=c(2,5,8),labels = c("Poisson", "Polylogaritmica", "Watts Strogatz"),tck=FALSE)
polygon(c(0.5,3.5,3.5,0.5),c(-1,-1,15,15),col=rgb(.9,0.9,0.9,.5),border=NA)
polygon(c(6.5,9.9,9.5,6.5),c(-1,-1,15,15),col=rgb(.8,0.8,0.8,.5),border=NA)
lines(md50000,col=colores[1],pch=19,cex=1.2,t="o")  
points(md10000,col=colores[2],pch=19,cex=1.1,t="o")
points(md1000,col=colores[3],pch=19,cex=1,t="o")
points(md500,col=colores[4],pch=19,cex=.9,t="o")
points(md100,col=colores[5],pch=19,cex=.8,t="o")
legend("topright",inset = c(-0.26, 0),pch=19,xpd = TRUE, c("50000","10000","1000","500","100"),
       col= colores,lwd=2,cex=.8,bty="n",title = "Núm. Nodos")




#Eigen-centralidad
par(mar = c(5, 5, 4, 6))
plot(1:9,ylim=c(3.5,6),xlab = "", ylab = "", main="Eigen-Centralidad",xaxt="n", t="n")  #no graficamos ahora
axis(1,at=c(2,5,8),labels = c("Poisson", "Polylogaritmica", "Watts Strogatz"),tck=FALSE)
polygon(c(0.5,3.5,3.5,0.5),c(0,0,10,10),col=rgb(.9,0.9,0.9,.5),border=NA)
polygon(c(6.5,9.9,9.5,6.5),c(0,0,10,10),col=rgb(.8,0.8,0.8,.5),border=NA)
lines(e50000,col=colores[1],pch=19,cex=1.2,t="o")  
points(e10000,col=colores[2],pch=19,cex=1.1,t="o")
points(e1000,col=colores[3],pch=19,cex=1,t="o")
points(e500,col=colores[4],pch=19,cex=.9,t="o")
points(e100,col=colores[5],pch=19,cex=.8,t="o")
legend("topright",inset = c(-0.26, 0),pch=19,xpd = TRUE, c("50000","10000","1000","500","100"),
       col= colores,lwd=2,cex=.8,bty="n",title = "Núm. Nodos")


