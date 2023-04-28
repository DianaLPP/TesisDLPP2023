

rm(list=ls())
setwd("C:\\Users\\Diana\\Dropbox\\Diana\\CodigoR")
load("Brotes_SIR_Pois_Poly_WS.RData")
load("Brotes_SEIR_Pois_Poly_WS.RData")

###################################################################################


graf_infectados_SIR<-function(nodos){
  if (nodos==100) {
    x<-400
    y<-8
    ylim<-c(0,y)
    xlim<-c(0,x)
  }
  else if (nodos==500) {
    x<-600
    y<-20
    ylim<-c(0,y)
    xlim<-c(0,x)
  }
  else if (nodos==1000) {
    x<-500
    y<-30
    ylim<-c(0,y)
    xlim<-c(0,x)
  }
  else if (nodos==10000) {
    x<-700
    y<-200
    ylim<-c(0,y)
    xlim<-c(0,x)
  }
  else if (nodos==50000) {
    x<-900
    y<-950
    ylim<-c(0,y)
    xlim<-c(0,x)
  }
  else print("Solo se aceptan 100,500,1000,10000,50000 nodos")
 
  mat_valores_hora<-matrix(NA,x,45)
  ii<-1  
   
 for (t in 1:3){
   if (t==1) {
     red<-c("pois")
     col=rgb(0, 0, 1, 0.3) #azul
   }
   else if (t==2) {
     red<-c("poly")
     col=rgb(1, 0, 0, 0.3) #rojo
   }
   else {
     red<-c("ws")
     col=rgb(0, 0.8, 0, 0.3) #verde
   }
   
    for (r in 1:3){
      for (k in 1:5){
        aa<-get(paste("sim_SIR",red,nodos,r,k,sep="_"))$indexes$InfTime
        aa<-sort(aa[!is.na(aa)])
        interv_SIR<-seq(0,floor(max(aa)+1)+1,by=1)
        pmediointerv_SIR<-(interv_SIR[-1]+interv_SIR[-length(interv_SIR)])/2
        newinfdiscr_SIR<-table(cut(c(aa,pmediointerv_SIR),interv_SIR,label=FALSE,include.lowest = TRUE))-1
        datos<-data.frame(pmediointerv_SIR,newinfdiscr_SIR) 
        if (r==1 && k==1 && t==1) matplot(datos$Var1,datos$Freq,col=col,main=paste("Brotes SIR"," de ",nodos," nodos",sep = ""),cex.main=2, pch=1,cex=0.5,ylim=ylim,xlim=xlim,xlab = "Horas",ylab = "Infectados",cex.lab=1.3, axes=FALSE)
        else points(datos$Var1,datos$Freq,col=col, cex=0.5) 
          axis(1,at=seq(0, x, 10))
          axis(2,at=seq(0, y, 1))
        mat_valores_hora[,ii]<-c(datos$Freq[1:nrow(datos)],rep(0,(x-nrow(datos)))) # Promedio general para los 3 tipos de redes
        ii<-ii+1
      } #for k
    } #for r
  } #for t
  legend("topright", c("Poisson","Poly-log","Watts-S","Promedio"),col= c(rgb(0, 0, 1, 0.3),rgb(1, 0, 0, 0.3),rgb(0, 0.8, 0, 0.3),"black"),lwd=2,cex=1.3,bty="n")
  no_casos<-rowSums(!is.na(mat_valores_hora))
  promedio<-rowSums(mat_valores_hora,na.rm=TRUE)/no_casos
  lines(1:x,promedio,col="black",lwd=2)
} #funcion

#Aqui se ilustran las 5 graficas SIR que se generan respecto a los nodos
pdf("Brotes juntos NewInf, 5 graficas SIR.pdf")
graf_infectados_SIR(100)
graf_infectados_SIR(500)
graf_infectados_SIR(1000)
graf_infectados_SIR(10000)
graf_infectados_SIR(50000)
dev.off()




###################################################################################


graf_infectados_SEIR<-function(nodos){
  if (nodos==100) {
    x<-650
    y<-6
    ylim<-c(0,y)
    xlim<-c(0,x)
  }
  else if (nodos==500) {
    x<-750
    y<-14
    ylim<-c(0,y)
    xlim<-c(0,x)
  }
  else if (nodos==1000) {
    x<-800
    y<-18
    ylim<-c(0,y)
    xlim<-c(0,x)
  }
  else if (nodos==10000) {
    x<-1250
    y<-110
    ylim<-c(0,y)
    xlim<-c(0,x)
  }
  else if (nodos==50000) {
    x<-1400
    y<-420
    ylim<-c(0,y)
    xlim<-c(0,x)
  }
  else print("Solo se aceptan 100,500,1000,10000,50000 nodos")
  
  mat_valores_hora<-matrix(NA,x,45)
  ii<-1  
  
  for (t in 1:3){
    if (t==1) {
      red<-c("pois")
      col=rgb(0, 0, 1, 0.3) #azul
    }
    else if (t==2) {
      red<-c("poly")
      col=rgb(1, 0, 0, 0.3) #rojo
    }
    else {
      red<-c("ws")
      col=rgb(0, 0.8, 0, 0.3) #verde
    }
    
    for (r in 1:3){
      for (k in 1:5){
        aa<-get(paste("sim_SEIR",red,nodos,r,k,sep="_"))$indexes$InfTime
        aa<-sort(aa[!is.na(aa)])
        interv_SEIR<-seq(0,floor(max(aa)+1)+1,by=1)
        pmediointerv_SEIR<-(interv_SEIR[-1]+interv_SEIR[-length(interv_SEIR)])/2
        newinfdiscr_SEIR<-table(cut(c(aa,pmediointerv_SEIR),interv_SEIR,label=FALSE,include.lowest = TRUE))-1
        datos<-data.frame(pmediointerv_SEIR,newinfdiscr_SEIR) 
        if (r==1 && k==1 && t==1) matplot(datos$Var1,datos$Freq,col=col,main=paste("Brotes SEIR"," de ",nodos," nodos",sep = ""),cex.main=2, pch=1,cex=0.5,ylim=ylim,xlim=xlim,xlab = "Horas",ylab = "Infectados",cex.lab=1.3, axes=FALSE)
        else points(datos$Var1,datos$Freq,col=col, cex=0.5) 
        axis(1,at=seq(0, x, 10))
        axis(2,at=seq(0, y, 1))
        mat_valores_hora[,ii]<-c(datos$Freq[1:nrow(datos)],rep(0,(x-nrow(datos)))) # Promedio general para los 3 tipos de redes
        ii<-ii+1
      } #for k
    } #for r
  } #for t
  legend("topright", c("Poisson","Poly-log","Watts-S","Promedio"),col= c(rgb(0, 0, 1, 0.3),rgb(1, 0, 0, 0.3),rgb(0, 0.8, 0, 0.3),"black"),lwd=2,cex=1.3,bty="n")
  no_casos<-rowSums(!is.na(mat_valores_hora))
  promedio<-rowSums(mat_valores_hora,na.rm=TRUE)/no_casos
  lines(1:x,promedio,col="black",lwd=2)
} #funcion

##Aqui se ilustran las 5 graficas SIR que se generan respecto a los nodos
pdf("Brotes juntos NewInf, 5 graficas SEIR.pdf")
graf_infectados_SEIR(100)
graf_infectados_SEIR(500)
graf_infectados_SEIR(1000)
graf_infectados_SEIR(10000)
graf_infectados_SEIR(50000)
dev.off()






















######################################################################################3




###promedio por cada hora de todos los brotes, por cada distribucion por separado manualmente
#se cambia tipo de red y numero de nodos (SIR)

graf_infectados_SIR(50000)

red<-"pois"
nodos<-50000
mat_valores_hora<-matrix(NA,525,15)
ii<-1
for (r in 1:3){
  for (k in 1:5){
    aa<-get(paste("sim_SIR",red,nodos,r,k,sep="_"))$indexes$InfTime
    aa<-sort(aa[!is.na(aa)])
    interv_SIR<-seq(0,floor(max(aa)+1)+1,by=1)
    pmediointerv_SIR<-(interv_SIR[-1]+interv_SIR[-length(interv_SIR)])/2
    newinfdiscr_SIR<-table(cut(c(aa,pmediointerv_SIR),interv_SIR,label=FALSE,include.lowest = TRUE))-1
    datos<-data.frame(pmediointerv_SIR,newinfdiscr_SIR) 
    mat_valores_hora[,ii]<-datos$Freq[1:525]
    ii<-ii+1
  }
}

no_casos<-rowSums(!is.na(mat_valores_hora))
promedio<-rowSums(mat_valores_hora,na.rm=TRUE)/no_casos

lines(1:525,promedio,col="yellow",lwd=2)

