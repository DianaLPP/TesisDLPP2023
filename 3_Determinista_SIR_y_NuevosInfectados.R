
setwd("C:\\Users\\Diana\\Dropbox\\Diana\\CodigoR")
source("04_funciones_Redes_LR.R")
load("Brotes_SIR_Pois_Poly_WS.RData")

library(deSolve)
betareal<-0.03
gamareal<-0.01


#### Solver del modelo SIR determinista
X_theta_SIR<-function(theta_SIR,t=tiempos){
  SIRmod <- function(t, x, theta_SIR) { 
    with(as.list(c(theta_SIR, x)), 
         { 
           ds <- -bet*s*i 
           di <- bet*s*i -gam*i
           dr <- gam*i
           die<- bet*s*i
           res <- c(ds, di, dr, die) 
           list(res) }
    ) }
  ## Solver
  out <- lsoda(X_ini_SIR, t, SIRmod, theta_SIR)
  return(out)
}


########################################################################################


sim_SIR_pois<-list(sim_SIR_pois_100_1_1, sim_SIR_pois_100_1_2, sim_SIR_pois_100_1_3, sim_SIR_pois_100_1_4, sim_SIR_pois_100_1_5, 
                   sim_SIR_pois_100_2_1, sim_SIR_pois_100_2_2, sim_SIR_pois_100_2_3, sim_SIR_pois_100_2_4, sim_SIR_pois_100_2_5,
                   sim_SIR_pois_100_3_1, sim_SIR_pois_100_3_2, sim_SIR_pois_100_3_3, sim_SIR_pois_100_3_4, sim_SIR_pois_100_3_5,
                   sim_SIR_pois_500_1_1, sim_SIR_pois_500_1_2, sim_SIR_pois_500_1_3, sim_SIR_pois_500_1_4, sim_SIR_pois_500_1_5, 
                   sim_SIR_pois_500_2_1, sim_SIR_pois_500_2_2, sim_SIR_pois_500_2_3, sim_SIR_pois_500_2_4, sim_SIR_pois_500_2_5,
                   sim_SIR_pois_500_3_1, sim_SIR_pois_500_3_2, sim_SIR_pois_500_3_3, sim_SIR_pois_500_3_4, sim_SIR_pois_500_3_5,
                   sim_SIR_pois_1000_1_1, sim_SIR_pois_1000_1_2, sim_SIR_pois_1000_1_3, sim_SIR_pois_1000_1_4, sim_SIR_pois_1000_1_5, 
                   sim_SIR_pois_1000_2_1, sim_SIR_pois_1000_2_2, sim_SIR_pois_1000_2_3, sim_SIR_pois_1000_2_4, sim_SIR_pois_1000_2_5,
                   sim_SIR_pois_1000_3_1, sim_SIR_pois_1000_3_2, sim_SIR_pois_1000_3_3, sim_SIR_pois_1000_3_4, sim_SIR_pois_1000_3_5,
                   sim_SIR_pois_10000_1_1, sim_SIR_pois_10000_1_2, sim_SIR_pois_10000_1_3, sim_SIR_pois_10000_1_4, sim_SIR_pois_10000_1_5, 
                   sim_SIR_pois_10000_2_1, sim_SIR_pois_10000_2_2, sim_SIR_pois_10000_2_3, sim_SIR_pois_10000_2_4, sim_SIR_pois_10000_2_5,
                   sim_SIR_pois_10000_3_1, sim_SIR_pois_10000_3_2, sim_SIR_pois_10000_3_3, sim_SIR_pois_10000_3_4, sim_SIR_pois_10000_3_5,
                   sim_SIR_pois_50000_1_1, sim_SIR_pois_50000_1_2, sim_SIR_pois_50000_1_3, sim_SIR_pois_50000_1_4, sim_SIR_pois_50000_1_5, 
                   sim_SIR_pois_50000_2_1, sim_SIR_pois_50000_2_2, sim_SIR_pois_50000_2_3, sim_SIR_pois_50000_2_4, sim_SIR_pois_50000_2_5,
                   sim_SIR_pois_50000_3_1, sim_SIR_pois_50000_3_2, sim_SIR_pois_50000_3_3, sim_SIR_pois_50000_3_4, sim_SIR_pois_50000_3_5)

for(j in 1:75){
  aa<-sim_SIR_pois[[j]]$indexes$InfTime
  aa<-sort(aa[!is.na(aa)])
  interv_SIR<-seq(0,floor(max(aa)+1)+24,by=24)
  pmediointerv_SIR<-(interv_SIR[-1]+interv_SIR[-length(interv_SIR)])/2
  newinfdiscr_SIR<-table(cut(c(aa,pmediointerv_SIR),interv_SIR,label=FALSE,include.lowest = TRUE))-1
  datos<-data.frame(pmediointerv_SIR,newinfdiscr_SIR)
  #plot(datos$pmediointerv_SIR,datos$Freq,t="o",pch=19,ylab="Infectados",xlab="Tiempo")
  
  #EJEMPLO SIR ODE:
  if (j<=15) N<-100
  else if (j<=30) N<-500
  else if (j<=45) N<-1000
  else if (j<=60) N<-10000
  else N<-50000
  
  theta_SIR<-c(bet=betareal,gam=gamareal)
  i=2/N
  X_ini_SIR=c(s = 1-i, i = i, r = 0, ie=i) 
  
  r<-X_theta_SIR(theta_SIR=theta_SIR,t=seq(0,1000,length=1000))
  
  interv_SIR<-seq(0,max(r[,1])+24,by=24)
  pmediointerv_SIR<-(interv_SIR[-1]+interv_SIR[-length(interv_SIR)])/2
  clase<-cut(r[,1],interv_SIR,label=FALSE,include.lowest = TRUE)
  maxclas<-function(i){max(which(clase==i))}
  indices<-sapply(1:max(clase),FUN=maxclas)
  iemas0<-c(0,r[,5])
  newinfdiscr_SIR<-iemas0[c(1,indices+1)]
  newinfdiscr_SIR<-newinfdiscr_SIR[-1]-newinfdiscr_SIR[-length(newinfdiscr_SIR)]
  
  #setwd("d:\\Users\\Diana\\Documents\\Tesis CIMAT\\Cluster\\Poisson\\EpidyODE")
  #pdf(paste("pois_SIR", j, "pdf", sep = "."))
  plot(datos$pmediointerv_SIR,datos$Freq,t="o",pch=19,ylab="Infectados",xlab="Tiempo")
  lines(pmediointerv_SIR,newinfdiscr_SIR*N,t="o",pch=19,col="red",ylab="Infectados",xlab="Tiempo")
  legend("topright", c("M. epidemiológico","M. Determinista"),col= c("black", "red"),lwd=2)
  #dev.off()
}


########################################################################################


sim_SIR_poly<-list(sim_SIR_poly_100_1_1, sim_SIR_poly_100_1_2, sim_SIR_poly_100_1_3, sim_SIR_poly_100_1_4, sim_SIR_poly_100_1_5, 
                   sim_SIR_poly_100_2_1, sim_SIR_poly_100_2_2, sim_SIR_poly_100_2_3, sim_SIR_poly_100_2_4, sim_SIR_poly_100_2_5,
                   sim_SIR_poly_100_3_1, sim_SIR_poly_100_3_2, sim_SIR_poly_100_3_3, sim_SIR_poly_100_3_4, sim_SIR_poly_100_3_5,
                   sim_SIR_poly_500_1_1, sim_SIR_poly_500_1_2, sim_SIR_poly_500_1_3, sim_SIR_poly_500_1_4, sim_SIR_poly_500_1_5, 
                   sim_SIR_poly_500_2_1, sim_SIR_poly_500_2_2, sim_SIR_poly_500_2_3, sim_SIR_poly_500_2_4, sim_SIR_poly_500_2_5,
                   sim_SIR_poly_500_3_1, sim_SIR_poly_500_3_2, sim_SIR_poly_500_3_3, sim_SIR_poly_500_3_4, sim_SIR_poly_500_3_5,
                   sim_SIR_poly_1000_1_1, sim_SIR_poly_1000_1_2, sim_SIR_poly_1000_1_3, sim_SIR_poly_1000_1_4, sim_SIR_poly_1000_1_5, 
                   sim_SIR_poly_1000_2_1, sim_SIR_poly_1000_2_2, sim_SIR_poly_1000_2_3, sim_SIR_poly_1000_2_4, sim_SIR_poly_1000_2_5,
                   sim_SIR_poly_1000_3_1, sim_SIR_poly_1000_3_2, sim_SIR_poly_1000_3_3, sim_SIR_poly_1000_3_4, sim_SIR_poly_1000_3_5,
                   sim_SIR_poly_10000_1_1, sim_SIR_poly_10000_1_2, sim_SIR_poly_10000_1_3, sim_SIR_poly_10000_1_4, sim_SIR_poly_10000_1_5, 
                   sim_SIR_poly_10000_2_1, sim_SIR_poly_10000_2_2, sim_SIR_poly_10000_2_3, sim_SIR_poly_10000_2_4, sim_SIR_poly_10000_2_5,
                   sim_SIR_poly_10000_3_1, sim_SIR_poly_10000_3_2, sim_SIR_poly_10000_3_3, sim_SIR_poly_10000_3_4, sim_SIR_poly_10000_3_5,
                   sim_SIR_poly_50000_1_1, sim_SIR_poly_50000_1_2, sim_SIR_poly_50000_1_3, sim_SIR_poly_50000_1_4, sim_SIR_poly_50000_1_5, 
                   sim_SIR_poly_50000_2_1, sim_SIR_poly_50000_2_2, sim_SIR_poly_50000_2_3, sim_SIR_poly_50000_2_4, sim_SIR_poly_50000_2_5,
                   sim_SIR_poly_50000_3_1, sim_SIR_poly_50000_3_2, sim_SIR_poly_50000_3_3, sim_SIR_poly_50000_3_4, sim_SIR_poly_50000_3_5)

Nodos<-c(100,500,1000,10000,50000)


for(j in 1:75){
  #  plot(sim_SIR_poly[[j]]$Ctime.hist, sim_SIR_poly[[j]]$Istatus.hist, t="l",ylab="Infectados",xlab="Tiempo en que ocurre un evento",lwd=2)
  aa<-sim_SIR_poly[[j]]$indexes$InfTime
  aa<-sort(aa[!is.na(aa)])
  #aa
  #hist(aa)
  interv_SIR<-seq(0,floor(max(aa)+1)+24,by=24)
  pmediointerv_SIR<-(interv_SIR[-1]+interv_SIR[-length(interv_SIR)])/2
  newinfdiscr_SIR<-table(cut(c(aa,pmediointerv_SIR),interv_SIR,label=FALSE,include.lowest = TRUE))-1
  datos<-data.frame(pmediointerv_SIR,newinfdiscr_SIR)
  #plot(datos$pmediointerv_SIR,datos$Freq,t="o",pch=19,ylab="Infectados",xlab="Tiempo")
  
  #EJEMPLO SIR ODE:
  if (j<=15) N<-100
  else if (j<=30) N<-500
  else if (j<=45) N<-1000
  else if (j<=60) N<-10000
  else N<-50000
  
  theta_SIR<-c(bet=betareal,gam=gamareal)
  i=2/N
  X_ini_SIR=c(s = 1-i, i = i, r = 0, ie=i)  
  
  #5% de individuos infectados, 95% son susceptibles
  r<-X_theta_SIR(theta_SIR=theta_SIR,t=seq(0,1000,length=1000))
  
  interv_SIR<-seq(0,max(r[,1])+24,by=24)
  pmediointerv_SIR<-(interv_SIR[-1]+interv_SIR[-length(interv_SIR)])/2
  clase<-cut(r[,1],interv_SIR,label=FALSE,include.lowest = TRUE)
  maxclas<-function(i){max(which(clase==i))}
  indices<-sapply(1:max(clase),FUN=maxclas)
  iemas0<-c(0,r[,5])
  newinfdiscr_SIR<-iemas0[c(1,indices+1)]
  newinfdiscr_SIR<-newinfdiscr_SIR[-1]-newinfdiscr_SIR[-length(newinfdiscr_SIR)]
  
  #setwd("d:\\Users\\Diana\\Documents\\Tesis CIMAT\\Cluster\\Poly\\EpidyODE")
  #pdf(paste("poly_SIR", j, "pdf", sep = "."))
  plot(datos$pmediointerv_SIR,datos$Freq,t="o",pch=19,ylab="Infectados",xlab="Tiempo")
  lines(pmediointerv_SIR,newinfdiscr_SIR*N,t="o",pch=19,col="blue",ylab="Infectados",xlab="Tiempo")
  legend("topright", c("M. epidemiológico","M. Determinista"),col= c("black", "blue"),lwd=2)
  #dev.off()
}


#####################################################################



sim_SIR_ws<-list(sim_SIR_ws_100_1_1, sim_SIR_ws_100_1_2, sim_SIR_ws_100_1_3, sim_SIR_ws_100_1_4, sim_SIR_ws_100_1_5, 
                 sim_SIR_ws_100_2_1, sim_SIR_ws_100_2_2, sim_SIR_ws_100_2_3, sim_SIR_ws_100_2_4, sim_SIR_ws_100_2_5,
                 sim_SIR_ws_100_3_1, sim_SIR_ws_100_3_2, sim_SIR_ws_100_3_3, sim_SIR_ws_100_3_4, sim_SIR_ws_100_3_5,
                 sim_SIR_ws_500_1_1, sim_SIR_ws_500_1_2, sim_SIR_ws_500_1_3, sim_SIR_ws_500_1_4, sim_SIR_ws_500_1_5, 
                 sim_SIR_ws_500_2_1, sim_SIR_ws_500_2_2, sim_SIR_ws_500_2_3, sim_SIR_ws_500_2_4, sim_SIR_ws_500_2_5,
                 sim_SIR_ws_500_3_1, sim_SIR_ws_500_3_2, sim_SIR_ws_500_3_3, sim_SIR_ws_500_3_4, sim_SIR_ws_500_3_5,
                 sim_SIR_ws_1000_1_1, sim_SIR_ws_1000_1_2, sim_SIR_ws_1000_1_3, sim_SIR_ws_1000_1_4, sim_SIR_ws_1000_1_5, 
                 sim_SIR_ws_1000_2_1, sim_SIR_ws_1000_2_2, sim_SIR_ws_1000_2_3, sim_SIR_ws_1000_2_4, sim_SIR_ws_1000_2_5,
                 sim_SIR_ws_1000_3_1, sim_SIR_ws_1000_3_2, sim_SIR_ws_1000_3_3, sim_SIR_ws_1000_3_4, sim_SIR_ws_1000_3_5,
                 sim_SIR_ws_10000_1_1, sim_SIR_ws_10000_1_2, sim_SIR_ws_10000_1_3, sim_SIR_ws_10000_1_4, sim_SIR_ws_10000_1_5, 
                 sim_SIR_ws_10000_2_1, sim_SIR_ws_10000_2_2, sim_SIR_ws_10000_2_3, sim_SIR_ws_10000_2_4, sim_SIR_ws_10000_2_5,
                 sim_SIR_ws_10000_3_1, sim_SIR_ws_10000_3_2, sim_SIR_ws_10000_3_3, sim_SIR_ws_10000_3_4, sim_SIR_ws_10000_3_5,
                 sim_SIR_ws_50000_1_1, sim_SIR_ws_50000_1_2, sim_SIR_ws_50000_1_3, sim_SIR_ws_50000_1_4, sim_SIR_ws_50000_1_5, 
                 sim_SIR_ws_50000_2_1, sim_SIR_ws_50000_2_2, sim_SIR_ws_50000_2_3, sim_SIR_ws_50000_2_4, sim_SIR_ws_50000_2_5,
                 sim_SIR_ws_50000_3_1, sim_SIR_ws_50000_3_2, sim_SIR_ws_50000_3_3, sim_SIR_ws_50000_3_4, sim_SIR_ws_50000_3_5)

for(j in 1:75){
  aa<-sim_SIR_ws[[j]]$indexes$InfTime
  aa<-sort(aa[!is.na(aa)])
  interv_SIR<-seq(0,floor(max(aa)+1)+24,by=24)
  pmediointerv_SIR<-(interv_SIR[-1]+interv_SIR[-length(interv_SIR)])/2
  newinfdiscr_SIR<-table(cut(c(aa,pmediointerv_SIR),interv_SIR,label=FALSE,include.lowest = TRUE))-1
  datos<-data.frame(pmediointerv_SIR,newinfdiscr_SIR)
  #plot(datos$pmediointerv_SIR,datos$Freq,t="o",pch=19,ylab="Infectados",xlab="Tiempo")
  
  #EJEMPLO SIR ODE:
  if (j<=15) N<-100
  else if (j<=30) N<-500
  else if (j<=45) N<-1000
  else if (j<=60) N<-10000
  else N<-50000
  
  theta_SIR<-c(bet=betareal,gam=gamareal)
  i=2/N
  X_ini_SIR=c(s = 1-i, i = i, r = 0, ie=i) 
  
  r<-X_theta_SIR(theta_SIR=theta_SIR,t=seq(0,1000,length=1000))
  
  interv_SIR<-seq(0,max(r[,1])+24,by=24)
  pmediointerv_SIR<-(interv_SIR[-1]+interv_SIR[-length(interv_SIR)])/2
  clase<-cut(r[,1],interv_SIR,label=FALSE,include.lowest = TRUE)
  maxclas<-function(i){max(which(clase==i))}
  indices<-sapply(1:max(clase),FUN=maxclas)
  iemas0<-c(0,r[,5])
  newinfdiscr_SIR<-iemas0[c(1,indices+1)]
  newinfdiscr_SIR<-newinfdiscr_SIR[-1]-newinfdiscr_SIR[-length(newinfdiscr_SIR)]
  
  #setwd("d:\\Users\\Diana\\Documents\\Tesis CIMAT\\Cluster\\WS\\EpidyODE")
  #pdf(paste("ws_SIR", j, "pdf", sep = "."))
  plot(datos$pmediointerv_SIR,datos$Freq,t="o",pch=19,ylab="Infectados",xlab="Tiempo")
  lines(pmediointerv_SIR,newinfdiscr_SIR*N,t="o",pch=19,col="green",ylab="Infectados",xlab="Tiempo")
  legend("topright", c("M. epidemiológico","M. Determinista"),col= c("black", "green"),lwd=2)
  #dev.off()
}

