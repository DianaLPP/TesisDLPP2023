setwd("C:\\Users\\Diana\\Dropbox\\Diana\\CodigoR_extra")
load("Brotes_SIR_Pois_Poly_WS.RData")
load("pois_5difnodos_15sim.RData")
load("poly_5difnodos_15sim.RData")
load("WS_5difnodos_15sim.RData")
load("Cualidades de red, extra.RData")

crear_datas_extra<-function(m){

      for (n in 1:5){
        if (n==1) nodos<-100
        else if (n==2) nodos<-500
        else if (n==3) nodos<-1000
        else if (n==4) nodos<-10000
        else nodos<-50000
        
        ######################################## Poisson #########################################3
        assign(paste("data",nodos,"_pois_SIR", sep=""),"a")
        assign(paste("datagte",nodos,"_pois_SIR", sep=""),"b")
        
        for (r in 1:3){
          if (r==1){
            assign(paste("red_pois_",nodos, sep = ""),c(paste("sim_SIR_pois_",nodos,"_1_1", sep = ""), paste("sim_SIR_pois_",nodos,"_1_2", sep = ""), paste("sim_SIR_pois_",nodos,"_1_3", sep = ""), paste("sim_SIR_pois_",nodos,"_1_4", sep = ""), paste("sim_SIR_pois_",nodos,"_1_5", sep = "")))
            assign("diam",paste("diam_pois_",nodos,"_1", sep = ""))
            assign("med",paste("med_pois_",nodos,"_1", sep = ""))
            assign("var",paste("var_pois_",nodos,"_1", sep = ""))
            assign("coef",paste("coef_pois_",nodos,"_1", sep = ""))
            assign("assor",paste("assor_pois_",nodos,"_1", sep = ""))
            assign("mdist",paste("mdist_pois_",nodos,"_1", sep = ""))
            assign("eigen",paste("eigen_pois_",nodos,"_1", sep = ""))
          } 
          
          if (r==2){
            assign(paste("red_pois_",nodos, sep = ""),c(paste("sim_SIR_pois_",nodos,"_2_1", sep = ""), paste("sim_SIR_pois_",nodos,"_2_2", sep = ""), paste("sim_SIR_pois_",nodos,"_2_3", sep = ""), paste("sim_SIR_pois_",nodos,"_2_4", sep = ""), paste("sim_SIR_pois_",nodos,"_2_5", sep = "")))
            assign("diam",paste("diam_pois_",nodos,"_2", sep = ""))
            assign("med",paste("med_pois_",nodos,"_2", sep = ""))
            assign("var",paste("var_pois_",nodos,"_2", sep = ""))
            assign("coef",paste("coef_pois_",nodos,"_2", sep = ""))
            assign("assor",paste("assor_pois_",nodos,"_2", sep = ""))
            assign("mdist",paste("mdist_pois_",nodos,"_2", sep = ""))
            assign("eigen",paste("eigen_pois_",nodos,"_2", sep = ""))
          }
          
          if (r==3){
            assign(paste("red_pois_",nodos, sep = ""),c(paste("sim_SIR_pois_",nodos,"_3_1", sep = ""), paste("sim_SIR_pois_",nodos,"_3_2", sep = ""), paste("sim_SIR_pois_",nodos,"_3_3", sep = ""), paste("sim_SIR_pois_",nodos,"_3_4", sep = ""), paste("sim_SIR_pois_",nodos,"_3_5", sep = "")))
            assign("diam",paste("diam_pois_",nodos,"_3", sep = ""))
            assign("med",paste("med_pois_",nodos,"_3", sep = ""))
            assign("var",paste("var_pois_",nodos,"_3", sep = ""))
            assign("coef",paste("coef_pois_",nodos,"_3", sep = ""))
            assign("assor",paste("assor_pois_",nodos,"_3", sep = ""))
            assign("mdist",paste("mdist_pois_",nodos,"_3", sep = ""))
            assign("eigen",paste("eigen_pois_",nodos,"_3", sep = ""))
          }
          
          for (q in 1:5){
            aa<-get(get(paste("red_pois_",nodos, sep = ""))[[q]])$indexes$InfTime
            aa<-sort(aa[!is.na(aa)])
            interv_SIR<-seq(0,floor(max(aa)+1)+1,by=1)
            pmediointerv_SIR<-(interv_SIR[-1]+interv_SIR[-length(interv_SIR)])/2
            newinfdiscr_SIR<-table(cut(c(aa,pmediointerv_SIR),interv_SIR,label=FALSE,include.lowest = TRUE))-1
            datos<-data.frame(pmediointerv_SIR,newinfdiscr_SIR)
            
            if (length(datos$Freq)>(m+10)){  
              
              y<-datos$Freq    #time_series <- rnorm(nodos)
              u<-length(y)
              tmp <- rep(NA, u)
              lagged_data <- data.frame(V1 = tmp)
              for(i in seq_len(m)) {     #seq_len(max_lag) secuencia de 1 hasta max_lag
                lagged_data[, i] <- c(rep(NA, i), head(y[-1], u - i))
              }
              lagged_data <- lagged_data[complete.cases(lagged_data),]
              
              assign(paste("data",nodos,"pSIR", r, q, sep=""),data.frame(yinf = y[(m+2):(u+1)], lagged_data, x1=get(diam),x2=get(med),x3=get(var),x4=get(coef),x5=get(assor),x6=get(mdist),x7=get(eigen)))
              
              filas<-nrow(get(paste("data",nodos,"pSIR", r, q, sep="")))
              #numero de suceptibles al incio de cada intervalo (de una hora)
              suc<-nodos-cumsum(newinfdiscr_SIR) 
              yp<-((get(paste("data",nodos,"pSIR", r, q, sep="")))[,1])/(suc[(m+1):length(suc)])
              W1<-suc[(m+1):length(suc)]
              
              assign(paste("datagte",nodos,"pSIR", r, q, sep=""),(data.frame(yinfp=yp,yinf = y[(m+2):(u+1)], lagged_data, x1=get(diam),x2=get(med),x3=get(var),x4=get(coef),x5=get(assor),x6=get(mdist),x7=get(eigen),red=rep("Pois",filas),sim=rep(r,filas),brote=rep(q,filas),long=(1:filas),W1=W1))[-filas,])
              
              if (get(paste("data",nodos,"_pois_SIR", sep=""))=="a") assign(paste("data",nodos,"_pois_SIR", sep=""),rbind(get(paste("data",nodos,"pSIR", r, q, sep=""))))
              else assign(paste("data",nodos,"_pois_SIR", sep=""),rbind(get(paste("data",nodos,"_pois_SIR", sep="")),get(paste("data",nodos,"pSIR", r, q, sep=""))))
              
              if (get(paste("datagte",nodos,"_pois_SIR", sep=""))=="b") assign(paste("datagte",nodos,"_pois_SIR", sep=""),rbind(get(paste("datagte",nodos,"pSIR", r, q, sep=""))))
              else assign(paste("datagte",nodos,"_pois_SIR", sep=""),rbind(get(paste("datagte",nodos,"_pois_SIR", sep="")),get(paste("datagte",nodos,"pSIR", r, q, sep=""))))
            } #if
          } #for i
        }  #for r
        ################################################################################ 
        ######################################## Polylog #########################################3
      
        assign(paste("data",nodos,"_poly_SIR", sep=""),"a")
        assign(paste("datagte",nodos,"_poly_SIR", sep=""),"b")
        
        for (r in 1:3){
          if (r==1){
            assign(paste("red_poly_",nodos, sep = ""),c(paste("sim_SIR_poly_",nodos,"_1_1", sep = ""), paste("sim_SIR_poly_",nodos,"_1_2", sep = ""), paste("sim_SIR_poly_",nodos,"_1_3", sep = ""), paste("sim_SIR_poly_",nodos,"_1_4", sep = ""), paste("sim_SIR_poly_",nodos,"_1_5", sep = "")))
            assign("diam",paste("diam_poly_",nodos,"_1", sep = ""))
            assign("med",paste("med_poly_",nodos,"_1", sep = ""))
            assign("var",paste("var_poly_",nodos,"_1", sep = ""))
            assign("coef",paste("coef_poly_",nodos,"_1", sep = ""))
            assign("assor",paste("assor_poly_",nodos,"_1", sep = ""))
            assign("mdist",paste("mdist_poly_",nodos,"_1", sep = ""))
            assign("eigen",paste("eigen_poly_",nodos,"_1", sep = ""))
          } 
          
          if (r==2){
            assign(paste("red_poly_",nodos, sep = ""),c(paste("sim_SIR_poly_",nodos,"_2_1", sep = ""), paste("sim_SIR_poly_",nodos,"_2_2", sep = ""), paste("sim_SIR_poly_",nodos,"_2_3", sep = ""), paste("sim_SIR_poly_",nodos,"_2_4", sep = ""), paste("sim_SIR_poly_",nodos,"_2_5", sep = "")))
            assign("diam",paste("diam_poly_",nodos,"_2", sep = ""))
            assign("med",paste("med_poly_",nodos,"_2", sep = ""))
            assign("var",paste("var_poly_",nodos,"_2", sep = ""))
            assign("coef",paste("coef_poly_",nodos,"_2", sep = ""))
            assign("assor",paste("assor_poly_",nodos,"_2", sep = ""))
            assign("mdist",paste("mdist_poly_",nodos,"_2", sep = ""))
            assign("eigen",paste("eigen_poly_",nodos,"_2", sep = ""))
          }
          
          if (r==3){
            assign(paste("red_poly_",nodos, sep = ""),c(paste("sim_SIR_poly_",nodos,"_3_1", sep = ""), paste("sim_SIR_poly_",nodos,"_3_2", sep = ""), paste("sim_SIR_poly_",nodos,"_3_3", sep = ""), paste("sim_SIR_poly_",nodos,"_3_4", sep = ""), paste("sim_SIR_poly_",nodos,"_3_5", sep = "")))
            assign("diam",paste("diam_poly_",nodos,"_3", sep = ""))
            assign("med",paste("med_poly_",nodos,"_3", sep = ""))
            assign("var",paste("var_poly_",nodos,"_3", sep = ""))
            assign("coef",paste("coef_poly_",nodos,"_3", sep = ""))
            assign("assor",paste("assor_poly_",nodos,"_3", sep = ""))
            assign("mdist",paste("mdist_poly_",nodos,"_3", sep = ""))
            assign("eigen",paste("eigen_poly_",nodos,"_3", sep = ""))
          }
          
          for (q in 1:5){
            aa<-get(get(paste("red_poly_",nodos, sep = ""))[[q]])$indexes$InfTime
            aa<-sort(aa[!is.na(aa)])
            interv_SIR<-seq(0,floor(max(aa)+1)+1,by=1)
            pmediointerv_SIR<-(interv_SIR[-1]+interv_SIR[-length(interv_SIR)])/2
            newinfdiscr_SIR<-table(cut(c(aa,pmediointerv_SIR),interv_SIR,label=FALSE,include.lowest = TRUE))-1
            datos<-data.frame(pmediointerv_SIR,newinfdiscr_SIR)
            
            if (length(datos$Freq)>(m+10)){  
              
              y<-datos$Freq    #time_series <- rnorm(nodos)
              u<-length(y)
              tmp <- rep(NA, u)
              lagged_data <- data.frame(V1 = tmp)
              for(i in seq_len(m)) {     #seq_len(max_lag) secuencia de 1 hasta max_lag
                lagged_data[, i] <- c(rep(NA, i), head(y[-1], u - i))
              }
              lagged_data <- lagged_data[complete.cases(lagged_data),]
              
              assign(paste("data",nodos,"pySIR", r, q, sep=""),data.frame(yinf = y[(m+2):(u+1)], lagged_data, x1=get(diam),x2=get(med),x3=get(var),x4=get(coef),x5=get(assor),x6=get(mdist),x7=get(eigen)))
              
              filas<-nrow(get(paste("data",nodos,"pySIR", r, q, sep="")))
              #numero de suceptibles al incio de cada intervalo (de una hora)
              suc<-nodos-cumsum(newinfdiscr_SIR) 
              yp<-((get(paste("data",nodos,"pySIR", r, q, sep="")))[,1])/(suc[(m+1):length(suc)])
              W1<-suc[(m+1):length(suc)]
              
              assign(paste("datagte",nodos,"pySIR", r, q, sep=""),(data.frame(yinfp=yp,yinf = y[(m+2):(u+1)], lagged_data, x1=get(diam),x2=get(med),x3=get(var),x4=get(coef),x5=get(assor),x6=get(mdist),x7=get(eigen),red=rep("Poly",filas),sim=rep(r,filas),brote=rep(q,filas),long=(1:filas),W1=W1))[-filas,])
              
              if (get(paste("data",nodos,"_poly_SIR", sep=""))=="a") assign(paste("data",nodos,"_poly_SIR", sep=""),rbind(get(paste("data",nodos,"pySIR", r, q, sep=""))))
              else assign(paste("data",nodos,"_poly_SIR", sep=""),rbind(get(paste("data",nodos,"_poly_SIR", sep="")),get(paste("data",nodos,"pySIR", r, q, sep=""))))
              
              if (get(paste("datagte",nodos,"_poly_SIR", sep=""))=="b") assign(paste("datagte",nodos,"_poly_SIR", sep=""),rbind(get(paste("datagte",nodos,"pySIR", r, q, sep=""))))
              else assign(paste("datagte",nodos,"_poly_SIR", sep=""),rbind(get(paste("datagte",nodos,"_poly_SIR", sep="")),get(paste("datagte",nodos,"pySIR", r, q, sep=""))))
            } #if  
          } #for i
        }  #for r
        ################################################################################ 
        ######################################## Watts Strogatz #########################################3
      
        assign(paste("data",nodos,"_ws_SIR", sep=""),"a")
        assign(paste("datagte",nodos,"_ws_SIR", sep=""),"b")
        
        for (r in 1:3){
          if (r==1){
            assign(paste("red_ws_",nodos, sep = ""),c(paste("sim_SIR_ws_",nodos,"_1_1", sep = ""), paste("sim_SIR_ws_",nodos,"_1_2", sep = ""), paste("sim_SIR_ws_",nodos,"_1_3", sep = ""), paste("sim_SIR_ws_",nodos,"_1_4", sep = ""), paste("sim_SIR_ws_",nodos,"_1_5", sep = "")))
            assign("diam",paste("diam_ws_",nodos,"_1", sep = ""))
            assign("med",paste("med_ws_",nodos,"_1", sep = ""))
            assign("var",paste("var_ws_",nodos,"_1", sep = ""))
            assign("coef",paste("coef_ws_",nodos,"_1", sep = ""))
            assign("assor",paste("assor_ws_",nodos,"_1", sep = ""))
            assign("mdist",paste("mdist_ws_",nodos,"_1", sep = ""))
            assign("eigen",paste("eigen_ws_",nodos,"_1", sep = ""))
          } 
          
          if (r==2){
            assign(paste("red_ws_",nodos, sep = ""),c(paste("sim_SIR_ws_",nodos,"_2_1", sep = ""), paste("sim_SIR_ws_",nodos,"_2_2", sep = ""), paste("sim_SIR_ws_",nodos,"_2_3", sep = ""), paste("sim_SIR_ws_",nodos,"_2_4", sep = ""), paste("sim_SIR_ws_",nodos,"_2_5", sep = "")))
            assign("diam",paste("diam_ws_",nodos,"_2", sep = ""))
            assign("med",paste("med_ws_",nodos,"_2", sep = ""))
            assign("var",paste("var_ws_",nodos,"_2", sep = ""))
            assign("coef",paste("coef_ws_",nodos,"_2", sep = ""))
            assign("assor",paste("assor_ws_",nodos,"_2", sep = ""))
            assign("mdist",paste("mdist_ws_",nodos,"_2", sep = ""))
            assign("eigen",paste("eigen_ws_",nodos,"_2", sep = ""))
          }
          
          if (r==3){
            assign(paste("red_ws_",nodos, sep = ""),c(paste("sim_SIR_ws_",nodos,"_3_1", sep = ""), paste("sim_SIR_ws_",nodos,"_3_2", sep = ""), paste("sim_SIR_ws_",nodos,"_3_3", sep = ""), paste("sim_SIR_ws_",nodos,"_3_4", sep = ""), paste("sim_SIR_ws_",nodos,"_3_5", sep = "")))
            assign("diam",paste("diam_ws_",nodos,"_3", sep = ""))
            assign("med",paste("med_ws_",nodos,"_3", sep = ""))
            assign("var",paste("var_ws_",nodos,"_3", sep = ""))
            assign("coef",paste("coef_ws_",nodos,"_3", sep = ""))
            assign("assor",paste("assor_ws_",nodos,"_3", sep = ""))
            assign("mdist",paste("mdist_ws_",nodos,"_3", sep = ""))
            assign("eigen",paste("eigen_ws_",nodos,"_3", sep = ""))
          }
          
          for (q in 1:5){
            aa<-get(get(paste("red_ws_",nodos, sep = ""))[[q]])$indexes$InfTime
            aa<-sort(aa[!is.na(aa)])
            interv_SIR<-seq(0,floor(max(aa)+1)+1,by=1)
            pmediointerv_SIR<-(interv_SIR[-1]+interv_SIR[-length(interv_SIR)])/2
            newinfdiscr_SIR<-table(cut(c(aa,pmediointerv_SIR),interv_SIR,label=FALSE,include.lowest = TRUE))-1
            datos<-data.frame(pmediointerv_SIR,newinfdiscr_SIR)
            
            if (length(datos$Freq)>(m+10)){  
              
              y<-datos$Freq    #time_series <- rnorm(nodos)
              u<-length(y)
              tmp <- rep(NA, u)
              lagged_data <- data.frame(V1 = tmp)
              for(i in seq_len(m)) {     #seq_len(max_lag) secuencia de 1 hasta max_lag
                lagged_data[, i] <- c(rep(NA, i), head(y[-1], u - i))
              }
              lagged_data <- lagged_data[complete.cases(lagged_data),]
              
              assign(paste("data",nodos,"wsSIR", r, q, sep=""),data.frame(yinf = y[(m+2):(u+1)], lagged_data, x1=get(diam),x2=get(med),x3=get(var),x4=get(coef),x5=get(assor),x6=get(mdist),x7=get(eigen)))
              
              filas<-nrow(get(paste("data",nodos,"wsSIR", r, q, sep="")))
              #numero de suceptibles al incio de cada intervalo (de una hora)
              suc<-nodos-cumsum(newinfdiscr_SIR) 
              yp<-((get(paste("data",nodos,"wsSIR", r, q, sep="")))[,1])/(suc[(m+1):length(suc)])
              W1<-suc[(m+1):length(suc)]
              
              assign(paste("datagte",nodos,"wsSIR", r, q, sep=""),(data.frame(yinfp=yp,yinf = y[(m+2):(u+1)], lagged_data, x1=get(diam),x2=get(med),x3=get(var),x4=get(coef),x5=get(assor),x6=get(mdist),x7=get(eigen),red=rep("WS",filas),sim=rep(r,filas),brote=rep(q,filas),long=(1:filas),W1=W1))[-filas,])
              
              if (get(paste("data",nodos,"_ws_SIR", sep=""))=="a") assign(paste("data",nodos,"_ws_SIR", sep=""),rbind(get(paste("data",nodos,"wsSIR", r, q, sep=""))))
              else assign(paste("data",nodos,"_ws_SIR", sep=""),rbind(get(paste("data",nodos,"_ws_SIR", sep="")),get(paste("data",nodos,"wsSIR", r, q, sep=""))))
              
              if (get(paste("datagte",nodos,"_ws_SIR", sep=""))=="b") assign(paste("datagte",nodos,"_ws_SIR", sep=""),rbind(get(paste("datagte",nodos,"wsSIR", r, q, sep=""))))
              else assign(paste("datagte",nodos,"_ws_SIR", sep=""),rbind(get(paste("datagte",nodos,"_ws_SIR", sep="")),get(paste("datagte",nodos,"wsSIR", r, q, sep=""))))
            } #if
          } #for i
        }  #for r
        ################################################################################ 
      } #for n
      
      data_SIR_100<-rbind(data100_pois_SIR,data100_poly_SIR,data100_ws_SIR)
      data_SIR_500<-rbind(data500_pois_SIR,data500_poly_SIR,data500_ws_SIR)
      data_SIR_1000<-rbind(data1000_pois_SIR,data1000_poly_SIR,data1000_ws_SIR)
      data_SIR_10000<-rbind(data10000_pois_SIR,data10000_poly_SIR,data10000_ws_SIR)
      data_SIR_50000<-rbind(data50000_pois_SIR,data50000_poly_SIR,data50000_ws_SIR)
      
      datagte_SIR_100<-rbind(datagte100_pois_SIR,datagte100_poly_SIR,datagte100_ws_SIR)
      datagte_SIR_500<-rbind(datagte500_pois_SIR,datagte500_poly_SIR,datagte500_ws_SIR)
      datagte_SIR_1000<-rbind(datagte1000_pois_SIR,datagte1000_poly_SIR,datagte1000_ws_SIR)
      datagte_SIR_10000<-rbind(datagte10000_pois_SIR,datagte10000_poly_SIR,datagte10000_ws_SIR)
      datagte_SIR_50000<-rbind(datagte50000_pois_SIR,datagte50000_poly_SIR,datagte50000_ws_SIR)
      
      save(datagte100_pois_SIR,datagte100_poly_SIR,datagte100_ws_SIR,
           datagte100pSIR11,datagte100pSIR12,datagte100pSIR13,datagte100pSIR14,
           datagte100pSIR21,datagte100pSIR22,datagte100pSIR23,datagte100pSIR24,datagte100pSIR25,
           datagte100pSIR31,datagte100pSIR32,datagte100pSIR33,datagte100pSIR34,datagte100pSIR35,
           datagte100pySIR12,datagte100pySIR14,datagte100pySIR15,
           datagte100pySIR23,datagte100pySIR24,datagte100pySIR25,
           datagte100pySIR31,datagte100pySIR34,
           datagte100wsSIR12,datagte100wsSIR13,datagte100wsSIR14,datagte100wsSIR15,
           datagte100wsSIR21,datagte100wsSIR22,datagte100wsSIR23,datagte100wsSIR24,datagte100wsSIR25,
           datagte100wsSIR31,datagte100wsSIR33,datagte100wsSIR35,
           datagte500_pois_SIR,datagte500_poly_SIR,datagte500_ws_SIR,
           datagte500pSIR11,datagte500pSIR12,datagte500pSIR13,datagte500pSIR15,
           datagte500pSIR21,datagte500pSIR22,datagte500pSIR23,datagte500pSIR24,datagte500pSIR25,
           datagte500pSIR31,datagte500pSIR32,datagte500pSIR33,datagte500pSIR34,datagte500pSIR35, 
           datagte500pySIR12,datagte500pySIR13,datagte500pySIR14,datagte500pySIR15,
           datagte500pySIR21,datagte500pySIR22,datagte500pySIR23,datagte500pySIR24,datagte500pySIR25,
           datagte500pySIR31,datagte500pySIR32,datagte500pySIR33,datagte500pySIR34,datagte500pySIR35, 
           datagte500wsSIR11,datagte500wsSIR12,datagte500wsSIR13,datagte500wsSIR14,datagte500wsSIR15,
           datagte500wsSIR21,datagte500wsSIR22,datagte500wsSIR23,datagte500wsSIR24,datagte500wsSIR25,
           datagte500wsSIR31,datagte500wsSIR32,datagte500wsSIR33,datagte500wsSIR34,datagte500wsSIR35, 
           datagte1000_pois_SIR,datagte1000_poly_SIR,datagte1000_ws_SIR,
           datagte1000pSIR11,datagte1000pSIR12,datagte1000pSIR15,
           datagte1000pSIR21,datagte1000pSIR22,datagte1000pSIR23,datagte1000pSIR24,datagte1000pSIR25,
           datagte1000pSIR31,datagte1000pSIR33,datagte1000pSIR34,datagte1000pSIR35,  
           datagte1000pySIR11,datagte1000pySIR12,datagte1000pySIR13,datagte1000pySIR14,datagte1000pySIR15,
           datagte1000pySIR21,datagte1000pySIR22,datagte1000pySIR23,datagte1000pySIR24,datagte1000pySIR25,
           datagte1000pySIR31,datagte1000pySIR32,datagte1000pySIR33,datagte1000pySIR34,datagte1000pySIR35,  
           datagte1000wsSIR11,datagte1000wsSIR12,datagte1000wsSIR13,datagte1000wsSIR14,datagte1000wsSIR15,
           datagte1000wsSIR21,datagte1000wsSIR22,datagte1000wsSIR23,datagte1000wsSIR24,datagte1000wsSIR25,
           datagte1000wsSIR31,datagte1000wsSIR32,datagte1000wsSIR33,datagte1000wsSIR34,datagte1000wsSIR35, 
           datagte10000_pois_SIR,datagte10000_poly_SIR,datagte10000_ws_SIR,
           datagte10000pSIR12,datagte10000pSIR13,datagte10000pSIR14,datagte10000pSIR15,
           datagte10000pSIR21,datagte10000pSIR22,datagte10000pSIR23,datagte10000pSIR24,datagte10000pSIR25,
           datagte10000pSIR31,datagte10000pSIR32,datagte10000pSIR33,datagte10000pSIR34,datagte10000pSIR35,  
           datagte10000pySIR11,datagte10000pySIR12,datagte10000pySIR13,datagte10000pySIR14,datagte10000pySIR15,
           datagte10000pySIR21,datagte10000pySIR22,datagte10000pySIR23,datagte10000pySIR24,datagte10000pySIR25,
           datagte10000pySIR31,datagte10000pySIR32,datagte10000pySIR33,datagte10000pySIR34,datagte10000pySIR35,  
           datagte10000wsSIR11,datagte10000wsSIR12,datagte10000wsSIR13,datagte10000wsSIR14,datagte10000wsSIR15,
           datagte10000wsSIR21,datagte10000wsSIR22,datagte10000wsSIR23,datagte10000wsSIR24,datagte10000wsSIR25,
           datagte10000wsSIR31,datagte10000wsSIR32,datagte10000wsSIR33,datagte10000wsSIR34,datagte10000wsSIR35,  
           datagte50000_pois_SIR,datagte50000_poly_SIR,datagte50000_ws_SIR,
           datagte50000pSIR11,datagte50000pSIR12,datagte50000pSIR13,datagte50000pSIR14,datagte50000pSIR15,
           datagte50000pSIR21,datagte50000pSIR22,datagte50000pSIR23,datagte50000pSIR24,
           datagte50000pSIR31,datagte50000pSIR32,datagte50000pSIR33,datagte50000pSIR34,datagte50000pSIR35, 
           datagte50000pySIR11,datagte50000pySIR12,datagte50000pySIR14,
           datagte50000pySIR21,datagte50000pySIR22,datagte50000pySIR23,datagte50000pySIR25,
           datagte50000pySIR31,datagte50000pySIR32,datagte50000pySIR33,datagte50000pySIR34,datagte50000pySIR35, 
           datagte50000wsSIR11,datagte50000wsSIR12,datagte50000wsSIR13,datagte50000wsSIR14,datagte50000wsSIR15,
           datagte50000wsSIR21,datagte50000wsSIR22,datagte50000wsSIR23,datagte50000wsSIR24,datagte50000wsSIR25,
           datagte50000wsSIR31,datagte50000wsSIR33,datagte50000wsSIR34,datagte50000wsSIR35, 
           datagte_SIR_100,datagte_SIR_500,datagte_SIR_1000,datagte_SIR_10000,datagte_SIR_50000,
           file=paste("Todas_las_datagtes_SIR_","100_extra",".RData",sep = ""))
      
      
} #fin de la función
 
crear_datas_extra(100)

load("Todas_las_datagtes_SIR_40_extra.RData")



setwd("C:\\Users\\Diana\\Dropbox\\Diana\\CodigoR")
load("Todas_las_datagtes_SIR_40.RData")
setwd("C:\\Users\\Diana\\Dropbox\\Diana\\CodigoR_extra")
load("datas_extras.RData")




all(datagte_SIR_100_extra[,-c(m+8,m+9)]==datagte_SIR_100) #comprueba si efectivamente son iguales 

data_mod<-datagte_SIR_50000
data_mod[,1]<-data_mod[,2]/data_mod[,m+12]  #y infp tiene la %
ajuste1p<-glm(yinfp~.,data = data_mod[,-c(1,(m+8),(m+9),(m+10),(m+11))],family = binomial (link = "logit"))
ajuste1<-glm(yinf~.,data = data_mod[,-c(1,(m+8),(m+9),(m+10),(m+11))],family = poisson (link = "log"))



data_mod<-datagte_SIR_50000_extra
data_mod[,1]<-data_mod[,2]/data_mod[,m+14]  #y infp tiene la %
ajuste2p<-glm(yinfp~.,data = data_mod[,-c(1,(m+8),(m+9),(m+10),(m+11),(m+12),(m+13))],family = binomial (link = "logit"))
ajuste2<-glm(yinf~.,data = data_mod[,-c(1,(m+8),(m+9),(m+10),(m+11),(m+12),(m+13))],family = poisson (link = "log"))






save(datagte500_pois_SIR,datagte500_poly_SIR,datagte500_ws_SIR,
     datagte500pSIR11,datagte500pSIR12,datagte500pSIR13,datagte500pSIR15,
     datagte500pSIR21,datagte500pSIR22,datagte500pSIR23,datagte500pSIR24,datagte500pSIR25,
     datagte500pSIR31,datagte500pSIR32,datagte500pSIR33,datagte500pSIR34,datagte500pSIR35, 
     datagte500pySIR12,datagte500pySIR13,datagte500pySIR14,datagte500pySIR15,
     datagte500pySIR21,datagte500pySIR22,datagte500pySIR23,datagte500pySIR24,datagte500pySIR25,
     datagte500pySIR31,datagte500pySIR32,datagte500pySIR33,datagte500pySIR34,datagte500pySIR35, 
     datagte500wsSIR11,datagte500wsSIR12,datagte500wsSIR13,datagte500wsSIR14,datagte500wsSIR15,
     datagte500wsSIR21,datagte500wsSIR22,datagte500wsSIR23,datagte500wsSIR24,datagte500wsSIR25,
     datagte500wsSIR31,datagte500wsSIR32,datagte500wsSIR33,datagte500wsSIR34,datagte500wsSIR35, 
     datagte1000_pois_SIR,datagte1000_poly_SIR,datagte1000_ws_SIR,
     datagte1000pSIR11,datagte1000pSIR12,datagte1000pSIR14,datagte1000pSIR15,
     datagte1000pSIR21,datagte1000pSIR22,datagte1000pSIR23,datagte1000pSIR24,datagte1000pSIR25,
     datagte1000pSIR31,datagte1000pSIR33,datagte1000pSIR34,datagte1000pSIR35,  
     datagte1000pySIR11,datagte1000pySIR12,datagte1000pySIR13,datagte1000pySIR14,datagte1000pySIR15,
     datagte1000pySIR21,datagte1000pySIR22,datagte1000pySIR23,datagte1000pySIR24,datagte1000pySIR25,
     datagte1000pySIR31,datagte1000pySIR32,datagte1000pySIR33,datagte1000pySIR34,datagte1000pySIR35,  
     datagte1000wsSIR11,datagte1000wsSIR12,datagte1000wsSIR13,datagte1000wsSIR14,datagte1000wsSIR15,
     datagte1000wsSIR21,datagte1000wsSIR22,datagte1000wsSIR23,datagte1000wsSIR24,datagte1000wsSIR25,
     datagte1000wsSIR31,datagte1000wsSIR32,datagte1000wsSIR33,datagte1000wsSIR34,datagte1000wsSIR35, 
     datagte10000_pois_SIR,datagte10000_poly_SIR,datagte10000_ws_SIR,
     datagte10000pSIR12,datagte10000pSIR13,datagte10000pSIR14,datagte10000pSIR15,
     datagte10000pSIR21,datagte10000pSIR22,datagte10000pSIR23,datagte10000pSIR24,datagte10000pSIR25,
     datagte10000pSIR31,datagte10000pSIR32,datagte10000pSIR33,datagte10000pSIR34,datagte10000pSIR35,  
     datagte10000pySIR11,datagte10000pySIR12,datagte10000pySIR13,datagte10000pySIR14,datagte10000pySIR15,
     datagte10000pySIR21,datagte10000pySIR22,datagte10000pySIR23,datagte10000pySIR24,datagte10000pySIR25,
     datagte10000pySIR31,datagte10000pySIR32,datagte10000pySIR33,datagte10000pySIR34,datagte10000pySIR35,  
     datagte10000wsSIR11,datagte10000wsSIR12,datagte10000wsSIR13,datagte10000wsSIR14,datagte10000wsSIR15,
     datagte10000wsSIR21,datagte10000wsSIR22,datagte10000wsSIR23,datagte10000wsSIR24,datagte10000wsSIR25,
     datagte10000wsSIR31,datagte10000wsSIR32,datagte10000wsSIR33,datagte10000wsSIR34,datagte10000wsSIR35,  
     datagte50000_pois_SIR,datagte50000_poly_SIR,datagte50000_ws_SIR,
     datagte50000pSIR11,datagte50000pSIR12,datagte50000pSIR13,datagte50000pSIR14,datagte50000pSIR15,
     datagte50000pSIR21,datagte50000pSIR22,datagte50000pSIR23,datagte50000pSIR24,datagte50000pSIR25,
     datagte50000pSIR31,datagte50000pSIR32,datagte50000pSIR33,datagte50000pSIR34,datagte50000pSIR35, 
     datagte50000pySIR11,datagte50000pySIR12,datagte50000pySIR14,
     datagte50000pySIR21,datagte50000pySIR22,datagte50000pySIR23,datagte50000pySIR25,
     datagte50000pySIR31,datagte50000pySIR32,datagte50000pySIR33,datagte50000pySIR34,datagte50000pySIR35, 
     datagte50000wsSIR11,datagte50000wsSIR12,datagte50000wsSIR13,datagte50000wsSIR14,datagte50000wsSIR15,
     datagte50000wsSIR21,datagte50000wsSIR22,datagte50000wsSIR23,datagte50000wsSIR24,datagte50000wsSIR25,
     datagte50000wsSIR31,datagte50000wsSIR33,datagte50000wsSIR34,datagte50000wsSIR35, 
     datagte_SIR_100,datagte_SIR_500,datagte_SIR_1000,datagte_SIR_10000,datagte_SIR_50000,
     file=paste("Menos100_las_datagtes_SIR_",m,".RData",sep = ""))


save(datagte100_pois_SIR,datagte100_poly_SIR,datagte100_ws_SIR,
     datagte100pSIR11,datagte100pSIR12,datagte100pSIR13,datagte100pSIR14,
     datagte100pSIR21,datagte100pSIR22,datagte100pSIR23,datagte100pSIR24,datagte100pSIR25,
     datagte100pSIR31,datagte100pSIR32,datagte100pSIR33,datagte100pSIR34,datagte100pSIR35,
     datagte100pySIR12,datagte100pySIR14,datagte100pySIR15,
     datagte100pySIR23,datagte100pySIR24,datagte100pySIR25,
     datagte100pySIR31,datagte100pySIR33,datagte100pySIR34,
     datagte100wsSIR11,datagte100wsSIR12,datagte100wsSIR13,datagte100wsSIR14,datagte100wsSIR15,
     datagte100wsSIR21,datagte100wsSIR22,datagte100wsSIR23,datagte100wsSIR24,datagte100wsSIR25,
     datagte100wsSIR31,datagte100wsSIR32,datagte100wsSIR33,datagte100wsSIR34,datagte100wsSIR35,
     datagte500_pois_SIR,datagte500_poly_SIR,datagte500_ws_SIR,
     datagte500pSIR11,datagte500pSIR12,datagte500pSIR13,datagte500pSIR15,
     datagte500pSIR21,datagte500pSIR22,datagte500pSIR23,datagte500pSIR24,datagte500pSIR25,
     datagte500pSIR31,datagte500pSIR32,datagte500pSIR33,datagte500pSIR34,datagte500pSIR35, 
     datagte500pySIR12,datagte500pySIR13,datagte500pySIR14,datagte500pySIR15,
     datagte500pySIR21,datagte500pySIR22,datagte500pySIR23,datagte500pySIR24,datagte500pySIR25,
     datagte500pySIR31,datagte500pySIR32,datagte500pySIR33,datagte500pySIR34,datagte500pySIR35, 
     datagte500wsSIR11,datagte500wsSIR12,datagte500wsSIR13,datagte500wsSIR14,datagte500wsSIR15,
     datagte500wsSIR21,datagte500wsSIR22,datagte500wsSIR23,datagte500wsSIR24,datagte500wsSIR25,
     datagte500wsSIR31,datagte500wsSIR32,datagte500wsSIR33,datagte500wsSIR34,datagte500wsSIR35, 
     datagte1000_pois_SIR,datagte1000_poly_SIR,datagte1000_ws_SIR,
     datagte1000pSIR11,datagte1000pSIR12,datagte1000pSIR14,datagte1000pSIR15,
     datagte1000pSIR21,datagte1000pSIR22,datagte1000pSIR23,datagte1000pSIR24,datagte1000pSIR25,
     datagte1000pSIR31,datagte1000pSIR33,datagte1000pSIR34,datagte1000pSIR35,  
     datagte1000pySIR11,datagte1000pySIR12,datagte1000pySIR13,datagte1000pySIR14,datagte1000pySIR15,
     datagte1000pySIR21,datagte1000pySIR22,datagte1000pySIR23,datagte1000pySIR24,datagte1000pySIR25,
     datagte1000pySIR31,datagte1000pySIR32,datagte1000pySIR33,datagte1000pySIR34,datagte1000pySIR35,  
     datagte1000wsSIR11,datagte1000wsSIR12,datagte1000wsSIR13,datagte1000wsSIR14,datagte1000wsSIR15,
     datagte1000wsSIR21,datagte1000wsSIR22,datagte1000wsSIR23,datagte1000wsSIR24,datagte1000wsSIR25,
     datagte1000wsSIR31,datagte1000wsSIR32,datagte1000wsSIR33,datagte1000wsSIR34,datagte1000wsSIR35, 
     datagte10000_pois_SIR,datagte10000_poly_SIR,datagte10000_ws_SIR,
     datagte10000pSIR12,datagte10000pSIR13,datagte10000pSIR14,datagte10000pSIR15,
     datagte10000pSIR21,datagte10000pSIR22,datagte10000pSIR23,datagte10000pSIR24,datagte10000pSIR25,
     datagte10000pSIR31,datagte10000pSIR32,datagte10000pSIR33,datagte10000pSIR34,datagte10000pSIR35,  
     datagte10000pySIR11,datagte10000pySIR12,datagte10000pySIR13,datagte10000pySIR14,datagte10000pySIR15,
     datagte10000pySIR21,datagte10000pySIR22,datagte10000pySIR23,datagte10000pySIR24,datagte10000pySIR25,
     datagte10000pySIR31,datagte10000pySIR32,datagte10000pySIR33,datagte10000pySIR34,datagte10000pySIR35,  
     datagte10000wsSIR11,datagte10000wsSIR12,datagte10000wsSIR13,datagte10000wsSIR14,datagte10000wsSIR15,
     datagte10000wsSIR21,datagte10000wsSIR22,datagte10000wsSIR23,datagte10000wsSIR24,datagte10000wsSIR25,
     datagte10000wsSIR31,datagte10000wsSIR32,datagte10000wsSIR33,datagte10000wsSIR34,datagte10000wsSIR35,  
     datagte50000_pois_SIR,datagte50000_poly_SIR,datagte50000_ws_SIR,
     datagte50000pSIR11,datagte50000pSIR12,datagte50000pSIR13,datagte50000pSIR14,datagte50000pSIR15,
     datagte50000pSIR21,datagte50000pSIR22,datagte50000pSIR23,datagte50000pSIR24,datagte50000pSIR25,
     datagte50000pSIR31,datagte50000pSIR32,datagte50000pSIR33,datagte50000pSIR34,datagte50000pSIR35, 
     datagte50000pySIR11,datagte50000pySIR12,datagte50000pySIR14,
     datagte50000pySIR21,datagte50000pySIR22,datagte50000pySIR23,datagte50000pySIR25,
     datagte50000pySIR31,datagte50000pySIR32,datagte50000pySIR33,datagte50000pySIR34,datagte50000pySIR35, 
     datagte50000wsSIR11,datagte50000wsSIR12,datagte50000wsSIR13,datagte50000wsSIR14,datagte50000wsSIR15,
     datagte50000wsSIR21,datagte50000wsSIR22,datagte50000wsSIR23,datagte50000wsSIR24,datagte50000wsSIR25,
     datagte50000wsSIR31,datagte50000wsSIR33,datagte50000wsSIR34,datagte50000wsSIR35, 
     datagte_SIR_100,datagte_SIR_500,datagte_SIR_1000,datagte_SIR_10000,datagte_SIR_50000,
     file="Todas_las_datagtes_SIR_20_extra.RData")

save(datagte100_pois_SIR,datagte100_poly_SIR,datagte100_ws_SIR,
     datagte100pSIR11,datagte100pSIR12,datagte100pSIR13,datagte100pSIR14,
     datagte100pSIR21,datagte100pSIR22,datagte100pSIR23,datagte100pSIR24,datagte100pSIR25,
     datagte100pSIR31,datagte100pSIR32,datagte100pSIR33,datagte100pSIR34,datagte100pSIR35,
     datagte100pySIR12,datagte100pySIR14,datagte100pySIR15,
     datagte100pySIR23,datagte100pySIR24,datagte100pySIR25,
     datagte100pySIR31,datagte100pySIR33,datagte100pySIR34,
     datagte100wsSIR11,datagte100wsSIR12,datagte100wsSIR13,datagte100wsSIR14,datagte100wsSIR15,
     datagte100wsSIR21,datagte100wsSIR22,datagte100wsSIR23,datagte100wsSIR24,datagte100wsSIR25,
     datagte100wsSIR31,datagte100wsSIR32,datagte100wsSIR33,datagte100wsSIR34,datagte100wsSIR35,
     datagte500_pois_SIR,datagte500_poly_SIR,datagte500_ws_SIR,
     datagte500pSIR11,datagte500pSIR12,datagte500pSIR13,datagte500pSIR15,
     datagte500pSIR21,datagte500pSIR22,datagte500pSIR23,datagte500pSIR24,datagte500pSIR25,
     datagte500pSIR31,datagte500pSIR32,datagte500pSIR33,datagte500pSIR34,datagte500pSIR35, 
     datagte500pySIR12,datagte500pySIR13,datagte500pySIR14,datagte500pySIR15,
     datagte500pySIR21,datagte500pySIR22,datagte500pySIR23,datagte500pySIR24,datagte500pySIR25,
     datagte500pySIR31,datagte500pySIR32,datagte500pySIR33,datagte500pySIR34,datagte500pySIR35, 
     datagte500wsSIR11,datagte500wsSIR12,datagte500wsSIR13,datagte500wsSIR14,datagte500wsSIR15,
     datagte500wsSIR21,datagte500wsSIR22,datagte500wsSIR23,datagte500wsSIR24,datagte500wsSIR25,
     datagte500wsSIR31,datagte500wsSIR32,datagte500wsSIR33,datagte500wsSIR34,datagte500wsSIR35, 
     datagte1000_pois_SIR,datagte1000_poly_SIR,datagte1000_ws_SIR,
     datagte1000pSIR11,datagte1000pSIR12,datagte1000pSIR14,datagte1000pSIR15,
     datagte1000pSIR21,datagte1000pSIR22,datagte1000pSIR23,datagte1000pSIR24,datagte1000pSIR25,
     datagte1000pSIR31,datagte1000pSIR33,datagte1000pSIR34,datagte1000pSIR35,  
     datagte1000pySIR11,datagte1000pySIR12,datagte1000pySIR13,datagte1000pySIR14,datagte1000pySIR15,
     datagte1000pySIR21,datagte1000pySIR22,datagte1000pySIR23,datagte1000pySIR24,datagte1000pySIR25,
     datagte1000pySIR31,datagte1000pySIR32,datagte1000pySIR33,datagte1000pySIR34,datagte1000pySIR35,  
     datagte1000wsSIR11,datagte1000wsSIR12,datagte1000wsSIR13,datagte1000wsSIR14,datagte1000wsSIR15,
     datagte1000wsSIR21,datagte1000wsSIR22,datagte1000wsSIR23,datagte1000wsSIR24,datagte1000wsSIR25,
     datagte1000wsSIR31,datagte1000wsSIR32,datagte1000wsSIR33,datagte1000wsSIR34,datagte1000wsSIR35, 
     datagte10000_pois_SIR,datagte10000_poly_SIR,datagte10000_ws_SIR,
     datagte10000pSIR12,datagte10000pSIR13,datagte10000pSIR14,datagte10000pSIR15,
     datagte10000pSIR21,datagte10000pSIR22,datagte10000pSIR23,datagte10000pSIR24,datagte10000pSIR25,
     datagte10000pSIR31,datagte10000pSIR32,datagte10000pSIR33,datagte10000pSIR34,datagte10000pSIR35,  
     datagte10000pySIR11,datagte10000pySIR12,datagte10000pySIR13,datagte10000pySIR14,datagte10000pySIR15,
     datagte10000pySIR21,datagte10000pySIR22,datagte10000pySIR23,datagte10000pySIR24,datagte10000pySIR25,
     datagte10000pySIR31,datagte10000pySIR32,datagte10000pySIR33,datagte10000pySIR34,datagte10000pySIR35,  
     datagte10000wsSIR11,datagte10000wsSIR12,datagte10000wsSIR13,datagte10000wsSIR14,datagte10000wsSIR15,
     datagte10000wsSIR21,datagte10000wsSIR22,datagte10000wsSIR23,datagte10000wsSIR24,datagte10000wsSIR25,
     datagte10000wsSIR31,datagte10000wsSIR32,datagte10000wsSIR33,datagte10000wsSIR34,datagte10000wsSIR35,  
     datagte50000_pois_SIR,datagte50000_poly_SIR,datagte50000_ws_SIR,
     datagte50000pSIR11,datagte50000pSIR12,datagte50000pSIR13,datagte50000pSIR14,datagte50000pSIR15,
     datagte50000pSIR21,datagte50000pSIR22,datagte50000pSIR23,datagte50000pSIR24,datagte50000pSIR25,
     datagte50000pSIR31,datagte50000pSIR32,datagte50000pSIR33,datagte50000pSIR34,datagte50000pSIR35, 
     datagte50000pySIR11,datagte50000pySIR12,datagte50000pySIR14,
     datagte50000pySIR21,datagte50000pySIR22,datagte50000pySIR23,datagte50000pySIR25,
     datagte50000pySIR31,datagte50000pySIR32,datagte50000pySIR33,datagte50000pySIR34,datagte50000pySIR35, 
     datagte50000wsSIR11,datagte50000wsSIR12,datagte50000wsSIR13,datagte50000wsSIR14,datagte50000wsSIR15,
     datagte50000wsSIR21,datagte50000wsSIR22,datagte50000wsSIR23,datagte50000wsSIR24,datagte50000wsSIR25,
     datagte50000wsSIR31,datagte50000wsSIR33,datagte50000wsSIR34,datagte50000wsSIR35, 
     datagte_SIR_100,datagte_SIR_500,datagte_SIR_1000,datagte_SIR_10000,datagte_SIR_50000,
     file=paste("Todas_las_datagtes_SIR_","40_extra",".RData",sep = ""))



save(datagte100_pois_SIR,datagte100_poly_SIR,datagte100_ws_SIR,
     datagte100pSIR11,datagte100pSIR12,datagte100pSIR13,datagte100pSIR14,
     datagte100pSIR21,datagte100pSIR22,datagte100pSIR23,datagte100pSIR24,datagte100pSIR25,
     datagte100pSIR31,datagte100pSIR32,datagte100pSIR33,datagte100pSIR34,datagte100pSIR35,
     datagte100pySIR12,datagte100pySIR14,datagte100pySIR15,
     datagte100pySIR23,datagte100pySIR24,datagte100pySIR25,
     datagte100pySIR31,datagte100pySIR34,
     datagte100wsSIR12,datagte100wsSIR13,datagte100wsSIR14,datagte100wsSIR15,
     datagte100wsSIR21,datagte100wsSIR22,datagte100wsSIR23,datagte100wsSIR24,datagte100wsSIR25,
     datagte100wsSIR31,datagte100wsSIR33,datagte100wsSIR35,
     datagte500_pois_SIR,datagte500_poly_SIR,datagte500_ws_SIR,
     datagte500pSIR11,datagte500pSIR12,datagte500pSIR13,datagte500pSIR15,
     datagte500pSIR21,datagte500pSIR22,datagte500pSIR23,datagte500pSIR24,datagte500pSIR25,
     datagte500pSIR31,datagte500pSIR32,datagte500pSIR33,datagte500pSIR34,datagte500pSIR35, 
     datagte500pySIR12,datagte500pySIR13,datagte500pySIR14,datagte500pySIR15,
     datagte500pySIR21,datagte500pySIR22,datagte500pySIR23,datagte500pySIR24,datagte500pySIR25,
     datagte500pySIR31,datagte500pySIR32,datagte500pySIR33,datagte500pySIR34,datagte500pySIR35, 
     datagte500wsSIR11,datagte500wsSIR12,datagte500wsSIR13,datagte500wsSIR14,datagte500wsSIR15,
     datagte500wsSIR21,datagte500wsSIR22,datagte500wsSIR23,datagte500wsSIR24,datagte500wsSIR25,
     datagte500wsSIR31,datagte500wsSIR32,datagte500wsSIR33,datagte500wsSIR34,datagte500wsSIR35, 
     datagte1000_pois_SIR,datagte1000_poly_SIR,datagte1000_ws_SIR,
     datagte1000pSIR11,datagte1000pSIR12,datagte1000pSIR15,
     datagte1000pSIR21,datagte1000pSIR22,datagte1000pSIR23,datagte1000pSIR24,datagte1000pSIR25,
     datagte1000pSIR31,datagte1000pSIR33,datagte1000pSIR34,datagte1000pSIR35,  
     datagte1000pySIR11,datagte1000pySIR12,datagte1000pySIR13,datagte1000pySIR14,datagte1000pySIR15,
     datagte1000pySIR21,datagte1000pySIR22,datagte1000pySIR23,datagte1000pySIR24,datagte1000pySIR25,
     datagte1000pySIR31,datagte1000pySIR32,datagte1000pySIR33,datagte1000pySIR34,datagte1000pySIR35,  
     datagte1000wsSIR11,datagte1000wsSIR12,datagte1000wsSIR13,datagte1000wsSIR14,datagte1000wsSIR15,
     datagte1000wsSIR21,datagte1000wsSIR22,datagte1000wsSIR23,datagte1000wsSIR24,datagte1000wsSIR25,
     datagte1000wsSIR31,datagte1000wsSIR32,datagte1000wsSIR33,datagte1000wsSIR34,datagte1000wsSIR35, 
     datagte10000_pois_SIR,datagte10000_poly_SIR,datagte10000_ws_SIR,
     datagte10000pSIR12,datagte10000pSIR13,datagte10000pSIR14,datagte10000pSIR15,
     datagte10000pSIR21,datagte10000pSIR22,datagte10000pSIR23,datagte10000pSIR24,datagte10000pSIR25,
     datagte10000pSIR31,datagte10000pSIR32,datagte10000pSIR33,datagte10000pSIR34,datagte10000pSIR35,  
     datagte10000pySIR11,datagte10000pySIR12,datagte10000pySIR13,datagte10000pySIR14,datagte10000pySIR15,
     datagte10000pySIR21,datagte10000pySIR22,datagte10000pySIR23,datagte10000pySIR24,datagte10000pySIR25,
     datagte10000pySIR31,datagte10000pySIR32,datagte10000pySIR33,datagte10000pySIR34,datagte10000pySIR35,  
     datagte10000wsSIR11,datagte10000wsSIR12,datagte10000wsSIR13,datagte10000wsSIR14,datagte10000wsSIR15,
     datagte10000wsSIR21,datagte10000wsSIR22,datagte10000wsSIR23,datagte10000wsSIR24,datagte10000wsSIR25,
     datagte10000wsSIR31,datagte10000wsSIR32,datagte10000wsSIR33,datagte10000wsSIR34,datagte10000wsSIR35,  
     datagte50000_pois_SIR,datagte50000_poly_SIR,datagte50000_ws_SIR,
     datagte50000pSIR11,datagte50000pSIR12,datagte50000pSIR13,datagte50000pSIR14,datagte50000pSIR15,
     datagte50000pSIR21,datagte50000pSIR22,datagte50000pSIR23,datagte50000pSIR24,
     datagte50000pSIR31,datagte50000pSIR32,datagte50000pSIR33,datagte50000pSIR34,datagte50000pSIR35, 
     datagte50000pySIR11,datagte50000pySIR12,datagte50000pySIR14,
     datagte50000pySIR21,datagte50000pySIR22,datagte50000pySIR23,datagte50000pySIR25,
     datagte50000pySIR31,datagte50000pySIR32,datagte50000pySIR33,datagte50000pySIR34,datagte50000pySIR35, 
     datagte50000wsSIR11,datagte50000wsSIR12,datagte50000wsSIR13,datagte50000wsSIR14,datagte50000wsSIR15,
     datagte50000wsSIR21,datagte50000wsSIR22,datagte50000wsSIR23,datagte50000wsSIR24,datagte50000wsSIR25,
     datagte50000wsSIR31,datagte50000wsSIR33,datagte50000wsSIR34,datagte50000wsSIR35, 
     datagte_SIR_100,datagte_SIR_500,datagte_SIR_1000,datagte_SIR_10000,datagte_SIR_50000,
     file=paste("Todas_las_datagtes_SIR_","100_extra",".RData",sep = ""))
