
# creacion de datas SEIR con las covariables extras


setwd("C:\\Users\\Diana\\Dropbox\\Diana\\CodigoR")
load("Brotes_SEIR_Pois_Poly_WS.RData")
load("pois_5difnodos_15sim.RData")
load("poly_5difnodos_15sim.RData")
load("WS_5difnodos_15sim.RData")
load("Cualidades de red, extra.RData")



m<-100



#crear_datas<-function(m){

      for (n in 1:5){
        if (n==1) nodos<-100
        else if (n==2) nodos<-500
        else if (n==3) nodos<-1000
        else if (n==4) nodos<-10000
        else nodos<-50000
        
        ######################################## Poisson #########################################3
        assign(paste("data",nodos,"_pois_SEIR", sep=""),"a")
        assign(paste("datagte",nodos,"_pois_SEIR", sep=""),"b")
        
        for (r in 1:3){
          if (r==1){
            assign(paste("red_pois_",nodos, sep = ""),c(paste("sim_SEIR_pois_",nodos,"_1_1", sep = ""), paste("sim_SEIR_pois_",nodos,"_1_2", sep = ""), paste("sim_SEIR_pois_",nodos,"_1_3", sep = ""), paste("sim_SEIR_pois_",nodos,"_1_4", sep = ""), paste("sim_SEIR_pois_",nodos,"_1_5", sep = "")))
            assign("diam",paste("diam_pois_",nodos,"_1", sep = ""))
            assign("med",paste("med_pois_",nodos,"_1", sep = ""))
            assign("var",paste("var_pois_",nodos,"_1", sep = ""))
            assign("coef",paste("coef_pois_",nodos,"_1", sep = ""))
            assign("assor",paste("assor_pois_",nodos,"_1", sep = ""))
            assign("mdist",paste("mdist_pois_",nodos,"_1", sep = ""))
            assign("eigen",paste("eigen_pois_",nodos,"_1", sep = ""))
          } 
          
          if (r==2){
            assign(paste("red_pois_",nodos, sep = ""),c(paste("sim_SEIR_pois_",nodos,"_2_1", sep = ""), paste("sim_SEIR_pois_",nodos,"_2_2", sep = ""), paste("sim_SEIR_pois_",nodos,"_2_3", sep = ""), paste("sim_SEIR_pois_",nodos,"_2_4", sep = ""), paste("sim_SEIR_pois_",nodos,"_2_5", sep = "")))
            assign("diam",paste("diam_pois_",nodos,"_2", sep = ""))
            assign("med",paste("med_pois_",nodos,"_2", sep = ""))
            assign("var",paste("var_pois_",nodos,"_2", sep = ""))
            assign("coef",paste("coef_pois_",nodos,"_2", sep = ""))
            assign("assor",paste("assor_pois_",nodos,"_2", sep = ""))
            assign("mdist",paste("mdist_pois_",nodos,"_2", sep = ""))
            assign("eigen",paste("eigen_pois_",nodos,"_2", sep = ""))
          }
          
          if (r==3){
            assign(paste("red_pois_",nodos, sep = ""),c(paste("sim_SEIR_pois_",nodos,"_3_1", sep = ""), paste("sim_SEIR_pois_",nodos,"_3_2", sep = ""), paste("sim_SEIR_pois_",nodos,"_3_3", sep = ""), paste("sim_SEIR_pois_",nodos,"_3_4", sep = ""), paste("sim_SEIR_pois_",nodos,"_3_5", sep = "")))
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
            interv_SEIR<-seq(0,floor(max(aa)+1)+1,by=1)
            pmediointerv_SEIR<-(interv_SEIR[-1]+interv_SEIR[-length(interv_SEIR)])/2
            newinfdiscr_SEIR<-table(cut(c(aa,pmediointerv_SEIR),interv_SEIR,label=FALSE,include.lowest = TRUE))-1
            datos<-data.frame(pmediointerv_SEIR,newinfdiscr_SEIR)
            
            if (length(datos$Freq)>(m+10)){  
              
              y<-datos$Freq    #time_series <- rnorm(nodos)
              u<-length(y)
              tmp <- rep(NA, u)
              lagged_data <- data.frame(V1 = tmp)
              for(i in seq_len(m)) {     #seq_len(max_lag) secuencia de 1 hasta max_lag
                lagged_data[, i] <- c(rep(NA, i), head(y[-1], u - i))
              }
              lagged_data <- lagged_data[complete.cases(lagged_data),]
              
              assign(paste("data",nodos,"pSEIR", r, q, sep=""),data.frame(yinf = y[(m+2):(u+1)], lagged_data, x1=get(diam),x2=get(med),x3=get(var),x4=get(coef),x5=get(assor)))
              
              filas<-nrow(get(paste("data",nodos,"pSEIR", r, q, sep="")))
              #numero de suceptibles al incio de cada intervalo (de una hora)
              suc<-nodos-cumsum(newinfdiscr_SEIR) 
              yp<-((get(paste("data",nodos,"pSEIR", r, q, sep="")))[,1])/(suc[(m+1):length(suc)])
              W1<-suc[(m+1):length(suc)]
              
              assign(paste("datagte",nodos,"pSEIR", r, q, sep=""),(data.frame(yinfp=yp,yinf = y[(m+2):(u+1)], lagged_data, x1=get(diam),x2=get(med),x3=get(var),x4=get(coef),x5=get(assor),x6=get(mdist),x7=get(eigen),red=rep("Pois",filas),sim=rep(r,filas),brote=rep(q,filas),long=(1:filas),W1=W1))[-filas,])
              
              if (get(paste("data",nodos,"_pois_SEIR", sep=""))=="a") assign(paste("data",nodos,"_pois_SEIR", sep=""),rbind(get(paste("data",nodos,"pSEIR", r, q, sep=""))))
              else assign(paste("data",nodos,"_pois_SEIR", sep=""),rbind(get(paste("data",nodos,"_pois_SEIR", sep="")),get(paste("data",nodos,"pSEIR", r, q, sep=""))))
              
              if (get(paste("datagte",nodos,"_pois_SEIR", sep=""))=="b") assign(paste("datagte",nodos,"_pois_SEIR", sep=""),rbind(get(paste("datagte",nodos,"pSEIR", r, q, sep=""))))
              else assign(paste("datagte",nodos,"_pois_SEIR", sep=""),rbind(get(paste("datagte",nodos,"_pois_SEIR", sep="")),get(paste("datagte",nodos,"pSEIR", r, q, sep=""))))
            } #if
          } #for i
        }  #for r
        ################################################################################ 
        ######################################## Polylog #########################################3
      
        assign(paste("data",nodos,"_poly_SEIR", sep=""),"a")
        assign(paste("datagte",nodos,"_poly_SEIR", sep=""),"b")
        
        for (r in 1:3){
          if (r==1){
            assign(paste("red_poly_",nodos, sep = ""),c(paste("sim_SEIR_poly_",nodos,"_1_1", sep = ""), paste("sim_SEIR_poly_",nodos,"_1_2", sep = ""), paste("sim_SEIR_poly_",nodos,"_1_3", sep = ""), paste("sim_SEIR_poly_",nodos,"_1_4", sep = ""), paste("sim_SEIR_poly_",nodos,"_1_5", sep = "")))
            assign("diam",paste("diam_poly_",nodos,"_1", sep = ""))
            assign("med",paste("med_poly_",nodos,"_1", sep = ""))
            assign("var",paste("var_poly_",nodos,"_1", sep = ""))
            assign("coef",paste("coef_poly_",nodos,"_1", sep = ""))
            assign("assor",paste("assor_poly_",nodos,"_1", sep = ""))
            assign("mdist",paste("mdist_poly_",nodos,"_1", sep = ""))
            assign("eigen",paste("eigen_poly_",nodos,"_1", sep = ""))
          } 
          
          if (r==2){
            assign(paste("red_poly_",nodos, sep = ""),c(paste("sim_SEIR_poly_",nodos,"_2_1", sep = ""), paste("sim_SEIR_poly_",nodos,"_2_2", sep = ""), paste("sim_SEIR_poly_",nodos,"_2_3", sep = ""), paste("sim_SEIR_poly_",nodos,"_2_4", sep = ""), paste("sim_SEIR_poly_",nodos,"_2_5", sep = "")))
            assign("diam",paste("diam_poly_",nodos,"_2", sep = ""))
            assign("med",paste("med_poly_",nodos,"_2", sep = ""))
            assign("var",paste("var_poly_",nodos,"_2", sep = ""))
            assign("coef",paste("coef_poly_",nodos,"_2", sep = ""))
            assign("assor",paste("assor_poly_",nodos,"_2", sep = ""))
            assign("mdist",paste("mdist_poly_",nodos,"_2", sep = ""))
            assign("eigen",paste("eigen_poly_",nodos,"_2", sep = ""))
          }
          
          if (r==3){
            assign(paste("red_poly_",nodos, sep = ""),c(paste("sim_SEIR_poly_",nodos,"_3_1", sep = ""), paste("sim_SEIR_poly_",nodos,"_3_2", sep = ""), paste("sim_SEIR_poly_",nodos,"_3_3", sep = ""), paste("sim_SEIR_poly_",nodos,"_3_4", sep = ""), paste("sim_SEIR_poly_",nodos,"_3_5", sep = "")))
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
            interv_SEIR<-seq(0,floor(max(aa)+1)+1,by=1)
            pmediointerv_SEIR<-(interv_SEIR[-1]+interv_SEIR[-length(interv_SEIR)])/2
            newinfdiscr_SEIR<-table(cut(c(aa,pmediointerv_SEIR),interv_SEIR,label=FALSE,include.lowest = TRUE))-1
            datos<-data.frame(pmediointerv_SEIR,newinfdiscr_SEIR)
            
            if (length(datos$Freq)>(m+10)){  
              
              y<-datos$Freq    #time_series <- rnorm(nodos)
              u<-length(y)
              tmp <- rep(NA, u)
              lagged_data <- data.frame(V1 = tmp)
              for(i in seq_len(m)) {     #seq_len(max_lag) secuencia de 1 hasta max_lag
                lagged_data[, i] <- c(rep(NA, i), head(y[-1], u - i))
              }
              lagged_data <- lagged_data[complete.cases(lagged_data),]
              
              assign(paste("data",nodos,"pySEIR", r, q, sep=""),data.frame(yinf = y[(m+2):(u+1)], lagged_data, x1=get(diam),x2=get(med),x3=get(var),x4=get(coef),x5=get(assor)))
              
              filas<-nrow(get(paste("data",nodos,"pySEIR", r, q, sep="")))
              #numero de suceptibles al incio de cada intervalo (de una hora)
              suc<-nodos-cumsum(newinfdiscr_SEIR) 
              yp<-((get(paste("data",nodos,"pySEIR", r, q, sep="")))[,1])/(suc[(m+1):length(suc)])
              W1<-suc[(m+1):length(suc)]
              
              assign(paste("datagte",nodos,"pySEIR", r, q, sep=""),(data.frame(yinfp=yp,yinf = y[(m+2):(u+1)], lagged_data, x1=get(diam),x2=get(med),x3=get(var),x4=get(coef),x5=get(assor),x6=get(mdist),x7=get(eigen),red=rep("Poly",filas),sim=rep(r,filas),brote=rep(q,filas),long=(1:filas),W1=W1))[-filas,])
              
              if (get(paste("data",nodos,"_poly_SEIR", sep=""))=="a") assign(paste("data",nodos,"_poly_SEIR", sep=""),rbind(get(paste("data",nodos,"pySEIR", r, q, sep=""))))
              else assign(paste("data",nodos,"_poly_SEIR", sep=""),rbind(get(paste("data",nodos,"_poly_SEIR", sep="")),get(paste("data",nodos,"pySEIR", r, q, sep=""))))
              
              if (get(paste("datagte",nodos,"_poly_SEIR", sep=""))=="b") assign(paste("datagte",nodos,"_poly_SEIR", sep=""),rbind(get(paste("datagte",nodos,"pySEIR", r, q, sep=""))))
              else assign(paste("datagte",nodos,"_poly_SEIR", sep=""),rbind(get(paste("datagte",nodos,"_poly_SEIR", sep="")),get(paste("datagte",nodos,"pySEIR", r, q, sep=""))))
            } #if  
          } #for i
        }  #for r
        ################################################################################ 
        ######################################## Watts Strogatz #########################################3
      
        assign(paste("data",nodos,"_ws_SEIR", sep=""),"a")
        assign(paste("datagte",nodos,"_ws_SEIR", sep=""),"b")
        
        for (r in 1:3){
          if (r==1){
            assign(paste("red_ws_",nodos, sep = ""),c(paste("sim_SEIR_ws_",nodos,"_1_1", sep = ""), paste("sim_SEIR_ws_",nodos,"_1_2", sep = ""), paste("sim_SEIR_ws_",nodos,"_1_3", sep = ""), paste("sim_SEIR_ws_",nodos,"_1_4", sep = ""), paste("sim_SEIR_ws_",nodos,"_1_5", sep = "")))
            assign("diam",paste("diam_ws_",nodos,"_1", sep = ""))
            assign("med",paste("med_ws_",nodos,"_1", sep = ""))
            assign("var",paste("var_ws_",nodos,"_1", sep = ""))
            assign("coef",paste("coef_ws_",nodos,"_1", sep = ""))
            assign("assor",paste("assor_ws_",nodos,"_1", sep = ""))
            assign("mdist",paste("mdist_ws_",nodos,"_1", sep = ""))
            assign("eigen",paste("eigen_ws_",nodos,"_1", sep = ""))
          } 
          
          if (r==2){
            assign(paste("red_ws_",nodos, sep = ""),c(paste("sim_SEIR_ws_",nodos,"_2_1", sep = ""), paste("sim_SEIR_ws_",nodos,"_2_2", sep = ""), paste("sim_SEIR_ws_",nodos,"_2_3", sep = ""), paste("sim_SEIR_ws_",nodos,"_2_4", sep = ""), paste("sim_SEIR_ws_",nodos,"_2_5", sep = "")))
            assign("diam",paste("diam_ws_",nodos,"_2", sep = ""))
            assign("med",paste("med_ws_",nodos,"_2", sep = ""))
            assign("var",paste("var_ws_",nodos,"_2", sep = ""))
            assign("coef",paste("coef_ws_",nodos,"_2", sep = ""))
            assign("assor",paste("assor_ws_",nodos,"_2", sep = ""))
            assign("mdist",paste("mdist_ws_",nodos,"_2", sep = ""))
            assign("eigen",paste("eigen_ws_",nodos,"_2", sep = ""))
          }
          
          if (r==3){
            assign(paste("red_ws_",nodos, sep = ""),c(paste("sim_SEIR_ws_",nodos,"_3_1", sep = ""), paste("sim_SEIR_ws_",nodos,"_3_2", sep = ""), paste("sim_SEIR_ws_",nodos,"_3_3", sep = ""), paste("sim_SEIR_ws_",nodos,"_3_4", sep = ""), paste("sim_SEIR_ws_",nodos,"_3_5", sep = "")))
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
            interv_SEIR<-seq(0,floor(max(aa)+1)+1,by=1)
            pmediointerv_SEIR<-(interv_SEIR[-1]+interv_SEIR[-length(interv_SEIR)])/2
            newinfdiscr_SEIR<-table(cut(c(aa,pmediointerv_SEIR),interv_SEIR,label=FALSE,include.lowest = TRUE))-1
            datos<-data.frame(pmediointerv_SEIR,newinfdiscr_SEIR)
            
            if (length(datos$Freq)>(m+10)){  
              
              y<-datos$Freq    #time_series <- rnorm(nodos)
              u<-length(y)
              tmp <- rep(NA, u)
              lagged_data <- data.frame(V1 = tmp)
              for(i in seq_len(m)) {     #seq_len(max_lag) secuencia de 1 hasta max_lag
                lagged_data[, i] <- c(rep(NA, i), head(y[-1], u - i))
              }
              lagged_data <- lagged_data[complete.cases(lagged_data),]
              
              assign(paste("data",nodos,"wsSEIR", r, q, sep=""),data.frame(yinf = y[(m+2):(u+1)], lagged_data, x1=get(diam),x2=get(med),x3=get(var),x4=get(coef),x5=get(assor)))
              
              filas<-nrow(get(paste("data",nodos,"wsSEIR", r, q, sep="")))
              #numero de suceptibles al incio de cada intervalo (de una hora)
              suc<-nodos-cumsum(newinfdiscr_SEIR) 
              yp<-((get(paste("data",nodos,"wsSEIR", r, q, sep="")))[,1])/(suc[(m+1):length(suc)])
              W1<-suc[(m+1):length(suc)]
              
              assign(paste("datagte",nodos,"wsSEIR", r, q, sep=""),(data.frame(yinfp=yp,yinf = y[(m+2):(u+1)], lagged_data, x1=get(diam),x2=get(med),x3=get(var),x4=get(coef),x5=get(assor),x6=get(mdist),x7=get(eigen),red=rep("WS",filas),sim=rep(r,filas),brote=rep(q,filas),long=(1:filas),W1=W1))[-filas,])
              
              if (get(paste("data",nodos,"_ws_SEIR", sep=""))=="a") assign(paste("data",nodos,"_ws_SEIR", sep=""),rbind(get(paste("data",nodos,"wsSEIR", r, q, sep=""))))
              else assign(paste("data",nodos,"_ws_SEIR", sep=""),rbind(get(paste("data",nodos,"_ws_SEIR", sep="")),get(paste("data",nodos,"wsSEIR", r, q, sep=""))))
              
              if (get(paste("datagte",nodos,"_ws_SEIR", sep=""))=="b") assign(paste("datagte",nodos,"_ws_SEIR", sep=""),rbind(get(paste("datagte",nodos,"wsSEIR", r, q, sep=""))))
              else assign(paste("datagte",nodos,"_ws_SEIR", sep=""),rbind(get(paste("datagte",nodos,"_ws_SEIR", sep="")),get(paste("datagte",nodos,"wsSEIR", r, q, sep=""))))
            } #if
          } #for i
        }  #for r
        ################################################################################ 
      } #for n
      
      data_SEIR_100<-rbind(data100_pois_SEIR,data100_poly_SEIR,data100_ws_SEIR)
      data_SEIR_500<-rbind(data500_pois_SEIR,data500_poly_SEIR,data500_ws_SEIR)
      data_SEIR_1000<-rbind(data1000_pois_SEIR,data1000_poly_SEIR,data1000_ws_SEIR)
      data_SEIR_10000<-rbind(data10000_pois_SEIR,data10000_poly_SEIR,data10000_ws_SEIR)
      data_SEIR_50000<-rbind(data50000_pois_SEIR,data50000_poly_SEIR,data50000_ws_SEIR)
      
      datagte_SEIR_100<-rbind(datagte100_pois_SEIR,datagte100_poly_SEIR,datagte100_ws_SEIR)
      datagte_SEIR_500<-rbind(datagte500_pois_SEIR,datagte500_poly_SEIR,datagte500_ws_SEIR)
      datagte_SEIR_1000<-rbind(datagte1000_pois_SEIR,datagte1000_poly_SEIR,datagte1000_ws_SEIR)
      datagte_SEIR_10000<-rbind(datagte10000_pois_SEIR,datagte10000_poly_SEIR,datagte10000_ws_SEIR)
      datagte_SEIR_50000<-rbind(datagte50000_pois_SEIR,datagte50000_poly_SEIR,datagte50000_ws_SEIR)

      
#}










save(datagte100_pois_SEIR,datagte100_poly_SEIR,datagte100_ws_SEIR,
     datagte100pSEIR11,datagte100pSEIR12,datagte100pSEIR13,datagte100pSEIR14,datagte100pSEIR15,
     datagte100pSEIR21,datagte100pSEIR22,datagte100pSEIR23,datagte100pSEIR24,datagte100pSEIR25,
     datagte100pSEIR32,datagte100pSEIR33,datagte100pSEIR34,datagte100pSEIR35,
     datagte100pySEIR11,datagte100pySEIR12,datagte100pySEIR13,datagte100pySEIR14,datagte100pySEIR15,
     datagte100pySEIR21,datagte100pySEIR22,datagte100pySEIR23,datagte100pySEIR24,datagte100pySEIR25,
     datagte100pySEIR31,datagte100pySEIR32,datagte100pySEIR33,datagte100pySEIR34,datagte100pySEIR35,
     datagte100wsSEIR11,datagte100wsSEIR12,datagte100wsSEIR13,datagte100wsSEIR14,datagte100wsSEIR15,
     datagte100wsSEIR21,datagte100wsSEIR22,datagte100wsSEIR23,datagte100wsSEIR24,datagte100wsSEIR25,
     datagte100wsSEIR31,datagte100wsSEIR32,datagte100wsSEIR33,datagte100wsSEIR34,datagte100wsSEIR35,
     datagte500_pois_SEIR,datagte500_poly_SEIR,datagte500_ws_SEIR,
     datagte500pSEIR11,datagte500pSEIR12,datagte500pSEIR13,datagte500pSEIR15,
     datagte500pSEIR21,datagte500pSEIR22,datagte500pSEIR23,datagte500pSEIR24,datagte500pSEIR25,
     datagte500pSEIR31,datagte500pSEIR33,datagte500pSEIR34,
     datagte500pySEIR11,datagte500pySEIR12,datagte500pySEIR13,datagte500pySEIR14,datagte500pySEIR15,
     datagte500pySEIR21,datagte500pySEIR22,datagte500pySEIR23,datagte500pySEIR24,datagte500pySEIR25,
     datagte500pySEIR31,datagte500pySEIR32,datagte500pySEIR33,datagte500pySEIR34,datagte500pySEIR35, 
     datagte500wsSEIR12,datagte500wsSEIR13,datagte500wsSEIR14,datagte500wsSEIR15,
     datagte500wsSEIR21,datagte500wsSEIR22,datagte500wsSEIR23,datagte500wsSEIR24,datagte500wsSEIR25,
     datagte500wsSEIR31,datagte500wsSEIR32,datagte500wsSEIR33,datagte500wsSEIR34,datagte500wsSEIR35, 
     datagte1000_pois_SEIR,datagte1000_poly_SEIR,datagte1000_ws_SEIR,
     datagte1000pSEIR11,datagte1000pSEIR12,datagte1000pSEIR13,datagte1000pSEIR14,datagte1000pSEIR15,
     datagte1000pSEIR21,datagte1000pSEIR22,datagte1000pSEIR23,datagte1000pSEIR24,datagte1000pSEIR25,
     datagte1000pSEIR31,datagte1000pSEIR32,datagte1000pSEIR33,datagte1000pSEIR34,datagte1000pSEIR35,  
     datagte1000pySEIR11,datagte1000pySEIR12,datagte1000pySEIR13,datagte1000pySEIR14,datagte1000pySEIR15,
     datagte1000pySEIR21,datagte1000pySEIR22,datagte1000pySEIR23,datagte1000pySEIR24,datagte1000pySEIR25,
     datagte1000pySEIR31,datagte1000pySEIR32,datagte1000pySEIR33,datagte1000pySEIR34,datagte1000pySEIR35,  
     datagte1000wsSEIR11,datagte1000wsSEIR12,datagte1000wsSEIR13,datagte1000wsSEIR14,datagte1000wsSEIR15,
     datagte1000wsSEIR21,datagte1000wsSEIR22,datagte1000wsSEIR23,datagte1000wsSEIR24,datagte1000wsSEIR25,
     datagte1000wsSEIR31,datagte1000wsSEIR32,datagte1000wsSEIR33,datagte1000wsSEIR34,datagte1000wsSEIR35, 
     datagte10000_pois_SEIR,datagte10000_poly_SEIR,datagte10000_ws_SEIR,
     datagte10000pSEIR11,datagte10000pSEIR12,datagte10000pSEIR13,datagte10000pSEIR14,datagte10000pSEIR15,
     datagte10000pSEIR21,datagte10000pSEIR22,datagte10000pSEIR24,datagte10000pSEIR25,
     datagte10000pSEIR31,datagte10000pSEIR32,datagte10000pSEIR33,datagte10000pSEIR34,datagte10000pSEIR35,  
     datagte10000pySEIR11,datagte10000pySEIR12,datagte10000pySEIR13,datagte10000pySEIR14,datagte10000pySEIR15,
     datagte10000pySEIR21,datagte10000pySEIR22,datagte10000pySEIR23,datagte10000pySEIR24,datagte10000pySEIR25,
     datagte10000pySEIR31,datagte10000pySEIR33,datagte10000pySEIR34,datagte10000pySEIR35,  
     datagte10000wsSEIR11,datagte10000wsSEIR12,datagte10000wsSEIR13,datagte10000wsSEIR14,datagte10000wsSEIR15,
     datagte10000wsSEIR21,datagte10000wsSEIR22,datagte10000wsSEIR23,datagte10000wsSEIR24,datagte10000wsSEIR25,
     datagte10000wsSEIR31,datagte10000wsSEIR32,datagte10000wsSEIR33,datagte10000wsSEIR34,datagte10000wsSEIR35,  
     datagte50000_pois_SEIR,datagte50000_poly_SEIR,datagte50000_ws_SEIR,
     datagte50000pSEIR12,datagte50000pSEIR13,datagte50000pSEIR14,datagte50000pSEIR15,
     datagte50000pSEIR21,datagte50000pSEIR22,datagte50000pSEIR23,datagte50000pSEIR25,
     datagte50000pSEIR31,datagte50000pSEIR32,datagte50000pSEIR33,datagte50000pSEIR34,datagte50000pSEIR35, 
     datagte50000pySEIR11,datagte50000pySEIR12,datagte50000pySEIR13,datagte50000pySEIR14,datagte50000pySEIR15,
     datagte50000pySEIR21,datagte50000pySEIR22,datagte50000pySEIR23,datagte50000pySEIR24,datagte50000pySEIR25,
     datagte50000pySEIR31,datagte50000pySEIR33,datagte50000pySEIR34,datagte50000pySEIR35, 
     datagte50000wsSEIR11,datagte50000wsSEIR12,datagte50000wsSEIR13,datagte50000wsSEIR14,datagte50000wsSEIR15,
     datagte50000wsSEIR21,datagte50000wsSEIR22,datagte50000wsSEIR23,datagte50000wsSEIR24,datagte50000wsSEIR25,
     datagte50000wsSEIR31,datagte50000wsSEIR32,datagte50000wsSEIR33,datagte50000wsSEIR34,datagte50000wsSEIR35, 
     datagte_SEIR_100,datagte_SEIR_500,datagte_SEIR_1000,datagte_SEIR_10000,datagte_SEIR_50000,
     file="Todas_las_datagtes_SEIR_20_extra.RData")


save(datagte100_pois_SEIR,datagte100_poly_SEIR,datagte100_ws_SEIR,
     datagte100pSEIR11,datagte100pSEIR12,datagte100pSEIR13,datagte100pSEIR14,datagte100pSEIR15,
     datagte100pSEIR21,datagte100pSEIR22,datagte100pSEIR23,datagte100pSEIR24,datagte100pSEIR25,
     datagte100pSEIR32,datagte100pSEIR33,datagte100pSEIR34,datagte100pSEIR35,
     datagte100pySEIR11,datagte100pySEIR12,datagte100pySEIR13,datagte100pySEIR14,datagte100pySEIR15,
     datagte100pySEIR21,datagte100pySEIR22,datagte100pySEIR23,datagte100pySEIR24,datagte100pySEIR25,
     datagte100pySEIR31,datagte100pySEIR32,datagte100pySEIR33,datagte100pySEIR34,datagte100pySEIR35,
     datagte100wsSEIR11,datagte100wsSEIR12,datagte100wsSEIR13,datagte100wsSEIR14,datagte100wsSEIR15,
     datagte100wsSEIR21,datagte100wsSEIR22,datagte100wsSEIR23,datagte100wsSEIR24,datagte100wsSEIR25,
     datagte100wsSEIR31,datagte100wsSEIR32,datagte100wsSEIR33,datagte100wsSEIR34,datagte100wsSEIR35,
     datagte500_pois_SEIR,datagte500_poly_SEIR,datagte500_ws_SEIR,
     datagte500pSEIR11,datagte500pSEIR12,datagte500pSEIR13,datagte500pSEIR14,datagte500pSEIR15,
     datagte500pSEIR21,datagte500pSEIR22,datagte500pSEIR23,datagte500pSEIR24,datagte500pSEIR25,
     datagte500pSEIR31,datagte500pSEIR33,datagte500pSEIR34,
     datagte500pySEIR12,datagte500pySEIR13,datagte500pySEIR14,datagte500pySEIR15,
     datagte500pySEIR21,datagte500pySEIR22,datagte500pySEIR23,datagte500pySEIR24,datagte500pySEIR25,
     datagte500pySEIR31,datagte500pySEIR32,datagte500pySEIR33,datagte500pySEIR34,datagte500pySEIR35, 
     datagte500wsSEIR12,datagte500wsSEIR13,datagte500wsSEIR14,datagte500wsSEIR15,
     datagte500wsSEIR21,datagte500wsSEIR22,datagte500wsSEIR23,datagte500wsSEIR24,datagte500wsSEIR25,
     datagte500wsSEIR31,datagte500wsSEIR32,datagte500wsSEIR33,datagte500wsSEIR34,datagte500wsSEIR35, 
     datagte1000_pois_SEIR,datagte1000_poly_SEIR,datagte1000_ws_SEIR,
     datagte1000pSEIR11,datagte1000pSEIR12,datagte1000pSEIR14,datagte1000pSEIR15,
     datagte1000pSEIR21,datagte1000pSEIR22,datagte1000pSEIR23,datagte1000pSEIR24,datagte1000pSEIR25,
     datagte1000pSEIR31,datagte1000pSEIR32,datagte1000pSEIR33,datagte1000pSEIR34,datagte1000pSEIR35,  
     datagte1000pySEIR11,datagte1000pySEIR12,datagte1000pySEIR13,datagte1000pySEIR14,datagte1000pySEIR15,
     datagte1000pySEIR21,datagte1000pySEIR22,datagte1000pySEIR23,datagte1000pySEIR24,datagte1000pySEIR25,
     datagte1000pySEIR31,datagte1000pySEIR32,datagte1000pySEIR33,datagte1000pySEIR34,datagte1000pySEIR35,  
     datagte1000wsSEIR11,datagte1000wsSEIR12,datagte1000wsSEIR13,datagte1000wsSEIR14,datagte1000wsSEIR15,
     datagte1000wsSEIR21,datagte1000wsSEIR22,datagte1000wsSEIR23,datagte1000wsSEIR24,datagte1000wsSEIR25,
     datagte1000wsSEIR31,datagte1000wsSEIR32,datagte1000wsSEIR33,datagte1000wsSEIR34,datagte1000wsSEIR35, 
     datagte10000_pois_SEIR,datagte10000_poly_SEIR,datagte10000_ws_SEIR,
     datagte10000pSEIR11,datagte10000pSEIR12,datagte10000pSEIR13,datagte10000pSEIR14,datagte10000pSEIR15,
     datagte10000pSEIR21,datagte10000pSEIR22,datagte10000pSEIR24,datagte10000pSEIR25,
     datagte10000pSEIR31,datagte10000pSEIR32,datagte10000pSEIR33,datagte10000pSEIR34,datagte10000pSEIR35,  
     datagte10000pySEIR11,datagte10000pySEIR12,datagte10000pySEIR13,datagte10000pySEIR14,datagte10000pySEIR15,
     datagte10000pySEIR21,datagte10000pySEIR22,datagte10000pySEIR23,datagte10000pySEIR24,datagte10000pySEIR25,
     datagte10000pySEIR31,datagte10000pySEIR33,datagte10000pySEIR34,datagte10000pySEIR35,  
     datagte10000wsSEIR11,datagte10000wsSEIR12,datagte10000wsSEIR13,datagte10000wsSEIR14,datagte10000wsSEIR15,
     datagte10000wsSEIR21,datagte10000wsSEIR22,datagte10000wsSEIR23,datagte10000wsSEIR24,datagte10000wsSEIR25,
     datagte10000wsSEIR31,datagte10000wsSEIR32,datagte10000wsSEIR33,datagte10000wsSEIR34,datagte10000wsSEIR35,  
     datagte50000_pois_SEIR,datagte50000_poly_SEIR,datagte50000_ws_SEIR,
     datagte50000pSEIR12,datagte50000pSEIR13,datagte50000pSEIR14,datagte50000pSEIR15,
     datagte50000pSEIR21,datagte50000pSEIR22,datagte50000pSEIR23,datagte50000pSEIR25,
     datagte50000pSEIR31,datagte50000pSEIR32,datagte50000pSEIR33,datagte50000pSEIR34,datagte50000pSEIR35, 
     datagte50000pySEIR11,datagte50000pySEIR13,datagte50000pySEIR14,datagte50000pySEIR15,
     datagte50000pySEIR21,datagte50000pySEIR22,datagte50000pySEIR23,datagte50000pySEIR24,datagte50000pySEIR25,
     datagte50000pySEIR31,datagte50000pySEIR33,datagte50000pySEIR34,datagte50000pySEIR35, 
     datagte50000wsSEIR11,datagte50000wsSEIR12,datagte50000wsSEIR13,datagte50000wsSEIR14,datagte50000wsSEIR15,
     datagte50000wsSEIR21,datagte50000wsSEIR22,datagte50000wsSEIR23,datagte50000wsSEIR24,datagte50000wsSEIR25,
     datagte50000wsSEIR31,datagte50000wsSEIR32,datagte50000wsSEIR33,datagte50000wsSEIR34,datagte50000wsSEIR35, 
     datagte_SEIR_100,datagte_SEIR_500,datagte_SEIR_1000,datagte_SEIR_10000,datagte_SEIR_50000,
     file="Todas_las_datagtes_SEIR_40_extra.RData")



save(datagte100_pois_SEIR,datagte100_poly_SEIR,datagte100_ws_SEIR,
     datagte100pSEIR11,datagte100pSEIR12,datagte100pSEIR13,datagte100pSEIR14,datagte100pSEIR15,
     datagte100pSEIR21,datagte100pSEIR22,datagte100pSEIR23,datagte100pSEIR24,datagte100pSEIR25,
     datagte100pSEIR32,datagte100pSEIR33,datagte100pSEIR34,datagte100pSEIR35,
     datagte100pySEIR11,datagte100pySEIR12,datagte100pySEIR13,datagte100pySEIR14,datagte100pySEIR15,
     datagte100pySEIR21,datagte100pySEIR23,datagte100pySEIR24,datagte100pySEIR25,
     datagte100pySEIR31,datagte100pySEIR33,datagte100pySEIR34,datagte100pySEIR35,
     datagte100wsSEIR11,datagte100wsSEIR12,datagte100wsSEIR13,datagte100wsSEIR14,datagte100wsSEIR15,
     datagte100wsSEIR21,datagte100wsSEIR22,datagte100wsSEIR23,datagte100wsSEIR24,datagte100wsSEIR25,
     datagte100wsSEIR31,datagte100wsSEIR32,datagte100wsSEIR33,datagte100wsSEIR34,datagte100wsSEIR35,
     datagte500_pois_SEIR,datagte500_poly_SEIR,datagte500_ws_SEIR,
     datagte500pSEIR11,datagte500pSEIR12,datagte500pSEIR13,datagte500pSEIR15,
     datagte500pSEIR21,datagte500pSEIR22,datagte500pSEIR23,datagte500pSEIR24,datagte500pSEIR25,
     datagte500pSEIR31,datagte500pSEIR33,datagte500pSEIR34,
     datagte500pySEIR11,datagte500pySEIR12,datagte500pySEIR13,datagte500pySEIR14,datagte500pySEIR15,
     datagte500pySEIR21,datagte500pySEIR22,datagte500pySEIR23,datagte500pySEIR24,datagte500pySEIR25,
     datagte500pySEIR31,datagte500pySEIR32,datagte500pySEIR33,datagte500pySEIR34,datagte500pySEIR35, 
     datagte500wsSEIR12,datagte500wsSEIR13,datagte500wsSEIR14,datagte500wsSEIR15,
     datagte500wsSEIR21,datagte500wsSEIR22,datagte500wsSEIR23,datagte500wsSEIR24,datagte500wsSEIR25,
     datagte500wsSEIR31,datagte500wsSEIR32,datagte500wsSEIR33,datagte500wsSEIR34,datagte500wsSEIR35, 
     datagte1000_pois_SEIR,datagte1000_poly_SEIR,datagte1000_ws_SEIR,
     datagte1000pSEIR11,datagte1000pSEIR12,datagte1000pSEIR13,datagte1000pSEIR14,datagte1000pSEIR15,
     datagte1000pSEIR21,datagte1000pSEIR22,datagte1000pSEIR23,datagte1000pSEIR24,datagte1000pSEIR25,
     datagte1000pSEIR31,datagte1000pSEIR32,datagte1000pSEIR33,datagte1000pSEIR34,datagte1000pSEIR35,  
     datagte1000pySEIR11,datagte1000pySEIR12,datagte1000pySEIR13,datagte1000pySEIR14,datagte1000pySEIR15,
     datagte1000pySEIR21,datagte1000pySEIR22,datagte1000pySEIR23,datagte1000pySEIR24,datagte1000pySEIR25,
     datagte1000pySEIR31,datagte1000pySEIR32,datagte1000pySEIR34,datagte1000pySEIR35,  
     datagte1000wsSEIR11,datagte1000wsSEIR12,datagte1000wsSEIR13,datagte1000wsSEIR14,datagte1000wsSEIR15,
     datagte1000wsSEIR21,datagte1000wsSEIR22,datagte1000wsSEIR23,datagte1000wsSEIR24,datagte1000wsSEIR25,
     datagte1000wsSEIR31,datagte1000wsSEIR32,datagte1000wsSEIR33,datagte1000wsSEIR34,datagte1000wsSEIR35, 
     datagte10000_pois_SEIR,datagte10000_poly_SEIR,datagte10000_ws_SEIR,
     datagte10000pSEIR11,datagte10000pSEIR12,datagte10000pSEIR13,datagte10000pSEIR14,datagte10000pSEIR15,
     datagte10000pSEIR21,datagte10000pSEIR22,datagte10000pSEIR24,datagte10000pSEIR25,
     datagte10000pSEIR31,datagte10000pSEIR32,datagte10000pSEIR33,datagte10000pSEIR34,datagte10000pSEIR35,  
     datagte10000pySEIR11,datagte10000pySEIR12,datagte10000pySEIR13,datagte10000pySEIR14,datagte10000pySEIR15,
     datagte10000pySEIR21,datagte10000pySEIR22,datagte10000pySEIR23,datagte10000pySEIR24,datagte10000pySEIR25,
     datagte10000pySEIR31,datagte10000pySEIR33,datagte10000pySEIR34,datagte10000pySEIR35,  
     datagte10000wsSEIR11,datagte10000wsSEIR12,datagte10000wsSEIR13,datagte10000wsSEIR14,datagte10000wsSEIR15,
     datagte10000wsSEIR21,datagte10000wsSEIR22,datagte10000wsSEIR23,datagte10000wsSEIR24,datagte10000wsSEIR25,
     datagte10000wsSEIR31,datagte10000wsSEIR32,datagte10000wsSEIR33,datagte10000wsSEIR34,datagte10000wsSEIR35,  
     datagte50000_pois_SEIR,datagte50000_poly_SEIR,datagte50000_ws_SEIR,
     datagte50000pSEIR12,datagte50000pSEIR13,datagte50000pSEIR14,datagte50000pSEIR15,
     datagte50000pSEIR21,datagte50000pSEIR22,datagte50000pSEIR23,datagte50000pSEIR25,
     datagte50000pSEIR31,datagte50000pSEIR32,datagte50000pSEIR33,datagte50000pSEIR34,datagte50000pSEIR35, 
     datagte50000pySEIR11,datagte50000pySEIR13,datagte50000pySEIR14,datagte50000pySEIR15,
     datagte50000pySEIR21,datagte50000pySEIR22,datagte50000pySEIR23,datagte50000pySEIR24,datagte50000pySEIR25,
     datagte50000pySEIR31,datagte50000pySEIR33,datagte50000pySEIR34,datagte50000pySEIR35, 
     datagte50000wsSEIR11,datagte50000wsSEIR12,datagte50000wsSEIR13,datagte50000wsSEIR14,datagte50000wsSEIR15,
     datagte50000wsSEIR21,datagte50000wsSEIR22,datagte50000wsSEIR23,datagte50000wsSEIR24,datagte50000wsSEIR25,
     datagte50000wsSEIR31,datagte50000wsSEIR32,datagte50000wsSEIR33,datagte50000wsSEIR34,datagte50000wsSEIR35, 
     datagte_SEIR_100,datagte_SEIR_500,datagte_SEIR_1000,datagte_SEIR_10000,datagte_SEIR_50000,
     file="Todas_las_datagtes_SEIR_100_extra.RData")


# este solo es para el SIR 
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
     file="Todas_las_datagtes_SIR_100_extra.RData")
