
#############################################################################################

setwd("C:\\Users\\Diana\\Dropbox\\Diana\\CodigoR") #Es la fuente donde estan los códigos  
source("04_funciones_Redes_LR.R") #Este cógido me lo paso la Dra. Leticia

#############################################################################################


#Creamos una red poisson(3)

redes<-c("pois_100_1", "pois_100_2", "pois_100_3", "pois_500_1", "pois_500_2", "pois_500_3", "pois_1000_1", "pois_1000_2", "pois_1000_3", "pois_10000_1", "pois_10000_2", "pois_10000_3", "pois_50000_1", "pois_50000_2", "pois_50000_3")    
nodos<-c(100, 500, 1000, 10000, 50000)


for(i in 1:15){
  if (i<4) assign(redes[i], local.network(nodos[1], distrib="pois", param=3, degree=NULL, take.p=0.01))
  else if (i<7) assign(redes[i], local.network(nodos[2], distrib="pois", param=3, degree=NULL, take.p=0.01))
  else if (i<10) assign(redes[i], local.network(nodos[3], distrib="pois", param=3, degree=NULL, take.p=0.01))
  else if (i<13) assign(redes[i], local.network(nodos[4], distrib="pois", param=3, degree=NULL, take.p=0.01))
  else assign(redes[i], local.network(nodos[5], distrib="pois", param=3, degree=NULL, take.p=0.01))
  #edges es una matriz de dos columnas que especifica que arcos existen entre dos vertices
  #degree es el grado de cada vertice que se ha utilizado para crear la red
  #degree.left son los grados que no se pudieron asignar en la red construida. Lo ideal es que sea un vector de ceros
  
  assign(paste("grafo", redes[i], sep = "_"), graph_from_edgelist(get(redes[i])$edges))
  
  #Cualidades de la red 
  assign(paste("diam", redes[i], sep = "_"),diameter(get(paste("grafo", redes[i], sep = "_")))) #Diametro 
  assign(paste("med", redes[i], sep = "_"),mean(get(redes[i])$degree)) #Media de grados 
  assign(paste("var", redes[i], sep = "_"),var(get(redes[i])$degree)) #Varianza de grados 
  #transitivity(net_pois,type="local") #Coeficiente de agrupamiento
  assign(paste("coef", redes[i], sep = "_"),transitivity(get(paste("grafo", redes[i], sep = "_")),type="global")) #Coeficiente de agrupamiento medio de un red 
  assign(paste("assor", redes[i], sep = "_"),assortativity_degree(get(paste("grafo", redes[i], sep = "_")))) #assortativity
}  #Termina el for

save(pois_100_1, pois_100_2, pois_100_3, pois_500_1, pois_500_2, pois_500_3, pois_1000_1, pois_1000_2, pois_1000_3, pois_10000_1, pois_10000_2, pois_10000_3, pois_50000_1, pois_50000_2, pois_50000_3,
     grafo_pois_100_1, grafo_pois_100_2, grafo_pois_100_3, grafo_pois_500_1, grafo_pois_500_2, grafo_pois_500_3, grafo_pois_1000_1, grafo_pois_1000_2, grafo_pois_1000_3,grafo_pois_10000_1, grafo_pois_10000_2, grafo_pois_10000_3, grafo_pois_50000_1, grafo_pois_50000_2, grafo_pois_50000_3,
     diam_pois_100_1, diam_pois_100_2, diam_pois_100_3, diam_pois_500_1, diam_pois_500_2, diam_pois_500_3, diam_pois_1000_1, diam_pois_1000_2, diam_pois_1000_3, diam_pois_10000_1, diam_pois_10000_2, diam_pois_10000_3, diam_pois_50000_1, diam_pois_50000_2, diam_pois_50000_3,
     med_pois_100_1, med_pois_100_2, med_pois_100_3, med_pois_500_1, med_pois_500_2, med_pois_500_3, med_pois_1000_1, med_pois_1000_2, med_pois_1000_3, med_pois_10000_1, med_pois_10000_2, med_pois_10000_3, med_pois_50000_1, med_pois_50000_2, med_pois_50000_3,
     var_pois_100_1, var_pois_100_2, var_pois_100_3, var_pois_500_1, var_pois_500_2, var_pois_500_3, var_pois_1000_1, var_pois_1000_2, var_pois_1000_3, var_pois_10000_1, var_pois_10000_2, var_pois_10000_3, var_pois_50000_1, var_pois_50000_2, var_pois_50000_3,
     coef_pois_100_1, coef_pois_100_2, coef_pois_100_3, coef_pois_500_1, coef_pois_500_2, coef_pois_500_3, coef_pois_1000_1, coef_pois_1000_2, coef_pois_1000_3, coef_pois_10000_1, coef_pois_10000_2, coef_pois_10000_3, coef_pois_50000_1, coef_pois_50000_2, coef_pois_50000_3,
     assor_pois_100_1, assor_pois_100_2, assor_pois_100_3, assor_pois_500_1, assor_pois_500_2, assor_pois_500_3, assor_pois_1000_1, assor_pois_1000_2, assor_pois_1000_3, assor_pois_10000_1, assor_pois_10000_2, assor_pois_10000_3, assor_pois_50000_1, assor_pois_50000_2, assor_pois_50000_3,
     file="pois_5difnodos_15sim.RData")
save(redes, nodos, file="Redes y Nodos Poisson.RData")


#############################################################################################


#creamos una red con grados distribuidos polylogaritmica

redespoly<-c("poly_100_1", "poly_100_2", "poly_100_3", "poly_500_1", "poly_500_2", "poly_500_3", "poly_1000_1", "poly_1000_2", "poly_1000_3", "poly_10000_1", "poly_10000_2", "poly_10000_3", "poly_50000_1", "poly_50000_2", "poly_50000_3")    
nodos<-c(100, 500, 1000, 10000, 50000)


for(p in 1:15){
  if (p<4) assign(redespoly[p], local.network(nodos[1], "poly.log",c(0.1,2.055)))
  else if (p<7) assign(redespoly[p], local.network(nodos[2], "poly.log",c(0.1,2.055)))
  else if (p<10) assign(redespoly[p], local.network(nodos[3], "poly.log",c(0.1,2.055)))
  else if (p<13) assign(redespoly[p], local.network(nodos[4], "poly.log",c(0.1,2.055)))
  else assign(redespoly[p], local.network(nodos[5], "poly.log",c(0.1,2.055)))
  
  assign(paste("grafo", redespoly[p], sep = "_"), graph_from_edgelist(get(redespoly[p])$edges))  
  
  #Cualidades de la red 
  assign(paste("diam", redespoly[p], sep = "_"),diameter(get(paste("grafo", redespoly[p], sep = "_")))) #Diametro 
  assign(paste("med", redespoly[p], sep = "_"),mean(get(redespoly[p])$degree)) #Media de grados 
  assign(paste("var", redespoly[p], sep = "_"),var(get(redespoly[p])$degree)) #Varianza de grados 
  assign(paste("coef", redespoly[p], sep = "_"),transitivity(get(paste("grafo", redespoly[p], sep = "_")),type="global")) #Coeficiente de agrupamiento medio de un red 
  assign(paste("assor", redespoly[p], sep = "_"),assortativity_degree(get(paste("grafo", redespoly[p], sep = "_")))) #assortativity
}  #Termina el for p

save(poly_100_1, poly_100_2, poly_100_3, poly_500_1, poly_500_2, poly_500_3, poly_1000_1, poly_1000_2, poly_1000_3, poly_10000_1, poly_10000_2, poly_10000_3, poly_50000_1, poly_50000_2, poly_50000_3,
     grafo_poly_100_1, grafo_poly_100_2, grafo_poly_100_3, grafo_poly_500_1, grafo_poly_500_2, grafo_poly_500_3, grafo_poly_1000_1, grafo_poly_1000_2, grafo_poly_1000_3,grafo_poly_10000_1, grafo_poly_10000_2, grafo_poly_10000_3, grafo_poly_50000_1, grafo_poly_50000_2, grafo_poly_50000_3,
     diam_poly_100_1, diam_poly_100_2, diam_poly_100_3, diam_poly_500_1, diam_poly_500_2, diam_poly_500_3, diam_poly_1000_1, diam_poly_1000_2, diam_poly_1000_3, diam_poly_10000_1, diam_poly_10000_2, diam_poly_10000_3, diam_poly_50000_1, diam_poly_50000_2, diam_poly_50000_3,
     med_poly_100_1, med_poly_100_2, med_poly_100_3, med_poly_500_1, med_poly_500_2, med_poly_500_3, med_poly_1000_1, med_poly_1000_2, med_poly_1000_3, med_poly_10000_1, med_poly_10000_2, med_poly_10000_3, med_poly_50000_1, med_poly_50000_2, med_poly_50000_3,
     var_poly_100_1, var_poly_100_2, var_poly_100_3, var_poly_500_1, var_poly_500_2, var_poly_500_3, var_poly_1000_1, var_poly_1000_2, var_poly_1000_3, var_poly_10000_1, var_poly_10000_2, var_poly_10000_3, var_poly_50000_1, var_poly_50000_2, var_poly_50000_3,
     coef_poly_100_1, coef_poly_100_2, coef_poly_100_3, coef_poly_500_1, coef_poly_500_2, coef_poly_500_3, coef_poly_1000_1, coef_poly_1000_2, coef_poly_1000_3, coef_poly_10000_1, coef_poly_10000_2, coef_poly_10000_3, coef_poly_50000_1, coef_poly_50000_2, coef_poly_50000_3,
     assor_poly_100_1, assor_poly_100_2, assor_poly_100_3, assor_poly_500_1, assor_poly_500_2, assor_poly_500_3, assor_poly_1000_1, assor_poly_1000_2, assor_poly_1000_3, assor_poly_10000_1, assor_poly_10000_2, assor_poly_10000_3, assor_poly_50000_1, assor_poly_50000_2, assor_poly_50000_3,
     file="poly_5difnodos_15sim.RData")
save(redespoly, nodos, file="Redes y Nodos Poly.RData")


#############################################################################################

#creamos una red con grados distribuidos watts strongtz

redesw<-c("ws_100_1", "ws_100_2", "ws_100_3", "ws_500_1", "ws_500_2", "ws_500_3", "ws_1000_1", "ws_1000_2", "ws_1000_3", "ws_10000_1", "ws_10000_2", "ws_10000_3", "ws_50000_1", "ws_50000_2", "ws_50000_3")    
redesws<-c("ws_100_1", "ws_100_2", "ws_100_3", "ws_500_1", "ws_500_2", "ws_500_3", "ws_1000_1", "ws_1000_2", "ws_1000_3", "ws_10000_1", "ws_10000_2", "ws_10000_3", "ws_50000_1", "ws_50000_2", "ws_50000_3")    
nodos<-c(100, 500, 1000, 10000, 50000)


for(w in 1:9){
  if (w<4) assign(redesw[w],  watts.strogatz.game(1,nodos[1], 2, 0.5))
  else if (w<7) assign(redesw[w], watts.strogatz.game(1,nodos[2], 2, 0.5))
  else if (w<10) assign(redesw[w], watts.strogatz.game(1,nodos[3], 2, 0.5))
  else if (w<13) assign(redesw[w], watts.strogatz.game(1,nodos[4], 2, 0.5))
  else assign(redesw[w], watts.strogatz.game(1,nodos[5], 2, 0.5))
  
  ####### En esta parte se hace la transformación para tener el mismo formato que las Poisson y Poly-log
  ed<-get.edgelist(get(redesw[w]))
  ed<-ed[order(ed[,1],ed[,2]),] 
  if (w<4) assign(redesws[w],list(edges=ed,degree=degree(get(redesw[w])),degree.left=NULL,n=nodos[1]))
  else if (w<7) assign(redesws[w],list(edges=ed,degree=degree(get(redesw[w])),degree.left=NULL,n=nodos[2]))
  else if (w<10) assign(redesws[w],list(edges=ed,degree=degree(get(redesw[w])),degree.left=NULL,n=nodos[3]))
  else if (w<13) assign(redesws[w],list(edges=ed,degree=degree(get(redesw[w])),degree.left=NULL,n=nodos[4]))
  else assign(redesws[w],list(edges=ed,degree=degree(get(redesw[w])),degree.left=NULL,n=nodos[5]))
  #######
  
  assign(paste("grafo", redesws[w], sep = "_"),graph_from_edgelist(get(redesw[w])$edges))
  
  #Cualidades de la red 
  assign(paste("diam", redesws[w], sep = "_"),diameter(get(paste("grafo", redesws[w], sep = "_")))) #Diametro 
  assign(paste("med", redesws[w], sep = "_"),mean(get(redesws[w])$degree)) #Media de grados 
  assign(paste("var", redesws[w], sep = "_"),var(get(redesws[w])$degree)) #Varianza de grados 
  assign(paste("coef", redesws[w], sep = "_"),transitivity(get(paste("grafo", redesws[w], sep = "_")),type="global")) #Coeficiente de agrupamiento medio de un red 
  assign(paste("assor", redesws[w], sep = "_"),assortativity_degree(get(paste("grafo", redesws[w], sep = "_")))) #assortativity
}#Termina el for p

save(ws_100_1, ws_100_2, ws_100_3, ws_500_1, ws_500_2, ws_500_3, ws_1000_1, ws_1000_2, ws_1000_3, ws_10000_1, ws_10000_2, ws_10000_3, ws_50000_1, ws_50000_2, ws_50000_3,
     grafo_ws_100_1, grafo_ws_100_2, grafo_ws_100_3, grafo_ws_500_1, grafo_ws_500_2, grafo_ws_500_3, grafo_ws_1000_1, grafo_ws_1000_2, grafo_ws_1000_3,grafo_ws_10000_1, grafo_ws_10000_2, grafo_ws_10000_3, grafo_ws_50000_1, grafo_ws_50000_2, grafo_ws_50000_3,
     diam_ws_100_1, diam_ws_100_2, diam_ws_100_3, diam_ws_500_1, diam_ws_500_2, diam_ws_500_3, diam_ws_1000_1, diam_ws_1000_2, diam_ws_1000_3, diam_ws_10000_1, diam_ws_10000_2, diam_ws_10000_3, diam_ws_50000_1, diam_ws_50000_2, diam_ws_50000_3,
     med_ws_100_1, med_ws_100_2, med_ws_100_3, med_ws_500_1, med_ws_500_2, med_ws_500_3, med_ws_1000_1, med_ws_1000_2, med_ws_1000_3, med_ws_10000_1, med_ws_10000_2, med_ws_10000_3, med_ws_50000_1, med_ws_50000_2, med_ws_50000_3,
     var_ws_100_1, var_ws_100_2, var_ws_100_3, var_ws_500_1, var_ws_500_2, var_ws_500_3, var_ws_1000_1, var_ws_1000_2, var_ws_1000_3, var_ws_10000_1, var_ws_10000_2, var_ws_10000_3, var_ws_50000_1, var_ws_50000_2, var_ws_50000_3,
     coef_ws_100_1, coef_ws_100_2, coef_ws_100_3, coef_ws_500_1, coef_ws_500_2, coef_ws_500_3, coef_ws_1000_1, coef_ws_1000_2, coef_ws_1000_3, coef_ws_10000_1, coef_ws_10000_2, coef_ws_10000_3, coef_ws_50000_1, coef_ws_50000_2, coef_ws_50000_3,
     assor_ws_100_1, assor_ws_100_2, assor_ws_100_3, assor_ws_500_1, assor_ws_500_2, assor_ws_500_3, assor_ws_1000_1, assor_ws_1000_2, assor_ws_1000_3, assor_ws_10000_1, assor_ws_10000_2, assor_ws_10000_3, assor_ws_50000_1, assor_ws_50000_2, assor_ws_50000_3,
     file="WS_5difnodos_15sim.RData")
save(redesws, nodos, file="Redes y Nodos WS.RData")




