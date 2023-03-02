
#############################################################################################

setwd("C:\\Users\\Diana\\Dropbox\\Diana\\CodigoR")
source("04_funciones_Redes_LR.R")
load("pois_5difnodos_15sim.RData")
load("Redes y Nodos Poisson.RData")
load("poly_5difnodos_15sim.RData")
load("Redes y Nodos Poly.RData")
load("WS_5difnodos_15sim.RData")
load("Redes y Nodos WS.RData")

ini_inf<-2       #InIciamos con 2 infectados
#parametros fijos
betreal<-0.03
gamreal<-0.01
lambdareal<-0.05
t_max<-10000 # Tiempo máximo de observación
t_repo<-24 # Tiempo entre los reportes (en este caso, en horas)

#############################################################################################


######################################################
# Simulando 5 brotes SIR y SEIR en la red Poisson

for(i in 1:15){
  for(j in 1:5) {
    x<-i*10+j
    set.seed(x)
    assign(paste("sim_SIR", redes[i], j ,sep = "_"),epidemic.sim(get(redes[i]),obs.time=t_max,ini.infected=0, ini.infective=ini_inf,BETA=betreal,distrib.inf="exp", GAMMA=gamreal))
    assign(paste("sim_SEIR", redes[i], j, sep = "_"), epidemic.sim(get(redes[i]), obs.time=t_max, ini.infected=0, ini.infective=ini_inf,BETA=betreal,distrib.lat="exp", LAMBDA=lambdareal, distrib.inf="exp", GAMMA=gamreal))
  } #Termina for j
} #Termina for i

save(sim_SIR_pois_100_1_1, sim_SIR_pois_100_1_2, sim_SIR_pois_100_1_3, sim_SIR_pois_100_1_4, sim_SIR_pois_100_1_5, 
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
     sim_SIR_pois_50000_3_1, sim_SIR_pois_50000_3_2, sim_SIR_pois_50000_3_3, sim_SIR_pois_50000_3_4, sim_SIR_pois_50000_3_5,
     file="2sim_SIR_pois.RData")

save(sim_SEIR_pois_100_1_1, sim_SEIR_pois_100_1_2, sim_SEIR_pois_100_1_3, sim_SEIR_pois_100_1_4, sim_SEIR_pois_100_1_5, 
     sim_SEIR_pois_100_2_1, sim_SEIR_pois_100_2_2, sim_SEIR_pois_100_2_3, sim_SEIR_pois_100_2_4, sim_SEIR_pois_100_2_5,
     sim_SEIR_pois_100_3_1, sim_SEIR_pois_100_3_2, sim_SEIR_pois_100_3_3, sim_SEIR_pois_100_3_4, sim_SEIR_pois_100_3_5,
     sim_SEIR_pois_500_1_1, sim_SEIR_pois_500_1_2, sim_SEIR_pois_500_1_3, sim_SEIR_pois_500_1_4, sim_SEIR_pois_500_1_5, 
     sim_SEIR_pois_500_2_1, sim_SEIR_pois_500_2_2, sim_SEIR_pois_500_2_3, sim_SEIR_pois_500_2_4, sim_SEIR_pois_500_2_5,
     sim_SEIR_pois_500_3_1, sim_SEIR_pois_500_3_2, sim_SEIR_pois_500_3_3, sim_SEIR_pois_500_3_4, sim_SEIR_pois_500_3_5,
     sim_SEIR_pois_1000_1_1, sim_SEIR_pois_1000_1_2, sim_SEIR_pois_1000_1_3, sim_SEIR_pois_1000_1_4, sim_SEIR_pois_1000_1_5, 
     sim_SEIR_pois_1000_2_1, sim_SEIR_pois_1000_2_2, sim_SEIR_pois_1000_2_3, sim_SEIR_pois_1000_2_4, sim_SEIR_pois_1000_2_5,
     sim_SEIR_pois_1000_3_1, sim_SEIR_pois_1000_3_2, sim_SEIR_pois_1000_3_3, sim_SEIR_pois_1000_3_4, sim_SEIR_pois_1000_3_5,
     sim_SEIR_pois_10000_1_1, sim_SEIR_pois_10000_1_2, sim_SEIR_pois_10000_1_3, sim_SEIR_pois_10000_1_4, sim_SEIR_pois_10000_1_5, 
     sim_SEIR_pois_10000_2_1, sim_SEIR_pois_10000_2_2, sim_SEIR_pois_10000_2_3, sim_SEIR_pois_10000_2_4, sim_SEIR_pois_10000_2_5,
     sim_SEIR_pois_10000_3_1, sim_SEIR_pois_10000_3_2, sim_SEIR_pois_10000_3_3, sim_SEIR_pois_10000_3_4, sim_SEIR_pois_10000_3_5,
     sim_SEIR_pois_50000_1_1, sim_SEIR_pois_50000_1_2, sim_SEIR_pois_50000_1_3, sim_SEIR_pois_50000_1_4, sim_SEIR_pois_50000_1_5, 
     sim_SEIR_pois_50000_2_1, sim_SEIR_pois_50000_2_2, sim_SEIR_pois_50000_2_3, sim_SEIR_pois_50000_2_4, sim_SEIR_pois_50000_2_5,
     sim_SEIR_pois_50000_3_1, sim_SEIR_pois_50000_3_2, sim_SEIR_pois_50000_3_3, sim_SEIR_pois_50000_3_4, sim_SEIR_pois_50000_3_5,
     file="2sim_SEIR_pois.RData")
#############################################################################################

######################################################
# Simulando 5 brotes SIR y SEIR en la red Poly  

for(i in 1:15){
  for(j in 1:5) {
    x<-200+i*10+j
    set.seed(x)
    assign(paste("sim_SIR", redespoly[i], j ,sep = "_"),epidemic.sim(get(redespoly[i]),obs.time=t_max,ini.infected=0, ini.infective=ini_inf,BETA=betreal,distrib.inf="exp", GAMMA=gamreal))
    assign(paste("sim_SEIR", redespoly[i], j, sep = "_"), epidemic.sim(get(redespoly[i]), obs.time=t_max, ini.infected=0, ini.infective=ini_inf,BETA=betreal,distrib.lat="exp", LAMBDA=lambdareal, distrib.inf="exp", GAMMA=gamreal))
  } #Termina for j
} #Termina for i

save(sim_SIR_poly_100_1_1, sim_SIR_poly_100_1_2, sim_SIR_poly_100_1_3, sim_SIR_poly_100_1_4, sim_SIR_poly_100_1_5, 
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
     sim_SIR_poly_50000_3_1, sim_SIR_poly_50000_3_2, sim_SIR_poly_50000_3_3, sim_SIR_poly_50000_3_4, sim_SIR_poly_50000_3_5,
     file="2sim_SIR_poly.RData")

save(sim_SEIR_poly_100_1_1, sim_SEIR_poly_100_1_2, sim_SEIR_poly_100_1_3, sim_SEIR_poly_100_1_4, sim_SEIR_poly_100_1_5, 
     sim_SEIR_poly_100_2_1, sim_SEIR_poly_100_2_2, sim_SEIR_poly_100_2_3, sim_SEIR_poly_100_2_4, sim_SEIR_poly_100_2_5,
     sim_SEIR_poly_100_3_1, sim_SEIR_poly_100_3_2, sim_SEIR_poly_100_3_3, sim_SEIR_poly_100_3_4, sim_SEIR_poly_100_3_5,
     sim_SEIR_poly_500_1_1, sim_SEIR_poly_500_1_2, sim_SEIR_poly_500_1_3, sim_SEIR_poly_500_1_4, sim_SEIR_poly_500_1_5, 
     sim_SEIR_poly_500_2_1, sim_SEIR_poly_500_2_2, sim_SEIR_poly_500_2_3, sim_SEIR_poly_500_2_4, sim_SEIR_poly_500_2_5,
     sim_SEIR_poly_500_3_1, sim_SEIR_poly_500_3_2, sim_SEIR_poly_500_3_3, sim_SEIR_poly_500_3_4, sim_SEIR_poly_500_3_5,
     sim_SEIR_poly_1000_1_1, sim_SEIR_poly_1000_1_2, sim_SEIR_poly_1000_1_3, sim_SEIR_poly_1000_1_4, sim_SEIR_poly_1000_1_5, 
     sim_SEIR_poly_1000_2_1, sim_SEIR_poly_1000_2_2, sim_SEIR_poly_1000_2_3, sim_SEIR_poly_1000_2_4, sim_SEIR_poly_1000_2_5,
     sim_SEIR_poly_1000_3_1, sim_SEIR_poly_1000_3_2, sim_SEIR_poly_1000_3_3, sim_SEIR_poly_1000_3_4, sim_SEIR_poly_1000_3_5,
     sim_SEIR_poly_10000_1_1, sim_SEIR_poly_10000_1_2, sim_SEIR_poly_10000_1_3, sim_SEIR_poly_10000_1_4, sim_SEIR_poly_10000_1_5, 
     sim_SEIR_poly_10000_2_1, sim_SEIR_poly_10000_2_2, sim_SEIR_poly_10000_2_3, sim_SEIR_poly_10000_2_4, sim_SEIR_poly_10000_2_5,
     sim_SEIR_poly_10000_3_1, sim_SEIR_poly_10000_3_2, sim_SEIR_poly_10000_3_3, sim_SEIR_poly_10000_3_4, sim_SEIR_poly_10000_3_5,
     sim_SEIR_poly_50000_1_1, sim_SEIR_poly_50000_1_2, sim_SEIR_poly_50000_1_3, sim_SEIR_poly_50000_1_4, sim_SEIR_poly_50000_1_5, 
     sim_SEIR_poly_50000_2_1, sim_SEIR_poly_50000_2_2, sim_SEIR_poly_50000_2_3, sim_SEIR_poly_50000_2_4, sim_SEIR_poly_50000_2_5,
     sim_SEIR_poly_50000_3_1, sim_SEIR_poly_50000_3_2, sim_SEIR_poly_50000_3_3, sim_SEIR_poly_50000_3_4, sim_SEIR_poly_50000_3_5,
     file="2sim_SEIR_poly.RData")
#############################################################################################

######################################################
# Simulando 5 brotes SIR y SEIR en la red W-S

for(i in 1:15){
  for(j in 1:5) {
    x<-400+i*10+j
    set.seed(x)
    assign(paste("sim_SIR", redesws[i], j ,sep = "_"),epidemic.sim(get(redesws[i]),obs.time=t_max,ini.infected=0, ini.infective=ini_inf,BETA=betreal,distrib.inf="exp", GAMMA=gamreal))
    assign(paste("sim_SEIR", redesws[i], j, sep = "_"), epidemic.sim(get(redesws[i]), obs.time=t_max, ini.infected=0, ini.infective=ini_inf,BETA=betreal,distrib.lat="exp", LAMBDA=lambdareal, distrib.inf="exp", GAMMA=gamreal))
  } #Termina for j
} #Termina for i

save(sim_SIR_ws_100_1_1, sim_SIR_ws_100_1_2, sim_SIR_ws_100_1_3, sim_SIR_ws_100_1_4, sim_SIR_ws_100_1_5, 
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
     sim_SIR_ws_50000_3_1, sim_SIR_ws_50000_3_2, sim_SIR_ws_50000_3_3, sim_SIR_ws_50000_3_4, sim_SIR_ws_50000_3_5,
     file="1sim_SIR_WS.RData")

save(sim_SEIR_ws_100_1_1, sim_SEIR_ws_100_1_2, sim_SEIR_ws_100_1_3, sim_SEIR_ws_100_1_4, sim_SEIR_ws_100_1_5, 
     sim_SEIR_ws_100_2_1, sim_SEIR_ws_100_2_2, sim_SEIR_ws_100_2_3, sim_SEIR_ws_100_2_4, sim_SEIR_ws_100_2_5,
     sim_SEIR_ws_100_3_1, sim_SEIR_ws_100_3_2, sim_SEIR_ws_100_3_3, sim_SEIR_ws_100_3_4, sim_SEIR_ws_100_3_5,
     sim_SEIR_ws_500_1_1, sim_SEIR_ws_500_1_2, sim_SEIR_ws_500_1_3, sim_SEIR_ws_500_1_4, sim_SEIR_ws_500_1_5, 
     sim_SEIR_ws_500_2_1, sim_SEIR_ws_500_2_2, sim_SEIR_ws_500_2_3, sim_SEIR_ws_500_2_4, sim_SEIR_ws_500_2_5,
     sim_SEIR_ws_500_3_1, sim_SEIR_ws_500_3_2, sim_SEIR_ws_500_3_3, sim_SEIR_ws_500_3_4, sim_SEIR_ws_500_3_5,
     sim_SEIR_ws_1000_1_1, sim_SEIR_ws_1000_1_2, sim_SEIR_ws_1000_1_3, sim_SEIR_ws_1000_1_4, sim_SEIR_ws_1000_1_5, 
     sim_SEIR_ws_1000_2_1, sim_SEIR_ws_1000_2_2, sim_SEIR_ws_1000_2_3, sim_SEIR_ws_1000_2_4, sim_SEIR_ws_1000_2_5,
     sim_SEIR_ws_1000_3_1, sim_SEIR_ws_1000_3_2, sim_SEIR_ws_1000_3_3, sim_SEIR_ws_1000_3_4, sim_SEIR_ws_1000_3_5,
     sim_SEIR_ws_10000_1_1, sim_SEIR_ws_10000_1_2, sim_SEIR_ws_10000_1_3, sim_SEIR_ws_10000_1_4, sim_SEIR_ws_10000_1_5, 
     sim_SEIR_ws_10000_2_1, sim_SEIR_ws_10000_2_2, sim_SEIR_ws_10000_2_3, sim_SEIR_ws_10000_2_4, sim_SEIR_ws_10000_2_5,
     sim_SEIR_ws_10000_3_1, sim_SEIR_ws_10000_3_2, sim_SEIR_ws_10000_3_3, sim_SEIR_ws_10000_3_4, sim_SEIR_ws_10000_3_5,
     sim_SEIR_ws_50000_1_1, sim_SEIR_ws_50000_1_2, sim_SEIR_ws_50000_1_3, sim_SEIR_ws_50000_1_4, sim_SEIR_ws_50000_1_5, 
     sim_SEIR_ws_50000_2_1, sim_SEIR_ws_50000_2_2, sim_SEIR_ws_50000_2_3, sim_SEIR_ws_50000_2_4, sim_SEIR_ws_50000_2_5,
     sim_SEIR_ws_50000_3_1, sim_SEIR_ws_50000_3_2, sim_SEIR_ws_50000_3_3, sim_SEIR_ws_50000_3_4, sim_SEIR_ws_50000_3_5,
     file="1sim_SEIR_WS.RData")
#############################################################################################


###### Aqui se guardadan todos los brotes juntos
###### Brotes_SIR_Pois_Poly_WS.RData   contiene   2sim_SIR_pois.RData, 2sim_SIR_poly.RData, 1sim_SIR_WS.RData
###### Brotes_SEIR_Pois_Poly_WS.RData  contiene   2sim_SEIR_pois.RData, 2sim_SEIR_poly.RData, 1sim_SEIR_WS.RData 