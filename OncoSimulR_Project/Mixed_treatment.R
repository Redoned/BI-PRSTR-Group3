# Instalar #
'
if (!require("devtools"))
  install.packages("devtools")
library(devtools)
install_github("rdiaz02/OncoSimul/OncoSimulR", ref = "freq-dep-fitness")

library(OncoSimulR)
'

# Healthy     Sensitive   Resistance1  Resistance2   BothRes
a <- 0.85;    b <- 0.5;   c <- 0.5;    d <- 0.5;     e <- 0.5   # Healthy
f <- 0.9;     g <- 1.25;  h <- 0.8;    i <- 0.8;     j <- 0.5   # Sensitive
k <- 0.88;    l <- -0.3;  m <- 0.75;   n <- 0.75;    o <- 0.5   # Resistance1
p <- 0.8;     q <- -0.3;  r <- 0.75;   s <- 0.75;    t <- 0.8   # Resistance2
u <- 0.7;     v <- -0.5;  w <- 0.4;    x <- 0.4;     y <- 0.9   # Both

drug_eff<-0.001
wt_fitness<-paste0(a,"*f_+",b,"*f_S+",c,"*f_S_R+",d,"*f_S_T+", e,"*f_S_R_T")
ss_fitness<-paste0(f,"*f_+",g,"*f_S+",h,"*f_S_R+",i,"*f_S_T+", j,"*f_S_R_T")
r1_fitness<-paste0(k,"*f_+",l,"*f_S+",m,"*f_S_R+",n,"*f_S_T+", o,"*f_S_R_T")
r2_fitness<-paste0(p,"*f_+",q,"*f_S+",r,"*f_S_R+",s,"*f_S_T+", t,"*f_S_R_T")
rr_fitness<-paste0(u,"*f_+",v,"*f_S+",w,"*f_S_R+",x,"*f_S_T+", y,"*f_S_R_T")

players <- data.frame(Genotype = c("WT","S","R","T","R,T","S,R","S,T","S,R,T"),
                      Fitness = c(wt_fitness, #WT
                                  paste0("if (T>50) ", drug_eff*10, "*(",ss_fitness, ")",";  
                                           else ", ss_fitness, ";"), #SS
                                  "0","0","0",
                                  paste0("if (T>50) ", drug_eff, "*(",r1_fitness, ")",";  
                                           else ", r1_fitness, ";"), #R1
                                  paste0("if (T>50) ", drug_eff, "*(",r2_fitness, ")",";  
                                           else ", r2_fitness, ";"), #R2
                                  paste0("if (T>11000) ", drug_eff, "*(",rr_fitness, ")",";  
                                           else ", rr_fitness, ";")), #RR
                      stringsAsFactors = FALSE)

game <- allFitnessEffects(genotFitness = players,
                          frequencyDependentFitness = TRUE,
                          frequencyType = "rel")

gamesimul <- oncoSimulIndiv(game,
                                 model = "McFL",
                                 onlyCancer = FALSE,
                                 finalTime = 500,
                                 mu = 0.01,
                                 initSize = 5000,
                                 keepPhylog = FALSE,
                                 seed = NULL)

plot(gamesimul, show = "genotypes", type = "line",
     col=c("red","green","blue","black","yellow","pink"), ylim = c(20, 50000))


players <- data.frame(Genotype = c("WT","S","R","T","R,T","S,R","S,T","S,R,T"),
                      Fitness = c(wt_fitness, #WT
                                  paste0("if (T>50) ", drug_eff*10, "*(",ss_fitness, ")",";  
                                           else ", ss_fitness, ";"), #SS
                                  "0","0","0",
                                  paste0("if (T>50) ", drug_eff, "*(",r1_fitness, ")",";  
                                           else ", r1_fitness, ";"), #R1
                                  paste0("if (T>50) ", drug_eff, "*(",r2_fitness, ")",";  
                                           else ", r2_fitness, ";"), #R2
                                  paste0("if (T>11000) ", drug_eff, "*(",rr_fitness, ")",";  
                                           else ", rr_fitness, ";")), #RR
                      stringsAsFactors = FALSE)

game <- allFitnessEffects(genotFitness = players,
                          frequencyDependentFitness = TRUE,
                          frequencyType = "rel")

gamesimul <- oncoSimulIndiv(game,
                            model = "McFL",
                            onlyCancer = FALSE,
                            finalTime = 500,
                            mu = 0.01,
                            initSize = 5000,
                            keepPhylog = FALSE,
                            seed = NULL)

plot(gamesimul, show = "genotypes", type = "line",
     col=c("red","green","blue","black","yellow","pink"), ylim = c(20, 50000))
