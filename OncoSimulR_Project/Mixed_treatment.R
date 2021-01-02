# Instalar la última versión (≥2.99.8)
'
if (!require("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("OncoSimulR", version = "3.11")

if (!require("devtools"))
  install.packages("devtools")
library(devtools)
install_github("rdiaz02/OncoSimul/OncoSimulR", ref = "freq-dep-fitness")
'
library(OncoSimulR)

# Healthy   Sensitive   Resistance1  Resistance2   BothRes
a <- 1;     b <- 0.5;   c <- 0.5;    d <- 0.5;     e <- 0.5   # Healthy
f <- 1.25;  g <- 0.7;   h <- 1;      i <- 1;       j <- 1.1   # Sensitive
k <- 1.1;   l <- -0.2;  m <- 0.8;    n <- 0.8;     o <- 1     # Resistance1 (R)
p <- 1.1;   q <- -0.2;  r <- 0.8;    s <- 0.8;     t <- 1     # Resistance2 (T)
u <- 1.05;   v <- -0.3;  w <- 0.4;   x <- 0.4;     y <- 0.9   # BothRes
 

wt_fitness<-paste0(a,"*f_+",b,"*f_S+",c,"*f_S_R+",d,"*f_S_T+", e,"*f_S_R_T")
ss_fitness<-paste0(f,"*f_+",g,"*f_S+",h,"*f_S_R+",i,"*f_S_T+", j,"*f_S_R_T")
r1_fitness<-paste0(k,"*f_+",l,"*f_S+",m,"*f_S_R+",n,"*f_S_T+", o,"*f_S_R_T")
r2_fitness<-paste0(p,"*f_+",q,"*f_S+",r,"*f_S_R+",s,"*f_S_T+", t,"*f_S_R_T")
rr_fitness<-paste0(u,"*f_+",v,"*f_S+",w,"*f_S_R+",x,"*f_S_T+", y,"*f_S_R_T")

# FIRST SIMULATION: SINGLE-DOSE MIXED TREATMENT 
drug1_eff<-0.01
drug2_eff<-0.01

players <- data.frame(Genotype = c("WT","S","R","T","R,T","S,R","S,T","S,R,T"),
                      Fitness = c(wt_fitness, #WT
                                  paste0("if (T>50) ", drug1_eff*drug2_eff, "*(",ss_fitness, ")",";  
                                           else ", ss_fitness, ";"), #SS
                                  "0","0","0",
                                  paste0("if (T>50) ", drug1_eff, "*(",r1_fitness, ")",";  
                                           else ", r1_fitness, ";"), #R1
                                  paste0("if (T>50) ", drug2_eff, "*(",r2_fitness, ")",";  
                                           else ", r2_fitness, ";"), #R2
                                  rr_fitness), #RR
                      stringsAsFactors = FALSE)

game <- allFitnessEffects(genotFitness = players,
                          frequencyDependentFitness = TRUE,
                          frequencyType = "rel")

SDMT <- oncoSimulIndiv(game,
                                 model = "McFL",
                                 onlyCancer = FALSE,
                                 finalTime = 300,
                                 mu = 0.01,
                                 initSize = 5000,
                                 keepPhylog = FALSE,
                                 seed = NULL)
plot(SDMT, main="SINGLE-DOSE MIXED TREATMENT", show = "genotypes", type = "line",
     col=c("red","green","blue","black","yellow","pink"), ylim = c(20, 8000))
# SECOND SIMULATION: SINGLE-DOSE TREATMENT 1 (R is resistant)
drug1_eff<-0.1
drug2_eff<-1

players <- data.frame(Genotype = c("WT","S","R","T","R,T","S,R","S,T","S,R,T"),
                      Fitness = c(wt_fitness, #WT
                                  paste0("if (T>50) ", drug1_eff*drug2_eff, "*(",ss_fitness, ")",";  
                                           else ", ss_fitness, ";"), #SS
                                  "0","0","0",
                                  paste0("if (T>50) ", drug1_eff, "*(",r1_fitness, ")",";  
                                           else ", r1_fitness, ";"), #R1
                                  paste0("if (T>50) ", drug2_eff, "*(",r2_fitness, ")",";  
                                           else ", r2_fitness, ";"), #R2
                                  rr_fitness), #RR
                      stringsAsFactors = FALSE)

game <- allFitnessEffects(genotFitness = players,
                          frequencyDependentFitness = TRUE,
                          frequencyType = "rel")
SD1T <- oncoSimulIndiv(game,
                       model = "McFL",
                       onlyCancer = FALSE,
                       finalTime = 300,
                       mu = 0.01,
                       initSize = 5000,
                       keepPhylog = FALSE,
                       seed = NULL)
plot(SD1T, main="SINGLE-DOSE TREATMENT 1", show = "genotypes", type = "line",
     col=c("red","green","blue","black","yellow","pink"), ylim = c(20, 8000))
# THIRD SIMULATION: SINGLE-DOSE TREATMENT 2 (T is resistant)
drug1_eff<-1
drug2_eff<-0.1

players <- data.frame(Genotype = c("WT","S","R","T","R,T","S,R","S,T","S,R,T"),
                      Fitness = c(wt_fitness, #WT
                                  paste0("if (T>50) ", drug1_eff*drug2_eff, "*(",ss_fitness, ")",";  
                                           else ", ss_fitness, ";"), #SS
                                  "0","0","0",
                                  paste0("if (T>50) ", drug1_eff, "*(",r1_fitness, ")",";  
                                           else ", r1_fitness, ";"), #R1
                                  paste0("if (T>50) ", drug2_eff, "*(",r2_fitness, ")",";  
                                           else ", r2_fitness, ";"), #R2
                                  rr_fitness), #RR
                      stringsAsFactors = FALSE)

game <- allFitnessEffects(genotFitness = players,
                          frequencyDependentFitness = TRUE,
                          frequencyType = "rel")

SD2T <- oncoSimulIndiv(game,
                       model = "McFL",
                       onlyCancer = FALSE,
                       finalTime = 300,
                       mu = 0.01,
                       initSize = 5000,
                       keepPhylog = FALSE,
                       seed = NULL)
plot(SD2T, main="SINGLE-DOSE TREATMENT 2", show = "genotypes", type = "line",
     col=c("red","green","blue","black","yellow","pink"), ylim = c(20, 8000))

# FOURTH SIMULATION: ADJUSTED-DOSE MIXED TREATMENT 
drug1_eff<-0.01
drug2_eff<-0.01

players <- data.frame(Genotype = c("WT","S","R","T","R,T","S,R","S,T","S,R,T"),
                      Fitness = c(wt_fitness, #WT
                                  paste0("if (T>50) (cos(T)/1.5) *", drug1_eff*drug2_eff*10, "*(",ss_fitness, ")",";  
                                           else ", ss_fitness, ";"), #SS
                                  "0","0","0",
                                  paste0("if (T>50) (cos(T)/1.5) *", drug1_eff, "*(",r1_fitness, ")",";  
                                           else ", r1_fitness, ";"), #R1
                                  paste0("if (T>50) (cos(T)/1.5) *", drug2_eff, "*(",r2_fitness, ")",";  
                                           else ", r2_fitness, ";"), #R2
                                  rr_fitness), #RR
                      stringsAsFactors = FALSE)

game <- allFitnessEffects(genotFitness = players,
                          frequencyDependentFitness = TRUE,
                          frequencyType = "rel")
#plotFitnessLandscape(evalAllGenotypes(game, spPopSizes = c(
#  WT=5000,S=1000,R=500,T=500,"R, T"=2000, "S, R"=2000, "S, T"=2000, "S, R, T"=100)))

ADMT <- oncoSimulIndiv(game,
                       model = "McFL",
                       onlyCancer = FALSE,
                       finalTime = 300,
                       mu = 0.01,
                       initSize = 5000,
                       keepPhylog = FALSE,
                       seed = NULL)
plot(ADMT, main="ADJUSTED-DOSE MIXED TREATMENT", show = "genotypes", type = "line",
     col=c("red","green","blue","black","yellow","pink"), ylim = c(20, 8000))
