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

# Oligo     Degrader      S-degrader     Opport      Predator
a <- 1.1;     b <- 0.2;   c <- 0.2;     d <- -0.2;    e <- -0.3   # Oligo (WT)
f <- 0.3;     g <- 0.3;   h <- 0.1;     i <- -0.1;     j <- -0.1   # Degrader (D)
k <- 0.8;     l <- 0.5;   m <- 0.7;     n <- -0.1;     o <- -0.1   # S-degrader (S)
p <- 0.1;   q <- 0.5;   r <- 0.5;    s <- 0.1;     t <- -0.5   # Opportunist (O)
u <- 0.8;     v <- 0.3;   w <- 0.4;     x <- 0.2;     y <- 0.01    # Predator (P)
 

oligo_fitness<-paste0(a,"*f_+",b,"*f_D+",c,"*f_D_S+",d,"*f_O+", e,"*f_O_P")
degrd_fitness<-paste0(f,"*f_+",g,"*f_D+",h,"*f_D_S+",i,"*f_O+", j,"*f_O_P")
sdegr_fitness<-paste0(k,"*f_+",l,"*f_D+",m,"*f_D_S+",n,"*f_O+", o,"*f_O_P")
oppor_fitness<-paste0(p,"*f_+",q,"*f_D+",r,"*f_D_S+",s,"*f_O+", t,"*f_O_P")
preda_fitness<-paste0(u,"*f_+",v,"*f_D+",w,"*f_D_S+",x,"*f_O+", y,"*f_O_P")

bloom=2

poset <- data.frame(parent = c("Root", "Root", "D", "O"),
                 child = c("D", "O", "S", "P"),
                 s = c(0.05, 0.02, 0.1, 0.02),
                 sh = c(rep(-60, 4)), #no dev allowed 
                 typeDep = "MN"
                 )

GFt = data.frame(
  Genotype = c("WT","D","S","O","P","D,S","O,P"),
  Fitness = c(oligo_fitness, #WT
              paste0("if (T>20 & T<50) ", bloom, "*(",degrd_fitness, ")",";  
                                           else 0.05*(", degrd_fitness, ");"), #D
              "0", #S
              paste0("if (T>20 & T<50) ", bloom, "*(",oppor_fitness, ")",";  
                                           else 0.05*(", oppor_fitness, ");"), #O
              "0", #P
              paste0("if (T>20 & T<55) ", bloom, "*(",sdegr_fitness, ")",";  
                                           else 0.05*(", sdegr_fitness, ");"), #DS
              preda_fitness)) #OP

tryal<- allFitnessEffects(poset,
                          genotFitness = GFt,
                          frequencyDependentFitness = TRUE,
                          frequencyType = "rel")

letse<- oncoSimulIndiv(
  tryal,
  model = "McFLD",
  onlyCancer = FALSE,
  finalTime = 100,
  mu = 0.1,
  initSize = 5000,
  keepPhylog = FALSE,
  seed = NULL
)

plot(letse, main="SINGLE-DOSE MIXED TREATMENT", show = "genotypes", type = "line",
     col=c("red","green","blue","black","yellow","pink"), xlim=c(-100,100), )

# FIRST SIMULATION: SINGLE-DOSE MIXED TREATMENT 

players <- data.frame(Genotype = c("WT","D","S","O","P","D,S","O,P"),
                      Fitness = c(wt_fitness, #WT
                                  paste0("if (T>50) ", bloom, "*(",ss_fitness, ")",";  
                                           else ", ss_fitness, ";"), #SS
                                  paste0("if (T>50) ", bloom, "*(",r1_fitness, ")",";  
                                           else ", r1_fitness, ";"), #R1
                                  paste0("if (T>50) ", bloom, "*(",r2_fitness, ")",";  
                                           else ", r2_fitness, ";"), #R2
                                  rr_fitness,rr_fitness,rr_fitness), #RR
                      stringsAsFactors = FALSE)

game <- allFitnessEffects(genotFitness = players,
                          frequencyDependentFitness = TRUE,
                          frequencyType = "rel")

SDMT <- oncoSimulIndiv(game,
                                 model = "McFLD",
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
                                  paste0("if (T>50) (cos(T)/1.5) *", drug1_eff*drug2_eff*10000, "*(",ss_fitness, ")",";  
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
