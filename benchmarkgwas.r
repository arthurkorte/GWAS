suppressPackageStartupMessages(require("doMC"))
setwd("/storage/evolgen/jan/Gwas")

#source("t2029_par.r")
source("t2029_gwas.r")

my_dorm <- read.csv("Dormancy.csv")
my_dorm <- my_dorm[,c(2,8)]

time_sc <- Sys.time()
results_sc <- t2029_gwas(my_dorm, u.max = 1, mc= FALSE)
jane_sc <- Sys.time() - time_sc

time_mc <- Sys.time()
results_mc <- t2029_gwas(my_dorm,u.max = 1, mult.core = TRUE)
jane_mc <- Sys.time() - time_mc

print(jane_mc)
print(jane_sc)

save(results_mc,jane_mc, file = 'bench_mc.rda')
save(results_sc,jane_sc, file= 'bench_sc.rda')


load("X_2029_1.rda")
X.rda <- X
X.rds <- readRDS("X_2029_1.rds")



time <- Sys.time()
load("X_2029_1.rda")
jane <- Sys.time() - time 




