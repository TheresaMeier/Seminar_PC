# Seminar Modelling under Dependence
# Simulation Study

# required libraries
library(pcalg)
library(ParallelPC)
library(parallel)
library(ggplot2)
library(pROC)
library(scales)
library(ggh4x)
set.seed(1)

#################### varying dimensional setting ###############################
## parameter setting for simulating data
n <- 1000  # number of samples
ens <- 5  # expected neighborhood size (p-1)*s = ens
p.vars <- c(10, 100, 250, 500)   # number of nodes

# tables for evaluation
runtime_p <- as.data.frame(matrix(data = 0, nrow = 3*length(p.vars), ncol = 3))
colnames(runtime_p) <- c("alg","p", "runtime")
runtime_p$alg <- rep(c("original", "stable", "parallel"), each = length(p.vars))
runtime_p$p <- rep(p.vars, 3)

eval_p <- as.data.frame(matrix(data = 0, nrow = 3*length(p.vars), ncol = 9))
colnames(eval_p) = c("alg", "p","SHD", "TPR", "FPR", "TDR", "citests", 
                     "num_edges_true", "num_edges_pc")
eval_p$alg <- rep(c("original", "stable", "parallel"), each = length(p.vars))
eval_p$p <- rep(p.vars, 3)

j = 1
num <- 10    # number of simulated graphs to take the mean in the end

while(j < num + 1) {
  for (dim in p.vars) {
    s = ens / (dim-1)          # sparsity parameter (the lower the more sparse) 
    g <- randomDAG(dim,s)                  # true DAG
    d <- rmvDAG(n,g, errDist = "normal")   # generate random samples from true DAG
    
    print(dim)
    for(alg in c("original", "stable", "parallel")) {
      print(alg)
      if(alg == "parallel"){
        start_time <- Sys.time()
        pc.fit <- pc_parallel(list(C=cor(d), n = n), indepTest = gaussCItest,
                              p = dim, alpha = 0.01, skel.method = "parallel", 
                              num.cores = 4)
        end_time <- Sys.time()
      }
      else {
        start_time <- Sys.time()
        pc.fit <- pc_stable(list(C=cor(d), n = n), indepTest = gaussCItest,
                            p = dim, alpha = 0.01, skel.method = alg)
        end_time <- Sys.time()
      }  
      print(difftime(end_time, start_time, units = "secs"))
      runtime_p[runtime_p$alg == alg & runtime_p$p == dim, ]$runtime <- 
        runtime_p[runtime_p$alg == alg & runtime_p$p == dim, ]$runtime + 
        difftime(end_time, start_time, units = "secs")
      
      eval_p[eval_p$alg == alg & eval_p$p == dim,]$SHD <- 
        eval_p[eval_p$alg == alg & eval_p$p == dim,]$SHD + shd(g,pc.fit)
      eval_p[eval_p$alg == alg & eval_p$p == dim,]$TPR <- 
        eval_p[eval_p$alg == alg & eval_p$p == dim,]$TPR + compareGraphs(g, pc.fit@graph)[1]
      eval_p[eval_p$alg == alg & eval_p$p == dim,]$FPR <- 
        eval_p[eval_p$alg == alg & eval_p$p == dim,]$FPR + compareGraphs(g, pc.fit@graph)[2]
      eval_p[eval_p$alg == alg & eval_p$p == dim,]$TDR <- 
        eval_p[eval_p$alg == alg & eval_p$p == dim,]$TDR + compareGraphs(g, pc.fit@graph)[3]
      eval_p[eval_p$alg == alg & eval_p$p == dim,]$citests <- 
        eval_p[eval_p$alg == alg & eval_p$p == dim,]$citests + sum(pc.fit@n.edgetests)
      eval_p[eval_p$alg == alg & eval_p$p == dim,]$num_edges_true <- 
        eval_p[eval_p$alg == alg & eval_p$p == dim,]$num_edges_true + length(g@edgeData@data)
      eval_p[eval_p$alg == alg & eval_p$p == dim,]$num_edges_pc <- 
        eval_p[eval_p$alg == alg & eval_p$p == dim,]$num_edges_pc + length(pc.fit@graph@edgeData@data)
    }
  }
  
  j = j + 1
  print(j)
}

eval_p[,c(3:9)] <- eval_p[,c(3:9)] / num
runtime_p$runtime <- runtime_p$runtime / num
runtime_p$runtime_citests <- runtime_p$runtime / eval_p$citests

write.table(runtime_p, "runtime_p.txt")
write.table(eval_p, "eval_p.txt")

############################ varying sample size ###############################
## parameter setting for simulating data
p.vars <- c(10,100)   # number of random variables
ens <- 5  # expected neighborhood size (p-1)*s = ens
n.vars <- c(50,100,500,1000,10000,20000)     # different sample sizes

# tables for evaluation
runtime_n <- as.data.frame(matrix(data = 0,nrow = 2 * 3 *length(n.vars), ncol = 4))
colnames(runtime_n) <- c("alg","n", "p", "runtime")
runtime_n$alg <- rep(c("original", "stable", "parallel"), each = 2*length(n.vars))
runtime_n$n <- rep(n.vars, length(n.vars))
runtime_n$p <- rep(p.vars, each = length(n.vars))

eval_n <- as.data.frame(matrix(data = 0, nrow = 2*3*length(n.vars), ncol = 10))
colnames(eval_n) = c("alg", "n","p","SHD", "TPR", "FPR", "TDR", "citests", 
                     "num_edges_true", "num_edges_pc")
eval_n$alg <- rep(c("original", "stable", "parallel"), each = 2*length(n.vars))
eval_n$n <- rep(n.vars, length(n.vars))
eval_n$p <- rep(p.vars, each = length(n.vars))
j=1
num = 10

while(j < num + 1){
  for (p in p.vars){
    s = ens / (p-1)       # sparsity parameter (the lower the more sparse) 
    g <- randomDAG(p,s)   # true DAG
    print(j)
    for (size in n.vars){
      d <- rmvDAG(size,g, errDist = "normal")  # generate random samples from true DAG
      print(size)
      for (alg in c("original", "stable", "parallel")) {
        print(alg)
        if (alg == "parallel")
        {
          start_time <- Sys.time()
          pc.fit <- pc_parallel(list(C=cor(d), n = size), indepTest = gaussCItest, p = p,
                                alpha = 0.01, skel.method = "parallel", num.cores = 6)
          end_time <- Sys.time()
        }
        else
        {
          start_time <- Sys.time()
          pc.fit <- pc_stable(list(C=cor(d), n = size), indepTest = gaussCItest, p = p,
                              alpha = 0.01, skel.method = alg)
          end_time <- Sys.time()
        }
        print(difftime(end_time, start_time, units = "secs"))
        runtime_n[runtime_n$alg == alg & runtime_n$n == size & eval_n$p == p,]$runtime <- 
          runtime_n[runtime_n$alg == alg & runtime_n$n == size & eval_n$p == p,]$runtime + difftime(end_time, start_time, units = "secs")
        eval_n[eval_n$alg == alg & eval_n$n == size & eval_n$p == p,]$SHD <- 
          eval_n[eval_n$alg == alg & eval_n$n == size & eval_n$p == p,]$SHD + shd(g, pc.fit)
        eval_n[eval_n$alg == alg & eval_n$n == size & eval_n$p == p,]$TPR <- 
          eval_n[eval_n$alg == alg & eval_n$n == size & eval_n$p == p,]$TPR + compareGraphs(g, pc.fit@graph)[1]
        eval_n[eval_n$alg == alg & eval_n$n == size & eval_n$p == p,]$FPR <-
          eval_n[eval_n$alg == alg & eval_n$n == size & eval_n$p == p,]$FPR + compareGraphs(g, pc.fit@graph)[2]
        eval_n[eval_n$alg == alg & eval_n$n == size & eval_n$p == p,]$TDR <- eval_n[eval_n$alg == alg & eval_n$n == size & eval_n$p == p,]$TDR + compareGraphs(g, pc.fit@graph)[3]
        eval_n[eval_n$alg == alg & eval_n$n == size & eval_n$p == p,]$num_edges_true <- 
          eval_n[eval_n$alg == alg & eval_n$n == size & eval_n$p == p,]$num_edges_true + length(g@edgeData@data)
        eval_n[eval_n$alg == alg & eval_n$n == size & eval_n$p == p,]$num_edges_pc <- 
          eval_n[eval_n$alg == alg & eval_n$n == size & eval_n$p == p,]$num_edges_pc + length(pc.fit@graph@edgeData@data)
        eval_n[eval_n$alg == alg & eval_n$n == size & eval_n$p == p,]$citests <- 
          eval_n[eval_n$alg == alg & eval_n$n == size & eval_n$p == p,]$citests + sum(pc.fit@n.edgetests)
      }
    }
  }
    j = j + 1
}

runtime_n$runtime <- runtime_n$runtime / num
eval_n[,c(4:10)] <- eval_n[,c(4:10)] /num

write.table(eval_n, "eval_n.txt")
write.table(runtime_n, "runtime_n.txt")

##################### different distributional assumptions #####################
p.vars <- c(10, 100)    # number of random variables
n.vars <- c(100, 1000)  # number of samples
ens <- 5  # sparsity parameter (the lower the more sparse)
dist.vals = c("normal", "cauchy", "t4")

eval_dist<- as.data.frame(matrix(data = 0, nrow = 3*4, ncol = 7))
colnames(eval_dist) = c("dist","p", "n", "SHD", "TPR", "FPR", "TDR")
eval_dist$dist <- rep(dist.vals, each = 4)
eval_dist$p <- rep(p.vars, each = 2)
eval_dist$n <- rep(n.vars, 3)

j = 1
num = 10
while (j < num + 1){
  for (p in p.vars) {
    s = ens / (p-1)          # sparsity parameter (the lower the more sparse) 
    g <- randomDAG(p,s)      # true DAG
    
    for (size in n.vars) {
      for (dist in dist.vals) {
        d <- rmvDAG(size,g, errDist = dist)   # generate random samples from true DAG
        
        pc.fit <- pc_parallel(list(C=cor(d), n = size), indepTest = gaussCItest, p = p,
                          alpha = 0.01, skel.method = "parallel", num.cores = 4)
        
        eval_dist[eval_dist$dist == dist & eval_dist$p == p & eval_dist$n == size,]$SHD <- 
          eval_dist[eval_dist$dist == dist & eval_dist$p == p & eval_dist$n == size,]$SHD + shd(g, pc.fit)
        eval_dist[eval_dist$dist == dist & eval_dist$p == p & eval_dist$n == size,]$TPR <- 
          eval_dist[eval_dist$dist == dist & eval_dist$p == p & eval_dist$n == size,]$TPR + compareGraphs(g, pc.fit@graph)[1]
        eval_dist[eval_dist$dist == dist & eval_dist$p == p & eval_dist$n == size,]$FPR <- 
          eval_dist[eval_dist$dist == dist & eval_dist$p == p & eval_dist$n == size,]$FPR + compareGraphs(g, pc.fit@graph)[2]
        eval_dist[eval_dist$dist == dist & eval_dist$p == p & eval_dist$n == size,]$TDR <- 
          eval_dist[eval_dist$dist == dist & eval_dist$p == p & eval_dist$n == size,]$TDR + compareGraphs(g, pc.fit@graph)[3]
        }
    }
  }
  j = j + 1
}

eval_dist[,c(4:7)] = eval_dist[,c(4:7)]/num
eval_dist$szenario <- rep(c(1,2,3,4), 3)
eval_dist$dist <- ifelse(eval_dist$dist == "t4", "t (k=4)", eval_dist$dist)

write.table(eval_dist, "eval_dist.txt")

##################### different number of cores (parallel) #####################
p <- 100   # number of random variables
n <- 1000  # number of samples
ens <- 30   # sparsity parameter (the lower the more sparse)
cores.vars = c(2:6)

# Evaluation table
runtime_core <- as.data.frame(matrix(data = 0, nrow = length(cores.vars)+1, ncol = 3))
colnames(runtime_core) <- c("cores", "runtime", "citests")
runtime_core$cores <-c(1,cores.vars)

j=1
num = 10
while (j < num + 1){
  s = ens / (p-1)         # sparsity parameter (the lower the more sparse) 
  g <- randomDAG(p,s)     # true DAG
  d <- rmvDAG(n,g)        # generate random samples from true DAG
  
  start_time <- Sys.time()
  pc.parallel <- pc_stable(list(C=cor(d), n = n), indepTest = gaussCItest, p = p,
                             alpha = 0.01, skel.method = "stable")
  end_time <- Sys.time()
  
  runtime_core[runtime_core$cores == 1,]$runtime <- 
    runtime_core[runtime_core$cores == 1,]$runtime + difftime(end_time, start_time, units = "secs")
  runtime_core[runtime_core$cores == 1,]$citests <- 
    runtime_core[runtime_core$cores == 1,]$citests + sum(pc.parallel@n.edgetests)
  
  for (core in cores.vars){
    start_time <- Sys.time()
    pc.parallel <- pc_parallel(list(C=cor(d), n = n), indepTest = gaussCItest, p = p,
                      alpha = 0.01, skel.method = "parallel", num.cores = core)
    end_time <- Sys.time()
    
    runtime_core[runtime_core$cores == core,]$runtime <- 
      runtime_core[runtime_core$cores == core,]$runtime + difftime(end_time, start_time, units = "secs")
    runtime_core[runtime_core$cores == core,]$citests <- 
      runtime_core[runtime_core$cores == core,]$citests + sum(pc.parallel@n.edgetests)
  }
  
  j = j + 1
}

runtime_core$runtime <- runtime_core$runtime/num
runtime_core$citests <- runtime_core$citests/num
runtime_core$runtime_citests <- runtime_core$runtime / runtime_core$citests

write.table(runtime_core, "runtime_core.txt")

############################ different orderings ##############################
p <- 6     # number of random variables
n <- 20    # number of samples
s <- 0.6   # sparsity parameter (the lower the more sparse)
ens = 3

set.seed(3)
g <- randomDAG(p,s)     # true DAG
d <- rmvDAG(n,g)        # generate random samples from true DAG

plots_original_skeleton = list()
plots_stable_skeleton = list()

for (i in c(1:5)) {
  ord = sample(c(1:p), p, replace = FALSE)
  d_perm <- d[,ord]
  pc.original = skeleton(list(C=cor(d_perm), n = n), indepTest = gaussCItest, labels = as.character(ord),
                          alpha = 0.05, method = "original")
  pc.stable = skeleton(list(C=cor(d_perm), n = n), indepTest = gaussCItest, labels = as.character(ord),
                        alpha = 0.05, method = "stable")
  
  plots_original_skeleton[i] = pc.original@graph
  plots_stable_skeleton[i] = pc.stable@graph
}

par(mfrow=c(2,2))
plot(plots_original_skeleton[[2]])
plot(plots_original_skeleton[[3]])
plot(plots_original_skeleton[[4]])
plot(plots_original_skeleton[[5]])

plot(plots_stable_skeleton[[2]])
plot(plots_stable_skeleton[[3]])
plot(plots_stable_skeleton[[4]])
plot(plots_stable_skeleton[[5]])

# plots recreated in LaTeX

p <- 50   # number of random variables
n <- 100  # number of samples
ens <- 3  # expected neighborhood size
s = ens / (p-1)          # sparsity parameter 

g <- randomDAG(p,s)      # true DAG
d <- rmvDAG(n,g)         # generate random samples from true DAG

num_edges_original = c()
num_edges_stable = c()

for (i in c(1:10)) {
  ord = sample(c(1:p), p, replace = FALSE)
  d_perm <- d[,ord]
  pc.original = skeleton(list(C=cor(d_perm), n = n), indepTest = gaussCItest, labels = as.character(ord),
                         alpha = 0.05, method = "original")
  pc.stable = skeleton(list(C=cor(d_perm), n = n), indepTest = gaussCItest, labels = as.character(ord),
                       alpha = 0.05, method = "stable")
  num_edges_original[i] = length(pc.original@graph@edgeData@data)
  num_edges_stable[i] = length(pc.stable@graph@edgeData@data)
  
}

t(cbind(num_edges_original, num_edges_stable))


############################ nonparanormal data ################################
set.seed(1)
## parameter setting for simulating data
n <- 1000  # number of samples
ens <- 5   # expected neighborhood size (p-1)*s = ens
p <- 100   # number of nodes

s = ens / (p -1)    # sparsity parameter (the lower the more sparse) (p-1)*s = expected neighborhood size

# tables for evaluation

eval <- as.data.frame(matrix(data = 0, nrow = 3*2, ncol = 9))
colnames(eval) = c("alg", "data", "SHD", "TPR", "FPR", "TDR", "citests", "num_edges_true", "num_edges_pc")
eval$alg <- rep(c("parallel", "ranktau", "rankspearman"), each = 2)
eval$data <- rep(c("d_norm", "d_cop"), 3)

j = 1
num <- 10    # number of simulated graphs to take the mean in the end

while(j < num + 1) {
  
  g <- randomDAG(p,s)                        # true DAG
  d_norm <- rmvDAG(n,g, errDist = "normal")  # generate random samples from true DAG
  
  # transform marginals to F(1,1)-distribution
  d_cop = qf(pnorm(d_norm),1,1)
  
  for (data in c("d_norm", "d_cop")){
    
    if(data == "d_norm"){dat = d_norm}
    if(data == "d_cop"){dat = d_cop}
    
    for (alg in c("parallel", "ranktau", "rankspearman")){
      if (alg == "parallel") {suffStat = list(C=cor(dat), n = n)}
      if (alg == "ranktau") {suffStat = list(C = sin(cor(dat, method = "kendall") * pi/2), n = n)}
      if (alg == "rankspearman") {suffStat = list(C = 2 * sin(cor(dat, method = "spearman") * pi/6), n = n)}
      
      pc.fit <- pc_parallel(suffStat = suffStat, indepTest = gaussCItest, p = p,
                            alpha = 0.01, skel.method = "parallel", num.cores = 4, 
                            u2pd = "retry")
      
      eval[eval$alg == alg & eval$data == data,]$SHD <- eval[eval$alg == alg & eval$data == data,]$SHD  + shd(g,pc.fit)
      eval[eval$alg == alg & eval$data == data,]$TPR <- eval[eval$alg == alg & eval$data == data,]$TPR + compareGraphs(g, pc.fit@graph)[1]
      eval[eval$alg == alg & eval$data == data,]$FPR <- eval[eval$alg == alg & eval$data == data,]$FPR + compareGraphs(g, pc.fit@graph)[2]
      eval[eval$alg == alg & eval$data == data,]$TDR <- eval[eval$alg == alg & eval$data == data,]$TDR + compareGraphs(g, pc.fit@graph)[3]
      eval[eval$alg == alg & eval$data == data,]$citests <- eval[eval$alg == alg & eval$data == data,]$citests + sum(pc.fit@n.edgetests)
      eval[eval$alg == alg & eval$data == data,]$num_edges_true <- eval[eval$alg == alg & eval$data == data,]$num_edges_true + length(g@edgeData@data)
      eval[eval$alg == alg & eval$data == data,]$num_edges_pc <- eval[eval$alg == alg & eval$data == data,]$num_edges_pc + length(pc.fit@graph@edgeData@data)
      
      
    }
  }
  
  j = j + 1
  print(j)
}

eval[,c(3:9)] <- eval[,c(3:9)] / num
write.table(eval, "/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Seminar/Plots Sim Study/tables/eval_RPC.txt")

# Plot densities of d_cop and d_norm
pdf(file="/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Seminar/Plots Sim Study/RPC/plot_densities.pdf")
plot(density(d_norm), ylim = c(0,1.5), main = "", xlab = "", lwd = 2, col = "red")
lines(density(d_cop[d_cop<15]), col = "blue", lwd = 2)
legend("topright", legend = c("Gaussian", "Nonparanormal"), col = c("red", "blue"),lty = 1, cex = 0.8)
dev.off()
