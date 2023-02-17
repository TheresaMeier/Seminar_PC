# Seminar Modeling under Dependence
# Visualization of the simulation study

colors = c("#F8766D", "#00BA38", "#619CFF","#009E73", "#D55E00")

# required libraries
library(ggplot2)
library(scales)
library(cowplot)

##################### varying dimensional setting ##############################
runtime_p <- read.delim("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Seminar/Plots Sim Study/tables/runtime_p.txt", sep = "")
eval_p <- read.delim("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Seminar/Plots Sim Study/tables/eval_p.txt", sep = "")

# run time
runtime_p$alg <- factor(runtime_p$alg, levels = c("original", "stable", "parallel"))
p1 = ggplot(data = runtime_p, mapping = aes(x = as.factor(p), y = runtime, col = alg)) +
  geom_point(size = 5) + ylab("runtime in s") + xlab("number of nodes") + labs(colour="PC algorithm") + 
  scale_color_manual(values=c("#F8766D","#56B4E9","#E69F00"))

p2 = ggplot(data = runtime_p, mapping = aes(x = as.factor(p), y = runtime_citests, col = alg)) +
  geom_point(size = 5) + ylab("runtime in s") + xlab("number of nodes") + labs(colour="PC algorithm") + 
  scale_color_manual(values=c("#F8766D","#56B4E9","#E69F00"))

p3 = ggplot(data = eval_p[eval_p$alg != "parallel",], 
            mapping = aes(x = as.factor(p), y = eval_p[eval_p$alg != "parallel",]$citests, col = alg)) +
  geom_point(size = 5) + ylab("number of CI tests") + xlab("number of nodes") + labs(colour="PC algorithm") +  
  scale_y_continuous(labels = unit_format(unit = "k", scale = 1e-3)) + 
  scale_color_manual(values=c("#F8766D","#56B4E9","#E69F00"))

# Evaluation via Structural Hamming Distance (SHD)
eval_p$alg <- factor(eval_p$alg, levels = c("original", "stable", "parallel"))

p4 = ggplot(data = eval_p[eval_p$alg != "parallel",], mapping = aes(x = as.factor(p), y = SHD, col = alg)) +
  geom_point(size = 5) + ylab("SHD") + xlab("number of nodes") + labs(colour="PC algorithm") + 
  scale_color_manual(values=c("#F8766D","#56B4E9","#E69F00"))

# Evaluation via TPR, FPR and TDR
p5 = ggplot(data = eval_p[eval_p$alg != "parallel",], mapping = aes(x = as.factor(p), y = TPR, col = alg)) +
  geom_point(size = 5) + ylab("TPR") + xlab("number of nodes") + labs(colour="PC algorithm") + 
  scale_color_manual(values=c("#F8766D","#56B4E9","#E69F00"))

p6 = ggplot(data = eval_p[eval_p$alg != "parallel",], mapping = aes(x = as.factor(p), y = FPR, col = alg)) +
  geom_point(size = 5) + ylab("FPR") + xlab("number of nodes") + labs(colour="PC algorithm") + 
  scale_color_manual(values=c("#F8766D","#56B4E9","#E69F00"))

p7 = ggplot(data = eval_p[eval_p$alg != "parallel",], mapping = aes(x = as.factor(p), y = TDR, col = alg)) +
  geom_point(size = 5) + ylab("TDR") + xlab("number of nodes") + labs(colour="PC algorithm")+ 
  scale_color_manual(values=c("#F8766D","#56B4E9","#E69F00")) 


############################ varying sample size ###############################
runtime_n <- read.delim("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Seminar/Plots Sim Study/tables/runtime_n.txt", sep = "")
eval_n <- read.delim("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Seminar/Plots Sim Study/tables/eval_n.txt", sep = "")

#run time
runtime_n$alg <- factor(runtime_n$alg, levels = c("original", "stable", "parallel"))
eval_n$alg <- factor(eval_n$alg, levels = c("original", "stable", "parallel"))

n1 = ggplot(data = runtime_n[runtime_n$p == 10,], mapping = aes(x = as.factor(n), y = runtime, col = alg)) +
  geom_point(size = 5) + ylab("runtime in s") + xlab("sample size") + labs(colour="PC algorithm") + 
  scale_color_manual(values=c("#F8766D","#56B4E9","#E69F00")) 

n2 = ggplot(data = runtime_n[runtime_n$p == 10,], mapping = aes(x = as.factor(n), y = runtime, col = alg)) +
  geom_point(size = 5) + ylab("runtime in s") + xlab("sample size") + labs(colour="PC algorithm") + 
  scale_color_manual(values=c("#F8766D","#56B4E9","#E69F00")) +ylim(c(0,30))

n3 = ggplot(data = runtime_n[runtime_n$p == 100,], mapping = aes(x = as.factor(n), y = runtime, col = alg)) +
  geom_point(size = 5) + ylab("runtime in s") + xlab("sample size") + labs(colour="PC algorithm") + 
  scale_color_manual(values=c("#F8766D","#56B4E9","#E69F00")) 

n4 = ggplot(data = eval_n[eval_n$alg != "parallel" & eval_n$p == 10,], mapping = aes(x = as.factor(n), y = citests, col = alg)) +
  geom_point(size = 5) + ylab("number of CI tests") + xlab("sample size") + labs(colour="PC algorithm") + 
  scale_color_manual(values=c("#F8766D","#56B4E9","#E69F00")) 

n5 = ggplot(data = eval_n[eval_n$alg != "parallel" & eval_n$p == 100,], mapping = aes(x = as.factor(n), y = citests, col = alg)) +
  geom_point(size = 5) + ylab("number of CI tests") + xlab("sample size") + labs(colour="PC algorithm") + 
  scale_color_manual(values=c("#F8766D","#56B4E9","#E69F00")) + scale_y_continuous(labels = unit_format(unit = "k", scale = 1e-3)) 

# Evaluation via Structural Hamming Distance (SHD)
n6 = ggplot(data = eval_n[eval_n$alg != "parallel" & eval_n$p == 10,], mapping = aes(x = as.factor(n), y = SHD, col = alg)) +
  geom_point(size = 5) + ylab("SHD") + xlab("sample size") + labs(colour="PC algorithm") + 
  scale_color_manual(values=c("#F8766D","#56B4E9","#E69F00")) 

n7 = ggplot(data = eval_n[eval_n$alg != "parallel" & eval_n$p == 100,], mapping = aes(x = as.factor(n), y = SHD, col = alg)) +
  geom_point(size = 5) + ylab("SHD") + xlab("sample size") + labs(colour="PC algorithm") + 
  scale_color_manual(values=c("#F8766D","#56B4E9","#E69F00"))

# Evaluation via TPR, FPR and TDR
n8 = ggplot(data = eval_n[eval_n$alg != "parallel"& eval_n$p == 10,], mapping = aes(x = as.factor(n), y = TPR, col = alg)) +
  geom_point(size = 5) + ylab("TPR") + xlab("sample size") + labs(colour="PC algorithm") + 
  scale_color_manual(values=c("#F8766D","#56B4E9","#E69F00")) 

n9 = ggplot(data = eval_n[eval_n$alg != "parallel"& eval_n$p == 100,], mapping = aes(x = as.factor(n), y = TPR, col = alg)) +
  geom_point(size = 5) + ylab("TPR") + xlab("sample size") + labs(colour="PC algorithm") + 
  scale_color_manual(values=c("#F8766D","#56B4E9","#E69F00")) 

n10 = ggplot(data = eval_n[eval_n$alg != "parallel"& eval_n$p == 10,], mapping = aes(x = as.factor(n), y = FPR, col = alg)) +
  geom_point(size = 5) + ylab("FPR") + xlab("sample size") + labs(colour="PC algorithm") + 
  scale_color_manual(values=c("#F8766D","#56B4E9","#E69F00"))

n11 = ggplot(data = eval_n[eval_n$alg != "parallel"& eval_n$p == 100,], mapping = aes(x = as.factor(n), y = FPR, col = alg)) +
  geom_point(size = 5) + ylab("FPR") + xlab("sample size") + labs(colour="PC algorithm") + 
  scale_color_manual(values=c("#F8766D","#56B4E9","#E69F00"))

n12 = ggplot(data = eval_n[eval_n$alg != "parallel"& eval_n$p == 10,], mapping = aes(x = as.factor(n), y = TDR, col = alg)) +
  geom_point(size = 5) + ylab("TDR") + xlab("sample size") + labs(colour="PC algorithm") + 
  scale_color_manual(values=c("#F8766D","#56B4E9","#E69F00")) 

n13 = ggplot(data = eval_n[eval_n$alg != "parallel"& eval_n$p == 100,], mapping = aes(x = as.factor(n), y = TDR, col = alg)) +
  geom_point(size = 5) + ylab("TDR") + xlab("sample size") + labs(colour="PC algorithm") + 
  scale_color_manual(values=c("#F8766D","#56B4E9","#E69F00"))

n14 = ggplot(data = eval_n[eval_n$alg != "parallel"& eval_n$p == 10,], mapping = aes(x = as.factor(n), y = num_edges_pc, col = alg)) +
  geom_point(size = 5) + ylab("number of edges") + xlab("sample size") + labs(colour="PC algorithm") + 
  scale_color_manual(values=c("#F8766D","#56B4E9","#E69F00")) + geom_hline(yintercept=eval_n[eval_n$alg == "original"& eval_n$p == 10,]$num_edges_true, linewidth=1, col = "black", linetype="dashed")

n15 = ggplot(data = eval_n[eval_n$alg != "parallel"& eval_n$p == 100,], mapping = aes(x = as.factor(n), y = num_edges_pc, col = alg)) +
  geom_point(size = 5) + ylab("number of edges") + xlab("sample size") + labs(colour="PC algorithm") + 
  scale_color_manual(values=c("#F8766D","#56B4E9","#E69F00")) + geom_hline(yintercept=eval_n[eval_n$alg == "original"& eval_n$p == 100,]$num_edges_true, linewidth=1, col = "black", linetype="dashed")


##################### different distributional assumptions #####################
eval_dist <- read.delim("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Seminar/Plots Sim Study/tables/eval_dist.txt", sep = "")
eval_dist$dist <- factor(eval_dist$dist, levels = c("normal", "cauchy", "t (k=4)"))

# plot of distributions
d1 = ggplot(data = data.frame(x = c(-3, 3)), aes(x)) +
  stat_function(fun = dnorm, n = 101, args = list(mean = 0, sd = 1), aes(colour = factor(c("normal"), levels = c("normal", "cauchy", "t (k=4)"))), linewidth = 1, show.legend = TRUE) + ylab("") +
  stat_function(fun = dcauchy, n = 101, args = list(location = 0, scale = 1), aes(colour = factor(c("cauchy"), levels = c("normal", "cauchy", "t (k=4)"))), linewidth = 1, show.legend = TRUE) +
  stat_function(fun = dt, n = 101, args = list(df = 4), aes(colour = factor(c("t (k=4)"), levels = c("normal", "cauchy", "t (k=4)"))), linewidth = 1, show.legend = TRUE) +
  scale_colour_manual("Distribution", values = c("#619CFF","#00BA38","#CC79A7"), breaks =c("normal", "cauchy", "t (k=4)"))

# Evaluation via Structural Hamming Distance (SHD)
d2 = ggplot(data = eval_dist[eval_dist$p == 10,], mapping = aes(x=as.factor(n), y=SHD, col = dist)) + 
  geom_point(size = 5)  + xlab("sample size") + labs(colour="Distribution") + 
  scale_color_manual(values=c("#619CFF","#00BA38","#CC79A7"))

d3 = ggplot(data = eval_dist[eval_dist$p == 100,], mapping = aes(x=as.factor(n), y=SHD, col = dist)) + 
  geom_point(size = 5)  + xlab("sample size") + labs(colour="Distribution") + 
  scale_color_manual(values=c("#619CFF","#00BA38","#CC79A7"))

# Evaluation via TPR, FPR and TDR
d4 = ggplot(data = eval_dist, mapping = aes(x=as.factor(szenario), y=TPR, col = dist)) + 
  geom_point(size = 5)  + xlab("") + labs(colour="Distribution") + 
  scale_color_manual(values=c("#619CFF","#00BA38","#CC79A7")) +
  scale_x_discrete(labels=c("1" = "p = 10\nn = 100", "2" = "p = 10\nn = 1000",
                            "3" = "p = 100\nn = 100", "4" = "p = 100\nn = 1000"))

d5 = ggplot(data = eval_dist, mapping = aes(x=as.factor(szenario), y=FPR, col = dist)) + 
  geom_point(size = 5)  + xlab("") + labs(colour="Distribution") + 
  scale_color_manual(values=c("#619CFF","#00BA38","#CC79A7")) +
  scale_x_discrete(labels=c("1" = "p = 10\nn = 100", "2" = "p = 10\nn = 1000",
                            "3" = "p = 100\nn = 100", "4" = "p = 100\nn = 1000"))

d6 = ggplot(data = eval_dist, mapping = aes(x=as.factor(szenario), y=TDR, col = dist)) + 
  geom_point(size = 5)  + xlab("") + labs(colour="Distribution") + 
  scale_color_manual(values=c("#619CFF","#00BA38","#CC79A7")) +
  scale_x_discrete(labels=c("1" = "p = 10\nn = 100", "2" = "p = 10\nn = 1000",
                            "3" = "p = 100\nn = 100", "4" = "p = 100\nn = 1000"))

##################### different number of cores (parallel) #####################
runtime_core <- read.delim("/home/theresa/Schreibtisch/Theresa/STUDIUM/Master Statistics and Data Science/Seminar/Plots Sim Study/tables/runtime_core.txt", sep = "")

# runtime
r1 = ggplot(data = runtime_core, aes(x = as.factor(cores), y = runtime)) +
  geom_point(size = 5, col = "#56B4E9") + ylab("runtime in s") + xlab("number of cores") +
  ylim(c(1.5,5))

########################### export plots #######################################

plots <- align_plots(p1, p2, p3, p4, p5, p6, p7, 
                     n1, n2, n3, n4, n5, n6, n7, n8, n9, n10, n11, n12, n13, n14, n15, 
                     d1, d2, d3, d4, d5, d6, 
                     r1, align="v", axis = "rbtl")
ggdraw(plots[[1]])
ggsave("runtime_p.jpg")
ggdraw(plots[[2]])
ggsave("runtime_citest_p.jpg")
ggdraw(plots[[3]])
ggsave("citests_p.jpg")
ggdraw(plots[[4]])
ggsave("SHD_p.jpg")
ggdraw(plots[[5]])
ggsave("TPR_p.jpg")
ggdraw(plots[[6]])
ggsave("FPR_p.jpg")
ggdraw(plots[[7]])
ggsave("TDR_p.jpg")
ggdraw(plots[[8]])
ggsave("runtime_n_10.jpg")
ggdraw(plots[[9]])
ggsave("runtime_n_10_scale.jpg")
ggdraw(plots[[10]])
ggsave("runtime_n_100.jpg")
ggdraw(plots[[11]])
ggsave("citests_n_10.jpg")
ggdraw(plots[[12]])
ggsave("citests_n_100.jpg")
ggdraw(plots[[13]])
ggsave("SHD_n_10.jpg")
ggdraw(plots[[14]])
ggsave("SHD_n_100.jpg")
ggdraw(plots[[15]])
ggsave("TPR_n_10.jpg")
ggdraw(plots[[16]])
ggsave("TPR_n_100.jpg")
ggdraw(plots[[17]])
ggsave("FPR_n_10.jpg")
ggdraw(plots[[18]])
ggsave("FPR_n_100.jpg")
ggdraw(plots[[19]])
ggsave("TDR_n_10.jpg")
ggdraw(plots[[20]])
ggsave("TDR_n_100.jpg")
ggdraw(plots[[21]])
ggsave("num_edges_n_10.jpg")
ggdraw(plots[[22]])
ggsave("num_edges_n_100.jpg")
ggdraw(plots[[23]])
ggsave("distributions_dist.jpg")
ggdraw(plots[[24]])
ggsave("SHD_dist_10.jpg")
ggdraw(plots[[25]])
ggsave("SHD_dist_100.jpg")
ggdraw(plots[[26]])
ggsave("TPR_dist.jpg")
ggdraw(plots[[27]])
ggsave("FPR_dist.jpg")
ggdraw(plots[[28]])
ggsave("TDR_dist.jpg")
ggdraw(plots[[29]])
ggsave("runtime_core.jpg")