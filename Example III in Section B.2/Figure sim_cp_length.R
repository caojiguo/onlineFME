#### Simultaneous coverage + length
rm(list=ls())
library(dplyr)
library(ggplot2)
library(gridExtra)
library(patchwork)

# =====================================================
# 1 path
# =====================================================
dir_prop <- "./multi_linear/result"
n_list <- c(10000,20000,30000,40000,50000)

# =====================================================
# 2 Simultaneous coverage (ONLY beta1)
# =====================================================
sim_list <- list()
length_box <- list()

for(n_val in n_list){
  
  files_prop <- list.files(dir_prop,
                           pattern=paste0("res_n", n_val, "_J5_rep.*\\.rds$"),
                           full.names=TRUE)
  R_prop <- length(files_prop)
  vec1 <- numeric(R_prop)
  
  for(r in seq_along(files_prop)){
    res <- readRDS(files_prop[r])
    vec1[r] <- res$cover_sim1
    
    length_box[[length(length_box)+1]] <- data.frame(
      n = n_val,
      method = "beta1",
      length = mean(res$length_sim1)
    )
  }
  
  # Simultaneous coverage beta1
  p1  <- mean(vec1)
  se1 <- sqrt(p1*(1-p1)/R_prop)
  
  sim_list[[length(sim_list)+1]] <- data.frame(
    n = n_val,
    method = "beta1",
    mean = p1,
    lower = p1 - 1.96*se1,
    upper = p1 + 1.96*se1
  )
}

df_sim <- bind_rows(sim_list)
df_length_box <- bind_rows(length_box)

df_sim$n_plot <- df_sim$n/10000
df_length_box$n_plot <- df_length_box$n/10000

df_sim$method <- factor(df_sim$method, levels=c("beta1"))
df_length_box$method <- factor(df_length_box$method, levels=c("beta1"))

# =====================================================
# 3 color
# =====================================================
cols <- c("beta1"="#FF6961")

# =====================================================
# 4 Coverage figure
# =====================================================
p1 <- ggplot(df_sim,
             aes(x=n_plot, y=mean, color=method, group=method)) +
  geom_line(aes(linetype=method), linewidth=0.9) +
  geom_point(aes(shape=method), size=2.5) +
  #geom_errorbar(aes(ymin=lower, ymax=upper), width=0, linewidth=0.5) +
  geom_hline(yintercept=0.95, linetype="dashed", color="grey40") +
  scale_color_manual(values=cols, labels=expression(beta[1])) +
  scale_linetype_manual(values=c("beta1"="solid"), labels=expression(beta[1])) +
  scale_shape_manual(values=c("beta1"=16), labels=expression(beta[1])) +
  labs(
    x = expression(plain("Sample Size ") * n * "(×10"^4*")"),
    y = "Simul. CP"
  ) +
  coord_cartesian(ylim=c(0.8,1)) +
  theme_bw(base_size=14) +
  theme(
    legend.position="none",   
    panel.grid.major=element_line(colour="#E5E5E5", linewidth=0.5),
    panel.grid.minor=element_blank(),
    axis.title=element_text(face="bold", size=14)
  )
# =====================================================
# 5 Length figure
# =====================================================
p2 <- ggplot(df_length_box,
             aes(x=factor(n_plot), y=length, fill=method)) +
  geom_boxplot(
    outlier.shape=NA,
    colour="black",
    linewidth=0.3,
    width=0.4
  ) +
  scale_fill_manual(values=cols) +
  labs(
    x = expression(plain("Sample Size ") * n * "(×10"^4*")"),
    y = "Wid. of Simul. CB"
  ) +
  coord_cartesian(ylim=c(0,0.3)) +
  theme_bw(base_size=14) +
  theme(
    legend.position="none",
    panel.grid.major=element_line(colour="#E5E5E5", linewidth=0.5),
    panel.grid.minor=element_blank(),
    axis.title=element_text(face="bold", size=14)
  )

# =====================================================
# 6 attatch figure
# =====================================================


final_plot <- p1 | p2
final_plot

ggsave(
  filename = "./1_sim_cp_length.pdf",
  plot = final_plot,
  width = 6,      
  height = 3,      
  device = "pdf"    
)










#### Simultaneous coverage + length
rm(list=ls())
library(dplyr)
library(ggplot2)
library(gridExtra)
library(patchwork)

# =====================================================
# 1 path
# =====================================================
dir_prop <- "./multi_linear/result"
n_list <- c(10000,20000,30000,40000,50000)

# =====================================================
# 2 Simultaneous coverage (ONLY beta2)
# =====================================================
sim_list <- list()
length_box <- list()

for(n_val in n_list){
  
  files_prop <- list.files(dir_prop,
                           pattern=paste0("res_n", n_val, "_J5_rep.*\\.rds$"),
                           full.names=TRUE)
  R_prop <- length(files_prop)
  vec2 <- numeric(R_prop)
  
  for(r in seq_along(files_prop)){
    res <- readRDS(files_prop[r])
    vec2[r] <- res$cover_sim2
    
    length_box[[length(length_box)+1]] <- data.frame(
      n = n_val,
      method = "beta2",
      length = mean(res$length_sim2)
    )
  }
  
  # Simultaneous coverage beta2
  p2_mean  <- mean(vec2)
  se2 <- sqrt(p2_mean*(1-p2_mean)/R_prop)
  
  sim_list[[length(sim_list)+1]] <- data.frame(
    n = n_val,
    method = "beta2",
    mean = p2_mean,
    lower = p2_mean - 1.96*se2,
    upper = p2_mean + 1.96*se2
  )
}

df_sim <- bind_rows(sim_list)
df_length_box <- bind_rows(length_box)

df_sim$n_plot <- df_sim$n/10000
df_length_box$n_plot <- df_length_box$n/10000

df_sim$method <- factor(df_sim$method, levels=c("beta2"))
df_length_box$method <- factor(df_length_box$method, levels=c("beta2"))

# =====================================================
# 3 color
# =====================================================
cols <- c("beta2"="#836953")

# =====================================================
# 4 Coverage figure
# =====================================================
p1 <- ggplot(df_sim,
             aes(x=n_plot, y=mean, color=method, group=method)) +
  geom_line(aes(linetype=method), linewidth=0.9) +
  geom_point(aes(shape=method), size=2.5) +
  #geom_errorbar(aes(ymin=lower, ymax=upper), width=0, linewidth=0.5) +
  geom_hline(yintercept=0.95, linetype="dashed", color="grey40") +
  scale_color_manual(values=cols, labels=expression(beta[2])) +
  scale_linetype_manual(values=c("beta2"="solid"), labels=expression(beta[2])) +
  scale_shape_manual(values=c("beta2"=16), labels=expression(beta[2])) +
  labs(
    x = expression(plain("Sample Size ") * n * "(×10"^4*")"),
    y = "Simul. CP"
  ) +
  coord_cartesian(ylim=c(0.8,1)) +
  theme_bw(base_size=14) +
  theme(
    legend.position="none",
    panel.grid.major=element_line(colour="#E5E5E5", linewidth=0.5),
    panel.grid.minor=element_blank(),
    axis.title=element_text(face="bold", size=14)
  )

# =====================================================
# 5 Length figure
# =====================================================
p2 <- ggplot(df_length_box,
             aes(x=factor(n_plot), y=length, fill=method)) +
  geom_boxplot(
    outlier.shape=NA,
    colour="black",
    linewidth=0.3,
    width=0.4
  ) +
  scale_fill_manual(values=cols) +
  labs(
    x = expression(plain("Sample Size ") * n * "(×10"^4*")"),
    y = "Wid. of Simul. CB"
  ) +
  coord_cartesian(ylim=c(0,0.3)) +
  theme_bw(base_size=14) +
  theme(
    legend.position="none",
    panel.grid.major=element_line(colour="#E5E5E5", linewidth=0.5),
    panel.grid.minor=element_blank(),
    axis.title=element_text(face="bold", size=14)
  )

# =====================================================
# 6 attatch figure
# =====================================================
final_plot <- p1 | p2
final_plot

ggsave(
  filename = "./2_sim_cp_length.pdf",
  plot = final_plot,
  width = 6,
  height = 3,
  device = "pdf"
)




