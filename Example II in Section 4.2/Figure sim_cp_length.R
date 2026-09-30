#### Simultaneous coverage + length 
rm(list=ls())

library(dplyr)
library(ggplot2)
library(gridExtra)
library(patchwork)

# =====================================================
# 1 choose beta
# =====================================================

beta_id <- 1   

# =====================================================
# 2 path and sample size
# =====================================================

dir_prop <- "./multi_logit/result"

n_list <- c(10000,20000,30000,40000,50000)

# =====================================================
# 3 save result
# =====================================================

sim_list <- list()
length_box <- list()

# =====================================================
# 4 rusult
# =====================================================

for(n_val in n_list){
  
  files_prop <- list.files(
    dir_prop,
    pattern=paste0("res_n", n_val, "_J5_rep.*\\.rds$"),
    full.names=TRUE
  )
  
  R_prop <- length(files_prop)
  
  vec <- numeric(R_prop)
  
  for(r in seq_along(files_prop)){
    
    res <- readRDS(files_prop[r])
    
    cover_name  <- paste0("cover_sim", beta_id)
    length_name <- paste0("length_sim", beta_id)
    
    vec[r] <- res[[cover_name]]
    
    length_box[[length(length_box)+1]] <- data.frame(
      n = n_val,
      method = paste0("beta",beta_id),
      length = mean(res[[length_name]])
    )
    
  }
  
  p  <- mean(vec)
  se <- sqrt(p*(1-p)/R_prop)
  
  sim_list[[length(sim_list)+1]] <- data.frame(
    n = n_val,
    method = paste0("beta",beta_id),
    mean = p,
    lower = p - 1.96*se,
    upper = p + 1.96*se
  )
}

# =====================================================
# 5 sort result
# =====================================================

df_sim <- bind_rows(sim_list)
df_length_box <- bind_rows(length_box)

df_sim$n_plot <- df_sim$n/10000
df_length_box$n_plot <- df_length_box$n/10000

beta_name <- paste0("beta",beta_id)

df_sim$method <- factor(df_sim$method, levels=c(beta_name))
df_length_box$method <- factor(df_length_box$method, levels=c(beta_name))

# =====================================================
# 6 color
# =====================================================

cols <- setNames("#FF6961", beta_name)

# =====================================================
# 7 CP figure
# =====================================================

p1 <- ggplot(df_sim,
             aes(x=n_plot, y=mean, color=method, group=method)) +
  
  geom_line(aes(linetype=method), linewidth=0.9) +
  
  geom_point(aes(shape=method), size=2.5) +
  
  #geom_errorbar(aes(ymin=lower, ymax=upper),width=0,linewidth=0.5) +
  
  geom_hline(yintercept=0.95,
             linetype="dashed",
             color="grey40") +
  
  scale_color_manual(
    values=cols,
    labels=bquote(beta[.(beta_id)])
  ) +
  
  scale_linetype_manual(
    values=setNames("solid", beta_name),
    labels=bquote(beta[.(beta_id)])
  ) +
  
  scale_shape_manual(
    values=setNames(16, beta_name),
    labels=bquote(beta[.(beta_id)])
  ) +
  
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
# 8 Length figure
# =====================================================

p2 <- ggplot(df_length_box,
             aes(x=factor(n_plot), y=length, fill=method)) +
  
  geom_boxplot(
    outlier.shape=NA,
    colour="black",
    linewidth=0.3,
    width=0.4
  ) +
  
  scale_fill_manual(
    values=cols,
    labels=bquote(beta[.(beta_id)])
  ) +
  
  labs(
    x = expression(plain("Sample Size ") * n * "(×10"^4*")"),
    y = "Wid. of Simul. CB"
  ) +
  
  coord_cartesian(ylim=c(0,1)) +
  
  theme_bw(base_size=14) +
  theme(
    legend.position="none",
    panel.grid.major=element_line(colour="#E5E5E5", linewidth=0.5),
    panel.grid.minor=element_blank(),
    axis.title=element_text(face="bold", size=14)
  )

# =====================================================
# 9 attach figure
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








#### Simultaneous coverage + length ( beta1 / beta2)
rm(list=ls())

library(dplyr)
library(ggplot2)
library(gridExtra)
library(patchwork)

# =====================================================
# 1 choose beta
# =====================================================

beta_id <- 2   

# =====================================================
# 2 path and sample size
# =====================================================

dir_prop <- "./multi_logit/result"

n_list <- c(10000,20000,30000,40000,50000)

# =====================================================
# 3 save result
# =====================================================

sim_list <- list()
length_box <- list()

# =====================================================
# 4 result
# =====================================================

for(n_val in n_list){
  
  files_prop <- list.files(
    dir_prop,
    pattern=paste0("res_n", n_val, "_J5_rep.*\\.rds$"),
    full.names=TRUE
  )
  
  R_prop <- length(files_prop)
  
  vec <- numeric(R_prop)
  
  for(r in seq_along(files_prop)){
    
    res <- readRDS(files_prop[r])
    
    cover_name  <- paste0("cover_sim", beta_id)
    length_name <- paste0("length_sim", beta_id)
    
    vec[r] <- res[[cover_name]]
    
    length_box[[length(length_box)+1]] <- data.frame(
      n = n_val,
      method = paste0("beta",beta_id),
      length = mean(res[[length_name]])
    )
    
  }
  
  p  <- mean(vec)
  se <- sqrt(p*(1-p)/R_prop)
  
  sim_list[[length(sim_list)+1]] <- data.frame(
    n = n_val,
    method = paste0("beta",beta_id),
    mean = p,
    lower = p - 1.96*se,
    upper = p + 1.96*se
  )
}

# =====================================================
# 5 sort result
# =====================================================

df_sim <- bind_rows(sim_list)
df_length_box <- bind_rows(length_box)

df_sim$n_plot <- df_sim$n/10000
df_length_box$n_plot <- df_length_box$n/10000

beta_name <- paste0("beta",beta_id)

df_sim$method <- factor(df_sim$method, levels=c(beta_name))
df_length_box$method <- factor(df_length_box$method, levels=c(beta_name))

# =====================================================
# 6 color
# =====================================================

cols <- setNames("#836953", beta_name)

# =====================================================
# 7 CP figure
# =====================================================

p1 <- ggplot(df_sim,
             aes(x=n_plot, y=mean, color=method, group=method)) +
  
  geom_line(aes(linetype=method), linewidth=0.9) +
  
  geom_point(aes(shape=method), size=2.5) +
  
  #geom_errorbar(aes(ymin=lower, ymax=upper),width=0,linewidth=0.5) +
  
  geom_hline(yintercept=0.95,
             linetype="dashed",
             color="grey40") +
  
  scale_color_manual(
    values=cols,
    labels=bquote(beta[.(beta_id)])
  ) +
  
  scale_linetype_manual(
    values=setNames("solid", beta_name),
    labels=bquote(beta[.(beta_id)])
  ) +
  
  scale_shape_manual(
    values=setNames(16, beta_name),
    labels=bquote(beta[.(beta_id)])
  ) +
  
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
# 8 Length figure
# =====================================================

p2 <- ggplot(df_length_box,
             aes(x=factor(n_plot), y=length, fill=method)) +
  
  geom_boxplot(
    outlier.shape=NA,
    colour="black",
    linewidth=0.3,
    width=0.4
  ) +
  
  scale_fill_manual(
    values=cols,
    labels=bquote(beta[.(beta_id)])
  ) +
  
  labs(
    x = expression(plain("Sample Size ") * n * "(×10"^4*")"),
    y = "Wid. of Simul. CB"
  ) +
  
  coord_cartesian(ylim=c(0,1)) +
  
  theme_bw(base_size=14) +
  theme(
    legend.position="none",
    panel.grid.major=element_line(colour="#E5E5E5", linewidth=0.5),
    panel.grid.minor=element_blank(),
    axis.title=element_text(face="bold", size=14)
  )

# =====================================================
# 9 attatch figure
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


















