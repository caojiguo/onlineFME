rm(list=ls())

library(gridExtra)
library(dplyr)
library(ggplot2)
library(patchwork)

# =====================================================
# path
# =====================================================

dir_prop <- "./proposed/result"
dir_rkhs <- "./uni_xie/result"

n_list <- c(10000,20000,30000,40000,50000)

lent <- 101
t_grid <- seq(0,1,length.out=lent)

beta_true <- sin(0.5*pi*t_grid)#exp(-t_grid)
# idx <- which(t_grid >= 0.025 & t_grid <= 0.975)
# idx <- which(t_grid >= 0.02 & t_grid <= 0.95)
idx <- which(t_grid >= 0.05 & t_grid <= 0.95)
# idx <- which(t_grid >= 0.1 & t_grid <= 0.9)

# =====================================================
# 
# =====================================================

summary_list <- list()
length_box   <- list()

# =====================================================
# loop
# =====================================================

for(n_val in n_list){
  
  cat("Processing n =",n_val,"\n")
  
  # ==========================
  # FME-Online
  # ==========================
  
  files_prop <- list.files(
    dir_prop,
    pattern=paste0("res_n",n_val,"_J5_rep.*\\.rds$"),
    full.names=TRUE
  )
  
  R_prop <- length(files_prop)
  vec_prop <- numeric(R_prop)
  
  for(r in seq_along(files_prop)){
    
    res <- readRDS(files_prop[r])
    
    vec_prop[r] <- res$cover_sim
    
    length_box[[length(length_box)+1]] <- data.frame(
      n = n_val,
      method = "FME-Online",
      length = mean(res$length_sim)
    )
  }
  
  p  <- mean(vec_prop)
  se <- sqrt(p*(1-p)/R_prop)
  
  summary_list[[length(summary_list)+1]] <-
    data.frame(
      n=n_val,
      method="FME-Online",
      mean=p,
      lower=p-1.96*se,
      upper=p+1.96*se
    )
  
  
  # ==========================
  # RKHS-BS
  # ==========================
  
  files_aos <- list.files(
    dir_rkhs,
    pattern=paste0("res_n",n_val,"_rep.*\\.rds$"),
    full.names=TRUE
  )
  
  R_aos <- length(files_aos)
  vec_aos <- numeric(R_aos)
  
  for(r in seq_along(files_aos)){
    
    res <- readRDS(files_aos[r])
    
    lower_sub <- res$beta_sim_lower[idx]
    upper_sub <- res$beta_sim_upper[idx]
    true_sub  <- beta_true[idx]
    
    covered <- all(true_sub >= lower_sub & true_sub <= upper_sub)
    
    vec_aos[r] <- as.numeric(covered)
    
    length_box[[length(length_box)+1]] <- data.frame(
      n = n_val,
      method = "RKHS-BS",
      length = mean(res$beta_sim_upper - res$beta_sim_lower)
    )
  }
  
  p  <- mean(vec_aos)
  se <- sqrt(p*(1-p)/R_aos)
  
  summary_list[[length(summary_list)+1]] <-
    data.frame(
      n=n_val,
      method="RKHS-BS",
      mean=p,
      lower=p-1.96*se,
      upper=p+1.96*se
    )
}

# =====================================================
# 
# =====================================================

df_summary <- bind_rows(summary_list)
df_length_box <- bind_rows(length_box)

df_summary$n_plot <- df_summary$n/10000
df_length_box$n_plot <- df_length_box$n/10000

df_summary$method <- factor(df_summary$method,
                            levels=c("FME-Online","RKHS-BS"))

df_length_box$method <- factor(df_length_box$method,
                               levels=c("FME-Online","RKHS-BS"))

# =====================================================
# 
# =====================================================

cols <- c(
  "FME-Online"="#FF6961",
  "RKHS-BS"="#836953"
)

# =====================================================
# 
# =====================================================

p1 <- ggplot(df_summary,
             aes(x=n_plot,y=mean,color=method,group=method)) +
  
  geom_line(aes(linetype=method),linewidth=0.9) +
  geom_point(aes(shape=method),size=2.5) +
  
  #geom_errorbar(aes(ymin=lower,ymax=upper),
  #              width=0,linewidth=0.5) +
  
  geom_hline(yintercept=0.95,
             linetype="dashed",
             color="grey40") +
  
  scale_color_manual(values=cols) +
  
  scale_linetype_manual(values=c(
    "FME-Online"="solid",
    "RKHS-BS"="dashed"
  )) +
  
  scale_shape_manual(values=c(
    "FME-Online"=16,
    "RKHS-BS"=17
  )) +
  
  labs(
    x = expression(plain("Sample Size ") * n * "(×10"^4*")"),
    y="Simul. CP",
    color="Method",
    shape="Method",
    linetype="Method"
  ) +
  
  coord_cartesian(ylim=c(0.8,1)) +
  
  theme_bw(base_size=14) +
  theme(
    legend.position="right",
    panel.grid.major=element_line(colour="#E5E5E5",linewidth=0.5),
    panel.grid.minor=element_blank(),
    axis.title=element_text(face="bold",size=14)
  )

# =====================================================
# 
# =====================================================

p2 <- ggplot(df_length_box,
             aes(x=factor(n_plot),
                 y=length,
                 fill=method)) +
  
  geom_boxplot(
    outlier.shape=NA,
    colour="black",
    linewidth=0.3,
    width=0.8
  ) +
  
  scale_fill_manual(values=cols) +
  
  labs(
    x = expression(plain("Sample Size ") * n * "(×10"^4*")"),
    y = "Wid. of Simul. CB"
  ) +
  
  coord_cartesian(ylim=c(0,0.4)) +
  
  theme_bw(base_size=14) +
  theme(
    legend.position="none",  
    panel.grid.major=element_line(colour="#E5E5E5",linewidth=0.5),
    panel.grid.minor=element_blank(),
    axis.title=element_text(face="bold",size=14)
  )

# =====================================================
# 
# =====================================================





final_plot <- p1 | p2  
final_plot


ggsave(
  filename = "./sim_cp_length.pdf",
  plot = final_plot,
  width = 7.3,       
  height = 3,       
  device = "pdf"    
)
