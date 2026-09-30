rm(list=ls())
library(gridExtra)
library(dplyr)
library(ggplot2)
library(patchwork)  

# =====================================================
# path
# =====================================================
dir_prop <- "./proposed/result"
dir_xie <- "./uni_xie/result"
n_list <- c(10000,20000,30000,40000,50000)

lent <- 101
tobs <- seq(0,1,length.out=lent)
beta_0 <- sin(0.5*pi*tobs)#exp(-tobs)

# =====================================================
# fixed t
# =====================================================
t_pick <- c(0.15,0.5,0.85)
# t_pick <- c(0.2,0.5,0.8)
# t_pick <- c(0.1,0.5,0.9)
idx_pick <- sapply(t_pick,function(tt) which.min(abs(tobs-tt)))

summary_list <- list()
length_list  <- list()
length_box <- list()  

# =====================================================
# loop
# =====================================================
for(n_val in n_list){
  
  cat("Processing n =",n_val,"\n")
  
  # FME-Online
  files <- list.files(
    dir_prop,
    pattern=paste0("res_n",n_val,"_J5_rep.*\\.rds$"),
    full.names=TRUE
  )
  
  R <- length(files)
  
  cover_mat <- matrix(0,R,lent)
  length_mat <- matrix(0,R,lent)
  
  for(r in seq_along(files)){
    res <- readRDS(files[r])
    beta_hat <- res$beta
    
    lower_pt <- beta_hat - res$length_pt/2
    upper_pt <- beta_hat + res$length_pt/2
    
    cover_mat[r,] <- (beta_0>=lower_pt & beta_0<=upper_pt)
    length_mat[r,] <- res$length_pt
    
 
    for(j in seq_along(idx_pick)){
      idx <- idx_pick[j]
      length_box[[length(length_box)+1]] <- data.frame(
        n = n_val,
        t = t_pick[j],
        method = "FME-Online",
        length = res$length_pt[idx]
      )
    }
  }
  
  for(j in seq_along(idx_pick)){
    idx <- idx_pick[j]
    vec <- cover_mat[,idx]
    p  <- mean(vec)
    se <- sqrt(p*(1-p)/R)
    
    summary_list[[length(summary_list)+1]] <-
      data.frame(
        n=n_val,
        t=t_pick[j],
        method="FME-Online",
        mean=p,
        lower=p-1.96*se,
        upper=p+1.96*se
      )
    
    length_list[[length(length_list)+1]] <-
      data.frame(
        n=n_val,
        t=t_pick[j],
        method="FME-Online",
        mean=mean(length_mat[,idx])
      )
  }
  
  # RKHS-BS
  files <- list.files(
    dir_xie,
    pattern=paste0("res_n",n_val,"_rep.*\\.rds$"),
    full.names=TRUE
  )
  
  R <- length(files)
  
  cover_sd  <- matrix(0,R,lent)
  length_sd <- matrix(0,R,lent)
  
  for(r in seq_along(files)){
    res <- readRDS(files[r])
    cover_sd[r,] <- res$count_sd
    length_sd[r,] <- res$beta_sd_upper - res$beta_sd_lower
    
    
    for(j in seq_along(idx_pick)){
      idx <- idx_pick[j]
      length_box[[length(length_box)+1]] <- data.frame(
        n = n_val,
        t = t_pick[j],
        method = "RKHS-BS",
        length = res$beta_sd_upper[idx] - res$beta_sd_lower[idx]
      )
    }
  }
  
  for(j in seq_along(idx_pick)){
    idx <- idx_pick[j]
    vec <- cover_sd[,idx]
    p  <- mean(vec)
    se <- sqrt(p*(1-p)/R)
    
    summary_list[[length(summary_list)+1]] <-
      data.frame(
        n=n_val,
        t=t_pick[j],
        method="RKHS-BS",
        mean=p,
        lower=p-1.96*se,
        upper=p+1.96*se
      )
    
    length_list[[length(length_list)+1]] <-
      data.frame(
        n=n_val,
        t=t_pick[j],
        method="RKHS-BS",
        mean=mean(length_sd[,idx])
      )
  }
}

# =====================================================
# 
# =====================================================
df_summary <- bind_rows(summary_list)
df_summary$type <- "Coverage"

df_length <- bind_rows(length_list)
df_length$type <- "Length"
df_length$lower <- NA
df_length$upper <- NA

df_plot <- bind_rows(df_summary,df_length)
df_plot$n_plot <- df_plot$n/10000
df_plot$method <- factor(df_plot$method, levels=c("FME-Online","RKHS-BS"))

df_length_box <- bind_rows(length_box)
df_length_box$n_plot <- df_length_box$n/10000
df_length_box$method <- factor(df_length_box$method, levels=c("FME-Online","RKHS-BS"))

# =====================================================

# =====================================================
cols <- c(
  "FME-Online"="#FF6961",
  "RKHS-BS"="#836953"
)

# =====================================================

# =====================================================


# Coverage 
p1 <- ggplot(subset(df_plot,type=="Coverage"),
             aes(x=n_plot, y=mean, color=method, group=method)) +
  geom_line(aes(linetype=method), linewidth=0.9) +
  geom_point(aes(shape=method), size=2.5) +
  #geom_errorbar(aes(ymin=lower, ymax=upper), width=0, linewidth=0.5) +
  geom_hline(yintercept=0.95, linetype="dashed", color="grey40") +
  facet_grid(. ~ t, labeller = labeller(t = function(x) paste0("t = ", x))) +
  scale_color_manual(values=cols) +
  scale_linetype_manual(values=c("FME-Online"="solid","RKHS-BS"="dashed")) +
  scale_shape_manual(values=c(16,17)) +
  labs(x = NULL, y="Point. CP") +   
  coord_cartesian(ylim = c(0.6, 1.0)) +  
  theme_bw(base_size=14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 12),
    strip.background = element_rect(fill = "#D9D9D9", colour = NA),
    strip.text = element_text(face = "bold", size = 13),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(colour = "#E5E5E5", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.2, "lines"),
    axis.text = element_text(size = 12),
    axis.title = element_text(face="bold", color="#4D4D4D", size = 14)
  )

# Length 
p2 <- ggplot(df_length_box,
             aes(x=factor(n_plot), y=length, fill=method)) +
  geom_boxplot(
    outlier.shape = NA,     
    colour = "black",       
    linewidth = 0.3,        
    width = 0.8             
  ) +
  facet_grid(. ~ t) +
  scale_fill_manual(values=cols) +
  labs(x = expression(plain("Sample Size ") * n * "(×10"^4*")"), 
       y = "Wid. of Point. CI")+
  coord_cartesian(ylim = c(0, 0.2)) +   
  theme_bw(base_size=14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 12),
    strip.text = element_blank(),   
    strip.background = element_rect(fill = "#D9D9D9", colour = NA),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_line(colour = "#E5E5E5", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.2, "lines"),
    axis.text = element_text(size = 12),
    axis.title = element_text(face="bold", color="#4D4D4D", size = 14)
  )




final_plot <- p1 / p2 


ggsave(
  filename = "./pt_cp_length.pdf",
  plot = final_plot,
  width = 9.3,       
  height = 6,       
  device = "pdf"    
)
