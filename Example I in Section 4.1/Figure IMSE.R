rm(list=ls())

library(dplyr)
library(ggplot2)
library(patchwork)

# =====================================================
# path
# =====================================================

dir_prop <- "./proposed/result"
dir_rkhs <- "./uni_xie/result"

# =====================================================
# setting
# =====================================================

n_list <- c(10000,20000,30000,40000,50000)

n_val <- 50000

lent <- 101
tobs <- seq(0,1,length.out=lent)

beta_true <- sin(0.5*pi*tobs)

t_pick <- c(0.15,0.5,0.85)
idx_pick <- sapply(t_pick,function(tt) which.min(abs(tobs-tt)))

cols <- c(
  "FME-Online"="#FF6961",
  "RKHS-BS"="#836953"
)

# =====================================================
# FME-Online curve
# =====================================================

files <- list.files(
  dir_prop,
  pattern=paste0("res_n",n_val,"_J5_rep.*\\.rds$"),
  full.names=TRUE
)

R <- length(files)

beta_mat  <- matrix(0,R,lent)
lower_mat <- matrix(0,R,lent)
upper_mat <- matrix(0,R,lent)

lower_sim_mat <- matrix(0,R,lent)
upper_sim_mat <- matrix(0,R,lent)

for(r in seq_along(files)){
  
  res <- readRDS(files[r])
  
  beta_hat <- res$beta
  
  beta_mat[r,]  <- beta_hat
  
  lower_mat[r,] <- beta_hat - res$length_pt/2
  upper_mat[r,] <- beta_hat + res$length_pt/2
  
  lower_sim_mat[r,] <- beta_hat - res$length_sim/2
  upper_sim_mat[r,] <- beta_hat + res$length_sim/2
}

beta_mean <- colMeans(beta_mat)
lower_mean <- colMeans(lower_mat)
upper_mean <- colMeans(upper_mat)

lower_sim_mean <- colMeans(lower_sim_mat)
upper_sim_mean <- colMeans(upper_sim_mat)

df_plot <- data.frame(
  t = tobs,
  truth = beta_true,
  estimate = beta_mean
)

df_ci <- data.frame(
  t = t_pick,
  mean  = beta_mean[idx_pick],
  lower = lower_mean[idx_pick],
  upper = upper_mean[idx_pick]
)

df_sim <- data.frame(
  t = tobs,
  lower = lower_sim_mean,
  upper = upper_sim_mean
)

p_fme <- ggplot(df_plot,aes(x=t)) +
  
  geom_ribbon(data=df_sim,
              aes(ymin=lower,ymax=upper),
              fill="grey80",alpha=0.5) +
  
  geom_line(aes(y=truth),
            linetype="dashed",
            linewidth=0.9) +
  
  geom_line(aes(y=estimate),
            color=cols["FME-Online"],
            linewidth=0.9) +
  
  geom_errorbar(data=df_ci,
                aes(x=t,ymin=lower,ymax=upper),
                width=0.03,
                linewidth=0.7,
                color=cols["FME-Online"]) +
  
  labs(x="t",
       y=expression(beta(t))) +
  
  theme_bw(base_size=14) +
  theme(
    legend.position="none",
    panel.grid.major=element_line(colour="#E5E5E5",linewidth=0.5),
    panel.grid.minor=element_blank(),
    axis.title=element_text(face="bold",size=14)
  )

# =====================================================
# RKHS-BS curve
# =====================================================

files <- list.files(
  dir_rkhs,
  pattern=paste0("res_n",n_val,"_rep.*\\.rds$"),
  full.names=TRUE
)

R <- length(files)

beta_mat  <- matrix(0,R,lent)
lower_mat <- matrix(0,R,lent)
upper_mat <- matrix(0,R,lent)

for(r in seq_along(files)){
  
  res <- readRDS(files[r])
  
  lower_mat[r,] <- res$beta_sd_lower
  upper_mat[r,] <- res$beta_sd_upper
  
  beta_mat[r,] <- (res$beta_sd_lower + res$beta_sd_upper)/2
}

beta_mean <- colMeans(beta_mat)
lower_mean <- colMeans(lower_mat)
upper_mean <- colMeans(upper_mat)

df_plot <- data.frame(
  t = tobs,
  truth = beta_true,
  estimate = beta_mean
)

df_ci <- data.frame(
  t = t_pick,
  mean  = beta_mean[idx_pick],
  lower = lower_mean[idx_pick],
  upper = upper_mean[idx_pick]
)

df_sim <- data.frame(
  t = tobs,
  lower = lower_mean,
  upper = upper_mean
)

p_rkhs <- ggplot(df_plot,aes(x=t)) +
  
  geom_ribbon(data=df_sim,
              aes(ymin=lower,ymax=upper),
              fill="grey80",alpha=0.5) +
  
  geom_line(aes(y=truth),
            linetype="dashed",
            linewidth=0.9) +
  
  geom_line(aes(y=estimate),
            color=cols["RKHS-BS"],
            linewidth=0.9) +
  
  geom_errorbar(data=df_ci,
                aes(x=t,ymin=lower,ymax=upper),
                width=0.03,
                linewidth=0.7,
                color=cols["RKHS-BS"]) +
  
  labs(x="t",
       y=expression(beta(t))) +
  
  theme_bw(base_size=14) +
  theme(
    legend.position="none",
    panel.grid.major=element_line(colour="#E5E5E5",linewidth=0.5),
    panel.grid.minor=element_blank(),
    axis.title=element_text(face="bold",size=14)
  )

# =====================================================
# IMSE
# =====================================================

imse_list <- list()

for(n_val in n_list){
  
  get_imse <- function(files){
    sapply(files,function(f) readRDS(f)$IMSE)
  }
  
  files_prop <- list.files(
    dir_prop,
    pattern=paste0("res_n",n_val,"_J5_rep.*\\.rds$"),
    full.names=TRUE
  )
  
  imse_vals <- get_imse(files_prop)
  
  imse_list[[length(imse_list)+1]] <-
    data.frame(n=n_val,
               method="FME-Online",
               mean=mean(imse_vals),
               sd=sd(imse_vals))
  
  files_rkhs <- list.files(
    dir_rkhs,
    pattern=paste0("res_n",n_val,"_rep.*\\.rds$"),
    full.names=TRUE
  )
  
  imse_vals <- get_imse(files_rkhs)
  
  imse_list[[length(imse_list)+1]] <-
    data.frame(n=n_val,
               method="RKHS-BS",
               mean=mean(imse_vals),
               sd=sd(imse_vals))
}

df_imse <- bind_rows(imse_list)

df_imse$n_plot <- df_imse$n/10000

df_imse$method <- factor(df_imse$method,
                         levels=c("FME-Online","RKHS-BS"))

df_imse$lower <- df_imse$mean - df_imse$sd
df_imse$upper <- df_imse$mean + df_imse$sd

# ★ IMSE ×1000
df_imse$mean  <- df_imse$mean  * 1000
df_imse$lower <- df_imse$lower * 1000
df_imse$upper <- df_imse$upper * 1000

p_imse <- ggplot(df_imse,
                 aes(x=n_plot,
                     y=mean,
                     color=method,
                     shape=method,
                     group=method)) +
  
  geom_line(aes(linetype=method),linewidth=0.9) +
  geom_point(size=2.5) +
  
  geom_errorbar(aes(ymin=lower,ymax=upper),
                width=0,
                linewidth=0.5) +
  
  scale_x_continuous(breaks=1:5) +
  
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
    y = expression("IMSE (×10"^{-3}*")"),
    color="Method",
    shape="Method",
    linetype="Method"
  ) +
  
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

library(patchwork)


p_curve <- (p_fme | p_rkhs | p_imse) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")

print(p_curve)



ggsave(
  filename = "./IMSE.pdf",
  plot = p_curve,
  width = 9,       
  height = 3.6,       
  device = "pdf"   
)
