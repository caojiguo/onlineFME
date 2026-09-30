rm(list=ls())

library(dplyr)
library(ggplot2)
library(patchwork)

# =====================================================
# 1 path
# =====================================================

dir_prop <- "./proposed/result"

n_list <- c(10000,20000,30000,40000,50000)
n_val <- 50000

# =====================================================
# 2 grid
# =====================================================

lent <- 101
tobs <- seq(0,1,length.out=lent)

beta1_true <- sin(0.5*pi*tobs)
beta2_true <- cos(0.5*pi*tobs)

t_pick <- c(0.2,0.5,0.8)
idx_pick <- sapply(t_pick,function(tt) which.min(abs(tobs-tt)))

# =====================================================
# 3. beta1 beta2
# =====================================================

files <- list.files(
  dir_prop,
  pattern=paste0("res_n",n_val,"_J5_rep.*\\.rds$"),
  full.names=TRUE
)

R <- length(files)

beta1_mat <- matrix(0,R,lent)
beta2_mat <- matrix(0,R,lent)

lower_pt1 <- matrix(0,R,lent)
upper_pt1 <- matrix(0,R,lent)

lower_pt2 <- matrix(0,R,lent)
upper_pt2 <- matrix(0,R,lent)

lower_sim1 <- matrix(0,R,lent)
upper_sim1 <- matrix(0,R,lent)

lower_sim2 <- matrix(0,R,lent)
upper_sim2 <- matrix(0,R,lent)

for(r in seq_along(files)){
  
  res <- readRDS(files[r])
  
  beta1 <- res$beta1
  beta2 <- res$beta2
  
  beta1_mat[r,] <- beta1
  beta2_mat[r,] <- beta2
  
  lower_pt1[r,] <- beta1 - res$length_pt1/2
  upper_pt1[r,] <- beta1 + res$length_pt1/2
  
  lower_pt2[r,] <- beta2 - res$length_pt2/2
  upper_pt2[r,] <- beta2 + res$length_pt2/2
  
  lower_sim1[r,] <- beta1 - res$length_sim1/2
  upper_sim1[r,] <- beta1 + res$length_sim1/2
  
  lower_sim2[r,] <- beta2 - res$length_sim2/2
  upper_sim2[r,] <- beta2 + res$length_sim2/2
}

# =====================================================
# 4 CP
# =====================================================

beta1_mean <- colMeans(beta1_mat)
beta2_mean <- colMeans(beta2_mat)

lower_pt1_mean <- colMeans(lower_pt1)
upper_pt1_mean <- colMeans(upper_pt1)

lower_pt2_mean <- colMeans(lower_pt2)
upper_pt2_mean <- colMeans(upper_pt2)

lower_sim1_mean <- colMeans(lower_sim1)
upper_sim1_mean <- colMeans(upper_sim1)

lower_sim2_mean <- colMeans(lower_sim2)
upper_sim2_mean <- colMeans(upper_sim2)

# =====================================================
# 5 beta1 figure
# =====================================================

df_plot1 <- data.frame(
  t=tobs,
  truth=beta1_true,
  estimate=beta1_mean
)

df_ci1 <- data.frame(
  t=t_pick,
  mean=beta1_mean[idx_pick],
  lower=lower_pt1_mean[idx_pick],
  upper=upper_pt1_mean[idx_pick]
)

df_sim1 <- data.frame(
  t=tobs,
  lower=lower_sim1_mean,
  upper=upper_sim1_mean
)
p_beta1 <- ggplot(df_plot1, aes(x = t)) +
  
  geom_ribbon(
    data = df_sim1,
    aes(ymin = lower, ymax = upper, fill = "Simul. CI"),
    alpha = 0.5
  ) +
  
  geom_line(
    aes(y = truth, linetype = "True"),
    linewidth = 0.8
  ) +
  
  geom_line(
    aes(y = estimate, color = "beta1"),
    linewidth = 0.8
  ) +
  
  geom_errorbar(
    data = df_ci1,
    aes(ymin = lower, ymax = upper, color = "Point. CI"),
    width = 0.03,
    linewidth = 0.7
  ) +
  
  scale_color_manual(
    values = c(
      "beta1" = "#FF6961",
      "Point. CI" = "#FF6961"
    ),
    labels = c(
      expression(beta[1](t)),
      "Point. CI"
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "Simul. CI" = "grey80"
    )
  ) +
  
  scale_linetype_manual(
    values = c(
      "True" = "dashed"
    )
  ) +
  
  labs(
    x = "t",
    y = expression(beta[1](t)),
    color = "",
    fill = "",
    linetype = ""
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_line(colour = "#E5E5E5", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold", size = 14)
  )

# =====================================================
# 6 beta2 figure
# =====================================================

df_plot2 <- data.frame(
  t=tobs,
  truth=beta2_true,
  estimate=beta2_mean
)

df_ci2 <- data.frame(
  t=t_pick,
  mean=beta2_mean[idx_pick],
  lower=lower_pt2_mean[idx_pick],
  upper=upper_pt2_mean[idx_pick]
)

df_sim2 <- data.frame(
  t=tobs,
  lower=lower_sim2_mean,
  upper=upper_sim2_mean
)
p_beta2 <- ggplot(df_plot2, aes(x = t)) +
  
  geom_ribbon(
    data = df_sim2,
    aes(ymin = lower, ymax = upper, fill = "Simul. CI"),
    alpha = 0.5
  ) +
  
  geom_line(
    aes(y = truth, linetype = "True"),
    linewidth = 0.8
  ) +
  
  geom_line(
    aes(y = estimate, color = "beta2"),
    linewidth = 0.8
  ) +
  
  geom_errorbar(
    data = df_ci2,
    aes(ymin = lower, ymax = upper, color = "Point. CI"),
    width = 0.03,
    linewidth = 0.7
  ) +
  
  scale_color_manual(
    values = c(
      "beta2" = "#836953",
      "Point. CI" = "#836953"
    ),
    labels = c(
      expression(beta[2](t)),
      "Point. CI"
    )
  ) +
  
  scale_fill_manual(
    values = c(
      "Simul. CI" = "grey80"
    )
  ) +
  
  scale_linetype_manual(
    values = c(
      "True" = "dashed"
    )
  ) +
  
  labs(
    x = "t",
    y = expression(beta[2](t)),
    color = "",
    fill = "",
    linetype = ""
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_line(colour = "#E5E5E5", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold", size = 14)
  )

# =====================================================
# 7 IMSE
# =====================================================

imse_list <- list()

for(n_val in n_list){
  
  files <- list.files(
    dir_prop,
    pattern=paste0("res_n",n_val,"_J5_rep.*\\.rds$"),
    full.names=TRUE
  )
  
  imse1 <- sapply(files,function(f) readRDS(f)$IMSE1)
  imse2 <- sapply(files,function(f) readRDS(f)$IMSE2)
  
  imse_list[[length(imse_list)+1]] <-
    data.frame(n=n_val,
               beta="beta1",
               mean=mean(imse1),
               sd=sd(imse1))
  
  imse_list[[length(imse_list)+1]] <-
    data.frame(n=n_val,
               beta="beta2",
               mean=mean(imse2),
               sd=sd(imse2))
}

df_imse <- bind_rows(imse_list)

df_imse$n_plot <- df_imse$n/10000

# =====================================================
# 8 IMSE beta1
# =====================================================

p_imse1 <- ggplot(
  df_imse %>%
    filter(beta == "beta1") %>%
    mutate(
      lower = (mean - sd)*1e3,
      upper = (mean + sd)*1e3,
      mean = mean*1e3
    ),
  aes(x = n_plot, y = mean)
) +
  
  geom_line(color = "#FF6961", linewidth = 0.8) +
  geom_point(color = "#FF6961", size = 1.5) +
  
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0,
    linewidth = 0.8,
    color = "#FF6961"
  ) +
  
  labs(
    x = expression(plain("Sample Size ") * n * "(×10"^4 * ")"),
    y = expression(IMSE~"(×10"^{-3}*")")
  ) +
  
  theme_bw(base_size = 14) +
  theme(axis.title = element_text(face = "bold"))

# =====================================================
# 9 IMSE beta2
# =====================================================

p_imse2 <- ggplot(
  df_imse %>%
    filter(beta == "beta2") %>%
    mutate(
      lower = (mean - sd)*1e3,
      upper = (mean + sd)*1e3,
      mean = mean*1e3
    ),
  aes(x = n_plot, y = mean)
) +
  
  geom_line(color = "#836953", linewidth = 0.8) +
  geom_point(color = "#836953", size = 1.5) +
  
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    width = 0,
    linewidth = 0.8,
    color = "#836953"
  ) +
  
  labs(
    x = expression(plain("Sample Size ") * n * "(×10"^4 * ")"),
    y = expression(IMSE~"(×10"^{-3}*")")
  ) +
  
  theme_bw(base_size = 14) +
  theme(axis.title = element_text(face = "bold"))

# =====================================================
# 10 attatch figure
# =====================================================

fig_beta1 <- (p_beta1 | p_imse1) &
  theme(legend.position = "right")

fig_beta2 <- (p_beta2 | p_imse2) &
  theme(legend.position = "right")

fig_beta1
fig_beta2




ggsave(
  filename = "./IMSE_1.pdf",
  plot = fig_beta1,
  width = 7,       
  height = 3,       
  device = "pdf"    
)

ggsave(
  filename = "./IMSE_2.pdf",
  plot = fig_beta2,
  width = 7,       
  height = 3,       
  device = "pdf"    
)
