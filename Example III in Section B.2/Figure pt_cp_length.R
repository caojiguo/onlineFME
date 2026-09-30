rm(list=ls())
library(gridExtra)
library(dplyr)
library(ggplot2)
library(patchwork)

# =====================================================
# 1 path
# =====================================================
dir_prop <- "./multi_linear/result"
n_list <- c(10000,20000,30000,40000,50000)

lent <- 101
tobs <- seq(0,1,length.out=lent)

# =====================================================
# 2 choose t
# =====================================================
t_pick  <- c(0.15,0.5,0.85)
idx_pick <- sapply(t_pick,function(tt) which.min(abs(tobs-tt)))

summary_list <- list()
length_list  <- list()
length_box <- list()

# =====================================================
# 3 result beta1
# =====================================================
for(n_val in n_list){
  cat("Processing n =", n_val, "\n")
  
  files <- list.files(
    dir_prop,
    pattern=paste0("res_n", n_val, "_J5_rep.*\\.rds$"),
    full.names=TRUE
  )
  
  R <- length(files)
  
  cover_mat1 <- matrix(0,R,lent)
  length_mat1 <- matrix(0,R,lent)
  
  for(r in seq_along(files)){
    res <- readRDS(files[r])
    
    cover_mat1[r,] <- res$cover_pt1
    length_mat1[r,] <- res$length_pt1
    
    for(j in seq_along(idx_pick)){
      idx <- idx_pick[j]
      length_box[[length(length_box)+1]] <- data.frame(
        n = n_val,
        t = t_pick[j],
        Fbeta = "b1",
        length = res$length_pt1[idx]
      )
    }
  }
  
  # summary
  for(j in seq_along(idx_pick)){
    idx <- idx_pick[j]
    vec <- cover_mat1[,idx]
    p <- mean(vec)
    se <- sqrt(p*(1-p)/R)
    
    summary_list[[length(summary_list)+1]] <- data.frame(
      n = n_val,
      t = t_pick[j],
      Fbeta = "b1",
      mean = p,
      lower = p - 1.96*se,
      upper = p + 1.96*se
    )
    
    length_list[[length(length_list)+1]] <- data.frame(
      n = n_val,
      t = t_pick[j],
      Fbeta = "b1",
      mean = mean(length_mat1[,idx])
    )
  }
}

# =====================================================
# 4 sort result
# =====================================================
df_summary <- bind_rows(summary_list)
df_summary$type <- "Coverage"

df_length <- bind_rows(length_list)
df_length$type <- "Length"
df_length$lower <- NA
df_length$upper <- NA

df_plot <- bind_rows(df_summary, df_length)
df_plot$n_plot <- df_plot$n/10000
df_plot$Fbeta <- factor(df_plot$Fbeta, levels=c("b1"))

df_length_box <- bind_rows(length_box)
df_length_box$n_plot <- df_length_box$n/10000
df_length_box$Fbeta <- factor(df_length_box$Fbeta, levels=c("b1"))

# =====================================================
# 5 facet 
# =====================================================
df_plot$t_lab <- factor(
  df_plot$t,
  levels=c(0.15,0.5,0.85),
  labels=c("t = 0.15","t = 0.5","t = 0.85")
)

df_length_box$t_lab <- factor(
  df_length_box$t,
  levels=c(0.15,0.5,0.85),
  labels=c("t = 0.15","t = 0.5","t = 0.85")
)

# =====================================================
# 6 color
# =====================================================
cols <- c("b1"="#FF6961")

# =====================================================
# 7 theme
# =====================================================
base_theme <- theme_bw(base_size=14) +
  theme(
    legend.position="right",
    legend.title = element_text(size=13),
    legend.text  = element_text(size=12),
    strip.background = element_rect(fill="#D9D9D9", colour=NA),
    strip.text = element_text(face="bold", size=13),
    panel.border = element_rect(colour="black", fill=NA, linewidth=0.8),
    panel.background = element_rect(fill="white"),
    panel.grid.major = element_line(colour="#E5E5E5", linewidth=0.5),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.2, "lines"),
    axis.text = element_text(size=12),
    axis.title = element_text(face="bold", color="#4D4D4D", size=14)
  )

# =====================================================
# 8 CP
# =====================================================
p1 <- ggplot(
  subset(df_plot,type=="Coverage"),
  aes(x=n_plot,y=mean,color=Fbeta,group=Fbeta)
) +
  geom_line(aes(linetype=Fbeta),linewidth=0.9) +
  geom_point(aes(shape=Fbeta),size=2.5) +
  geom_hline(yintercept=0.95, linetype="dashed", color="grey40") +
  facet_grid(. ~ t_lab) +
  scale_color_manual(values=cols, labels=c(expression(beta[1]))) +
  scale_linetype_manual(values=c("b1"="solid"), labels=c(expression(beta[1]))) +
  scale_shape_manual(values=c("b1"=16), labels=c(expression(beta[1]))) +
  labs(x=NULL, y="Point. CP", color="F-beta", linetype="F-beta", shape="F-beta") +
  coord_cartesian(ylim=c(0.6,1.0)) +
  base_theme

# =====================================================
# 9 Length
# =====================================================
p2 <- ggplot(
  df_length_box,
  aes(x=factor(n_plot),y=length,fill=Fbeta)
) +
  geom_boxplot(outlier.shape=NA, colour="black", linewidth=0.3, width=0.4) +
  facet_grid(. ~ t_lab) +
  scale_fill_manual(values=cols, labels=c(expression(beta[1]))) +
  labs(x=expression(plain("Sample Size ") * n * "(×10"^4*")"), y="Wid. of Point. CI", fill="F-beta") +
  coord_cartesian(ylim=c(0,0.15)) +
  base_theme +
  theme(strip.text.x = element_blank())  

# =====================================================
# 10 attatch figure
# =====================================================
final_plot <- (p1 / p2) & theme(legend.position = "none")
final_plot

ggsave(
  filename = "./1_pt_cp_length.pdf",
  plot = final_plot,
  width = 8,       
  height = 6,       
  device = "pdf"    
)




rm(list=ls())
library(gridExtra)
library(dplyr)
library(ggplot2)
library(patchwork)

# =====================================================
# 1 path
# =====================================================
dir_prop <- "./multi_linear/result"
n_list <- c(10000,20000,30000,40000,50000)

lent <- 101
tobs <- seq(0,1,length.out=lent)

# =====================================================
# 2 choose t
# =====================================================
t_pick  <- c(0.15,0.5,0.85)
idx_pick <- sapply(t_pick,function(tt) which.min(abs(tobs-tt)))

summary_list <- list()
length_list  <- list()
length_box <- list()

# =====================================================
# 3 result beta2
# =====================================================
for(n_val in n_list){
  cat("Processing n =", n_val, "\n")
  
  files <- list.files(
    dir_prop,
    pattern=paste0("res_n", n_val, "_J5_rep.*\\.rds$"),
    full.names=TRUE
  )
  
  R <- length(files)
  
  cover_mat2 <- matrix(0,R,lent)
  length_mat2 <- matrix(0,R,lent)
  
  for(r in seq_along(files)){
    res <- readRDS(files[r])
    
    cover_mat2[r,] <- res$cover_pt2
    length_mat2[r,] <- res$length_pt2
    
    for(j in seq_along(idx_pick)){
      idx <- idx_pick[j]
      length_box[[length(length_box)+1]] <- data.frame(
        n = n_val,
        t = t_pick[j],
        Fbeta = "b2",
        length = res$length_pt2[idx]
      )
    }
  }
  
  # summary
  for(j in seq_along(idx_pick)){
    idx <- idx_pick[j]
    vec <- cover_mat2[,idx]
    p <- mean(vec)
    se <- sqrt(p*(1-p)/R)
    
    summary_list[[length(summary_list)+1]] <- data.frame(
      n = n_val,
      t = t_pick[j],
      Fbeta = "b2",
      mean = p,
      lower = p - 1.96*se,
      upper = p + 1.96*se
    )
    
    length_list[[length(length_list)+1]] <- data.frame(
      n = n_val,
      t = t_pick[j],
      Fbeta = "b2",
      mean = mean(length_mat2[,idx])
    )
  }
}

# =====================================================
# 4 sort result
# =====================================================
df_summary <- bind_rows(summary_list)
df_summary$type <- "Coverage"

df_length <- bind_rows(length_list)
df_length$type <- "Length"
df_length$lower <- NA
df_length$upper <- NA

df_plot <- bind_rows(df_summary, df_length)
df_plot$n_plot <- df_plot$n/10000
df_plot$Fbeta <- factor(df_plot$Fbeta, levels=c("b2"))

df_length_box <- bind_rows(length_box)
df_length_box$n_plot <- df_length_box$n/10000
df_length_box$Fbeta <- factor(df_length_box$Fbeta, levels=c("b2"))

# =====================================================
# 5 facet
# =====================================================
df_plot$t_lab <- factor(
  df_plot$t,
  levels=c(0.15,0.5,0.85),
  labels=c("t = 0.15","t = 0.5","t = 0.85")
)

df_length_box$t_lab <- factor(
  df_length_box$t,
  levels=c(0.15,0.5,0.85),
  labels=c("t = 0.15","t = 0.5","t = 0.85")
)

# =====================================================
# 6 color
# =====================================================
cols <- c("b2"="#836953")

# =====================================================
# 7 theme
# =====================================================
base_theme <- theme_bw(base_size=14) +
  theme(
    legend.position="right",
    legend.title = element_text(size=13),
    legend.text  = element_text(size=12),
    strip.background = element_rect(fill="#D9D9D9", colour=NA),
    strip.text = element_text(face="bold", size=13),
    panel.border = element_rect(colour="black", fill=NA, linewidth=0.8),
    panel.background = element_rect(fill="white"),
    panel.grid.major = element_line(colour="#E5E5E5", linewidth=0.5),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.2, "lines"),
    axis.text = element_text(size=12),
    axis.title = element_text(face="bold", color="#4D4D4D", size=14)
  )

# =====================================================
# 8 Coverage figure
# =====================================================
p1 <- ggplot(
  subset(df_plot,type=="Coverage"),
  aes(x=n_plot,y=mean,color=Fbeta,group=Fbeta)
) +
  geom_line(aes(linetype=Fbeta),linewidth=0.9) +
  geom_point(aes(shape=Fbeta),size=2.5) +
  geom_hline(yintercept=0.95, linetype="dashed", color="grey40") +
  facet_grid(. ~ t_lab) +
  scale_color_manual(values=cols, labels=c(expression(beta[2]))) +
  scale_linetype_manual(values=c("b2"="solid"), labels=c(expression(beta[2]))) +
  scale_shape_manual(values=c("b2"=16), labels=c(expression(beta[2]))) +
  labs(x=NULL, y="Point. CP", color="F-beta", linetype="F-beta", shape="F-beta") +
  coord_cartesian(ylim=c(0.6,1.0)) +
  base_theme

# =====================================================
# 9 Length figure
# =====================================================
p2 <- ggplot(
  df_length_box,
  aes(x=factor(n_plot),y=length,fill=Fbeta)
) +
  geom_boxplot(outlier.shape=NA, colour="black", linewidth=0.3, width=0.4) +
  facet_grid(. ~ t_lab) +
  scale_fill_manual(values=cols, labels=c(expression(beta[2]))) +
  labs(x=expression(plain("Sample Size ") * n * "(×10"^4*")"), y="Wid. of Point. CI", fill="F-beta") +
  coord_cartesian(ylim=c(0,0.15)) +
  base_theme +
  theme(strip.text.x = element_blank())

# =====================================================
# 10 attatch figure
# =====================================================
final_plot <- (p1 / p2) & theme(legend.position = "none")
final_plot

ggsave(
  filename = "./2_pt_cp_length.pdf",
  plot = final_plot,
  width = 8,
  height = 6,
  device = "pdf"
)








