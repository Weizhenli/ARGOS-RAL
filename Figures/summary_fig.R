## ----setup, include=FALSE--------------------------------------------------------------
knitr::opts_chunk$set(echo = FALSE)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(scales)
library(cowplot)
library(ggh4x)
library(stringr)
library(ggpubr)
library(scales)
library(gridExtra)
library(R.matlab)
library(tibble)
# library(hrbrthemes)
library(dplyr)
library(plot3D)
library(stringr);
library(ggplot2);
library(ggpubr);
library(scales);
library(gridExtra);
library(tidyr)
library(grid)
library(latex2exp);
## main -------------------------------
solution_theme <- theme_minimal(18)+
  theme(plot.margin = unit(c(5,0,5,0),'pt'),legend.box.margin= unit(c(0,0,0,-20),'pt'),
        plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 16),
        plot.tag = element_text(size = 28, face = "bold", vjust = -3))

titles_plot <- function(plot_title, plot_lab, plot_stitle, title_size=18,eq_bd_h=0.6,eq_bd_w=0.4, ratio=c(2,3)){
  # plot_stitle can be expression
  title <- textGrob(plot_title,x=0.7,gp=gpar(col="black",fontsize=title_size,fontface ='bold'))
  tag <- textGrob(plot_lab,x=0.05,gp=gpar(col="black",fontsize=title_size+4,fontface ='bold'))
  pde_eq <- grobTree(roundrectGrob(gp=gpar(fill="grey90",col="grey90"),width=eq_bd_w,height=eq_bd_h,r=unit(0.3, "snpc")),
                     textGrob(plot_stitle, x=0.5, gp=gpar(col="black", fontsize=title_size-4)))
  # return(arrangeGrob(grobTree(tag,title), pde_eq, nrow=2, heights=ratio))
  return(arrangeGrob(grobTree(tag,title), pde_eq, nrow=1, widths=ratio))
  # return(grobTree(tag,title, pde_eq))
}
# layout_matrix <- matrix(c(rep(1,10),rep(2,12),rep(3,1)),ncol=1)
ggplot_heights <- c(15,12)
rect1<- data.frame(xmin=-Inf, xmax=Inf, ymin=0.8, ymax=Inf)
# layout_matrix_all <- matrix(c(rep(1,15),rep(3,1)),ncol=1)
base_text_size <- 16
plot_theme <- theme(
  axis.line = element_line(colour = "black"),
  axis.ticks.length = unit(5, "pt"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  # panel.grid.major.y = element_line(color='grey85'),
  legend.key = element_blank(),
  axis.text = element_text(size = base_text_size-2),
  axis.title.x = element_text(size =base_text_size),
  axis.title.y = element_text(
    size = base_text_size,
    #angle = 90,
    vjust = 0.5
  ),
  plot.title = element_text(size = base_text_size+2),
  legend.position = "none",
  legend.text.align = 0,
  legend.title = element_text(size = base_text_size+2),
  legend.text = element_text(size = base_text_size),
  legend.key.size = unit(18,'pt'),
  plot.margin = unit(c(0,5,5,5),'pt')
)
legend_guides <- guides(color = guide_legend(override.aes = list(size = 3, lwd=1.2)))

color_sacle <- c(ramp.col(col = c("#3B4CC0", "white"),n=100),ramp.col(col = c("white", "#B40426"),n=100))
create_separators <- function(x, extra_x, y, extra_y, angle = 45, scale = 1, length = .1){
  add_y <-  length * sin(angle * pi/180) /2
  add_x <- length * cos(angle * pi/180) 
  list(x = x - add_x*scale, xend = x + add_x*scale + extra_x, 
       y = rep(y - add_y*scale - extra_y, length(x)), yend = rep(y + add_y*scale - extra_y/2, length(x)))
}
y_labels = y_breaks <- seq(0, 1, 0.25)
xstart <- 41.5
xend <- 42.5
extra_x <- 1
y_sep=0
myseg <- create_separators(c(xstart, xend), extra_x = 1, y = y_sep, extra_y = 0.1, angle = 75)
myseg2 <- create_separators(c(xstart+20, xend+20), extra_x = 1, y = y_sep, extra_y = 0.1, angle = 75)
ggplot_heights2 <- c(10,1)
turn_off_legend <- function(p_ggplot){
  p_ggplot+theme(legend.position = "none")
}
data_path <- "C:/Users/cfzh32/Documents/GitHub/Weizhen-Li/paper_test/"
# data_path <- "~/GitHub/Weizhen-Li/paper_test/"
# data_path <- "D:/GitHub/Weizhen-Li/paper_test/"
setwd(paste0(data_path,".."))
shapes <- c(21,22,24,25)
# sizes <- c(6,5,4,4)-2
sizes <- c(3.5,2.5,2,2)
alpha=alpha_value <- 0.6
stroke_value <- 0.5
lwd_size=0.7
legend_labs <- c('ARGOS-RAL',expression(STRidge ~~ (d[tol] == 0.2)), expression(STRidge ~~ (d[tol] == 2)), expression(STRidge ~~ (d[tol] == 10)))
load('paper_test/more_test2/pde_data1.RData')
colors <- c("#F1766D","#21F5D7", "#56ABED", "#6261FF")
## ----bur plot, fig.height=4, fig.width=8-----------------------------------------------
## Burgers plots --------------------------------
path.burgers <- ('data/burgers.mat')
burgers.mat <- readMat(path.burgers)
burgers.t <- as.numeric(burgers.mat[['t']])
burgers.x <- as.numeric(burgers.mat[['x']])
burgers.usol <- as.data.frame(Re(burgers.mat[['u']]))

gg_bur_data <- cbind.data.frame(x = rep(burgers.x,length(burgers.t)), t = rep(burgers.t,each=length(burgers.x)), u = gather(burgers.usol)[,-1])
color_sacle <- ramp.col(col = c("#fafafa", "#d7191c"),n=200)
bur_sol_plot <- ggplot(gg_bur_data, aes(x, t, fill= u)) + geom_tile()+
  scale_fill_gradientn(colors = color_sacle) +
  # scale_fill_distiller(palette = "Spectral") + 
  # labs(title='Burgers', tag='A', subtitle = TeX('$u_t= -uu_x + 0.1 u_{xx}$')) +
  # labs(title='Burgers', tag='A') +
  solution_theme+
  # annotate('label',x=0,y=12.1,label='                                    ',fill = 'grey90',label.size=NA,size=5.7)+
  # annotate('label',x=5,y=1,label=TeX('$u_t= -uu_x + 0.1 u_{xx}$'),fill = 'grey90',label.size=NA,size=5.7)+
  coord_cartesian(clip="off")

## bur noisy data plot ---------------------------
load(paste0(data_path, "bur_SG_noise_seed_100_samp_100_snr_ada_lasso_pareto_AIC.RData"))
load(paste0(data_path, "bur_SG_noise_seed_100_samp_100_snr_rudy_d_thred.RData"))
d_tols <- c(0.2,2,10)
names(bur_density_noise_rudy) <- paste0("rudy_", d_tols)
names(bur_density_noise_ada_lasso_pareto) <- 'our'
dens_snr_all <- c(bur_density_noise_ada_lasso_pareto, bur_density_noise_rudy)

snr_db_seq <- c(seq(0, 40, 2), Inf) # SNR_dB grid
eta <- 10^(-snr_db_seq/20)
snr_db_seq1 <- c(seq(0, 40, 2), 45)
snr_db_seq2 <- c(seq(0, 40, 2),  45) # c(seq(20, 40, 2), 42.5, 45)
snr_db_seq3 <- c(as.character(seq(0, 40, 2)), expression(infinity)) # c(as.character(seq(20, 40, 2)), "...", "Inf")
# snr_db_seq3[c(1,3,5,7,9,11,13,15,17,19)+1] <- ''
snr_db_seq3[-(3*(0:(length(snr_db_seq)/2))+1)] <- ''

bur_count_vector_snr <- data.frame(SNR_dB = 0, Correct = 0, Incorrect = 0, method0 = '0')
for(m_names in names(dens_snr_all)){
  temp <- data.frame(SNR_dB = snr_db_seq1,
                     Correct = colSums(dens_snr_all[[m_names]] == 1) / 100,
                     Incorrect = colSums(dens_snr_all[[m_names]] == 0) / 100,
                     method0 = m_names)
  
  bur_count_vector_snr <- rbind.data.frame(bur_count_vector_snr, temp)
}
bur_count_vector_snr <- bur_count_vector_snr[-1, ]
rownames(bur_count_vector_snr) <- seq(1, nrow(bur_count_vector_snr))
bur_noise_gather_df <- bur_count_vector_snr %>%
  gather("Condition", "Value",
         3:ncol(bur_count_vector_snr)-1)

method_label <- c("RAL", paste0("STRidge (d_tol=", d_tols,")"))
bur_count_vector_snr$method <- NA
for(i in seq_along(method_label)){
  bur_count_vector_snr$method[which(bur_count_vector_snr$method0==unique(bur_count_vector_snr$method0)[i])] <- method_label[i]
}

x_limits <- snr_db_seq
x_breaks <- pretty_breaks()(snr_db_seq)
# colours <- c("#a670b0",
#              "#a58d48")
# ggcol <- c("#F8766D","#00BA38","#619CFF")
# rudy_color <- colorRampPalette(c("#4bf5db","#6261ff"))(length(d_tols))
# colors <- c(ggcol[1],rudy_color)

bur_SNR_plot <- ggplot() + 
  geom_hline(yintercept = 0.8, lty=2)+
  # geom_rect(aes(xmin=rect1$xmin,xmax=rect1$xmax,ymin=rect1$ymin,ymax=rect1$ymax),fill="#9de0e6",alpha=0.01)+
  geom_rect(data=rect1,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),alpha=0.2,fill="#9de0e6")+
  scale_y_continuous(breaks = seq(0,1,0.2))+
  geom_line(data=bur_count_vector_snr, aes(x=SNR_dB, y=Correct,color=method),lwd=lwd_size) +
  geom_point(data=bur_count_vector_snr,alpha=alpha_value,stroke=stroke_value,aes(x=SNR_dB, y=Correct, shape=method,size=method,fill = method)) + 
  # ylab("Probability of correct identification") +  xlab("SNR(dB)")+
  ylab("Success rate") +  xlab("SNR(dB)")+
  # xlab(TeX("Noise magnitude ($\\eta$*sd(u))")) +
  scale_x_continuous(breaks = snr_db_seq2, labels = snr_db_seq3)+
  scale_size_manual(values=sizes, breaks=unique(bur_count_vector_snr$method),
                    labels = legend_labs)+
  scale_fill_manual(values=colors, breaks=unique(bur_count_vector_snr$method),
                    labels = legend_labs)+
  scale_shape_manual(values=shapes, breaks=unique(bur_count_vector_snr$method),
                     labels = legend_labs)+
  scale_color_manual(values=colors, breaks=unique(bur_count_vector_snr$method),
                     labels = legend_labs)+
  labs(fill = 'Methods', shape='Methods', size='Methods',color='Methods')+
  # scale_y_continuous(labels = y_labels,breaks = y_breaks,limits = c(NA, max(y_breaks))) +
  plot_theme + legend_guides +
  guides(x = guide_axis_truncated(
    trunc_lower = c(-Inf, xend + extra_x/2),
    trunc_upper = c(xstart + extra_x/2, Inf)
  )) +
  annotate("segment", 
           x = myseg$x, xend = myseg$xend,
           y = myseg$y + 0.05, yend = myseg$yend)  +
  coord_cartesian(clip = "off", ylim = c(-0.0005, NA))

## bur noiseless plot ----------------------
load(paste0(data_path,"bur_SG_noiseless_seed_10_success_rate_ada_lasso_pareto_AIC.RData"))
load(paste0(data_path,"bur_SG_noiseless_seed_10_success_rate_rudy_d_thred.RData"))
num <- (seq(2,4.2,0.2))

d_tols <- c(0.2,2,10)
names(success_rate_rudy_list) <- paste0("rudy_", d_tols)
success_rate_our <- bur_density_noiseless_ada_lasso
dens_n_all <- c(list(our=success_rate_our), lapply(success_rate_rudy_list, function(x) x))

bur_count_vector_n <- data.frame(n = 0, Correct = 0, Incorrect = 0, method0 = '0')
for(m_names in names(dens_n_all)){
  temp <- data.frame(n = num,
                     Correct = colSums(dens_n_all[[m_names]] == 1) / 100,
                     Incorrect = colSums(dens_n_all[[m_names]] == 0) / 100,
                     method0 = m_names)
  
  bur_count_vector_n <- rbind.data.frame(bur_count_vector_n, temp)
}
bur_count_vector_n <- bur_count_vector_n[-1, ]

rownames(bur_count_vector_n) <- seq(1, nrow(bur_count_vector_n))
noise_gather_df_all <- bur_count_vector_n %>%
  gather("Condition", "Value",
         3:ncol(bur_count_vector_n)-1)

method_label <- c("RAL", paste0("STRidge (d_tol=", d_tols,")"))
bur_count_vector_n$method <- NA
for(i in seq_along(method_label)){
  bur_count_vector_n$method[which(bur_count_vector_n$method0==unique(bur_count_vector_n$method0)[i])] <- method_label[i]
}

colours <- c("#a670b0",
             "#a58d48")
ggcol <- c("#F8766D","#00BA38","#619CFF")
rudy_color <- colorRampPalette(c("#4bf5db","#6261ff"))(length(d_tols))
colors <- c(ggcol[1],rudy_color)
colors <- c("#F1766D","#21F5D7", "#56ABED", "#6261FF")

xlables <- unname(latex2exp::TeX(paste0('$10^{',round((num),1),'}$')))
# xlables[c(1,3,5,7,9)+1] <- ''
xlables[-(3*(0:(length(num)/2))+1)] <- ''
second_xlables <- round(10^num/(256*101)*100,2)
second_xlables[-(3*(0:(length(num)/2))+1)] <- ''

bur_N_plot <- ggplot() + 
  geom_hline(yintercept = 0.8, lty=2)+
  geom_rect(aes(xmin=rect1$xmin,xmax=rect1$xmax,ymin=rect1$ymin,ymax=rect1$ymax),alpha=0.2,fill="#9de0e6")+
  scale_y_continuous(breaks = seq(0,1,0.2))+
  geom_line(data=bur_count_vector_n, aes(x=n, y=Correct,color=method),lwd=lwd_size) +
  geom_point(data=bur_count_vector_n,alpha=alpha_value,stroke=stroke_value,aes(x=n, y=Correct,shape=method,size=method,fill = method)) + 
  scale_x_continuous(breaks = num, labels = xlables,
                     sec.axis = dup_axis(name='Percentage (%)', labels=second_xlables))+ #geom_line() + 
  # labs(x="SNR(dB)", y="Success rate", title = TeX('Reaction-Diffusion, $\\alpha$=0.3, Inf $\\alpha$=0.5'))+
  # labs(x="N", y="Success rate")+
  labs(x="N", y=" ")+
  scale_size_manual(values=sizes, breaks=unique(bur_count_vector_n$method),
                    labels = legend_labs)+
  scale_fill_manual(values=colors, breaks=unique(bur_count_vector_n$method),
                    labels = legend_labs)+
  scale_shape_manual(values=shapes, breaks=unique(bur_count_vector_n$method),
                     labels = legend_labs)+
  scale_color_manual(values=colors, breaks=unique(bur_count_vector_snr$method),
                     labels = legend_labs)+
  labs(fill = 'Methods', shape='Methods', size='Methods',color='Methods')+
  # scale_y_continuous(labels = y_labels,breaks = y_breaks,limits = c(NA, max(y_breaks))) +
  plot_theme + legend_guides
# theme_classic(18)
## out burgers -----------------------
legend_n <- get_legend(bur_N_plot+theme(legend.position='bottom'))
bur_plot <- arrangeGrob(bur_SNR_plot+guides(color='none')+theme(legend.position = 'none'), 
                        bur_N_plot+guides(color='none')+theme(legend.position = 'none'),
                        nrow=1,ncol=2)
bur_plot1 <- arrangeGrob(bur_plot,legend_n,nrow=2, heights = ggplot_heights2)
bur_title <- titles_plot('Burgers','A',expression(u[t] == -uu[x]+0.1*u[xx]))
bur_plot2 <- arrangeGrob(bur_title, bur_sol_plot, bur_plot1, nrow=3,heights =c(2,7,10))

# ggsave(bur_plot2, filename = 'Figures/new_figs/bur.png', width = 8, height = 8)
# ggsave(bur_plot2, filename = 'Figures/new_figs/bur.pdf', width = 8, height = 8)

## ----ad plot, fig.height=4, fig.width=8------------------------------------------------
## ad plots --------------------------------
ad.x <- ad_list$x; ad.t <- ad_list$t; ad.usol <- as.data.frame(ad_list$sol)
gg_ad_data <- cbind.data.frame(x = rep(ad.x,length(ad.t)), t = rep(ad.t,each=length(ad.x)), c = gather(ad.usol)[,-1])
color_sacle <- ramp.col(col = c("#fafafa", "#d7191c"),n=200)
ad_sol_plot <- ggplot(gg_ad_data, aes(x, t, fill= c)) + geom_tile()+
  scale_fill_gradientn(colors = color_sacle) +
  # scale_fill_distiller(palette = "Spectral") + 
  # labs(title='adgers', tag='A', subtitle = TeX('$u_t= -uu_x + 0.1 u_{xx}$')) +
  # labs(title='adgers', tag='A') +
  solution_theme+
  # annotate('label',x=0,y=12.1,label='                                    ',fill = 'grey90',label.size=NA,size=5.7)+
  # annotate('label',x=5,y=1,label=TeX('$u_t= u_{xx} - u_{x}$'),fill = 'grey90',label.size=NA,size=5.7)+
  coord_cartesian(clip="off")

## ad noisy data plot ---------------------------
load(paste0(data_path, "more_test2/ad_SG_noise_seed_100_samp_100_snr_ada_lasso_pareto.RData"))
load(paste0(data_path, "more_test2/ad_SG_noise_seed_100_samp_100_snr_rudy_d_thred.RData"))
d_tols <- c(0.2,2,10)
names(ad_density_noise_rudy) <- paste0("rudy_", d_tols)
names(ad_density_noise_ada_lasso_pareto) <- 'our'
dens_snr_all <- c(ad_density_noise_ada_lasso_pareto, ad_density_noise_rudy)

snr_db_seq <- c(seq(0, 40, 2), Inf) # SNR_dB grid
eta <- 10^(-snr_db_seq/20)
snr_db_seq1 <- c(seq(0, 40, 2), 45)
snr_db_seq2 <- c(seq(0, 40, 2),  45) # c(seq(20, 40, 2), 42.5, 45)
snr_db_seq3 <- c(as.character(seq(0, 40, 2)), expression(infinity)) # c(as.character(seq(20, 40, 2)), "...", "Inf")
# snr_db_seq3[c(1,3,5,7,9,11,13,15,17,19)+1] <- ''
snr_db_seq3[-(3*(0:(length(snr_db_seq)/2))+1)] <- ''
snr_db_seq3[length(snr_db_seq3)] <- expression(infinity)

ad_count_vector_snr <- data.frame(SNR_dB = 0, Correct = 0, Incorrect = 0, method0 = '0')
for(m_names in names(dens_snr_all)){
  temp <- data.frame(SNR_dB = snr_db_seq1,
                     Correct = colSums(dens_snr_all[[m_names]] == 1) / 100,
                     Incorrect = colSums(dens_snr_all[[m_names]] == 0) / 100,
                     method0 = m_names)
  
  ad_count_vector_snr <- rbind.data.frame(ad_count_vector_snr, temp)
}
ad_count_vector_snr <- ad_count_vector_snr[-1, ]
rownames(ad_count_vector_snr) <- seq(1, nrow(ad_count_vector_snr))
ad_noise_gather_df <- ad_count_vector_snr %>%
  gather("Condition", "Value",
         3:ncol(ad_count_vector_snr)-1)

method_label <- c("RAL", paste0("STRidge (d_tol=", d_tols,")"))
ad_count_vector_snr$method <- NA
for(i in seq_along(method_label)){
  ad_count_vector_snr$method[which(ad_count_vector_snr$method0==unique(ad_count_vector_snr$method0)[i])] <- method_label[i]
}

x_limits <- snr_db_seq
x_breaks <- pretty_breaks()(snr_db_seq)
# colours <- c("#a670b0",
#              "#a58d48")
# ggcol <- c("#F8766D","#00BA38","#619CFF")
# rudy_color <- colorRampPalette(c("#4bf5db","#6261ff"))(length(d_tols))
# colors <- c(ggcol[1],rudy_color)

ad_SNR_plot <- ggplot() + 
  geom_hline(yintercept = 0.8, lty=2)+
  geom_rect(aes(xmin=rect1$xmin,xmax=rect1$xmax,ymin=rect1$ymin,ymax=rect1$ymax),alpha=0.2,fill="#9de0e6")+
  scale_y_continuous(breaks = seq(0,1,0.2))+
  geom_line(data=ad_count_vector_snr, aes(x=SNR_dB, y=Correct,color=method),lwd=lwd_size) +
  geom_point(data=ad_count_vector_snr,alpha=alpha_value,stroke=stroke_value,aes(x=SNR_dB, y=Correct,shape=method,size=method,fill = method)) + 
  # ylab("Probability of correct identification") +  xlab("SNR(dB)")+
  ylab("Success rate") +  xlab("SNR(dB)")+
  # xlab(TeX("Noise magnitude ($\\eta$*sd(u))")) +
  scale_x_continuous(breaks = snr_db_seq2, labels = snr_db_seq3)+
  scale_size_manual(values=sizes, breaks=unique(ad_count_vector_snr$method),
                    labels = legend_labs)+
  scale_fill_manual(values=colors, breaks=unique(ad_count_vector_snr$method),
                    labels = legend_labs)+
  scale_shape_manual(values=shapes, breaks=unique(ad_count_vector_snr$method),
                     labels = legend_labs)+
  scale_color_manual(values=colors, breaks=unique(ad_count_vector_snr$method),
                     labels = legend_labs)+
  labs(fill = 'Methods', shape='Methods', size='Methods',color='Methods')+
  # scale_y_continuous(labels = y_labels,breaks = y_breaks,limits = c(NA, max(y_breaks))) +
  plot_theme + legend_guides +
  guides(x = guide_axis_truncated(
    trunc_lower = c(-Inf, xend + extra_x/2),
    trunc_upper = c(xstart + extra_x/2, Inf)
  )) +
  annotate("segment", 
           x = myseg$x, xend = myseg$xend,
           y = myseg$y + 0.05, yend = myseg$yend)  +
  coord_cartesian(clip = "off", ylim = c(-0.0005, NA))

## ad noiseless plot ----------------------
load(paste0(data_path,"more_test2/ad_SG_noiseless_seed_10_success_rate_ada_lasso_pareto.RData"))
load(paste0(data_path,"more_test2/ad_SG_noiseless_seed_10_success_rate_rudy_d_thred.RData"))
num <- (seq(2,5.2,0.2))

d_tols <- c(0.2,2,10)
names(success_rate_rudy_list) <- paste0("rudy_", d_tols)
success_rate_our <- ad_density_noiseless_ada_lasso
dens_n_all <- c(list(our=success_rate_our), lapply(success_rate_rudy_list, function(x) x))

ad_count_vector_n <- data.frame(n = 0, Correct = 0, Incorrect = 0, method0 = '0')
for(m_names in names(dens_n_all)){
  temp <- data.frame(n = num,
                     Correct = colSums(dens_n_all[[m_names]] == 1) / 100,
                     Incorrect = colSums(dens_n_all[[m_names]] == 0) / 100,
                     method0 = m_names)
  
  ad_count_vector_n <- rbind.data.frame(ad_count_vector_n, temp)
}
ad_count_vector_n <- ad_count_vector_n[-1, ]

rownames(ad_count_vector_n) <- seq(1, nrow(ad_count_vector_n))
noise_gather_df_all <- ad_count_vector_n %>%
  gather("Condition", "Value",
         3:ncol(ad_count_vector_n)-1)

method_label <- c("RAL", paste0("STRidge (d_tol=", d_tols,")"))
ad_count_vector_n$method <- NA
for(i in seq_along(method_label)){
  ad_count_vector_n$method[which(ad_count_vector_n$method0==unique(ad_count_vector_n$method0)[i])] <- method_label[i]
}

# colours <- c("#a670b0",
#              "#a58d48")
# ggcol <- c("#F8766D","#00BA38","#619CFF")
# rudy_color <- colorRampPalette(c("#4bf5db","#6261ff"))(length(d_tols))
# colors <- c(ggcol[1],rudy_color)

xlables <- unname(latex2exp::TeX(paste0('$10^{',round((num),1),'}$')))
# xlables[c(1,3,5,7,9)+1] <- ''
xlables[-(3*(0:(length(num)/2))+1)] <- ''
second_xlables <- round(10^num/(201*1001)*100,2)
second_xlables[-(3*(0:(length(num)/2))+1)] <- ''

ad_N_plot <- ggplot() + 
  geom_hline(yintercept = 0.8, lty=2)+
  geom_rect(aes(xmin=rect1$xmin,xmax=rect1$xmax,ymin=rect1$ymin,ymax=rect1$ymax),alpha=0.2,fill="#9de0e6")+
  scale_y_continuous(breaks = seq(0,1,0.2))+
  geom_line(data=ad_count_vector_n, aes(x=n, y=Correct,color=method),lwd=lwd_size) +
  geom_point(data=ad_count_vector_n,alpha=alpha_value,stroke=stroke_value,aes(x=n, y=Correct,shape=method,size=method,fill = method)) + 
  scale_x_continuous(breaks = num, labels = xlables,
                     sec.axis = dup_axis(name='Percentage (%)', labels=second_xlables))+ 
  # labs(x="SNR(dB)", y="Success rate", title = TeX('Reaction-Diffusion, $\\alpha$=0.3, Inf $\\alpha$=0.5'))+
  # labs(x="N", y="Success rate")+
  labs(x="N", y=" ")+
  scale_size_manual(values=sizes, breaks=unique(ad_count_vector_n$method),
                    labels = legend_labs)+
  scale_fill_manual(values=colors, breaks=unique(ad_count_vector_n$method),
                    labels = legend_labs)+
  scale_shape_manual(values=shapes, breaks=unique(ad_count_vector_n$method),
                     labels = legend_labs)+
  scale_color_manual(values=colors, breaks=unique(ad_count_vector_n$method),
                     labels = legend_labs)+
  labs(fill = 'Methods', shape='Methods', size='Methods',color='Methods')+
  # scale_y_continuous(labels = y_labels,breaks = y_breaks,limits = c(NA, max(y_breaks))) +
  plot_theme + legend_guides
# theme_classic(18)
## out adv -----------------------
legend_n <- get_legend(ad_N_plot+theme(legend.position='bottom'))
ad_plot <- arrangeGrob(ad_SNR_plot+guides(color='none')+theme(legend.position = 'none'), 
                        ad_N_plot+guides(color='none')+theme(legend.position = 'none'),
                        nrow=1,ncol=2)
ad_plot1 <- arrangeGrob(ad_plot,legend_n, nrow=2, heights = ggplot_heights2)
ad_title <- titles_plot('Advection-diffusion','A',expression(c[t] == c[xx]-c[x]))
ad_plot2 <- arrangeGrob(ad_title, ad_sol_plot,ad_plot1,nrow=3,heights =c(2,7,10))

# ggsave(ad_plot2, filename = 'Figures/new_figs/ad.png', width = 8, height = 8)
# ggsave(ad_plot2, filename = 'Figures/new_figs/ad.pdf', width = 8, height = 8)

## ----NS plot , fig.height=4, fig.width=8-----------------------------------------------
# NS plots --------------------------
path = './data/ibpm15300.plt'
library(ggplot2);library(plot3D);library(RColorBrewer)
gre <- read.table(path,sep='',header=F,fill=T,skip=6,
                  col.names=c('x','y','u','v','Vorticity'))
n <- 449
m <- 199
mid<-mean(gre[,'Vorticity'])
dt = 0.2
dx = 0.02
dy = 0.02
# Cut out the portion of the data before the cylinder
xmin = 100*dx
xmax = 425*dx
ymin = 15*dy
ymax = 185*dy
sam_data <- cbind.data.frame(x=runif(500,xmin+0.03,xmax-0.03),y=runif(500,ymin+0.03,ymax-0.03)) # make the plots looks better

# The center of mass of the cylinder is at (1,2)
plot_data <- cbind.data.frame(x=gre[,'x']+1,y=gre[,'y']+2,w=gre[,'Vorticity'])
plot_data$w[which(plot_data$w > 3)] <- 3
plot_data$w[which(plot_data$w < -3)] <- -3

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
# color_sacle <- c(ramp.col(col = c("#2c7bb6", "white"),n=100),ramp.col(col = c("white", "#d7191c"),n=100))
color_sacle <- c(ramp.col(col = c("#3288BD", "white"),n=100),ramp.col(col = c("white", "#d7191c"),n=100))
NS_sol_plot <- ggplot()+geom_raster(data=plot_data,aes(x=x,y=y,fill=w))+
  # scale_fill_distiller(palette = "Spectral",limits = c(-3,3))+
  scale_fill_gradientn(colors=color_sacle,limits = c(-3,3),breaks=c(-3,0,3))+
  # labs(title='Navier-Stokes', tag='B', subtitle = TeX('$\\omega_t=0.01\\omega_{xx}+0.01\\omega_{yy}-u\\omega_x-v\\omega_y$')) + 
  geom_segment(mapping=aes(x=c(xmin,xmin,xmin,xmax),y=c(ymin,ymax,ymin,ymin),
                           xend=c(xmax,xmax,xmin,xmax),yend=c(ymin,ymax,ymax,ymax)),
               color=2)+
  solution_theme + guides(fill = guide_colourbar(nbin = 3))+
  # annotate('label',x=4.5,y=4.9,label='                                        ',
  #          fill = 'grey90',label.size=NA,size=9)+
  coord_cartesian(clip="off")
## NS noise data plot -------------------------
load(paste0(data_path,"NS_SG_noise_density_seed_10_snr_ada_lasso_pareto.RData"))
load(paste0(data_path,"NS_SG_noise_seed_10_rudy_d_thred.RData"))

dens_fun <- function(noise_coeff, true_terms = c("w_{xx}", "w_{yy}", "uw_{x}", "vw_{y}")){
  dens_noise <- sapply(seq_along(noise_coeff), function(i){
    match_index <- match(names(noise_coeff[[i]]), true_terms)
    if(length(match_index) == 4 & all(!is.na(match_index))){
      return(1)
    }else{
      return(0)
    }
  })
  return(dens_noise)
}

noise = snr_db_seq <- c(seq(10, 40, 2), Inf) # SNR_dB grid
eta <- 10^(-snr_db_seq/20)
snr_db_seq1 <- c(seq(10, 40, 2), 45)
snr_db_seq2 <- c(seq(10, 40, 2),  45) # c(seq(20, 40, 2), 42.5, 45)
snr_db_seq3 <- c(as.character(seq(10, 40, 2)), expression(infinity)) # c(as.character(seq(20, 40, 2)), "...", "Inf")
# snr_db_seq3[c(2,4,6,8,10,12,14,16,18,20)] <- ''
snr_db_seq3[-(3*(0:(length(snr_db_seq)/2))+1)] <- ''
# snr_db_seq3 <- snr_db_seq3[-c(1:6)]
snr_db_seq3[length(snr_db_seq3)] <- expression(infinity)

dens_our <- cbind(0,sapply(seq_along(noise)[-1], function(i) dens_fun(NS_noise_ada_lasso_pareto[[i]])))
dens_rudy_list <- lapply(NS_noise_rudy_list, function(x) cbind(0,sapply(seq_along(noise)[-1], function(i) dens_fun(x[[i]]))))
d_tols <- c(0.2,2,10)
names(dens_rudy_list) <- paste0("rudy_", d_tols)
dens_all_noise <- c(list(our = dens_our), dens_rudy_list)

count_vector_snr <- data.frame(SNR_dB = 0, Correct = 0, Incorrect = 0, method0 = '0')
for(m_names in names(dens_all_noise)){
  temp <- data.frame(SNR_dB = snr_db_seq1,
                     Correct = colSums(dens_all_noise[[m_names]] == 1) / 100,
                     Incorrect = colSums(dens_all_noise[[m_names]] == 0) / 100,
                     method0 = m_names)
  count_vector_snr <- rbind.data.frame(count_vector_snr, temp)
}
count_vector_snr <- count_vector_snr[-1, ]

method_label <- c("RAL", paste0("STRidge (d_tol=", d_tols,")"))
count_vector_snr$method <- NA
for(i in seq_along(method_label)){
  count_vector_snr$method[which(count_vector_snr$method0==unique(count_vector_snr$method0)[i])] <- method_label[i]
}

rownames(count_vector_snr) <- seq(1, nrow(count_vector_snr))
noise_gather_df <- count_vector_snr %>%
  gather("Condition", "Value",
         3:ncol(count_vector_snr)-1)

x_limits <- snr_db_seq
x_breaks <- pretty_breaks()(snr_db_seq)
# colours <- c("#a670b0",
#              "#a58d48")
# ggcol <- c("#F8766D","#00BA38","#619CFF")
# rudy_color <- colorRampPalette(c("#4bf5db","#6261ff"))(length(d_tols))
# colors <- c(ggcol[1],rudy_color)

NS_SNR_plot <- ggplot() + 
  geom_hline(yintercept = 0.8, lty=2)+
  geom_rect(aes(xmin=rect1$xmin,xmax=rect1$xmax,ymin=rect1$ymin,ymax=rect1$ymax),alpha=0.2,fill="#9de0e6")+
  scale_y_continuous(breaks = seq(0,1,0.2))+
  geom_line(data=count_vector_snr, aes(x=SNR_dB, y=Correct,color=method),lwd=lwd_size) +
  geom_point(data=count_vector_snr, alpha=alpha_value,stroke=stroke_value,aes(x=SNR_dB, y=Correct,shape=method,size=method,fill = method)) + 
  ylab("Success rate") +  xlab("SNR(dB)")+ 
  scale_x_continuous(breaks = snr_db_seq2, labels = snr_db_seq3)+
  scale_size_manual(values=sizes, breaks=unique(count_vector_snr$method), 
                      labels = legend_labs)+
  scale_fill_manual(values=colors, breaks=unique(count_vector_snr$method), 
                      labels = legend_labs)+
  scale_shape_manual(values=shapes, breaks=unique(count_vector_snr$method),
                     labels = legend_labs)+
  scale_color_manual(values=colors, breaks=unique(count_vector_snr$method),
                     labels = legend_labs)+
  labs(fill = 'Methods', shape='Methods', size='Methods',color='Methods')+
  plot_theme + legend_guides +
  guides(x = guide_axis_truncated(
    trunc_lower = c(-Inf, xend + extra_x/2),
    trunc_upper = c(xstart + extra_x/2, Inf)
  ), fill=guide_legend(nrow=2,byrow=TRUE)) +
  annotate("segment", 
           x = myseg$x, xend = myseg$xend,
           y = myseg$y + 0.05, yend = myseg$yend) +
  coord_cartesian(clip = "off", ylim = c(-0.0005, NA))#+
  # guides(colour = guide_legend(override.aes = list(size=5)))

## noiseless data plot ---------------------------
load(paste0(data_path,"NS_SG_noiseless_density_seed_10_ada_lasso_pareto_AIC.RData"))
load(paste0(data_path,"NS_SG_noiseless_density_seed_10_rudy_d_thred.RData"))

# num <- round(10^(seq(log10(200), log10(3000*50), length=12)))
num <- (seq(2, 5.0, 0.2))
dens_fun_NS <- function(noise_coeff, true_terms){
  dens_noise <- sapply(seq_along(noise_coeff), function(i){
    match_index <- match(names(noise_coeff[[i]]), true_terms)
    if(length(match_index) == 4 & all(!is.na(match_index))){
      return(1)
    }else{
      return(0)
    }
  })
  return(dens_noise)
}
true_term_NS = c("w_{xx}", "w_{yy}", "uw_{x}", "vw_{y}")
success_rate_our <- sapply(seq_along(num), function(i) dens_fun_NS(NS_outs_ada_lasso_pareto_success_rate[[i]], true_term_NS))
success_rate_rudy_list <- lapply(NS_outs_rudy_list, function(x) sapply(seq_along(num), function(i) dens_fun_NS(x[[i]], true_term_NS)))
d_tols <- c(0.2,2,10)
names(success_rate_rudy_list) <- paste0("rudy_", d_tols)
dens_all_noiseless <- c(list(our = success_rate_our), success_rate_rudy_list)

count_vector_n <- data.frame(n = 0, Correct = 0, Incorrect = 0, method0 = '0')
for(m_names in names(dens_all_noiseless)){
  temp <- data.frame(n = num,
                     Correct = colSums(dens_all_noiseless[[m_names]] == 1) / 100,
                     Incorrect = colSums(dens_all_noiseless[[m_names]] == 0) / 100,
                     method0 = m_names)
  
  count_vector_n <- rbind.data.frame(count_vector_n, temp)
}
count_vector_n <- count_vector_n[-1, ]
method_label <- c("RAL", paste0("STRidge (d_tol=", d_tols,")"))
count_vector_n$method <- NA
for(i in seq_along(method_label)){
  count_vector_n$method[which(count_vector_n$method0==unique(count_vector_n$method0)[i])] <- method_label[i]
}

rownames(count_vector_n) <- seq(1, nrow(count_vector_n))
noise_gather_df_all <- count_vector_n %>%
  gather("Condition", "Value",
         3:ncol(count_vector_n)-1)

# colours <- c("#a670b0",
#              "#a58d48")
# ggcol <- c("#F8766D","#00BA38","#619CFF")
# rudy_color <- colorRampPalette(c("#4bf5db","#6261ff"))(length(d_tols))
# colors <- c(ggcol[1],rudy_color)

xlables <- unname(latex2exp::TeX(paste0('$10^{',round((num),1),'}$')))
# xlables[c(2,3,5,6,8,9,11,12,14,15)] <- ''
xlables[-(3*(0:(length(num)/2))+1)] <- ''
second_xlables <- round(10^num/(499*199*151)*100,3)
second_xlables[-(3*(0:(length(num)/2))+1)] <- ''

NS_N_plot <- ggplot() + 
  geom_hline(yintercept = 0.8, lty=2)+
  geom_rect(aes(xmin=rect1$xmin,xmax=rect1$xmax,ymin=rect1$ymin,ymax=rect1$ymax),alpha=0.2,fill="#9de0e6")+
  scale_y_continuous(breaks = seq(0,1,0.2))+
  geom_line(data=count_vector_n, aes(x=n, y=Correct,color=method),lwd=lwd_size) +
  geom_point(data=count_vector_n, alpha=alpha_value,stroke=stroke_value,aes(x=n, y=Correct,shape=method,size=method,fill = method)) + 
  scale_x_continuous(breaks = num, labels = xlables,
                     sec.axis = dup_axis(name='Percentage (%)', labels=second_xlables))+ 
  scale_size_manual(values=sizes, breaks=unique(count_vector_n$method), 
                    labels = legend_labs)+
  scale_fill_manual(values=colors, breaks=unique(count_vector_n$method), 
                    labels = legend_labs)+
  scale_shape_manual(values=shapes, breaks=unique(count_vector_n$method), 
                     labels = legend_labs)+
  scale_color_manual(values=colors, breaks=unique(count_vector_n$method),
                     labels = legend_labs)+
  # labs(x="N", y="Success rate")+
  labs(x="N", y=" ")+
  labs(colour = 'Methods', shape='Methods',color='Methods')+
  plot_theme + legend_guides
## compose plots ------------------
legend_n <- get_legend(NS_SNR_plot+theme(legend.position='bottom', ))
                                         # legend.text = element_text(size=16),
                                         # legend.title = element_text(size=18)))
NS_plot <- arrangeGrob(NS_SNR_plot+guides(color='none')+theme(legend.position = 'none'), 
                     NS_N_plot+guides(color='none')+theme(legend.position = 'none'),
                     nrow=1,ncol=2)
NS_plot1 <- arrangeGrob(NS_plot,legend_n,nrow=2, heights = ggplot_heights2)
NS_title <- titles_plot('Navier-Stokes','C',expression(omega[t] == 0.01*omega[xx]+0.01*omega[yy]-u*omega[x]-v*omega[y]),eq_bd_w=0.6)
NS_plot2 <- arrangeGrob(NS_title, NS_sol_plot, NS_plot1, nrow=3, heights =c(2,7,10))
# ggsave(NS_plot2, filename = 'Figures/new_figs/NS.png', width = 8, height = 8)
# ggsave(NS_plot2, filename = 'Figures/new_figs/NS.pdf', width = 8, height = 8)

## ----cable plot, fig.height=4, fig.width=8---------------------------------------------
## cable plots --------------------------------
cable.x <- cable_list$x; cable.t <- cable_list$t; cable.usol <- as.data.frame(cable_list$sol)
gg_cable_data <- cbind.data.frame(x = rep(cable.x,length(cable.t)), t = rep(cable.t,each=length(cable.x)), u = gather(cable.usol)[,-1])
color_sacle <- ramp.col(col = c("#fafafa", "#d7191c"),n=200)
cable_sol_plot <- ggplot(gg_cable_data, aes(x, t, fill= u)) + geom_tile()+
  scale_fill_gradientn(colors = color_sacle) +
  solution_theme+
  coord_cartesian(clip="off")

## cable noisy data plot ---------------------------
load(paste0(data_path, "more_test2/cable_SG_noise_seed_100_samp_100_snr_ada_lasso_pareto.RData"))
load(paste0(data_path, "more_test2/cable_SG_noise_seed_100_samp_100_snr_rudy_d_thred.RData"))
d_tols <- c(0.2,2,10)
names(cable_density_noise_rudy) <- paste0("rudy_", d_tols)
names(cable_density_noise_ada_lasso_pareto) <- 'our'
dens_snr_all <- c(cable_density_noise_ada_lasso_pareto, cable_density_noise_rudy)

snr_db_seq <- c(seq(10, 40, 2), Inf) # SNR_dB grid
eta <- 10^(-snr_db_seq/20)
snr_db_seq1 <- c(seq(10, 40, 2), 45)
snr_db_seq2 <- c(seq(10, 40, 2),  45) # c(seq(20, 40, 2), 42.5, 45)
snr_db_seq3 <- c(as.character(seq(10, 40, 2)), expression(infinity)) # c(as.character(seq(20, 40, 2)), "...", "Inf")
# snr_db_seq3[c(1,3,5,7,9,11,13,15,17,19)+1] <- ''
snr_db_seq3[-(3*(0:(length(snr_db_seq)/2))+1)] <- ''
snr_db_seq3[length(snr_db_seq3)] <- expression(infinity)

cable_count_vector_snr <- data.frame(SNR_dB = 0, Correct = 0, Incorrect = 0, method0 = '0')
for(m_names in names(dens_snr_all)){
  temp <- data.frame(SNR_dB = snr_db_seq1,
                     Correct = colSums(dens_snr_all[[m_names]] == 1)[-c(1:5)] / 100,
                     Incorrect = colSums(dens_snr_all[[m_names]] == 0)[-c(1:5)] / 100,
                     method0 = m_names)
  
  cable_count_vector_snr <- rbind.data.frame(cable_count_vector_snr, temp)
}
cable_count_vector_snr <- cable_count_vector_snr[-1, ]
rownames(cable_count_vector_snr) <- seq(1, nrow(cable_count_vector_snr))
cable_noise_gather_df <- cable_count_vector_snr %>%
  gather("Condition", "Value",
         3:ncol(cable_count_vector_snr)-1)

method_label <- c("RAL", paste0("STRidge (d_tol=", d_tols,")"))
cable_count_vector_snr$method <- NA
for(i in seq_along(method_label)){
  cable_count_vector_snr$method[which(cable_count_vector_snr$method0==unique(cable_count_vector_snr$method0)[i])] <- method_label[i]
}

x_limits <- snr_db_seq
x_breaks <- pretty_breaks()(snr_db_seq)
# colours <- c("#a670b0",
#              "#a58d48")
# ggcol <- c("#F8766D","#00BA38","#619CFF")
# rudy_color <- colorRampPalette(c("#4bf5db","#6261ff"))(length(d_tols))
# colors <- c(ggcol[1],rudy_color)

cable_SNR_plot <- ggplot() + 
  geom_hline(yintercept = 0.8, lty=2)+
  geom_rect(aes(xmin=rect1$xmin,xmax=rect1$xmax,ymin=rect1$ymin,ymax=rect1$ymax),alpha=0.2,fill="#9de0e6")+
  scale_y_continuous(breaks = seq(0,1,0.2))+
  geom_line(data=cable_count_vector_snr, aes(x=SNR_dB, y=Correct,color=method),lwd=lwd_size) +
  geom_point(data=cable_count_vector_snr, alpha=alpha_value,stroke=stroke_value,aes(x=SNR_dB, y=Correct,shape=method,size=method,fill = method)) + 
  # ylab("Probability of correct identification") +  xlab("SNR(dB)")+
  ylab("Success rate") +  xlab("SNR(dB)")+
  # xlab(TeX("Noise magnitude ($\\eta$*sd(u))")) +
  scale_x_continuous(breaks = snr_db_seq2, labels = snr_db_seq3)+
  scale_size_manual(values=sizes, breaks=unique(cable_count_vector_snr$method), 
                    labels = legend_labs)+
  scale_fill_manual(values=colors, breaks=unique(cable_count_vector_snr$method), 
                    labels = legend_labs)+
  scale_shape_manual(values=shapes, breaks=unique(cable_count_vector_snr$method), 
                     labels = legend_labs)+
  scale_color_manual(values=colors, breaks=unique(cable_count_vector_snr$method),
                     labels = legend_labs)+
  labs(fill = 'Methods', shape='Methods', size='Methods',color='Methods')+
  # scale_y_continuous(labels = y_labels,breaks = y_breaks,limits = c(NA, max(y_breaks))) +
  plot_theme + legend_guides +
  guides(x = guide_axis_truncated(
    trunc_lower = c(-Inf, xend + extra_x/2),
    trunc_upper = c(xstart + extra_x/2, Inf)
  )) +
  annotate("segment", 
           x = myseg$x, xend = myseg$xend,
           y = myseg$y + 0.05, yend = myseg$yend)  +
  coord_cartesian(clip = "off", ylim = c(-0.0005, NA))

## cable noiseless plot ----------------------
load(paste0(data_path,"more_test2/cable_SG_noiseless_seed_10_success_rate_ada_lasso_pareto.RData"))
load(paste0(data_path,"more_test2/cable_SG_noiseless_seed_10_success_rate_rudy_d_thred.RData"))
num <- (seq(2,4.4,0.2))

d_tols <- c(0.2,2,10)
names(success_rate_rudy_list) <- paste0("rudy_", d_tols)
success_rate_our <- cable_density_noiseless_ada_lasso
dens_n_all <- c(list(our=success_rate_our), lapply(success_rate_rudy_list, function(x) x))

cable_count_vector_n <- data.frame(n = 0, Correct = 0, Incorrect = 0, method0 = '0')
for(m_names in names(dens_n_all)){
  temp <- data.frame(n = num,
                     Correct = colSums(dens_n_all[[m_names]] == 1) / 100,
                     Incorrect = colSums(dens_n_all[[m_names]] == 0) / 100,
                     method0 = m_names)
  
  cable_count_vector_n <- rbind.data.frame(cable_count_vector_n, temp)
}
cable_count_vector_n <- cable_count_vector_n[-1, ]

rownames(cable_count_vector_n) <- seq(1, nrow(cable_count_vector_n))
noise_gather_df_all <- cable_count_vector_n %>%
  gather("Condition", "Value",
         3:ncol(cable_count_vector_n)-1)

method_label <- c("RAL", paste0("STRidge (d_tol=", d_tols,")"))
cable_count_vector_n$method <- NA
for(i in seq_along(method_label)){
  cable_count_vector_n$method[which(cable_count_vector_n$method0==unique(cable_count_vector_n$method0)[i])] <- method_label[i]
}

# colours <- c("#a670b0",
#              "#a58d48")
# ggcol <- c("#F8766D","#00BA38","#619CFF")
# rudy_color <- colorRampPalette(c("#4bf5db","#6261ff"))(length(d_tols))
# colors <- c(ggcol[1],rudy_color)

xlables <- unname(latex2exp::TeX(paste0('$10^{',round((num),1),'}$')))
# xlables[c(1,3,5,7,9)+1] <- ''
xlables[-(3*(0:(length(num)/2))+1)] <- ''
second_xlables <- round(10^num/(81*501)*100,2)
second_xlables[-(3*(0:(length(num)/2))+1)] <- ''

cable_N_plot <- ggplot() + 
  geom_hline(yintercept = 0.8, lty=2)+
  geom_rect(aes(xmin=rect1$xmin,xmax=rect1$xmax,ymin=rect1$ymin,ymax=rect1$ymax),alpha=0.2,fill="#9de0e6")+
  scale_y_continuous(breaks = seq(0,1,0.2))+
  geom_line(data=cable_count_vector_n, aes(x=n, y=Correct,color=method),lwd=lwd_size) +
  geom_point(data=cable_count_vector_n, alpha=alpha_value,stroke=stroke_value,aes(x=n, y=Correct,shape=method,size=method,fill = method)) + 
  scale_x_continuous(breaks = num, labels = xlables,
                     sec.axis = dup_axis(name='Percentage (%)', labels=second_xlables))+ 
  # labs(x="SNR(dB)", y="Success rate", title = TeX('Reaction-Diffusion, $\\alpha$=0.3, Inf $\\alpha$=0.5'))+
  # labs(x="N", y="Success rate")+
  labs(x="N", y=" ")+
  scale_size_manual(values=sizes, breaks=unique(cable_count_vector_n$method), 
                    labels = legend_labs)+
  scale_fill_manual(values=colors, breaks=unique(cable_count_vector_n$method), 
                    labels = legend_labs)+
  scale_shape_manual(values=shapes, breaks=unique(cable_count_vector_n$method), 
                     labels = legend_labs)+
  scale_color_manual(values=colors, breaks=unique(cable_count_vector_n$method),
                     labels = legend_labs)+
  labs(fill = 'Methods', shape='Methods', size='Methods',color='Methods')+
  # scale_y_continuous(labels = y_labels,breaks = y_breaks,limits = c(NA, max(y_breaks))) +
  plot_theme + legend_guides 
# theme_classic(18)
## out cable -----------------------
legend_n <- get_legend(cable_N_plot+theme(legend.position='bottom'))
cable_plot <- arrangeGrob(cable_SNR_plot+guides(color='none')+theme(legend.position = 'none'), 
                        cable_N_plot+guides(color='none')+theme(legend.position = 'none'),
                        nrow=1,ncol=2)
cable_plot1 <- arrangeGrob(cable_plot,legend_n,nrow=2, heights = ggplot_heights2)
cable_title <- titles_plot('Cable','B',expression(u[t] == u[xx] - u))
cable_plot2 <- arrangeGrob(cable_title, cable_sol_plot, cable_plot1, nrow=3,heights =c(2,7,10))

# ggsave(cable_plot2, filename = 'Figures/new_figs/cable.png', width = 8, height = 8)
# ggsave(cable_plot2, filename = 'Figures/new_figs/cable.pdf', width = 8, height = 8)

## ----RD plot, fig.height=4, fig.width=8------------------------------------------------
# RD plots -------------------------------
RD.mat <- readMat('data/reaction_diffusion_t_1.mat')
RD.t <- as.numeric(RD.mat[['t']])
RD.x <- as.numeric(RD.mat[['x']])
RD.y <- as.numeric(RD.mat[['y']])
RD.u <- as.data.frame(Re(RD.mat[['u.t.1']]))
RD.v <- as.data.frame(Re(RD.mat[['v.t.1']]))

gg_RD_data <- cbind.data.frame(
  x = rep(RD.x,length(RD.y)), y = rep(RD.y,each=length(RD.x)), 
  u = gather(RD.u)[,-1], v = gather(RD.v)[,-1]) %>% 
  reshape2::melt(id=1:2, variable.name='sol')

RD_sol_plot <- ggplot(gg_RD_data, aes(x, y, fill= value))+
  geom_raster()+
  facet_wrap(~sol,ncol=2)+
  scale_fill_distiller(palette = "Spectral") + 
  # labs(title='Reaction-diffusion', tag='C', 
  #      subtitle = expression(atop(u[t] == 0.1*u[xx] + 0.1*u[yy] + u - u*v^{2} - u^{3} + v^{3} + u^{2}*v,
  #                                 v[t] == 0.1*v[xx] + 0.1*v[yy] + v - u*v^{2} - u^{3} - v^{3} - u^{2}*v))) +
  # labs(title='Reaction-diffusion', tag='C', subtitle = expression('0.1*[xx]'))+
  solution_theme + 
  theme(strip.text = element_text(size = 16),
        # strip.background = element_rect(fill = 'grey90',color='grey90')
        ) +
  # annotate('label',x=c(5,-5),y=22,label='        ',
  #          fill = 'grey90',label.size=NA,size=20)+
  coord_cartesian(clip="off")
## noise data plot -------------------------
load(paste0(data_path,"../RD/RD_SG_noise_density_seed_10_snr_ada_lasso_pareto_AIC.RData"))
load(paste0(data_path,"RD_SG_noise_density_seed_10_snr_rudy_d_thred.RData"))

dens_fun_RD <- function(noise_coeff, true_terms, which){
  dens_noise <- sapply(seq_along(noise_coeff), function(i){
    match_index_ut <- match(names(noise_coeff[[i]][[1]]), true_terms[[1]])
    match_index_vt <- match(names(noise_coeff[[i]][[2]]), true_terms[[2]])
    match_TF_all <- length(match_index_ut) == 7 & all(!is.na(match_index_ut)) & 
      length(match_index_vt) == 7 & all(!is.na(match_index_vt))
    match_TF_ut <- length(match_index_ut) == 7 & all(!is.na(match_index_ut))
    match_TF_vt <- length(match_index_vt) == 7 & all(!is.na(match_index_vt))
    if(which=="ut"){
      if(match_TF_ut){
        return(1)
      }else{
        return(0)
      }
    }else if(which=="vt"){
      if(match_TF_vt){
        return(1)
      }else{
        return(0)
      }
    }else{
      if(match_TF_all){
        return(1)
      }else{
        return(0)
      }
    }
  })
  return(dens_noise)
}
true_term_RD <- list(u_t = c("u_{xx}", "u_{yy}", "u", "v^3", "uv^2", "u^2v", "u^3"),
                     v_t = c("v_{xx}", "v_{yy}", "v", "v^3", "uv^2", "u^2v", "u^3"))

RD_noise_ada_lasso_pareto2 <- lapply(RD_noise_ada_lasso_pareto, function(x){
  lapply(x, function(x1){
    lapply(x1, function(i){
      i[which(abs(i) >= 0.08)]
    })
  })
})
dens_our <- sapply(1:12, function(i) dens_fun_RD(RD_noise_ada_lasso_pareto[[i]], true_term_RD, "all"))
dens_our2 <- sapply(1:12, function(i) dens_fun_RD(RD_noise_ada_lasso_pareto2[[i]], true_term_RD, "all"))
dens_rudy_list <- lapply(RD_noise_rudy_list, function(k) sapply(1:12, function(i) dens_fun_RD(k[[i]], true_term_RD, "all"))) 
d_tols <- c(0.2,2,10)
names(dens_rudy_list) <- paste0("rudy_", d_tols)
dens_all_noise <- c(list(our = dens_our), dens_rudy_list)

noise = snr_db_seq <- c(seq(20, 40, 2), Inf) # SNR_dB grid
eta <- 10^(-snr_db_seq/20)/2
snr_db_seq1 <- c(seq(20, 40, 2), 45)
snr_db_seq2 <- c(seq(20, 40, 2),  45) # c(seq(20, 40, 2), 42.5, 45)
snr_db_seq3 <- c(as.character(seq(20, 40, 2)), expression(infinity)) # c(as.character(seq(20, 40, 2)), "...", "Inf")
snr_db_seq3[c(1,3,5,7,9)+1] <- ''

count_vector_snr <- data.frame(SNR_dB = 0, Correct = 0, Incorrect = 0, method0 = '0')
for(m_names in names(dens_all_noise)){
  temp <- data.frame(SNR_dB = snr_db_seq1,
                     Correct = colSums(dens_all_noise[[m_names]] == 1) / 100,
                     Incorrect = colSums(dens_all_noise[[m_names]] == 0) / 100,
                     method0 = m_names)
  
  count_vector_snr <- rbind.data.frame(count_vector_snr, temp)
}
count_vector_snr <- count_vector_snr[-1, ]

rownames(count_vector_snr) <- seq(1, nrow(count_vector_snr))
noise_gather_df_all <- count_vector_snr %>%
  gather("Condition", "Value",
         3:ncol(count_vector_snr)-1)
method_label <- c("RAL", paste0("STRidge (d_tol=", d_tols,")"))
count_vector_snr$method <- NA
for(i in seq_along(method_label)){
  count_vector_snr$method[which(count_vector_snr$method0==unique(count_vector_snr$method0)[i])] <- method_label[i]
}

x_limits <- snr_db_seq1
x_breaks <- pretty_breaks()(snr_db_seq)
# colours <- c("#a670b0",
#              "#a58d48")
# ggcol <- c("#F8766D","#00BA38","#619CFF")
# rudy_color <- colorRampPalette(c("#4bf5db","#6261ff"))(length(d_tols))
# colors <- c(ggcol[1],rudy_color)

RD_SNR_plot <- ggplot() + 
  geom_hline(yintercept = 0.8, lty=2)+
  geom_rect(aes(xmin=rect1$xmin,xmax=rect1$xmax,ymin=rect1$ymin,ymax=rect1$ymax),alpha=0.2,fill="#9de0e6")+
  scale_y_continuous(breaks = seq(0,1,0.2))+
  geom_line(data=count_vector_snr, aes(x=SNR_dB, y=Correct,color=method),lwd=lwd_size) +
  geom_point(data=count_vector_snr, alpha=alpha_value,stroke=stroke_value,aes(x=SNR_dB, y=Correct,shape=method,size=method,fill = method)) + 
  labs(x="SNR(dB)", y="Success rate")+
  scale_x_continuous(breaks = snr_db_seq1, labels = snr_db_seq3)+
  scale_size_manual(values=sizes, breaks=unique(count_vector_snr$method), 
                    labels = legend_labs)+
  scale_fill_manual(values=colors, breaks=unique(count_vector_snr$method), 
                    labels = legend_labs)+
  scale_shape_manual(values=shapes, breaks=unique(count_vector_snr$method), 
                     labels = legend_labs)+
  scale_color_manual(values=colors, breaks=unique(count_vector_snr$method),
                     labels = legend_labs)+
  labs(fill = 'Methods', shape='Methods', size='Methods',color='Methods')+
  plot_theme + legend_guides +
  guides(x = guide_axis_truncated(
    trunc_lower = c(-Inf, xend + extra_x/2),
    trunc_upper = c(xstart + extra_x/2, Inf)
  )) +
  annotate("segment", 
           x = myseg$x, xend = myseg$xend,
           y = myseg$y + 0.05, yend = myseg$yend)  +
  coord_cartesian(clip = "off", ylim = c(-0.0005, NA))

## noiseless data ----------------------------
load(paste0(data_path,"RD_SG_noiseless_seed_10_ada_lasso_pareto.RData"))
# load("RD/RD_SG_noiseless_200_150000_seed10_MSA_lasso_pareto_AIC.RData")
load(paste0(data_path,"RD_SG_noiseless_seed_10_rudy_d_thred.RData"))
# num <- round(10^(seq(2.2, 5, 0.2)))
num <- seq(2.2, 5, 0.2)

RD_outs_ada_lasso_pareto_success_rate2 <- lapply(RD_outs_ada_lasso_pareto_success_rate, function(x){
  lapply(x, function(x1){
    lapply(x1, function(i){
      i[which(abs(i) >= 0.08)]
    })
  })
})
success_rate_our <- sapply(seq_along(num), function(i) dens_fun_RD(RD_outs_ada_lasso_pareto_success_rate[[i]], true_term_RD, "all"))
success_rate_our2 <- sapply(seq_along(num), function(i) dens_fun_RD(RD_outs_ada_lasso_pareto_success_rate2[[i]], true_term_RD, "all"))
success_rate_rudy_list <- lapply(RD_outs_rudy_list, function(k) sapply(seq_along(num), function(i) dens_fun_RD(k[[i]], true_term_RD, "all")))
d_tols <- c(0.2,2,10)
names(success_rate_rudy_list) <- paste0("rudy_", d_tols)
dens_all_noiseless <- c(list(our = success_rate_our), success_rate_rudy_list)

count_vector_n <- data.frame(n = 0, Correct = 0, Incorrect = 0, method0 = '0')
for(m_names in names(dens_all_noiseless)){
  temp <- data.frame(n = num,
                     Correct = colSums(dens_all_noiseless[[m_names]] == 1) / 100,
                     Incorrect = colSums(dens_all_noiseless[[m_names]] == 0) / 100,
                     method0 = m_names)
  
  count_vector_n <- rbind.data.frame(count_vector_n, temp)
}
count_vector_n <- count_vector_n[-1, ]

rownames(count_vector_n) <- seq(1, nrow(count_vector_n))
noise_gather_df_all <- count_vector_n %>%
  gather("Condition", "Value",
         3:ncol(count_vector_n)-1)
method_label <- c("RAL", paste0("STRidge (d_tol=", d_tols,")"))
count_vector_n$method <- NA
for(i in seq_along(method_label)){
  count_vector_n$method[which(count_vector_n$method0==unique(count_vector_n$method0)[i])] <- method_label[i]
}

x_limits <- num
x_breaks <- pretty_breaks()(num)
# colours <- c("#a670b0",
#              "#a58d48")
# ggcol <- c("#F8766D","#00BA38","#619CFF")
# rudy_color <- colorRampPalette(c("#4bf5db","#6261ff"))(length(d_tols))
# colors <- c(ggcol[1],rudy_color)

xlables <- unname(latex2exp::TeX(paste0('$10^{',round((num),1),'}$')))
# xlables[c(2,4,6,8,10,12,14)] <- ''
# xlables[-c(3*(0:5)+1)] <- ''
xlables[-(3*(0:(length(num)/2))+1)] <- ''
second_xlables <- round(10^num/(512*512*101)*100,3)
second_xlables[-(3*(0:(length(num)/2))+1)] <- ''

RD_N_plot <- ggplot() + 
  geom_hline(yintercept = 0.8, lty=2)+
  geom_rect(aes(xmin=rect1$xmin,xmax=rect1$xmax,ymin=rect1$ymin,ymax=rect1$ymax),alpha=0.2,fill="#9de0e6")+
  scale_y_continuous(breaks = seq(0,1,0.2))+
  geom_line(data=count_vector_n, aes(x=n, y=Correct,color=method),lwd=lwd_size) +
  geom_point(data=count_vector_n,alpha=alpha_value,stroke=stroke_value,aes(x=n, y=Correct,shape=method,size=method,fill = method)) + 
  scale_x_continuous(breaks = num, labels = xlables,
                     sec.axis = dup_axis(name='Percentage (%)', labels=second_xlables))+ 
  # labs(x="N", y="Success rate")+
  labs(x="N", y=" ")+
  scale_size_manual(values=sizes, breaks=unique(count_vector_n$method), 
                    labels = legend_labs)+
  scale_fill_manual(values=colors, breaks=unique(count_vector_n$method), 
                    labels = legend_labs)+
  scale_shape_manual(values=shapes, breaks=unique(count_vector_n$method), 
                     labels = legend_labs)+
  scale_color_manual(values=colors, breaks=unique(count_vector_n$method),
                     labels = legend_labs)+
  labs(fill = 'Methods', shape='Methods', size='Methods',color='Methods')+
  plot_theme + legend_guides
  
# source('./RD_SG_noise_plot.R')
## compose plots ---------------------
# layout_matrix <- matrix(c(1,2,1,3),ncol=2)
legend_n <- get_legend(RD_N_plot+theme(legend.position='bottom'))
RD_plot <- arrangeGrob(RD_SNR_plot+guides(color='none')+theme(legend.position = 'none'), 
                     RD_N_plot+guides(color='none')+theme(legend.position = 'none'),
                     nrow=1,ncol=2)
RD_plot1 <- arrangeGrob(RD_plot,legend_n,nrow=2, heights = ggplot_heights2)
RD_title <- titles_plot('Reaction-diffusion','D',
                        expression(atop(u[t] == 0.1*u[xx] + 0.1*u[yy] + u - u*v^{2} - u^{3} + v^{3} + u^{2}*v,
                                        v[t] == 0.1*v[xx] + 0.1*v[yy] + v - u*v^{2} - u^{3} - v^{3} - u^{2}*v)),
                        eq_bd_h=1,eq_bd_w=0.7,title_size=16)
RD_title1 <- titles_plot('Reaction-diffusion','D',
                        expression(u[t] == 0.1*u[xx] + 0.1*u[yy] + u - u*v^{2} - u^{3} + v^{3} + u^{2}*v ~~~~~~~~~~~
                                   v[t] == 0.1*v[xx] + 0.1*v[yy] + v - u*v^{2} - u^{3} - v^{3} - u^{2}*v),
                        eq_bd_h=1,eq_bd_w=0.6,title_size=16)

RD_plot2 <- arrangeGrob(RD_title, RD_sol_plot, RD_plot1, nrow=3,heights =c(2,7,10))
# ggsave(RD_plot2, filename = 'Figures/new_figs/RD.png', width = 8, height = 8)
# ggsave(RD_plot2, filename = 'Figures/new_figs/RD.pdf', width = 8, height = 8)

## ----kdv plot, fig.height=4, fig.width=8-----------------------------------------------
## kdv plot ---------------------
source("pde_solver_data/kdv_solver.R")
kdv.usol <- as.data.frame(kdv.usol)
gg_kdv_data <- cbind.data.frame(x = rep(kdv.x,length(kdv.t)), t = rep(kdv.t,each=length(kdv.x)), u = gather(kdv.usol)[,-1])
kdv_sol_plot <- ggplot(gg_kdv_data) + 
  geom_raster(aes(x, t, fill = u))+ 
  # labs(title='KdV', tag='D', subtitle = TeX('$u_t=6uu_x-u_{xxx}$')) + 
  # geom_contour(aes(x, t, z = u))+
  scale_fill_distiller(palette = "Spectral") + 
  solution_theme +
  # annotate('label',x=10,y=24.2,label='                                 ',fill = 'grey90',label.size=NA,size=5.7)+
  coord_cartesian(clip="off")

## kdv noisy plot -----------------------
load(paste0(data_path,"kdv_SG_noise_seed_100_samp_100_snr_ada_lasso_pareto.RData"))
load(paste0(data_path,"kdv_SG_noise_seed_100_samp_100_snr_rudy_d_thred.RData"))
d_tols <- c(0.2,2,10)
names(kdv_density_noise_rudy_list) <- paste0("rudy_", d_tols)
dens_snr_all <- c(list(our=kdv_density_noise_ada_lasso_pareto[[1]]), kdv_density_noise_rudy_list)

noise = snr_db_seq <- c(seq(40, 60, 2), Inf) # SNR_dB grid
eta <- 10^(-snr_db_seq/20)
snr_db_seq1 <- c(seq(40, 60, 2), 65)
snr_db_seq2 <- c(seq(40, 60, 2),  65) # c(seq(20, 40, 2), 42.5, 45)
snr_db_seq3 <- c(as.character(seq(40, 60, 2)), TeX("$\\infty$")) # c(as.character(seq(20, 40, 2)), "...", "Inf")
snr_db_seq3[c(1,3,5,7,9)+1] <- ''

kdv_count_vector_snr <- data.frame(SNR_dB = 0, Correct = 0, Incorrect = 0, method0 = '0')
for(m_names in names(dens_snr_all)){
  temp <- data.frame(SNR_dB = snr_db_seq1,
                     Correct = colSums(dens_snr_all[[m_names]] == 1) / 100,
                     Incorrect = colSums(dens_snr_all[[m_names]] == 0) / 100,
                     method0 = m_names)
  
  kdv_count_vector_snr <- rbind.data.frame(kdv_count_vector_snr, temp)
}
kdv_count_vector_snr <- kdv_count_vector_snr[-1, ]

rownames(kdv_count_vector_snr) <- seq(1, nrow(kdv_count_vector_snr))
noise_gather_df <- kdv_count_vector_snr %>%
  gather("Condition", "Value",
         3:ncol(kdv_count_vector_snr)-1)

method_label <- c("RAL", paste0("STRidge (d_tol=", d_tols,")"))
kdv_count_vector_snr$method <- NA
for(i in seq_along(method_label)){
  kdv_count_vector_snr$method[which(kdv_count_vector_snr$method0==unique(kdv_count_vector_snr$method0)[i])] <- method_label[i]
}

x_limits <- snr_db_seq
x_breaks <- pretty_breaks()(snr_db_seq)
# colours <- c("#a670b0",
#              "#a58d48")
# ggcol <- c("#F8766D","#00BA38","#619CFF")
# rudy_color <- colorRampPalette(c("#4bf5db","#6261ff"))(length(d_tols))
# colors <- c(ggcol[1],rudy_color)

kdv_SNR_plot <- ggplot() + 
  geom_hline(yintercept = 0.8, lty=2)+
  geom_rect(aes(xmin=rect1$xmin,xmax=rect1$xmax,ymin=rect1$ymin,ymax=rect1$ymax),alpha=0.2,fill="#9de0e6")+
  scale_y_continuous(breaks = seq(0,1,0.2))+
  geom_line(data=kdv_count_vector_snr, aes(x=SNR_dB, y=Correct,color=method),lwd=lwd_size) +
  geom_point(data=kdv_count_vector_snr, alpha=alpha_value,stroke=stroke_value,aes(x=SNR_dB, y=Correct,shape=method,size=method,fill = method)) + 
  ylab("Success rate") +  xlab("SNR(dB)")+ 
  scale_x_continuous(breaks = snr_db_seq2, labels = snr_db_seq3)+
  scale_size_manual(values=sizes, breaks=unique(kdv_count_vector_snr$method), 
                    labels = legend_labs)+
  scale_fill_manual(values=colors, breaks=unique(kdv_count_vector_snr$method), 
                    labels = legend_labs)+
  scale_shape_manual(values=shapes, breaks=unique(kdv_count_vector_snr$method), 
                     labels = legend_labs)+
  scale_color_manual(values=colors, breaks=unique(kdv_count_vector_snr$method),
                     labels = legend_labs)+
  labs(fill = 'Methods', shape='Methods', size='Methods',color='Methods')+
  guides(x = guide_axis_truncated(
    trunc_lower = c(-Inf, xend+20 + extra_x/2),
    trunc_upper = c(xstart+20 + extra_x/2, Inf)
  )) +
  plot_theme + legend_guides +
  annotate("segment", 
           x = myseg2$x, xend = myseg2$xend,
           y = myseg2$y + 0.05, yend = myseg2$yend)  +
  coord_cartesian(clip = "off", ylim = c(-0.0005, NA))

## kdv noiseless plot -----------------------
load(paste0(data_path,"kdv_SG_noiseless_seed_10_success_rate_pareto_AIC.RData"))
load(paste0(data_path,"kdv_SG_noiseless_seed_10_success_rate_rudy_d_thred.RData"))
# num <- round(10^(seq(2,4.8,0.2)))
num <- seq(2,4.8,0.2)
d_tols <- c(0.2,2,10)
names(success_rate_rudy_list) <- paste0("rudy_", d_tols)
dens_n_all <- c(list(our=kdv_density_noiseless_ada_lasso_pareto[[1]]), success_rate_rudy_list)

kdv_count_vector_n <- data.frame(n = 0, Correct = 0, Incorrect = 0, method0 = '0')
for(m_names in names(dens_n_all)){
  temp <- data.frame(n = num,
                     Correct = colSums(dens_n_all[[m_names]] == 1) / 100,
                     Incorrect = colSums(dens_n_all[[m_names]] == 0) / 100,
                     method0 = m_names)
  
  kdv_count_vector_n <- rbind.data.frame(kdv_count_vector_n, temp)
}
kdv_count_vector_n <- kdv_count_vector_n[-1, ]

rownames(kdv_count_vector_n) <- seq(1, nrow(kdv_count_vector_n))
noise_gather_df_all <- kdv_count_vector_n %>%
  gather("Condition", "Value",
         3:ncol(kdv_count_vector_n)-1)
method_label <- c("RAL", paste0("STRidge (d_tol=", d_tols,")"))
kdv_count_vector_n$method <- NA
for(i in seq_along(method_label)){
  kdv_count_vector_n$method[which(kdv_count_vector_n$method0==unique(kdv_count_vector_n$method0)[i])] <- method_label[i]
}

# colours <- c("#a670b0",
#              "#a58d48")
# ggcol <- c("#F8766D","#00BA38","#619CFF")
# rudy_color <- colorRampPalette(c("#4bf5db","#6261ff"))(length(d_tols))
# colors <- c(ggcol[1],rudy_color)

xlables <- unname(latex2exp::TeX(paste0('$10^{',round((num),1),'}$')))
# xlables[c(2,4,6,8,10,12,14)] <- ''
xlables[-(3*(0:(length(num)/2))+1)] <- ''
second_xlables <- round(10^num/(201*512)*100,2)
second_xlables[-(3*(0:(length(num)/2))+1)] <- ''

kdv_N_plot <- ggplot() + 
  geom_hline(yintercept = 0.8, lty=2)+
  geom_rect(aes(xmin=rect1$xmin,xmax=rect1$xmax,ymin=rect1$ymin,ymax=rect1$ymax),alpha=0.2,fill="#9de0e6")+
  scale_y_continuous(breaks = seq(0,1,0.2))+
  geom_line(data=kdv_count_vector_n, aes(x=n, y=Correct,color=method),lwd=lwd_size) +
  geom_point(data=kdv_count_vector_n,alpha=alpha_value,stroke=stroke_value,aes(x=n, y=Correct,shape=method,size=method,fill = method)) + 
  scale_x_continuous(breaks = num, labels = xlables,
                     sec.axis = dup_axis(name='Percentage (%)', labels=second_xlables))+ 
  # labs(x="N", y="Success rate")+
  labs(x="N", y=" ")+
  scale_size_manual(values=sizes, breaks=unique(kdv_count_vector_n$method), 
                    labels = legend_labs)+
  scale_fill_manual(values=colors, breaks=unique(kdv_count_vector_n$method), 
                    labels = legend_labs)+
  scale_shape_manual(values=shapes, breaks=unique(kdv_count_vector_n$method), 
                     labels = legend_labs)+
  scale_color_manual(values=colors, breaks=unique(kdv_count_vector_n$method),
                     labels = legend_labs)+
  labs(fill = 'Methods', shape='Methods', size='Methods',color='Methods')+
  plot_theme + legend_guides

## compose plots ---------------------
legend_n <- get_legend(kdv_SNR_plot+theme(legend.position='bottom'))
kdv_plot <- arrangeGrob(kdv_SNR_plot+guides(color='none')+theme(legend.position = 'none'), 
                     kdv_N_plot+guides(color='none')+theme(legend.position = 'none'),
                     nrow=1,ncol=2)
kdv_plot1 <- arrangeGrob(kdv_plot,legend_n,nrow=2, heights = ggplot_heights2)
# grid.arrange(kdv_plot2)
kdv_title <- titles_plot('KdV','B',expression(u[t] == 6*uu[x]-u[xxx]))
kdv_plot2 <- arrangeGrob(kdv_title, kdv_sol_plot,kdv_plot1,nrow=3,heights =c(2,7,10))
# ggsave(kdv_plot2, filename = 'Figures/new_figs/kdv.png', width = 8, height = 8)
# ggsave(kdv_plot2, filename = 'Figures/new_figs/kdv.pdf', width = 8, height = 8)

## ----qho plot, fig.height=4, fig.width=8-----------------------------------------------
## qho plots --------------------------------
qho.x <- qho_list$x; qho.t <- qho_list$t; qho.abs_usol <- as.data.frame(abs(t(qho_list$sol)))
gg_qho_data <- cbind.data.frame(x = rep(qho.x,length(qho.t)), t = rep(qho.t,each=length(qho.x)), u = gather(qho.abs_usol)[,-1])
color_sacle <- ramp.col(col = c("#fafafa", "#d7191c"),n=200)
qho_sol_plot <- ggplot(gg_qho_data, aes(x, t, fill= u)) + geom_tile()+
  scale_fill_gradientn(colors = color_sacle) +
  labs(fill='|u|') +
  solution_theme+
  coord_cartesian(clip="off")

# point_shapes_qho <- c(15,rep(16,3),17)
# plot_legends_qho <- c('ARGOS-RAL',TeX('STRidge ($d_{tol}$=0.2)'), TeX('STRidge ($d_{tol}$=2)'), TeX('STRidge ($d_{tol}$=10)'), 'Backward')
# plot_legends_qho <- c('ARGOS-RAL',TeX('STRidge ($d_{tol}$=0.2)'), TeX('STRidge ($d_{tol}$=2)'), TeX('STRidge ($d_{tol}$=10)'))
# colors_qho <- c(ggcol[1],rudy_color,'#40d667')
## qho noisy plot -----------------------
load(paste0(data_path,"qho_SG_noise_seed_10_samp_100_snr_ada_rudy_back.RData"))
d_tols <- c(0.2,2,10)
dens_snr_all <- list(our = qho_density_noise_ada_AIC, rudy_0.2 = qho_density_noise_rudy_0.2,
                     rudy_2 = qho_density_noise_rudy_2, rudy_10 = qho_density_noise_rudy_10)
# back = qho_density_noise_back)

noise = snr_db_seq <- c(seq(40, 60, 2), Inf) # SNR_dB grid
eta <- 10^(-snr_db_seq/20)
snr_db_seq1 <- c(seq(40, 60, 2), 65)
snr_db_seq2 <- c(seq(40, 60, 2),  65) # c(seq(20, 40, 2), 42.5, 45)
snr_db_seq3 <- c(as.character(seq(40, 60, 2)), TeX("$\\infty$")) # c(as.character(seq(20, 40, 2)), "...", "Inf")
snr_db_seq3[c(1,3,5,7,9)+1] <- ''

qho_count_vector_snr <- data.frame(SNR_dB = 0, Correct = 0, Incorrect = 0, method0 = '0')
for(m_names in names(dens_snr_all)){
  temp <- data.frame(SNR_dB = snr_db_seq1,
                     Correct = colSums(dens_snr_all[[m_names]] == 1) / 100,
                     Incorrect = colSums(dens_snr_all[[m_names]] == 0) / 100,
                     method0 = m_names)
  
  qho_count_vector_snr <- rbind.data.frame(qho_count_vector_snr, temp)
}
qho_count_vector_snr <- qho_count_vector_snr[-1, ]

rownames(qho_count_vector_snr) <- seq(1, nrow(qho_count_vector_snr))
noise_gather_df <- qho_count_vector_snr %>%
  gather("Condition", "Value",
         3:ncol(qho_count_vector_snr)-1)

# method_label <- c("RAL", paste0("STRidge (d_tol=", d_tols,")"), 'Backward')
method_label <- c("RAL", paste0("STRidge (d_tol=", d_tols,")"))
qho_count_vector_snr$method <- NA
for(i in seq_along(method_label)){
  qho_count_vector_snr$method[which(qho_count_vector_snr$method0==unique(qho_count_vector_snr$method0)[i])] <- method_label[i]
}

x_limits <- snr_db_seq
x_breaks <- pretty_breaks()(snr_db_seq)
# colours <- c("#a670b0",
#              "#a58d48")
# ggcol <- c("#F8766D","#00BA38","#619CFF")
# rudy_color <- colorRampPalette(c("#4bf5db","#6261ff"))(length(d_tols))
# colors <- c(ggcol[1],rudy_color)

qho_SNR_plot <- ggplot() + 
  geom_hline(yintercept = 0.8, lty=2)+
  geom_rect(aes(xmin=rect1$xmin,xmax=rect1$xmax,ymin=rect1$ymin,ymax=rect1$ymax),alpha=0.2,fill="#9de0e6")+
  scale_y_continuous(breaks = seq(0,1,0.2))+
  geom_line(data=qho_count_vector_snr, aes(x=SNR_dB, y=Correct,color=method),lwd=lwd_size) +
  geom_point(data=qho_count_vector_snr, alpha=alpha_value,stroke=stroke_value,aes(x=SNR_dB, y=Correct,shape=method,size=method,fill = method)) + 
  ylab("Success rate") +  xlab("SNR(dB)")+ 
  scale_x_continuous(breaks= snr_db_seq2, labels = snr_db_seq3)+
  scale_size_manual(values=sizes, breaks=unique(qho_count_vector_snr$method), 
                    labels = legend_labs)+
  scale_fill_manual(values=colors, breaks=unique(qho_count_vector_snr$method), 
                    labels = legend_labs)+
  scale_shape_manual(values=shapes, breaks=unique(qho_count_vector_snr$method), 
                     labels = legend_labs)+
  scale_color_manual(values=colors, breaks=unique(qho_count_vector_snr$method),
                     labels = legend_labs)+
  labs(fill = 'Methods', shape='Methods', size='Methods',color='Methods')+
  # labs(colour = 'Methods', shape='Methods', size='Methods',title = 'QHO noisy data', tag = 'a')+
  guides(x = guide_axis_truncated(
    trunc_lower = c(-Inf, xend + 20 + extra_x/2),
    trunc_upper = c(xstart +20 + extra_x/2, Inf)
  )) +
  plot_theme + legend_guides +
  annotate("segment", 
           x = myseg2$x, xend = myseg2$xend,
           y = myseg2$y + 0.05, yend = myseg2$yend)+
  # annotate("segment", 
  #          x = myseg2$x, xend = myseg2$xend,
  #          y = myseg2$y + 0.07, yend = myseg2$yend-0.021)+
  coord_cartesian(clip = "off", ylim = c(-0.0005, NA)) #+
  # guides(colour = guide_legend(override.aes = list(size=5)))

## qho noiseless plot -----------------------
load(paste0(data_path,"qho_SG_noiseless_seed_10_n_ada_rudy_back.RData"))
# num <- round(10^(seq(2.2,5.2,0.2)))
num <- seq(2.2,5.2,0.2)
d_tols <- c(0.2,2,10)
dens_n_all <- list(our = qho_density_noiseless_ada_AIC, rudy_0.2 = qho_density_noiseless_rudy_0.2, 
                   rudy_2 = qho_density_noiseless_rudy_2, rudy_10 = qho_density_noiseless_rudy_10)
# back = qho_density_noiseless_back)

qho_count_vector_n <- data.frame(n = 0, Correct = 0, Incorrect = 0, method0 = '0')
for(m_names in names(dens_n_all)){
  temp <- data.frame(n = num,
                     Correct = colSums(dens_n_all[[m_names]] == 1) / 100,
                     Incorrect = colSums(dens_n_all[[m_names]] == 0) / 100,
                     method0 = m_names)
  
  qho_count_vector_n <- rbind.data.frame(qho_count_vector_n, temp)
}
qho_count_vector_n <- qho_count_vector_n[-1, ]

rownames(qho_count_vector_n) <- seq(1, nrow(qho_count_vector_n))
noise_gather_df_all <- qho_count_vector_n %>%
  gather("Condition", "Value",
         3:ncol(qho_count_vector_n)-1)
# method_label <- c("RAL", paste0("STRidge (d_tol=", d_tols,")"),'Backward')
method_label <- c("RAL", paste0("STRidge (d_tol=", d_tols,")"))
qho_count_vector_n$method <- NA
for(i in seq_along(method_label)){
  qho_count_vector_n$method[which(qho_count_vector_n$method0==unique(qho_count_vector_n$method0)[i])] <- method_label[i]
}

# colours <- c("#a670b0",
#              "#a58d48")
# ggcol <- c("#F8766D","#00BA38","#619CFF")
# rudy_color <- colorRampPalette(c("#4bf5db","#6261ff"))(length(d_tols))
# colors <- c(ggcol[1],rudy_color)

xlables <- unname(latex2exp::TeX(paste0('$10^{',round((num),1),'}$')))
# xlables[seq(2,16,2)] <- ''
# xlables[-c(3*(0:5)+1)] <- ''
xlables[-(3*(0:(length(num)/2))+1)] <- ''
second_xlables <- round(10^num/(401*513)*100,2)
second_xlables[-(3*(0:(length(num)/2))+1)] <- ''

qho_N_plot <- ggplot() + 
  geom_hline(yintercept = 0.8, lty=2)+
  geom_rect(aes(xmin=rect1$xmin,xmax=rect1$xmax,ymin=rect1$ymin,ymax=rect1$ymax),alpha=0.2,fill="#9de0e6")+
  scale_y_continuous(breaks = seq(0,1,0.2))+
  geom_line(data=qho_count_vector_n, aes(x=n, y=Correct,color=method),lwd=lwd_size) +
  geom_point(data=qho_count_vector_n, alpha=alpha_value,stroke=stroke_value,aes(x=n, y=Correct,shape=method,size=method,fill = method)) + 
  scale_x_continuous(breaks = num, labels = xlables,
                     sec.axis = dup_axis(name='Percentage (%)', labels=second_xlables))+ 
  # labs(x="N", y="Success rate", title = 'QHO noiseless data', tag = 'b')+
  # labs(x="N", y="Success rate")+
  labs(x="N", y=" ")+
  scale_size_manual(values=sizes, breaks=unique(qho_count_vector_n$method), 
                    labels = legend_labs)+
  scale_fill_manual(values=colors, breaks=unique(qho_count_vector_n$method), 
                    labels = legend_labs)+
  scale_shape_manual(values=shapes, breaks=unique(qho_count_vector_n$method), 
                     labels = legend_labs)+
  scale_color_manual(values=colors, breaks=unique(qho_count_vector_n$method),
                     labels = legend_labs)+
  labs(fill = 'Methods', shape='Methods', size='Methods',color='Methods')+
  plot_theme + legend_guides
  # guides(colour = guide_legend(override.aes = list(size=5)))

## compose plots --------------------------
legend_n <- get_legend(qho_N_plot+theme(legend.position='bottom'))
qho_plot <- arrangeGrob(qho_SNR_plot+guides(color='none')+theme(legend.position = 'none'), 
                       qho_N_plot+guides(color='none')+theme(legend.position = 'none'),
                       nrow=1,ncol=2)
qho_plot1 <- arrangeGrob(qho_plot,legend_n, nrow=2, heights = ggplot_heights2)
qho_title <- titles_plot('Quantum harmonic oscillator','E',expression(u[t] == i*u[xx]/2 - x^2*i*u/2))
qho_plot2 <- arrangeGrob(qho_title, qho_sol_plot, qho_plot1, nrow=3,heights =c(2,7,10))
# ggsave(qho_plot2, filename = 'Figures/new_figs/qho.png', width = 8, height = 8)
# ggsave(qho_plot2, filename = 'Figures/new_figs/qho.pdf', width = 8, height = 8)

## ----trans plot, fig.height=4, fig.width=8---------------------------------------------
## trans plots --------------------------------
trans.x <- trans_list$x; trans.t <- trans_list$t; trans_usol <- as.data.frame(trans_list$sol)
gg_trans_data <- cbind.data.frame(x = rep(trans.x,length(trans.t)), t = rep(trans.t,each=length(trans.x)), u = gather(trans_usol)[,-1])
color_sacle <- ramp.col(col = c("#fafafa", "#d7191c"),n=200)
trans_sol_plot <- ggplot(gg_trans_data, aes(x, t, fill= u)) + geom_tile()+
  scale_fill_gradientn(colors = color_sacle) +
  solution_theme+
  coord_cartesian(clip="off")
## trans noisy data plot ---------------------------
load(paste0(data_path, "more_test2/trans_SG_noise_seed_100_samp_100_snr_ada_lasso_pareto.RData"))
load(paste0(data_path, "more_test2/trans_SG_noise_seed_100_samp_100_snr_rudy_d_thred.RData"))
d_tols <- c(0.2,2,10)
names(trans_density_noise_rudy) <- paste0("rudy_", d_tols)
names(trans_density_noise_ada_lasso_pareto) <- 'our'
dens_snr_all <- c(trans_density_noise_ada_lasso_pareto, trans_density_noise_rudy)
dens_snr_all <- lapply(dens_snr_all, function(x) x[,11:32])

snr_db_seq <- c(seq(20, 60, 2), Inf) # SNR_dB grid
eta <- 10^(-snr_db_seq/20)
snr_db_seq1 <- c(seq(20, 60, 2), 65)
snr_db_seq2 <- c(seq(20, 60, 2), 65) # c(seq(20, 40, 2), 42.5, 45)
snr_db_seq3 <- c(as.character(seq(20, 60, 2)), expression(infinity)) # c(as.character(seq(20, 40, 2)), "...", "Inf")
# snr_db_seq3[-c((5*(0:12)+1),32)] <- ''
snr_db_seq3[c(1,3,5,7,9,11,13,15,17,19)+1] <- ''

trans_count_vector_snr <- data.frame(SNR_dB = 0, Correct = 0, Incorrect = 0, method0 = '0')
for(m_names in names(dens_snr_all)){
  temp <- data.frame(SNR_dB = snr_db_seq1,
                     Correct = colSums(dens_snr_all[[m_names]] == 1) / 100,
                     Incorrect = colSums(dens_snr_all[[m_names]] == 0) / 100,
                     method0 = m_names)
  
  trans_count_vector_snr <- rbind.data.frame(trans_count_vector_snr, temp)
}
trans_count_vector_snr <- trans_count_vector_snr[-1, ]
rownames(trans_count_vector_snr) <- seq(1, nrow(trans_count_vector_snr))
trans_noise_gather_df <- trans_count_vector_snr %>%
  gather("Condition", "Value",
         3:ncol(trans_count_vector_snr)-1)

method_label <- c("RAL", paste0("STRidge (d_tol=", d_tols,")"))
trans_count_vector_snr$method <- NA
for(i in seq_along(method_label)){
  trans_count_vector_snr$method[which(trans_count_vector_snr$method0==unique(trans_count_vector_snr$method0)[i])] <- method_label[i]
}

x_limits <- snr_db_seq
x_breaks <- pretty_breaks()(snr_db_seq)
# colours <- c("#a670b0",
#              "#a58d48")
# ggcol <- c("#F8766D","#00BA38","#619CFF")
# rudy_color <- colorRampPalette(c("#4bf5db","#6261ff"))(length(d_tols))
# colors <- c(ggcol[1],rudy_color)

trans_SNR_plot <- ggplot() + 
  geom_hline(yintercept = 0.8, lty=2)+
  geom_rect(aes(xmin=rect1$xmin,xmax=rect1$xmax,ymin=rect1$ymin,ymax=rect1$ymax),alpha=0.2,fill="#9de0e6")+
  scale_y_continuous(breaks = seq(0,1,0.2))+
  geom_line(data=trans_count_vector_snr, aes(x=SNR_dB, y=Correct,color=method),lwd=lwd_size) +
  geom_point(data=trans_count_vector_snr, alpha=alpha_value,stroke=stroke_value,aes(x=SNR_dB, y=Correct,shape=method,size=method,fill = method)) + 
  # ylab("Probability of correct identification") +  xlab("SNR(dB)")+
  ylab("Success rate") +  xlab("SNR(dB)")+
  # xlab(TeX("Noise magnitude ($\\eta$*sd(u))")) +
  scale_x_continuous(breaks = snr_db_seq2, labels = snr_db_seq3)+
  scale_size_manual(values=sizes, breaks=unique(trans_count_vector_snr$method), 
                    labels = legend_labs)+
  scale_fill_manual(values=colors, breaks=unique(trans_count_vector_snr$method), 
                    labels = legend_labs)+
  scale_shape_manual(values=shapes, breaks=unique(trans_count_vector_snr$method), 
                     labels = legend_labs)+
  scale_color_manual(values=colors, breaks=unique(trans_count_vector_snr$method),
                     labels = legend_labs)+
  labs(fill = 'Methods', shape='Methods', size='Methods',color='Methods')+
  # scale_y_continuous(labels = y_labels,breaks = y_breaks,limits = c(NA, max(y_breaks))) +
  plot_theme + legend_guides +
  guides(x = guide_axis_truncated(
    trunc_lower = c(-Inf, xend +20+ extra_x/2),
    trunc_upper = c(xstart +20+ extra_x/2, Inf)
  )) +
  annotate("segment", 
           x = myseg2$x, xend = myseg2$xend,
           y = myseg2$y + 0.05, yend = myseg2$yend)  +
  coord_cartesian(clip = "off", ylim = c(-0.0005, NA))

## trans noiseless plot ----------------------
load(paste0(data_path,"more_test2/trans_SG_noiseless_seed_10_success_rate_ada_lasso_pareto.RData"))
load(paste0(data_path,"more_test2/trans_SG_noiseless_seed_10_success_rate_rudy_d_thred.RData"))
num <- (seq(2,5,0.2))

d_tols <- c(0.2,2,10)
names(success_rate_rudy_list) <- paste0("rudy_", d_tols)
success_rate_our <- trans_density_noiseless_ada_lasso
dens_n_all <- c(list(our=success_rate_our), lapply(success_rate_rudy_list, function(x) x))

trans_count_vector_n <- data.frame(n = 0, Correct = 0, Incorrect = 0, method0 = '0')
for(m_names in names(dens_n_all)){
  temp <- data.frame(n = num,
                     Correct = colSums(dens_n_all[[m_names]] == 1) / 100,
                     Incorrect = colSums(dens_n_all[[m_names]] == 0) / 100,
                     method0 = m_names)
  
  trans_count_vector_n <- rbind.data.frame(trans_count_vector_n, temp)
}
trans_count_vector_n <- trans_count_vector_n[-1, ]

rownames(trans_count_vector_n) <- seq(1, nrow(trans_count_vector_n))
noise_gather_df_all <- trans_count_vector_n %>%
  gather("Condition", "Value",
         3:ncol(trans_count_vector_n)-1)

method_label <- c("RAL", paste0("STRidge (d_tol=", d_tols,")"))
trans_count_vector_n$method <- NA
for(i in seq_along(method_label)){
  trans_count_vector_n$method[which(trans_count_vector_n$method0==unique(trans_count_vector_n$method0)[i])] <- method_label[i]
}

# colours <- c("#a670b0",
#              "#a58d48")
# ggcol <- c("#F8766D","#00BA38","#619CFF")
# rudy_color <- colorRampPalette(c("#4bf5db","#6261ff"))(length(d_tols))
# colors <- c(ggcol[1],rudy_color)

xlables <- unname(latex2exp::TeX(paste0('$10^{',round((num),1),'}$')))
# xlables[-c(3*(0:5)+1)] <- ''
xlables[-(3*(0:(length(num)/2))+1)] <- ''
second_xlables <- round(10^num/(601*201)*100,2)
second_xlables[-(3*(0:(length(num)/2))+1)] <- ''

trans_N_plot <- ggplot() + 
  geom_hline(yintercept = 0.8, lty=2)+
  geom_rect(aes(xmin=rect1$xmin,xmax=rect1$xmax,ymin=rect1$ymin,ymax=rect1$ymax),alpha=0.2,fill="#9de0e6")+
  scale_y_continuous(breaks = seq(0,1,0.2))+
  geom_line(data=trans_count_vector_n, aes(x=n, y=Correct,color=method),lwd=lwd_size) +
  geom_point(data=trans_count_vector_n, alpha=alpha_value,stroke=stroke_value,aes(x=n, y=Correct,shape=method,size=method,fill = method)) + 
  scale_x_continuous(breaks = num, labels = xlables,
                     sec.axis = dup_axis(name='Percentage (%)', labels=second_xlables))+ 
  # labs(x="SNR(dB)", y="Success rate", title = TeX('Reaction-Diffusion, $\\alpha$=0.3, Inf $\\alpha$=0.5'))+
  # labs(x="N", y="Success rate")+
  labs(x="N", y=" ")+
  scale_size_manual(values=sizes, breaks=unique(trans_count_vector_n$method), 
                    labels = legend_labs)+
  scale_fill_manual(values=colors, breaks=unique(trans_count_vector_n$method), 
                    labels = legend_labs)+
  scale_shape_manual(values=shapes, breaks=unique(trans_count_vector_n$method), 
                     labels = legend_labs)+
  scale_color_manual(values=colors, breaks=unique(trans_count_vector_n$method),
                     labels = legend_labs)+
  labs(fill = 'Methods', shape='Methods', size='Methods',color='Methods')+
  # scale_y_continuous(labels = y_labels,breaks = y_breaks,limits = c(NA, max(y_breaks))) +
  plot_theme + legend_guides+
  coord_cartesian(clip = "off", ylim = c(-0.0005, NA))
# theme_classic(18)
## trans out -----------------------
legend_n <- get_legend(trans_N_plot+theme(legend.position='bottom'))
trans_plot <- arrangeGrob(trans_SNR_plot+guides(color='none')+theme(legend.position = 'none'), 
                        trans_N_plot+guides(color='none')+theme(legend.position = 'none'),
                        nrow=1,ncol=2)
trans_plot1 <- arrangeGrob(trans_plot,legend_n,nrow=2, heights = ggplot_heights2)
trans_title <- titles_plot('Transport','C',expression(u[t] == 3*u[x]))
trans_plot2 <- arrangeGrob(trans_title, trans_sol_plot,trans_plot1,nrow=3,heights =c(2,7,10))
# ggsave(trans_plot2, filename = 'Figures/new_figs/trans.png', width = 8, height = 8)
# ggsave(trans_plot2, filename = 'Figures/new_figs/trans.pdf', width = 8, height = 8)

## ----heat plot, fig.height=4, fig.width=8----------------------------------------------
## heat plots --------------------------------
heat.x <- heat_list$x; heat.t <- heat_list$t; heat_usol <- as.data.frame(heat_list$sol)
gg_heat_data <- cbind.data.frame(x = rep(heat.x,length(heat.t)), t = rep(heat.t,each=length(heat.x)), u = gather(heat_usol)[,-1])
color_sacle <- ramp.col(col = c("#fafafa", "#d7191c"),n=200)
heat_sol_plot <- ggplot(gg_heat_data, aes(x, t, fill= u)) + geom_tile()+ geom_contour(mapping=aes(z=u),color='grey80')+
  scale_fill_gradientn(colors = color_sacle) +
  solution_theme+
  coord_cartesian(clip="off")
## heat noisy data plot ---------------------------
load(paste0(data_path, "more_test2/heat_SG_noise_seed_100_samp_100_snr_ada_lasso_pareto.RData"))
load(paste0(data_path, "more_test2/heat_SG_noise_seed_100_samp_100_snr_rudy_d_thred.RData"))
d_tols <- c(0.2,2,10)
names(heat_density_noise_rudy) <- paste0("rudy_", d_tols)
names(heat_density_noise_ada_lasso_pareto) <- 'our'
dens_snr_all <- c(heat_density_noise_ada_lasso_pareto, heat_density_noise_rudy)

snr_db_seq <- c(seq(40, 60, 2), Inf) # SNR_dB grid
eta <- 10^(-snr_db_seq/20)
snr_db_seq1 <- c(seq(40, 60, 2), 65)
snr_db_seq2 <- c(seq(40, 60, 2), 65) # c(seq(20, 40, 2), 42.5, 45)
snr_db_seq3 <- c(as.character(seq(40, 60, 2)), expression(infinity)) # c(as.character(seq(20, 40, 2)), "...", "Inf")
snr_db_seq3[c(1,3,5,7,9)+1] <- ''

heat_count_vector_snr <- data.frame(SNR_dB = 0, Correct = 0, Incorrect = 0, method0 = '0')
for(m_names in names(dens_snr_all)){
  temp <- data.frame(SNR_dB = snr_db_seq1,
                     Correct = colSums(dens_snr_all[[m_names]] == 1) / 100,
                     Incorrect = colSums(dens_snr_all[[m_names]] == 0) / 100,
                     method0 = m_names)
  
  heat_count_vector_snr <- rbind.data.frame(heat_count_vector_snr, temp)
}
heat_count_vector_snr <- heat_count_vector_snr[-1, ]
rownames(heat_count_vector_snr) <- seq(1, nrow(heat_count_vector_snr))
heat_noise_gather_df <- heat_count_vector_snr %>%
  gather("Condition", "Value",
         3:ncol(heat_count_vector_snr)-1)

method_label <- c("RAL", paste0("STRidge (d_tol=", d_tols,")"))
heat_count_vector_snr$method <- NA
for(i in seq_along(method_label)){
  heat_count_vector_snr$method[which(heat_count_vector_snr$method0==unique(heat_count_vector_snr$method0)[i])] <- method_label[i]
}

x_limits <- snr_db_seq
x_breaks <- pretty_breaks()(snr_db_seq)
# colours <- c("#a670b0",
#              "#a58d48")
# ggcol <- c("#F8766D","#00BA38","#619CFF")
# rudy_color <- colorRampPalette(c("#4bf5db","#6261ff"))(length(d_tols))
# colors <- c(ggcol[1],rudy_color)

heat_SNR_plot <- ggplot() + 
  geom_hline(yintercept = 0.8, lty=2)+
  geom_rect(aes(xmin=rect1$xmin,xmax=rect1$xmax,ymin=rect1$ymin,ymax=rect1$ymax),alpha=0.2,fill="#9de0e6")+
  scale_y_continuous(breaks = seq(0,1,0.2))+
  geom_line(data=heat_count_vector_snr, aes(x=SNR_dB, y=Correct,color=method),lwd=lwd_size) +
  geom_point(data=heat_count_vector_snr,alpha=alpha_value,stroke=stroke_value,aes(x=SNR_dB, y=Correct,shape=method,size=method,fill = method)) +
  # ylab("Probability of correct identification") +  xlab("SNR(dB)")+
  ylab("Success rate") +  xlab("SNR(dB)")+
  # xlab(TeX("Noise magnitude ($\\eta$*sd(u))")) +
  scale_x_continuous(breaks = snr_db_seq2, labels = snr_db_seq3)+
  scale_size_manual(values=sizes, breaks=unique(heat_count_vector_snr$method), 
                    labels = legend_labs)+
  scale_fill_manual(values=colors, breaks=unique(heat_count_vector_snr$method), 
                    labels = legend_labs)+
  scale_shape_manual(values=shapes, breaks=unique(heat_count_vector_snr$method), 
                     labels = legend_labs)+
  scale_color_manual(values=colors, breaks=unique(bur_count_vector_snr$method),
                     labels = legend_labs)+
  labs(fill = 'Methods', shape='Methods', size='Methods',color='Methods')+
  # scale_y_continuous(labels = y_labels,breaks = y_breaks,limits = c(NA, max(y_breaks))) +
  plot_theme + legend_guides +
  guides(x = guide_axis_truncated(
    trunc_lower = c(-Inf, xend +20+ extra_x/2),
    trunc_upper = c(xstart +20+ extra_x/2, Inf)
  )) +
  annotate("segment", 
           x = myseg2$x, xend = myseg2$xend,
           y = myseg2$y + 0.05, yend = myseg2$yend)  +
  coord_cartesian(clip = "off", ylim = c(-0.0005, NA))

## heat noiseless plot ----------------------
load(paste0(data_path,"more_test2/heat_SG_noiseless_seed_10_success_rate_ada_lasso_pareto.RData"))
load(paste0(data_path,"more_test2/heat_SG_noiseless_seed_10_success_rate_rudy_d_thred.RData"))
num <- (seq(2,4.8,0.2))

d_tols <- c(0.2,2,10)
names(success_rate_rudy_list) <- paste0("rudy_", d_tols)
success_rate_our <- heat_density_noiseless_ada_lasso
dens_n_all <- c(list(our=success_rate_our), lapply(success_rate_rudy_list, function(x) x))

heat_count_vector_n <- data.frame(n = 0, Correct = 0, Incorrect = 0, method0 = '0')
for(m_names in names(dens_n_all)){
  temp <- data.frame(n = num,
                     Correct = colSums(dens_n_all[[m_names]] == 1) / 100,
                     Incorrect = colSums(dens_n_all[[m_names]] == 0) / 100,
                     method0 = m_names)
  
  heat_count_vector_n <- rbind.data.frame(heat_count_vector_n, temp)
}
heat_count_vector_n <- heat_count_vector_n[-1, ]

rownames(heat_count_vector_n) <- seq(1, nrow(heat_count_vector_n))
noise_gather_df_all <- heat_count_vector_n %>%
  gather("Condition", "Value",
         3:ncol(heat_count_vector_n)-1)

method_label <- c("RAL", paste0("STRidge (d_tol=", d_tols,")"))
heat_count_vector_n$method <- NA
for(i in seq_along(method_label)){
  heat_count_vector_n$method[which(heat_count_vector_n$method0==unique(heat_count_vector_n$method0)[i])] <- method_label[i]
}

# colours <- c("#a670b0",
#              "#a58d48")
# ggcol <- c("#F8766D","#00BA38","#619CFF")
# rudy_color <- colorRampPalette(c("#4bf5db","#6261ff"))(length(d_tols))
# colors <- c(ggcol[1],rudy_color)

xlables <- unname(latex2exp::TeX(paste0('$10^{',round((num),1),'}$')))
# xlables[c(1,3,5,7,9,11,13)+1] <- ''
xlables[-(3*(0:(length(num)/2))+1)] <- ''
second_xlables <- round(10^num/(501*151)*100,2)
second_xlables[-(3*(0:(length(num)/2))+1)] <- ''

heat_N_plot <- ggplot() + 
  geom_hline(yintercept = 0.8, lty=2)+
  geom_rect(aes(xmin=rect1$xmin,xmax=rect1$xmax,ymin=rect1$ymin,ymax=rect1$ymax),alpha=0.2,fill="#9de0e6")+
  scale_y_continuous(breaks = seq(0,1,0.2))+
  geom_line(data=heat_count_vector_n, aes(x=n, y=Correct,color=method),lwd=lwd_size) +
  geom_point(data=heat_count_vector_n,alpha=alpha_value,stroke=stroke_value,aes(x=n, y=Correct,shape=method,size=method,fill = method)) +
  scale_x_continuous(breaks = num, labels = xlables,
                     sec.axis = dup_axis(name='Percentage (%)', labels=second_xlables))+ 
  # labs(x="SNR(dB)", y="Success rate", title = TeX('Reaction-Diffusion, $\\alpha$=0.3, Inf $\\alpha$=0.5'))+
  # labs(x="N", y="Success rate")+
  labs(x="N", y=" ")+
  scale_size_manual(values=sizes, breaks=unique(heat_count_vector_n$method), 
                    labels = legend_labs)+
  scale_fill_manual(values=colors, breaks=unique(heat_count_vector_n$method), 
                    labels = legend_labs)+
  scale_shape_manual(values=shapes, breaks=unique(heat_count_vector_n$method), 
                     labels = legend_labs)+
  scale_color_manual(values=colors, breaks=unique(bur_count_vector_snr$method),
                     labels = legend_labs)+
  labs(fill = 'Methods', shape='Methods', size='Methods',color='Methods')+
  # scale_y_continuous(labels = y_labels,breaks = y_breaks,limits = c(NA, max(y_breaks))) +
  plot_theme + legend_guides
# theme_classic(18)
## heat out -----------------------
legend_n <- get_legend(heat_N_plot+theme(legend.position='bottom'))
heat_plot <- arrangeGrob(heat_SNR_plot+guides(color='none')+theme(legend.position = 'none'), 
                        heat_N_plot+guides(color='none')+theme(legend.position = 'none'),
                        nrow=1,ncol=2)
heat_plot1 <- arrangeGrob(heat_plot,legend_n,nrow=2, heights = ggplot_heights2)
heat_title <- titles_plot('Diffusion','D',expression(u[t] == u[xx]))
heat_plot2 <- arrangeGrob(heat_title, heat_sol_plot,heat_plot1,nrow=3,heights =c(2,7,10))
# ggsave(heat_plot2, filename = 'Figures/new_figs/heat.png', width = 8, height = 8)
# ggsave(heat_plot2, filename = 'Figures/new_figs/heat.pdf', width = 8, height = 8)

## ----keller plot, fig.height=4, fig.width=8--------------------------------------------
## keller plot ------------------------
keller.u <- as.data.frame(keller_list$sol$u)
keller.v <- as.data.frame(keller_list$sol$v)
keller.x <- keller_list$x
keller.t <- keller_list$t

gg_keller_data <- cbind.data.frame(
  x = rep(keller.x,length(keller.t)), t = rep(keller.t,each=length(keller.x)), 
  u = gather(keller.u)[,-1], v = gather(keller.v)[,-1]) %>% 
  reshape2::melt(id=1:2, variable.name='sol')

color_sacle <- ramp.col(col = c("#fafafa", "#d7191c"),n=200)
# color_sacle <- c(ramp.col(col = c("#2c7bb6", "#FFFFBF"),n=100),ramp.col(col = c("#FFFFBF", "#d7191c"),n=100))
color_sacle <- c(ramp.col(col = c("#3288BD", "#FFFFBF"),n=100),ramp.col(col = c("#FFFFBF", "#d7191c"),n=100))

# RColorBrewer::brewer.pal(n = 11, name = "Spectral")

keller_sol_plot <- ggplot(gg_keller_data, aes(x, t, fill= value))+
  geom_raster()+#geom_contour(aes(z=value),color='grey90')+
  facet_wrap(~sol,ncol=2)+
  # scale_fill_distiller(palette = "Spectral") +
  scale_fill_gradientn(colors=color_sacle)+
  solution_theme + 
  theme(strip.text = element_text(size = 16),
        # strip.background = element_rect(fill = 'grey90',color='grey90')
  ) +
  coord_cartesian(clip="off")

## keller noiseless ------------------
load(paste0(data_path,"more_test2/keller_segel_FD_noiseless_ham8_ada_rudy_back.RData"))

# num <- 10^(seq(2.4,5,0.2)) 
num <- (seq(2.4,5,0.2)) 
d_tols <- c(0.2,2,10)
names(dens_rudy_list) <- paste0("rudy_", d_tols)
success_rate_our <- dens_our
dens_n_all <- c(list(our=success_rate_our), lapply(dens_rudy_list, function(x) x))

keller_count_vector_n <- data.frame(n = 0, Correct = 0, Incorrect = 0, method0 = '0')
for(m_names in names(dens_n_all)){
  temp <- data.frame(n = num,
                     Correct = colSums(dens_n_all[[m_names]] == 1) / 100,
                     Incorrect = colSums(dens_n_all[[m_names]] == 0) / 100,
                     method0 = m_names)
  
  keller_count_vector_n <- rbind.data.frame(keller_count_vector_n, temp)
}
keller_count_vector_n <- keller_count_vector_n[-1, ]

rownames(keller_count_vector_n) <- seq(1, nrow(keller_count_vector_n))
noise_gather_df_all <- keller_count_vector_n %>%
  gather("Condition", "Value",
         3:ncol(keller_count_vector_n)-1)

method_label <- c("RAL", paste0("STRidge (d_tol=", d_tols,")"))
keller_count_vector_n$method <- NA
for(i in seq_along(method_label)){
  keller_count_vector_n$method[which(keller_count_vector_n$method0==unique(keller_count_vector_n$method0)[i])] <- method_label[i]
}

# colours <- c("#a670b0",
#              "#a58d48")
# ggcol <- c("#F8766D","#00BA38","#619CFF")
# rudy_color <- colorRampPalette(c("#4bf5db","#6261ff"))(length(d_tols))
# colors <- c(ggcol[1],rudy_color)

xlables <- unname(latex2exp::TeX(paste0('$10^{',round((num),1),'}$')))
# xlables[-c(3*(0:5)+1)] <- ''
xlables[-(3*(0:(length(num)/2))+1)] <- ''
second_xlables <- round(10^num/(1001*10001)*100,2)
second_xlables[-(3*(0:(length(num)/2))+1)] <- ''

keller_N_FD_plot <- ggplot() + 
  geom_hline(yintercept = 0.8, lty=2)+
  geom_rect(aes(xmin=rect1$xmin,xmax=rect1$xmax,ymin=rect1$ymin,ymax=rect1$ymax),alpha=0.2,fill="#9de0e6")+
  scale_y_continuous(breaks = seq(0,1,0.2))+
  geom_line(data=keller_count_vector_n, aes(x=n, y=Correct,color=method),lwd=lwd_size) +
  geom_point(data=keller_count_vector_n, alpha=alpha_value,stroke=stroke_value,aes(x=n, y=Correct,shape=method,size=method,fill = method)) + 
  scale_x_continuous(breaks = num, labels = xlables,
                     sec.axis = dup_axis(name='Percentage (%)', labels=second_xlables))+ 
  # labs(x="SNR(dB)", y="Success rate", title = TeX('Reaction-Diffusion, $\\alpha$=0.3, Inf $\\alpha$=0.5'))+
  labs(title='Finite difference')+
  # labs(x="N", y="Success rate")+
  labs(x="N", y=" ")+
  scale_size_manual(values=sizes, breaks=unique(keller_count_vector_n$method), 
                    labels = legend_labs)+
  scale_fill_manual(values=colors, breaks=unique(keller_count_vector_n$method), 
                    labels = legend_labs)+
  scale_shape_manual(values=shapes, breaks=unique(keller_count_vector_n$method), 
                     labels = legend_labs)+
  scale_color_manual(values=colors, breaks=unique(keller_count_vector_n$method),
                     labels = legend_labs)+
  labs(fill = 'Methods', shape='Methods', size='Methods',color='Methods')+
  # scale_y_continuous(labels = y_labels,breaks = y_breaks,limits = c(NA, max(y_breaks))) +
  plot_theme + legend_guides +
  coord_cartesian(clip = "off", ylim = c(NA, 1))

keller_count_vector_n_SG <- keller_count_vector_n
keller_count_vector_n_SG$Correct <- 0
keller_count_vector_n_SG$Incorrect <- 1

keller_N_SG_plot <- ggplot() + 
  geom_hline(yintercept = 0.8, lty=2)+
  geom_rect(aes(xmin=rect1$xmin,xmax=rect1$xmax,ymin=rect1$ymin,ymax=rect1$ymax),alpha=0.2,fill="#9de0e6")+
  scale_y_continuous(breaks = seq(0,1,0.2))+
  geom_line(data=keller_count_vector_n_SG, aes(x=n, y=Correct,color=method),lwd=lwd_size) +
  geom_point(data=keller_count_vector_n_SG,alpha=alpha_value,stroke=stroke_value,aes(x=n, y=Correct,shape=method,size=method,fill = method)) + 
  scale_x_continuous(breaks = num, labels = xlables,
                     sec.axis = dup_axis(name='Percentage (%)', labels=second_xlables))+ 
  # labs(x="SNR(dB)", y="Success rate", title = TeX('Reaction-Diffusion, $\\alpha$=0.3, Inf $\\alpha$=0.5'))+
  labs(title='Savitzky-Golay filter')+
  # labs(x="N", y="Success rate")+
  labs(x="N", y=" ")+
  scale_size_manual(values=sizes, breaks=unique(keller_count_vector_n$method), 
                    labels = legend_labs)+
  scale_fill_manual(values=colors, breaks=unique(keller_count_vector_n$method), 
                    labels = legend_labs)+
  scale_shape_manual(values=shapes, breaks=unique(keller_count_vector_n$method), 
                     labels = legend_labs)+
  scale_color_manual(values=colors, breaks=unique(keller_count_vector_n$method),
                     labels = legend_labs)+
  labs(fill = 'Methods', shape='Methods', size='Methods',color='Methods')+
  # scale_y_continuous(labels = y_labels,breaks = y_breaks,limits = c(NA, max(y_breaks))) +
  plot_theme + legend_guides +
  coord_cartesian(clip = "off", ylim = c(NA, 1))
## compose plots ---------------------
# layout_matrix <- matrix(c(1,2,1,3),ncol=2)
legend_n <- get_legend(keller_N_FD_plot+theme(legend.position='bottom'))
keller_plot <- arrangeGrob(keller_N_FD_plot+guides(color='none')+theme(legend.position = 'none'), 
                           keller_N_SG_plot+guides(color='none')+theme(legend.position = 'none'),
                         nrow=1,ncol=2)
keller_plot1 <- arrangeGrob(keller_plot,legend_n,nrow=2, heights = ggplot_heights2)
keller_title <- titles_plot('Keller-Segel','D',
                            expression(atop(u[t] == 0.1*u[xx] - u[x]*v[x] - uv[xx],
                                            v[t] == v[xx] + u - v)),
                            eq_bd_h=0.9,eq_bd_w=0.7,title_size=16)
keller_plot2 <- arrangeGrob(keller_title, keller_sol_plot, keller_plot1,nrow=3,heights =c(2,7,10))
# ggsave(keller_plot2, filename = 'Figures/new_figs/keller.png', width = 8, height = 8)
# ggsave(keller_plot2, filename = 'Figures/new_figs/keller.pdf', width = 8, height = 8)

## put in main paper plots ---------------------------
# NS_plot3 <- arrangeGrob(NS_title, NS_sol_plot, NS_plot, nrow=3, heights =c(3,7,10))
# RD_plot3 <- arrangeGrob(RD_title, RD_sol_plot, RD_plot, nrow=3,heights =c(3,7,10))
# qho_plot3 <- arrangeGrob(qho_title, qho_sol_plot, qho_plot, nrow=3,heights =c(3,7,10))
# bur_plot3 <- arrangeGrob(bur_title, bur_sol_plot, bur_plot, nrow=3,heights =c(3,7,10))
# # keller_plot3 <- arrangeGrob(keller_title, keller_sol_plot, keller_plot,nrow=3,heights =c(3,7,10))
# # put_in_main0 <- arrangeGrob(NS_plot3,RD_plot3,qho_plot3,keller_plot3, nrow=2,ncol=2)
# put_in_main0 <- arrangeGrob(bur_plot3,NS_plot3,RD_plot3,qho_plot3, nrow=2,ncol=2)
# put_in_main <- arrangeGrob(put_in_main0,legend_n,nrow=2, heights = c(30,1))
# 
# 
# ggsave(put_in_main, filename = 'Figures/thesis_plot/out_all_out.png', width = 13, height = 13)
# ggsave(put_in_main, filename = 'Figures/thesis_plot/out_all_out.pdf', width = 13, height = 13)

## main paperplot -----------------------
NS_plot31 <- arrangeGrob(NS_title, arrangeGrob(NS_sol_plot, NS_plot, nrow=1, widths=c(1,2)), nrow=2, heights =c(3,12))
RD_plot31 <- arrangeGrob(RD_title, arrangeGrob(RD_sol_plot, RD_plot, nrow=1, widths=c(1,2)), nrow=2, heights =c(3,12))
qho_plot31 <- arrangeGrob(qho_title, arrangeGrob(qho_sol_plot, qho_plot, nrow=1, widths=c(1,2)), nrow=2, heights =c(3,12))
bur_plot31 <- arrangeGrob(bur_title, arrangeGrob(bur_sol_plot, bur_plot, nrow=1, widths=c(1,2)), nrow=2, heights =c(3,12))
cable_plot31 <- arrangeGrob(cable_title, arrangeGrob(cable_sol_plot, cable_plot, nrow=1, widths=c(1,2)), nrow=2, heights =c(3,12))

put_in_main01 <- arrangeGrob(bur_plot31,cable_plot31,NS_plot31,RD_plot31,qho_plot31, ncol=1)
put_in_main1 <- arrangeGrob(put_in_main01,legend_n,nrow=2, heights = c(30,1))

ggsave(put_in_main1, filename = 'Figures/new_figs/out_all1.png', width = 14, height = 16)
ggsave(put_in_main1, filename = 'Figures/new_figs/out_all1.pdf', width = 14, height = 16)
ggsave(put_in_main1, filename = 'Figures/new_figs/out_all1.svg', width = 14, height = 16)
## supporting doc plot ------------------------
ad_plot31 <- arrangeGrob(ad_title, arrangeGrob(ad_sol_plot, ad_plot, nrow=1, widths=c(1,2)), nrow=2, heights =c(3,12))
kdv_plot31 <- arrangeGrob(kdv_title, arrangeGrob(kdv_sol_plot, kdv_plot, nrow=1, widths=c(1,2)), nrow=2, heights =c(3,12))
trans_plot31 <- arrangeGrob(trans_title, arrangeGrob(trans_sol_plot, trans_plot, nrow=1, widths=c(1,2)), nrow=2, heights =c(3,12))
heat_plot31 <- arrangeGrob(heat_title, arrangeGrob(heat_sol_plot, heat_plot, nrow=1, widths=c(1,2)), nrow=2, heights =c(3,12))
# cable_plot31 <- arrangeGrob(cable_title, arrangeGrob(cable_sol_plot, cable_plot, nrow=1, widths=c(1,2)), nrow=2, heights =c(3,12))

put_in_main02 <- arrangeGrob(ad_plot31,kdv_plot31,trans_plot31,heat_plot31, ncol=1)
put_in_main2 <- arrangeGrob(put_in_main02,legend_n,nrow=2, heights = c(30,1))

ggsave(put_in_main2, filename = 'Figures/new_figs/out_all2.png', width = 14, height = 15)
ggsave(put_in_main2, filename = 'Figures/new_figs/out_all2.pdf', width = 14, height = 15)
ggsave(put_in_main2, filename = 'Figures/new_figs/out_all2.svg', width = 14, height = 15)
