#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
  library(svglite); library(ragg)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg)) dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/")) else getwd()
repo_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = FALSE)
root <- Sys.getenv("PAPER1_PROJECT_ROOT", unset = repo_root)
stats_dir <- file.path(root, "stats", "supplementary_sensitivity")
out_dir <- file.path(root, "nature_fig", "Supplementary_Fig_3")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
source(file.path(repo_root, "r_figures", "shared", "paper1_figure_style_config.R"))
theme_set(theme_paper1())

read_required <- function(path) {if (!file.exists(path)) stop("Missing: ", path); read.csv(path, stringsAsFactors=FALSE, check.names=FALSE)}
features <- c("S","source_sink_gini_gap","reciprocal_fraction","dense_fraction")
models <- c("role_reciprocity_hybrid","density_only","source_concentration", "overlap_synchronization","weight_strengthening","reciprocity_only")
model_labels <- c(role_reciprocity_hybrid="Role-reciprocity hybrid", density_only="Density only",
                  source_concentration="Source concentration", overlap_synchronization="Overlap synchronization",
                  weight_strengthening="Weight strengthening", reciprocity_only="Reciprocity only")
model_colours <- c(
  role_reciprocity_hybrid = "#D89000",
  density_only = "#3B6FB6",
  source_concentration = "#00897B",
  overlap_synchronization = "#B85C90",
  weight_strengthening = "#4F555A",
  reciprocity_only = "#7E6AB0"
)
model_shapes <- c(
  role_reciprocity_hybrid = 21,
  density_only = 22,
  source_concentration = 24,
  overlap_synchronization = 23,
  weight_strengthening = 25,
  reciprocity_only = 4
)
model_short_labels <- c(
  role_reciprocity_hybrid = "Hybrid",
  density_only = "Density",
  source_concentration = "Source",
  overlap_synchronization = "Overlap",
  weight_strengthening = "Weight",
  reciprocity_only = "Reciprocity"
)
figure_blue <- paper1_palette$pre
figure_orange <- paper1_palette$mid
rank_colours <- c(
  "1" = "#486F9F",
  "2" = "#6F91BD",
  "3" = "#69737D",
  "4" = "#808990",
  "5" = "#989FA5",
  "6" = "#B0B5B9"
)
rank_sizes <- c("1"=6.15,"2"=5.75,"3"=5.35,"4"=4.95,"5"=4.55,"6"=4.15)

emp_table <- read_required(file.path(root,"stats","phase15_elbow_empirical_vector.csv"))
emp <- emp_table$empirical_effect[match(features, emp_table$feature)]

fit_models <- function(summary_data) {
  out <- list(); cursor <- 0
  for (m in unique(summary_data$model)) {
    z <- summary_data %>% filter(model==m) %>% arrange(beta)
    base <- z %>% filter(beta==0)
    cand <- z %>% filter(beta>0)
    for (i in seq_len(nrow(cand))) {
      sim <- as.numeric(cand[i,features] - base[1,features])
      dist <- sqrt(sum((sim/(abs(emp)+.Machine$double.eps) - emp/(abs(emp)+.Machine$double.eps))^2))
      cursor <- cursor+1; out[[cursor]] <- data.frame(model=m,beta=cand$beta[i],distance=dist)
    }
  }
  bind_rows(out) %>% group_by(model) %>% slice_min(distance,n=1,with_ties=FALSE) %>% ungroup() %>% arrange(distance) %>% mutate(rank=row_number())
}

# (a) Network-size sensitivity
size_raw <- read_required(file.path(stats_dir,"supp_fig3_size_sensitivity_raw.csv"))
size_mean <- size_raw %>% group_by(n_nodes,model,beta) %>% summarise(across(all_of(features),mean),.groups="drop")
size_fit <- size_mean %>% group_by(n_nodes) %>% group_modify(~fit_models(.x)) %>% ungroup()
write.csv(size_fit,file.path(stats_dir,"supp_fig3_size_sensitivity_fit.csv"),row.names=FALSE)

p3a <- size_fit %>%
  mutate(model=factor(model,levels=rev(models),labels=rev(unname(model_short_labels[models]))),
         rank=factor(rank,levels=1:6)) %>%
  ggplot(aes(n_nodes,model))+
  geom_point(aes(size=rank,colour=rank),shape=21,stroke=.76,fill="white")+
  geom_text(aes(label=rank),size=FIG_HEATMAP_TEXT_PT/2.845276,fontface="bold",colour=paper1_palette$text_dark)+
  scale_size_manual(values=rank_sizes,guide="none")+
  scale_colour_manual(values=rank_colours,guide="none")+
  scale_x_continuous(breaks=c(30,45,60),expand=expansion(mult=c(.12,.12)))+
  labs(title=NULL,x="Simulated DMWA nodes",y=NULL)+
  theme(panel.grid.major.x=element_line(colour="#E6EAED",linewidth=.32),
        panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),
        axis.ticks.y=element_blank(),plot.margin=margin(5.5,6,5.5,5.5))

# (b) Replicate-count convergence
raw <- read_required(file.path(root,"simulation","competing_mechanism_raw.csv"))
k_grid <- c(10,20,40,60,80); B <- 200; set.seed(20260714)
boot_rows <- vector("list",length(k_grid)*B); cursor <- 0
for (k in k_grid) for (b in seq_len(B)) {
  ids <- sample(1:80,k,replace=TRUE); w <- as.data.frame(table(ids),stringsAsFactors=FALSE)
  w$replicate <- as.integer(as.character(w$ids)); w$weight <- as.numeric(w$Freq)
  sm <- raw %>% inner_join(w[,c("replicate","weight")],by="replicate") %>%
    group_by(model,beta) %>% summarise(across(all_of(features),~weighted.mean(.x,weight)),.groups="drop")
  f <- fit_models(sm); cursor <- cursor+1
  boot_rows[[cursor]] <- f %>% transmute(replicates=k,bootstrap=b,model,rank,distance)
}
boot_rank <- bind_rows(boot_rows)
rank_prob <- boot_rank %>% group_by(replicates,model) %>% summarise(probability_best=mean(rank==1),median_rank=median(rank),.groups="drop")
write.csv(rank_prob,file.path(stats_dir,"supp_fig3_replicate_convergence.csv"),row.names=FALSE)

p3b <- rank_prob %>%
  mutate(model=factor(model,levels=rev(models),labels=rev(unname(model_short_labels[models])))) %>%
  ggplot(aes(replicates,model,size=probability_best))+
  geom_point(shape=21,stroke=.72,colour=figure_blue,fill="white")+
  scale_size_continuous(name=expression(pi[m]^{(1)}),limits=c(0,1),breaks=c(.25,.50,.75),
                        labels=c("0.25","0.50","0.75"),range=c(1.25,6.2))+
  scale_x_continuous(breaks=k_grid,expand=expansion(mult=c(.08,.08)))+
  labs(title=NULL,x="Replicates per setting",y=NULL)+
  guides(size=guide_legend(direction="horizontal",title.position="left",label.position="bottom",
                           keywidth=grid::unit(12,"pt"),keyheight=grid::unit(9,"pt")))+
  theme(panel.grid.major.x=element_line(colour="#E6EAED",linewidth=.32),
        panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),
        axis.ticks.y=element_blank(),legend.position="bottom",legend.justification="left",
        legend.title=element_text(size=7.6),legend.text=element_text(size=7.6),
        plot.margin=margin(5.5,6,5.5,5.5))

# (c) Motif-sampling convergence
tri <- read_required(file.path(stats_dir,"supp_fig3_triad_sampling_convergence.csv"))
tri_long <- tri %>% select(model,replicate,sample_triples,abs_error_reciprocal,abs_error_dense) %>%
  pivot_longer(starts_with("abs_error"),names_to="metric",values_to="absolute_error") %>%
  mutate(metric=recode(metric,abs_error_reciprocal="Reciprocal fraction",abs_error_dense="Dense fraction"))
tri_sum <- tri_long %>% group_by(sample_triples,metric) %>%
  summarise(mean_error=mean(absolute_error),low=quantile(absolute_error,.025),high=quantile(absolute_error,.975),.groups="drop")
write.csv(tri_sum,file.path(stats_dir,"supp_fig3_triad_sampling_summary.csv"),row.names=FALSE)

p3c <- tri_sum %>%
  mutate(metric=factor(metric,levels=c("Reciprocal fraction","Dense fraction")),
         sample_label=factor(format(sample_triples,big.mark=",",scientific=FALSE,trim=TRUE),
                             levels=rev(c("250","500","1,000","2,500","5,000")))) %>%
  ggplot(aes(mean_error,sample_label,colour=metric,shape=metric))+
  geom_errorbar(aes(xmin=low,xmax=high),orientation="y",width=.17,linewidth=.55,
                position=position_dodge(width=.34))+
  geom_point(fill="white",size=2.5,stroke=.68,position=position_dodge(width=.34))+
  scale_colour_manual(values=c("Reciprocal fraction"=figure_blue,"Dense fraction"=figure_orange),
                      labels=c("Reciprocal fraction"="Reciprocal","Dense fraction"="Dense"),name=NULL)+
  scale_shape_manual(values=c("Reciprocal fraction"=21,"Dense fraction"=24),name=NULL)+
  scale_x_continuous(limits=c(0,.305),breaks=c(0,.1,.2,.3),expand=expansion(mult=c(0,.01)))+
  labs(title=NULL,x=expression(MAE[F](n[tri])),y=expression(n[tri]))+
  guides(colour=guide_legend(nrow=1,byrow=TRUE,override.aes=list(shape=c(21,24),fill="white"),
                             label.position="right"),shape="none")+
  theme(panel.grid.major.x=element_line(colour="#E6EAED",linewidth=.32),
        panel.grid.major.y=element_blank(),panel.grid.minor=element_blank(),
        axis.ticks.y=element_blank(),legend.position="bottom",legend.justification="left",
        legend.key.width=grid::unit(8,"pt"),legend.text=element_text(size=7.2),
        plot.margin=margin(5.5,5.5,5.5,5.5))

title_strip <- function(label) {
  wrap_elements(
    full = grid::textGrob(
      label,
      x = grid::unit(3.0, "mm"),
      y = grid::unit(0.50, "npc"),
      just = c("left", "center"),
      gp = grid::gpar(
        fontfamily = "Arial",
        fontsize = FIG_PANEL_PT,
        fontface = "bold",
        col = paper1_palette$text_dark
      )
    ),
    clip = FALSE
  )
}

panel3a <- title_strip(expression("(a) Ranking, "*r[m])) / p3a + plot_layout(heights=c(0.13,1))
panel3b <- title_strip("(b) Replicate stability") / p3b + plot_layout(heights=c(0.13,1))
panel3c <- title_strip("(c) Motif-sampling error") / p3c + plot_layout(heights=c(0.13,1))
fig <- wrap_plots(
  wrap_elements(full=panel3a,clip=FALSE),
  wrap_elements(full=panel3b,clip=FALSE),
  wrap_elements(full=panel3c,clip=FALSE),
  ncol=3,
  widths=c(1.00,1.08,1.12)
)
w <- 183/25.4; h <- 82/25.4; base <- file.path(out_dir,"Supplementary_Fig_3")
svglite::svglite(paste0(base,".svg"),width=w,height=h,bg="white"); print(fig); dev.off()
grDevices::cairo_pdf(paste0(base,".pdf"),width=w,height=h,family="Arial"); print(fig); dev.off()
ragg::agg_png(paste0(base,".png"),width=w,height=h,units="in",res=600,background="white"); print(fig); dev.off()
ragg::agg_tiff(paste0(base,".tiff"),width=w,height=h,units="in",res=600,background="white",compression="lzw"); print(fig); dev.off()
write.csv(data.frame(file=paste0(base,c(".svg",".pdf",".png",".tiff")),width_mm=183,height_mm=82,dpi=c(NA,NA,600,600)),
          file.path(out_dir,"Supplementary_Fig_3_export_manifest.csv"),row.names=FALSE)
cat("Supplementary Fig. 3 rendered.\n")
print(size_fit); print(rank_prob); print(tri_sum)
