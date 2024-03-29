#!/usr/bin/env Rscript

library(ggplot2)
library(gridExtra)
args <- commandArgs(trailingOnly=TRUE)
f_gs <- args[1]
real_dir <- args[2]
eval_dir <- args[3]
prefix <- args[4]
eval_tool <- args[5] # tail(strsplit(prefix, '/')[[1]], 1)

sourceDir <- this.path::here() # dirname(sys.frame(1)$ofile)
if (length(args) >= 6 && grepl('generate-Chinese', args[6], fixed = TRUE)) {
    source(paste(sourceDir, "/semantic2word-Chinese.r", sep=""))
    prefix <- paste0(prefix, "Chinese.")
} else {
    source(paste(sourceDir, "/semantic2word-English.r", sep=""))
    prefix <- paste0(prefix, "English.")
}

if (length(args) >= 6 && grepl('rmdup-sim', args[6], fixed = TRUE)) {
    sim.strategy <- GVAR_remove_duplicates_then_simulate
} else if (length(args) >= 6 && grepl('sim-rmdup', args[6], fixed = TRUE)) {
    sim.strategy <- GVAR_simulate_then_remove_duplicates
} else {
    sim.strategy <- '(MISSING-SIMULATION-STRATEGY)'
}

# plot mut
addline_format <- function(x,...){
  gsub('\\s','\n',x)
}

gs_df <- read.table(f_gs, header = F, sep = "\t")
#colnames(gs_df) <- c('mut', 'AF', 'Anno', 'CNT', 'AF_Class')
colnames(gs_df) <- c('mut', 'AF', 'gene', 'Anno', 'CNT', 'AF_Class', 'VAR_TYPE')
gs_df$AF <- gs_df$AF/5
gs_df$Class <- rep('GS_AF', nrow(gs_df))
gs_df$Batch <- rep('GS', nrow(gs_df))
df <- data.frame()
batches <- c(c(1:12), c(1:12))
for (dir in c(eval_dir, real_dir)){
  index <- 0
  for (f in list.files(dir, recursive = TRUE, pattern = '_simmut.info.txt')){
    index <- index + 1
    filename <- paste(dir, f, sep = '/')
    tmpdf <- read.table(filename, header = F, sep = '\t', col.names = c('mut','AF','Company', 'Class'), stringsAsFactors = FALSE)
    tmpdf$Batch <- rep(batches[index], nrow(tmpdf))
    df <- rbind(df, tmpdf)
  }
}

df$Batch <- factor(df$Batch, levels = c(1:12))

for (vartype in c('SNP', 'INDEL')) {
gs1df <- subset(gs_df, VAR_TYPE == vartype)

myplot <- function(n){
  start <- start_list[n]
  end   <- start_list[n] + 24
  if (end > nrow(gs1df)){ end <- nrow(gs1df) }
  df1 <- subset(df, df$mut %in% gs1df$mut[start:end])
  df1$mut <- factor(df1$mut, levels = gs1df$mut[start:end])
  cnt <- gs1df[start:end, c(1,4)]
  cnt$Company <- rep('ILM', nrow(cnt))
  info <- gs1df[start:end, c(1,1)]
  colnames(info) <- c("mut", 'pos')
  info$pos <- sapply(strsplit(info$mut, '-'), function(x){paste(x[1:2],collapse=':')})
  info$Company <- rep('IDT', nrow(info))
  p <- ggplot() +
       geom_point(data = df1, aes(x = mut, y = AF, color = Class, shape = Batch), size = 8, alpha = 0.7,
                  position = position_jitterdodge(dodge.width = 0.5), na.rm = TRUE) +
       geom_point(data = gs1df[start:end,], aes(x = mut, y = AF, color='GS'), shape = 95, size = 6) +
       facet_grid(Company~., scales = 'free') +
       scale_shape_manual(values = c(65:76)) +
       scale_color_manual(values = c('black', '#f8766d', '#00ba38'), labels = c('Reference truth', 'Real variant', 'Simulated variant')) +
       labs(title = paste(tools::toTitleCase(eval_tool), sim.strategy, sep = ' '),
            x = GVAR_Mutation, y = GVAR_Allele_Fraction) +
       scale_x_discrete(labels = addline_format(gs1df$Anno[start:end])) +
       guides(color = guide_legend(title="Allele fraction type")) +
       theme(panel.spacing = unit(10, "lines"),
             text = element_text(size = 16),
             axis.text = element_text(size = 16, color = "black"),
             axis.text.x = element_text(angle = 45,vjust = 0.75, hjust = 1),
             axis.title.x = element_text(vjust = 0), axis.ticks.x = element_blank(),
             plot.title = element_text(hjust = 0.5),
             panel.grid = element_blank()) + coord_cartesian(clip = "off") +
             geom_text(data = cnt, aes(x = mut, label = gs1df$CNT[start:end], y = -Inf), vjust = 1 ,size = 7) +
             geom_text(data = data.frame(tag = "CNT", Company = "ILM"), aes(x = Inf, label = tag, y = -Inf), vjust = 1.1, size = 7) +
             geom_text(data = info, aes(x = mut, label = pos, y = -Inf), vjust = 0.2, hjust = 1.1, size = 7, angle = 90) +
             geom_text(data = data.frame(tag = "", Company = "IDT"), aes(x = Inf, label = tag, y = -Inf), vjust = 1.9, size = 7) +
             geom_vline(xintercept = seq(1.5, length(unique(df1$mut)) - 0.5, 1), lwd = 1, linetype = 'dashed', colour = 'gray')
  return(p)
}

start_list <- seq(1,nrow(gs1df),25) # n mutations each plot

pl <- lapply(1:length(start_list),myplot)
ml <- marrangeGrob(pl, nrow = 1, ncol = 1, top = NULL)
ggsave(paste(prefix, 'ctDNA.mut.', eval_tool, '.', vartype, '.pdf', sep=''), ml, width = 32, height = 15)
}
# plot mse

f_real <- paste(real_dir, 'realdata.txt', sep = '/')
f_sim <- paste(eval_dir, 'simulation.txt', sep = '/')
realdf <- read.table(f_real, header = F, sep = "\t", row.names = 1)
simdf  <- read.table(f_sim, header = F ,sep = "\t" , row.names = 1)

# calculate Z-score, Mean, Variance
df <- data.frame()
df_var <- data.frame()
df_mean <- data.frame()

for (mut in rownames(realdf)){
x <- realdf[mut,1:12][!is.na(realdf[mut,1:12])]
y <- simdf[mut,1:12][!is.na(simdf[mut,1:12])]
if(length(x)!=0 && mean(x)!=0 && length(y)!=0 && mean(y)!=0){
zx <- (x-mean(x))/sd(x)
zy <- (y-mean(y))/sd(y)
zz <- rep('IDT',12)
df <- rbind(df, data.frame(x = zx, y = zy, Company=zz, mut=mut))
df_var <- rbind(df_var, data.frame(x = var(x), y = var(y), Company=zz, mut = mut))
df_mean <- rbind(df_mean, data.frame(x = mean(x), y =mean(y), Company=zz, mut=mut))
}
x <- realdf[mut,13:24][!is.na(realdf[mut,13:24])]
y <- simdf[mut,13:24][!is.na(simdf[mut,13:24])]
if(length(x)!=0 && mean(x)!=0 && length(y)!=0 && mean(y)!=0){
zx <- (x-mean(x))/sd(x)
zy <- (y-mean(y))/sd(y)
zz <- rep('ILM',12)
df <- rbind(df, data.frame(x = zx, y = zy, Company=zz, mut=mut))
df_var <- rbind(df_var, data.frame(x = var(x), y = var(y), Company = zz, mut = mut))
df_mean <- rbind(df_mean, data.frame(x = mean(x), y =mean(y), Company = zz, mut=mut))
}
}

# Add AF_level tag
af_level <- c()
for(s in df$mut){af_level <- append(af_level ,subset(gs_df, mut == s)[1,6])}
df$AF_level <- af_level
af_level <- c()
for(s in df_mean$mut){af_level <- append(af_level ,subset(gs_df, mut == s)[1,6])}
df_mean$AF_level <- af_level
af_level <- c()
for(s in df_var$mut){af_level <- append(af_level ,subset(gs_df, mut == s)[1,6])}
df_var$AF_level <- af_level

# Sort Z-score by Company and AF_level
IDT_l <- df[df$AF_level=="Low"&df$Company=="IDT",]
IDT_l$x <- sort(IDT_l$x)
IDT_l$y <- sort(IDT_l$y)
IDT_m <- df[df$AF_level=="Medium"&df$Company=="IDT",]
IDT_m$x <- sort(IDT_m$x)
IDT_m$y <- sort(IDT_m$y)
IDT_h <- df[df$AF_level=="High"&df$Company=="IDT",]
IDT_h$x <- sort(IDT_h$x)
IDT_h$y <- sort(IDT_h$y)
ILM_l <- df[df$AF_level=="Low"&df$Company=="ILM",]
ILM_l$x <- sort(ILM_l$x)
ILM_l$y <- sort(ILM_l$y)
ILM_m <- df[df$AF_level=="Medium"&df$Company=="ILM",]
ILM_m$x <- sort(ILM_m$x)
ILM_m$y <- sort(ILM_m$y)
ILM_h <- df[df$AF_level=="High"&df$Company=="ILM",]
ILM_h$x <- sort(ILM_h$x)
ILM_h$y <- sort(ILM_h$y)
df <- rbind(IDT_l, IDT_m, IDT_h, ILM_l, ILM_m, ILM_h)

df$AF_level <- factor(df$AF_level, levels = c('Low', 'Medium', 'High'))
#df_mean$AF_level <- factor(df_mean$AF_level, levels = c('Low', 'Medium', 'High'))
#df_var$AF_level <- factor(df_var$AF_level, levels = c('Low', 'Medium', 'High'))

#df$AF_level <- factor(df$AF_level, levels = c('Low', 'Medium', 'High'), labels = c(GVAR_Low, GVAR_Medium, GVAR_High))
df_mean$AF_level <- factor(df_mean$AF_level, levels = c('Low', 'Medium', 'High'), labels = c(GVAR_Low, GVAR_Medium, GVAR_High))
df_var$AF_level <- factor(df_var$AF_level, levels = c('Low', 'Medium', 'High'), labels = c(GVAR_Low, GVAR_Medium, GVAR_High))

df_mean$x <- log(df_mean$x, 10)
df_mean$y <- log(df_mean$y, 10)
df_var$x <- log(df_var$x, 10)
df_var$y <- log(df_var$y, 10)

# Calculate MSE
mse <- function(data){
  colnames(data) <- c("x", "y", "Company", 'mut', 'AF_level')
  ybar = mean(data$y)
  ss_residuals <- sum((data$y - data$x) ^ 2)
  mse <- format(ss_residuals / nrow(data), digits = 4)
  return(mse)
}

df_mean_IDT <- df_mean[df_mean$Company=="IDT",]
df_mean_ILM <- df_mean[df_mean$Company=="ILM",]
df_var_IDT <- df_var[df_var$Company=="IDT",]
df_var_ILM <- df_var[df_var$Company=="ILM",]

#df_qqplot$AF_level <- factor(df_qqplot$AF_level, levels = c('Low', 'Medium', 'High'), labels = c(GVAR_Low, GVAR_Medium, GVAR_High))
#df_mean$AF_level <- factor(df_mean$AF_level, levels = c('Low', 'Medium', 'High'), labels = c(GVAR_Low, GVAR_Medium, GVAR_High))
#df_var$AF_level <- factor(df_var$AF_level, levels = c('Low', 'Medium', 'High'), labels = c(GVAR_Low, GVAR_Medium, GVAR_High))

qq_anno_text <- data.frame(Company = rep(c('IDT', 'ILM'),3), AF_level = rep(c('Low', 'Medium', 'High'), c(2,2,2)),label = sprintf("italic(MSE)==%s",c(mse(IDT_l), mse(ILM_l),mse(IDT_m), mse(ILM_m), mse(IDT_h), mse(ILM_h))))
qq_anno_text$AF_level <- factor(qq_anno_text$AF_level, levels=c('Low', 'Medium', 'High'))
#qq_anno_text$AF_level <- factor(qq_anno_text$AF_level, levels=c('Low', 'Medium', 'High'), labels = c(GVAR_Low, GVAR_Medium, GVAR_High))
mean_anno_text <- data.frame(Company = c('IDT', 'ILM'), label = sprintf("italic(MSE)==%s",c(mse(df_mean_IDT), mse(df_mean_ILM))))
var_anno_text <- data.frame(Company = c('IDT', 'ILM'), label = sprintf("italic(MSE)==%s",c(mse(df_var_IDT), mse(df_var_ILM))))

# QQ-plot
int <- rep(0, 6)
slope <- rep(1, 6)
Company <- rep(c('IDT', 'ILM'),3)
AF_level <- rep(c('Low', 'Medium', 'High'), 2) # rep(c(GVAR_Low, GVAR_Medium, GVAR_High), 2) # rep(c('Low', 'Medium', 'High'), 2)
dl <- data.frame(int, slope, Company, AF_level)
p1 <- ggplot() + 
      geom_point(data = df, aes(x = x, y = y), size = 3, alpha = 0.2, shape = 4) +
      #facet_grid(factor(AF_level, levels = c('Low', 'Medium', 'High'))~Company) +
      facet_grid(factor(AF_level, levels = c('Low', 'Medium', 'High'), labels = c(GVAR_LowFA, GVAR_MediumFA, GVAR_HighFA))~Company) +
      labs(x = GVAR_Lbx_high_Allele_Fraction_Z_score,
           y = GVAR_Simulation_Allele_Fraction_Z_score,
           title = paste(tools::toTitleCase(eval_tool), ' ', GVAR_QQplot, ' ', sim.strategy, sep = ''),
           col = GVAR_AF_level, shape = GVAR_AF_level) +
      theme(text = element_text(size = 16),
            axis.text = element_text(size = 16, color = 'black'),
            plot.title = element_text(hjust = 0.5)) +
      geom_text(data = qq_anno_text, aes(x = -0.6, y = 2, label = label), parse = T, size = 7) +
      geom_abline(data = dl, aes(slope = slope, intercept = int), color = 'red', linetype = 1)
      #scale_x_continuous(sec.axis = sec_axis(~ . , name = GVAR_AF_level, breaks = NULL, labels = NULL))
ggsave(paste(prefix, 'ctDNA.qqplot.', eval_tool, '.pdf', sep = ''), p1, width = 8, height = 12)

# Variance scatter plot
p2 <- ggplot(data = df_var, aes(x = x, y = y, color = AF_level)) +
      geom_point(aes(shape = AF_level), size=3, alpha = 0.6) +
      geom_abline(aes(slope = 1, intercept = 0), color = 'black', linetype = 1) +
      geom_text(data = var_anno_text, aes(x = -5, y = -3, label = label), color = 'black', parse = T, size = 7) +
      facet_grid(.~Company) +
      labs(x = bquote(.(GVAR_LBX_HIGH) ~ log[10] ~ .(GVAR_Allele_Fraction_Variance)),
           y = bquote(.(GVAR_Simulation_space) ~ log[10] ~ .(GVAR_Allele_Fraction_Variance)),
           title = paste(tools::toTitleCase(eval_tool), GVAR_space_Variance, ' ', sim.strategy, sep = ""),
           col = GVAR_AF_level, shape = GVAR_AF_level) +
      theme(text = element_text(size = 16),
            axis.text = element_text(size = 16, color = "black"),
            plot.title = element_text(hjust = 0.5))
ggsave(paste(prefix, 'ctDNA.var.', eval_tool, '.pdf', sep = ''), p2, width = 8, height = 6)

# Mean scatter plot
p3 <- ggplot(data = df_mean, aes(x = x, y = y, color = AF_level)) +
      geom_point(aes(shape = AF_level), size = 3, alpha = 0.6) +
      geom_abline(aes(slope = 1, intercept = 0), color = 'black', linetype = 1) +
      geom_text(data = mean_anno_text, aes(x = -1.8, y = -1, label = label), color = 'black', parse = T, size = 7) +
      facet_grid(.~Company, scales = 'free') +
      labs(x = bquote(.(GVAR_LBX_HIGH) ~ log[10] ~ .(GVAR_Allele_Fraction_Mean)),
           y = bquote(.(GVAR_Simulation_space) ~ log[10] ~ .(GVAR_Allele_Fraction_Mean)),
           title = paste(tools::toTitleCase(eval_tool), " ", GVAR_mean, ' ', sim.strategy, sep = ""),
           col = GVAR_AF_level, shape = GVAR_AF_level) +
      theme(text = element_text(size = 16),
            axis.text = element_text(size = 16, color = "black"),
            plot.title = element_text(hjust = 0.5))
ggsave(paste(prefix, 'ctDNA.mean.', eval_tool, '.pdf', sep = ''), p3, width = 8, height = 6)
#p <- marrangeGrob(list(p1, p3, p2), nrow = 1, ncol = 1, top = NULL)
#ggsave(paste(prefix, '.stats.pdf', sep = ''), p, width = 12, height = 8)

