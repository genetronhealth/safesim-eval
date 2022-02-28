#!/usr/bin/env Rscript

library(ggplot2)
library(gridExtra)
args <- commandArgs(trailingOnly=TRUE)
f_gs <- args[1]
real_dir <- args[2]
eval_dir <- args[3]
prefix <- args[4]
eval_tool <- tail(strsplit(prefix, '/')[[1]], 1)

addline_format <- function(x,...){
  gsub('\\s','\n',x)
}

gs_df <- read.table(f_gs, header = F, sep = "\t")
colnames(gs_df) <- c('mut', 'value', 'Anno', 'CNT', 'AF_Class')
gs_df$Class <- rep('GS_AF', nrow(gs_df))
gs_df$variable <- rep('AF', nrow(gs_df))

df <- data.frame()
for (dir in c(eval_dir, real_dir)){
  for (f in list.files(dir, recursive = TRUE, pattern = '_simmut.info.txt')){
    filename <- paste(dir, f, sep = '/')
    tmpdf <- read.table(filename, header = F, sep = '\t', stringsAsFactors = FALSE)
    df <- rbind(df, tmpdf)
  }
}

variables <- rep(c('AF','DP'), c(nrow(df),nrow(df)))
fmt_df <- df[,c(1,2,4)]
colnames(fmt_df) <- c('mut', 'value', 'Class')
tmpdf <- df[,c(1,3,4)]
colnames(tmpdf) <- c('mut', 'value', 'Class')
tmpdf$value <- log(tmpdf$value, 10)
fmt_df <- rbind(fmt_df, tmpdf)
fmt_df$variable <- variables

myplot <- function(n){
  start <- start_list[n]
  end   <- start_list[n] + 47
  if (end > nrow(gs_df)){ end <- nrow(gs_df) }
  df1 <- subset(fmt_df, fmt_df$mut %in% gs_df$mut[start:end])
  df1$mut <- factor(df1$mut, levels = gs_df$mut[start:end])
  cnt <- gs_df[start:end, c(1,4)]
  cnt$variable <- rep('AF', nrow(cnt))
  info <- gs_df[start:end, c(1,1)]
  colnames(info) <- c("mut", 'pos')
  info$pos <- sapply(strsplit(info$mut, '-'), function(x){paste(x[1:2],collapse=':')})
  info$variable <- rep('DP', nrow(info))
  p <- ggplot() + 
       geom_point(data = df1, aes(x = mut, y = value, color = Class), 
                  size = 6, alpha = 0.7, position = position_jitterdodge(dodge.width = 0.5)) + 
       geom_point(data = gs_df[start:end,], aes(x = mut, y = value, color='GS_AF'), shape = 95, size = 6) + 
       facet_grid(factor(variable,levels=c('DP','AF'),labels=c('log10(DP)', 'AF'))~., scales = 'free_y') + 
       scale_color_manual(values = c('black', '#f8766d', '#00ba38')) + 
       labs(title = paste(eval_tool, "Simulation (WES)"), x = "Mutation", y = "Value") + 
       scale_x_discrete(labels = addline_format(gs_df$Anno[start:end])) + 
       theme(panel.spacing = unit(10, "lines"), text = element_text(size = 16), 
             axis.text = element_text(size = 16, color = "black"), 
             axis.text.x = element_text(angle = 45,vjust = 0.75, hjust = 1), 
             axis.title.x = element_text(vjust=0), axis.ticks.x = element_blank(), 
             plot.title = element_text(hjust = 0.5), panel.grid = element_blank()) + coord_cartesian(clip = "off") + 
       geom_text(data = cnt, aes(x = mut, label = CNT ,y = -Inf), vjust = 1 ,size = 7) + 
       geom_text(data = data.frame(tag = "CNT",variable = "AF"), aes(x = Inf, label = tag, y = -Inf), vjust = 1.1, size = 7) + 
       geom_text(data = info, aes(x = mut, label = pos ,y = -Inf), vjust = 0.2, hjust = 1.1, size = 6, angle = 90) + 
       geom_text(data = data.frame(tag = "MUT", variable = "DP"), aes(x = Inf, label = tag, y = -Inf), vjust = 1.9, size = 7) + 
       geom_vline(xintercept = seq(1.5, length(unique(df1$mut)) - 0.5, 1), lwd = 1, linetype = 'dashed', colour = 'gray')
  return(p)
}

start_list <- seq(1,nrow(gs_df), 48) # n mutations each plot

pl <- lapply(1:length(start_list), myplot)
ml <- marrangeGrob(pl, nrow = 1, ncol = 1, top = NULL)
ggsave(paste(eval_tool, '_WES.mut.pdf', sep = ''), ml, width = 32, height = 15)


# plot mse
f_real <- paste(real_dir, 'realdata.txt', sep = '/')
f_sim <- paste(eval_dir, 'simulation.txt', sep = '/')

realdf <- read.table(f_real, header = F, sep = "\t", row.names = 1)
simdf  <- read.table(f_sim, header = F ,sep = "\t" , row.names = 1)
df <- data.frame()
df_var <- data.frame()
df_mean <- data.frame()

for (mut in rownames(realdf)){
  x <- realdf[mut,][!is.na(realdf[mut,])]
  y <- simdf[mut,][!is.na(simdf[mut,])]
  if(mean(x)!=0 && mean(y)!=0){
    zx <- (x-mean(x))/sd(x)
    zy <- (y-mean(y))/sd(y)
    df <- rbind(df, data.frame(zx = zx, zy = zy, mut=mut))
    df_var <- rbind(df_var, data.frame(x = var(x), y = var(y), mut = mut))
    df_mean <- rbind(df_mean, data.frame(x = mean(x), y =mean(y), mut=mut))
  }
}

gs_df <- read.table(f_gs, header=F,sep = "\t")
af_level <- c()
for(s in df$mut){af_level <- append(af_level ,subset(gs_df, V1 == s)[1,5])}
df$AF_level <- af_level
df_l <- data.frame(x=sort(df[df$AF_level=="Low",1]),y=sort(df[df$AF_level=="Low",2]),AF_level=df[df$AF_level=="Low",4])
df_m <- data.frame(x=sort(df[df$AF_level=="Medium",1]),y=sort(df[df$AF_level=="Medium",2]),AF_level=df[df$AF_level=="Medium",4])
df_h <- data.frame(x=sort(df[df$AF_level=="High",1]),y=sort(df[df$AF_level=="High",2]),AF_level=df[df$AF_level=="High",4])
df_qqplot <- rbind(df_l,df_m,df_h)
af_level <- c()
for(s in df_mean$mut){af_level <- append(af_level ,subset(gs_df, V1 == s)[1,5])}
df_mean$AF_level <- af_level
af_level <- c()
for(s in df_var$mut){af_level <- append(af_level ,subset(gs_df, V1 == s)[1,5])}
df_var$AF_level <- af_level

df_qqplot$AF_level <- factor(df_qqplot$AF_level, levels = c('Low','Medium', 'High'))
df_mean$AF_level <- factor(df_mean$AF_level, levels = c('Low','Medium', 'High'))
df_var$AF_level <- factor(df_var$AF_level, levels = c('Low','Medium', 'High'))
df_mean$x <- log(df_mean$x, 10)
df_mean$y <- log(df_mean$y, 10)
df_var$x <- log(df_var$x, 10)
df_var$y <- log(df_var$y, 10)

mse <- function(data){
  colnames(data) <- c("x", "y")
  ybar = mean(data$y)
  ss_residuals <- sum((data$y - data$x) ^ 2)
  mse <- format(ss_residuals / nrow(data), digits = 4)
  return(mse)
}

mean_mse <- mse(df_mean)
var_mse <- mse(df_var)
anno_text <- data.frame(AF_level = c('Low', 'Medium', 'High'), label = sprintf("italic(y)==~italic(x)~','~italic(MSE)==%s",c(mse(df_l),mse(df_m),mse(df_h))))
anno_text$AF_level <- factor(anno_text$AF_level, levels = c('Low','Medium', 'High'))
p1 <- ggplot(data = df_qqplot, mapping = aes(x = x, y = y)) + 
      facet_grid(.~AF_level) + geom_point(size=3, alpha = 0.2, shape=4) + 
      geom_abline(aes(slope = 1, intercept = 0), color = 'red', linetype = 1) + 
      labs(x = "RealData Allele Fraction Z-Score" , y = "Simulation Allele Fraction Z-Score", title = paste(eval_tool, 'QQplot (WES)')) + 
      theme(text = element_text(size = 18), axis.text = element_text(size = 18, color = "black"), plot.title = element_text(hjust = 0.5)) + 
      geom_text(data = anno_text, mapping = aes(x = -0.6, y = 2, label = label), parse = T, size = 5)
ggsave(paste(prefix, '_WES.qqplot.pdf', sep = ''), p1, width = 12, height = 10)
p2 <- ggplot(data = df_var, mapping = aes(x = x, y = y, color = AF_level)) + 
      geom_point(size=3, alpha = 0.6) + 
      geom_abline(aes(slope = 1, intercept = 0), color = 'black', linetype = 1) + 
      labs(x = expression(paste("RealData ", log[10], "(Allele Fraction Variance)", sep='')) ,
           y = expression(paste("Simulation ", log[10], "(Allele Fraction Variance)", sep='')),
           title = paste(eval_tool, "Variance (WES)")) + 
      theme(text = element_text(size = 18), axis.text = element_text(size = 18, color = "black"), plot.title = element_text(hjust = 0.5)) + 
      annotate(geom = "text", x= -3.6, y= -1.5, label = paste("italic(y)==~italic(x)~','~italic(MSE)==", var_mse), size = 7, parse = T)
ggsave(paste(prefix, '_WES.var.pdf', sep = ''), p2, width = 12, height = 10)
p3 <- ggplot(data = df_mean, mapping = aes(x = x, y = y, color = AF_level)) + 
      geom_point(size=3, alpha = 0.6) + 
      geom_abline(aes(slope = 1, intercept = 0), color = 'black', linetype = 1) + 
      labs(x = expression(paste("RealData ", log[10], "(Allele Fraction Mean)", sep='')) , 
           y = expression(paste("Simulation ", log[10], "(Allele Fraction Mean)", sep = '')), 
           title = paste(eval_tool, "Mean (WES)")) + 
      theme(text = element_text(size = 18), axis.text = element_text(size = 18, color = "black"), plot.title = element_text(hjust = 0.5)) + 
      annotate(geom = "text", x = -1.2, y = -0.3, label = paste("italic(y)==~italic(x)~','~italic(MSE)==", mean_mse), size = 7, parse = T)
ggsave(paste(prefix, '_WES.mean.pdf', sep = ''), p3, width = 12, height = 10)
#p <- marrangeGrob(list(p1, p3, p2), nrow = 1, ncol = 1, top = NULL)
#ggsave(paste(eval_tool, '_WES.stats.pdf', sep = ''), p, width = 12, height = 10)

