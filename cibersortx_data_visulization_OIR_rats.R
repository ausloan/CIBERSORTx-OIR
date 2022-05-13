install.packages('corrplot')
install.packages('gmp')
install.packages('Rmpfr')
install.packages('PMCMRplus')
install.packages('pairwiseComparisons')
install.packages('ggbetweenstats')
install.packages('ggstatsplot')
install.packages('ggthemes')
install.packages('ggExtra')
install.packages('viridis')
library(RColorBrewer)
library(ggplot2)
library("limma")
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(viridis)
library(corrplot)
library(ggthemes)
library(ggplot2)
library(ggstatsplot)
library("ggstatsplot")
library("ggpubr")
library("ggExtra")

#cibersortX plots
filter <- read.table ("cibersortX/CIBERSORTx_Job43_Adjusted.txt",sep="\t",header = T)
filter$group = ifelse(as.numeric(str_sub(filter$Mixture,5,6))<4,"control","OIR")
table(filter$group)
filtered <- rbind(filter[filter$group=="control",],filter[filter$group=="OIR",])
filtered
dim(filtered)
tail(filtered)
filtered<-filtered[,c(1:23)]
write.table( filtered, "cibersortX/cibersortX_filtered.txt", sep="\t", quote = F, row.names = F, col.names = T )

#pheatmap
pp<-filtered
rownames(pp)<-pp$Mixture 
pp<-pp[,-1] 
pp <- t(scale(t(pp))) 
pp[pp < -2] <- -2
pp[pp > 2] <- 2
pp=t(pp) 
Group=c(rep("control",3),rep("OIR",3))   
names(Group)=colnames(pp)
Group=as.data.frame(Group)
pheatmap(pp, 
         annotation=Group, 
         cluster_cols =F,
         fontsize = 10,
         fontsize_row=14,
         fontsize_col=10, cluster_col = FALSE, border=FALSE)


#corpheatmap
cp<-filtered
cp
rownames(cp)<-cp$Mixture 
cp<-cp[,-1] 
class(cp)
pdf("corrHeatmap.pdf",height=13,width=13)    
corr<- cor(cp)
corrplot(corr,na.label ="NA",
         insig="blank",
         method = "color",
         tl.col="black",
         tl.srt=45,
         addCoef.col = "black",
         number.cex = 0.6,
         col=colorRampPalette(c("blue","white","red"))(50),)
dev.off()
cor(pp)
cor(cp)

#barplot
bk <- c(seq(0,0.2,by = 0.01),seq(0.21,0.85,by=0.01))
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
cibersort_barplot <- filtered %>%
  gather(key = Cell_type,value = Proportion,2:23)
colnames(cibersort_barplot) #[1] "Mixture"    "Cell_type"  "Proportion"
ggplot(cibersort_barplot,aes(Mixture,Proportion,fill = Cell_type)) + 
  geom_bar(position = "stack",stat = "identity") +
  labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") + theme_bw() +
  theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank()) + theme(legend.text=element_text(size=12)) +
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mypalette(23))

#barcharts
ggplot(cibersort_barplot,aes(Cell_type,Proportion,fill = Cell_type)) + 
  geom_boxplot(outlier.shape = 21,colour = "black") + theme_bw() + 
  labs(x = "Cell_Type", y = "Estimated Proportion") +
  theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank()) +
  scale_fill_manual(values = mypalette(23))  

ggplot(cibersort_barplot,aes(Cell_type,Proportion,fill = Cell_type)) + 
  geom_boxplot(outlier.shape = 21,colour = "black") + theme_bw() + 
  labs(x = "Cell_Type", y = "Estimated Proportion",main = "TME Cell composition") +
  theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank()) +
  scale_fill_manual(values = mypalette(23)) + theme_base() +theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 1))

cibersort_barplot
cibersort_barplot$Group = ifelse(as.numeric(str_sub(cibersort_barplot$Mixture,5,6))<4,"control","OIR")
cibersort_barplot$Group
ggplot(cibersort_barplot,aes(Cell_type,Proportion,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5, size=12))+
  scale_fill_manual(values = mypalette(22)[c(6,1)])+ stat_compare_means(aes(group = Group,label = ..p.signif..),method = "kruskal.test")

#barcharts with statistics
chayi<-filtered
chayi$group<-c(rep("control",3),rep("OIR",3))
colnames(chayi)
ggbetweenstats(data = chayi, 
               x = group, 
               y = Macrophages.M2,  
               plot.type = "violin", # type of plot ,"box", "violin", or "boxviolin"
               type = "parametric", # type of statistical test , p (parametric), np ( nonparametric), r(robust), bf (Bayes Factor).
               effsize.type = "unbiased", # type of effect size (unbiased = omega)
               partial = FALSE, # partial omega or omega?
               pairwise.comparisons = TRUE, # display results from pairwise comparisons
               pairwise.display = "significant", # display only significant pairwise comparisons
               pairwise.annotation = "p.value", # annotate the pairwise comparisons using p-values
               p.adjust.method = "fdr",   # adjust p-values for multiple tests using this method)
               title = "Macrophages.M2 expression",
               ylab = "log2 median-centered intensity",
               xlab = "group")

