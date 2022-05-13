install.packages('corrplot')
install.packages('gmp')
install.packages('Rmpfr')
install.packages('PMCMRplus')
install.packages('pairwiseComparisons')
install.packages('ggbetweenstats')
install.packages('ggstatsplot')
install.packages('ggthemes')
install.packages('ggExtra')
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
TPMs <- read.table ("data/hisat2_sw/raw_genename.txt",sep="\t",header = T)
dim(TPMs)
tail(sort(table(TPMs$Gene))) 
TPMs$mean<-apply(TPMs[,-1],1,median)
TPMs=TPMs[order(TPMs$Gene,TPMs$mean,decreasing = T),] 
TPMs=TPMs[!duplicated(TPMs$Gene),]
rownames(TPMs)<-TPMs$Gene
TPMs<-TPMs[,-1]
dim(TPMs) 
colnames(TPMs) 
TPMs<-TPMs[,-7]
TPMs[1:4,1:4]
TPMs<-as.data.frame(TPMs)
dim(TPMs)
save(TPMs,file = "cibersortX/TPMs.Rda")
load(file = "cibersortX/TPMs.Rda")
#data for cibersortX
rt=TPMs
rt=as.matrix(rt)
dimnames=list(rownames(rt),colnames(rt))
data=matrix(as.numeric(as.matrix(rt)),nrow=nrow(rt),dimnames=dimnames)
data=avereps(data) 
data=data[rowMeans(data)>0,] 
data=rbind(ID=colnames(data),data) 
write.table(data,file="cibersortX/data_for_cibersortX_1.txt",sep="\t",quote=F,col.names=F)  

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
install.packages('viridis')

pheatmap(pp, 
         annotation=Group, 
         color = colorRampPalette(c(plasma(3)))(40),
         cluster_cols =F,
         fontsize = 8,
         fontsize_row=7,
         fontsize_col=6)

###14.corpheatmap
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

###15.barplot
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
###16.barcharts
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
###17. barcharts with statistics
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

ggbetweenstats(data = chayi, 
               x = group, 
               y = Monocytes,  
               plot.type = "violin", # type of plot ,"box", "violin", or "boxviolin"
               type = "parametric", # type of statistical test , p (parametric), np ( nonparametric), r(robust), bf (Bayes Factor).
               effsize.type = "unbiased", # type of effect size (unbiased = omega)
               partial = FALSE, # partial omega or omega?
               pairwise.comparisons = TRUE, # display results from pairwise comparisons
               pairwise.display = "significant", # display only significant pairwise comparisons
               pairwise.annotation = "p.value", # annotate the pairwise comparisons using p-values
               p.adjust.method = "fdr",   # adjust p-values for multiple tests using this method)
               title = "Monocytes expression",
               ylab = "log2 median-centered intensity",
               xlab = "group")

#single gene's relation to all 22 immune cells
load(file = "cibersortX/TPMs.Rda")
TPMs <- as.data.frame(TPMs)
cor_shiyan1 = TPMs[,c(4:6)] 
dim(cor_shiyan1) #19641   471
library(stringr)
cor_shiyan1 <- as.data.frame(t(cor_shiyan1))
pre_filter<-read.table("cibersortX/CIBERSORTx_Job43_Adjusted.txt",sep="\t",header = T) #这个是免疫细胞表达量
dim(pre_filter) #[1] 514  26   
table(pre_filter$Mixture %in% rownames(cor_shiyan1)) #FALSE  65   TRUE  449
pre_filter <- pre_filter[pre_filter$Mixture %in% rownames(cor_shiyan1),] 
cor_shiyan1$Mixture <- rownames(cor_shiyan1) 
data<-merge(cor_shiyan1,pre_filter,by="Mixture")
rownames(data)<-data[,1]
data <- data[,-1] 
data
cor_data<- data.frame(matrix(NA,nrow(data),2))
dim(cor_data) #[1] 19666     2
rownames(cor_data)<-rownames(data)
cor_data$X1<-log2(data$`Ccl2`+1)
cor_data$X2<-data$Macrophages.M2  
colnames(cor_data)<-c("Ccl2","Macrophages.M2")
str(cor_data)
ggscatterstats(data = cor_data, 
               x = "Ccl2", 
               y = "Macrophages.M2",
               centrality.para = "mean", 
               margins = "both",
               xfill = "red", 
               yfill = "blue", 
               marginal.type = "densigram",
               title = "Relationship")