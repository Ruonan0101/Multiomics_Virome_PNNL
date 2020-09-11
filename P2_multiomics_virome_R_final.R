#######fig 1.c
setwd("~/Desktop/Project/P2_virome_multiomics/1_mapping_transcript_to_VC_contigs/")
d<- read.csv("DNA_V_transcript_mapped.csv", header = T)
colnames(d)<-c('treatment',"count")
library(ggplot2)
source("http://goo.gl/UUyEzD")
library(ggpubr)
#outlierKD(d, count)
ggplot(d, aes(x=treatment, y=count, color=treatment,add = "jitter",palette = "jco"))+ 
  geom_boxplot()+
  geom_point() + 
  theme_bw()+
  theme(axis.title.x=element_blank())+
  theme(text = element_text(size=12))+
  ylab('Number of active DNA viral contigs')+ 
  stat_compare_means(method = "t.test",label.x=1)

#######fig 1.d
setwd("~/Desktop/Project/P2_virome_multiomics/1_mapping_transcript_to_VC_contigs/")
d<- read.csv("NumberOfNormalizedNormalizedNormalized_transcriptCount_DNA_viralContigs.csv", header = T)
head (d)
colnames(d)<-c('TranscriptCount',"treatment")
library(ggplot2)
source("http://goo.gl/UUyEzD")
library(ggpubr)
#outlierKD(d, count)
ggplot(d, aes(x=treatment, y=TranscriptCount, color=treatment,add = "jitter",palette = "jco"))+ 
  geom_boxplot()+
  geom_point() + 
  theme_bw()+
  theme(axis.title.x=element_blank())+
  theme(text = element_text(size=13))+
  ylab('Total number of dna viral transcripts')+ 
  stat_compare_means(method = "t.test",label.x=1)

####fig2.a,b
setwd('/Users/wuru978/Desktop/Project/P2_virome_multiomics/Manuscript_dir//')
library(ggplot2)
library(DESeq2)
#abund_table<-read.csv("SignifantTest_transcriptCountMappedDNAViralContigs_DESeq2.csv",row.names=1,check.names=FALSE)
abund_table<-read.csv("DNA_all_transcript_mapped_normalized_416_final_ForDeseq2.csv",row.names=1,check.names=FALSE)
abund_table<-t(abund_table)
head(abund_table)
grouping_info<-data.frame(row.names=rownames(abund_table),t(as.data.frame(strsplit(rownames(abund_table),"_"))))
head(grouping_info[2])
countData = round(as(abund_table, "matrix"), digits = 0)
head(countData)
countData<-(t(countData+1)) 
dds <- DESeqDataSetFromMatrix(countData, grouping_info, as.formula(~ X2))
data_deseq_test = DESeq(dds, test="Wald", fitType="local")
res = results(data_deseq_test, cooksCutoff = FALSE)
head(res)
mcols(res, use.names=TRUE)
res_tax = cbind(as.data.frame(res), as.matrix(countData[rownames(res), ]), OTU = rownames(res))
sig = 0.05
fold = 0
plot.point.size = 2
label=T
tax.display = NULL
tax.aggregate = "OTU"
res_tax_sig = subset(res_tax, padj < sig & fold < abs(log2FoldChange))
res_tax_sig <- res_tax_sig[order(res_tax_sig$padj),]
res_tax$Significant <- ifelse(rownames(res_tax) %in% rownames(res_tax_sig) , "Yes", "No")
res_tax$Significant[is.na(res_tax$Significant)] <- "No"
p1 <- ggplot(data = res_tax, aes(x = baseMean, y = log2FoldChange, color = Significant)) +
  geom_point(size = plot.point.size) +
  scale_x_log10() +
  scale_color_manual(values=c("black", "red")) +
  labs(x = "Mean abundance", y = "Log2 fold change (wet relative to dry)")+theme_bw()
if(label == T){
  if (!is.null(tax.display)){
    rlab <- data.frame(res_tax, Display = apply(res_tax[,c(tax.display, tax.aggregate)], 1, paste, collapse="; "))
  } else {
    rlab <- data.frame(res_tax, Display = res_tax[,tax.aggregate])
  }
  p1 <- p1 + geom_text(data = subset(rlab, Significant == "Yes"),aes(label = Display), size = 2.5, vjust = 1.5)
}
p1
  


res_tax_sig_abund = cbind(as.data.frame(countData[rownames(res_tax_sig), ]), OTU = rownames(res_tax_sig), padj = res_tax[rownames(res_tax_sig),"padj"]) 
write.table(res_tax_sig_abund, file = "KS_wet_dry_normalized_normalized_sig_DNA_V.txt", sep = "\t", row.names = TRUE)
data<-log(abund_table)
data<-as.data.frame(data)

df<-NULL
for(i in res_tax[rownames(res_tax_sig),"OTU"]){
  tmp<-data.frame(data[,i],grouping_info[2],rep(paste(i," padj = ",round(res_tax[i,"padj"],5),sep=""),dim(data)[1]))
  if(is.null(df)){df<-tmp} else { df<-rbind(df,tmp)} 
}
colnames(df)<-c("Value","treatment","Taxa")
df

p<-ggplot(df,aes(x=treatment,y=Value,colour=treatment))+ylab("Log of the normalized counts")+facet_wrap( ~Taxa, scales="free")
p<-p+geom_boxplot()+geom_jitter()+theme_bw()
p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p

###fig2.c
library(dplyr)
library(broom)
library(tidyr)
library(ggplot2)
setwd("~/Desktop/Project/P2_virome_multiomics/Manuscript_dir/mapping_VC1_dirr/")
d <-read.csv("VC1_all_sample_mapping_profile_1_forR_removeCtrl.csv", header = TRUE)
densities.qtiles <-
  d %>%
  rename(Sepal.Length = x, dens = y) %>%
  ungroup() %>%
  group_by(Species, Sepal.Length) %>% 
  summarise(q05 = quantile(dens, 0.025),
            q50 = quantile(dens, 0.5),
            q95 = quantile(dens, 0.975),
            mean_d=mean(dens)) 
#head(d)
#densities.qtiles

ggplot(densities.qtiles, aes(Sepal.Length, mean_d)) +
  facet_wrap(~Species, nrow = 3) +
  theme_bw()+
  theme(strip.text = element_text(size =10, face = "bold"))+
  geom_ribbon(aes(ymin = q05, ymax = q95, fill='0.05-0.95 quantile'), alpha = 0.5) +
  scale_fill_manual('', values='grey50')+
  geom_line(size = 1, aes(color = Species)) + 
  labs(y = "Reads coverage")

###fig 3.b
setwd('/Users/wuru978/Desktop/Project/P2_virome_multiomics/Manuscript_dir/Manuscript_revision_dir/')
abund_table<-read.csv("fig2b_R_plot.csv")
#abund_table<-read.csv("transcript_peptide_wetDry_Count_heatmap.csv")
library(reshape2)
library(wesanderson)
library(ggplot2)
colnames(abund_table)<-c("Gene","Cat","Sample","Value")
pal <- wes_palette("Zissou1", 50, type = "continuous")
tail(abund_table)
level_order <- c('dry', 'wet')
ggplot(abund_table, aes(x=factor(Sample,level=level_order),y=Gene))+ 
  geom_tile(aes(fill = Value)) + 
  scale_fill_gradientn(colours = pal, na.value = 'white')+
  theme(legend.position = "right",axis.ticks = element_blank(),axis.text.y = element_text(size=8),axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5,size=8),strip.text = element_text(size = 8))+
  facet_grid(.~Cat,space="free", scales="free")+
  labs(fill = "log of the normalized counts")+theme_bw()
###fig 3.c
library(ggtree)
library(ggplot2)
setwd('/Users/wuru978/Desktop/Project/P2_virome_multiomics/chaperonin_dir/final_dir/final_final_dir/final_final_final_dir/')
info <- read.csv("1_GroEL_refseq_KS_MarineV_uniqClass_final_PC_12345_ID.csv")
tree <- read.tree("1_GroEL_refseq_KS_MarineV_uniqClass_final_PC_12345_rename.tree")
colnames(info)<-c("id1","group",'id3')
tail(info)
p<-ggtree(tree, layout='circular') %<+% info + 
  geom_tippoint(aes(color=group),size=0.5)+ 
  theme(legend.position="right", legend.text = element_text(size=10))+
  guides(color = guide_legend(override.aes = list(size = 2)))
p
p1<- msaplot(p,'1_GroEL_refseq_KS_MarineV_uniqClass_final_PC_12345_rename_cut_conserved_final.fasta', window=c(1,6), width=0.2)+ 
  theme(legend.position="bottom", legend.text = element_text(size=10))+
  guides(color = guide_legend(override.aes = list(size = 2)))
p1

###fig 4
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")
library(ggtree)
setwd('/Users/wuru978/Desktop/Project/P2_virome_multiomics/3_RNA_viruses/ggtree_dir/final_dir/thirt_dry_dir/')
info <- read.csv("/Users/wuru978/Desktop/Project/P2_virome_multiomics/3_RNA_viruses/ggtree_dir/re-rooted_info.csv")
tree <- read.tree("/Users/wuru978/Desktop/Project/P2_virome_multiomics/3_RNA_viruses/ggtree_dir/final_dir/thirt_dry_dir/210_Remove_QDH91182_BBI93117_QDH88671_mafft_aligned_renamed_rerooted_newick.txt")
#new_cols <- c(HCMC='black', Hue='purple2', KH='skyblue2')
colnames(info)<-c("id","taxa",'acc')
head(info)
p<-ggtree(tree, layout='circular') %<+% info + 
  geom_tippoint(aes(color=taxa),size=0.5) + 
  geom_tiplab2(aes(label=taxa), align=T, linetype=NA, 
               size=1.5, offset=2.5, hjust=0.5)
p
heatmapData=read.csv("/Users/wuru978/Desktop/Project/P2_virome_multiomics/3_RNA_viruses/ggtree_dir/final_dir/thirt_dry_dir/coverage_210_average_2treatment.csv", row.names=1)
gheatmap(p, log(heatmapData), offset = 3.2, color=NULL, 
         colnames_position="top", 
         colnames_angle=90, colnames_offset_y = 1, 
         hjust=0, font.size=2,width=0.2) + ylab("Time")

####fig 5.a
library(ggplot2)
library(pheatmap)
setwd('/Users/wuru978/Desktop/Project/P2_virome_multiomics/Manuscript_dir/')
data<-read.delim("RNA_by_cat_final.txt",header=T, row.names="taxa")
pheatmap(data)
cal_z_score <- function(x){
  (log(x+1))
}

data_norm <- t(apply(data, 1, cal_z_score))
pheatmap(data_norm)
my_hclust_gene <- hclust(dist(data_norm), method = "complete")
library(dendextend)
as.dendrogram(my_hclust_gene) %>%
  plot(horiz = TRUE)
my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = 4)
my_gene_col <- data.frame(cluster = ifelse(test = my_gene_col ==1, yes = "cluster 1", no = "cluster 2"))
my_sample_col <- data.frame(sample = rep(c("dry", "wet"), c(6,6)))
row.names(my_sample_col) <- colnames(data_norm)
my_sample_col
pheatmap(data_norm, annotation_row = my_gene_col, annotation_col = my_sample_col)
pheatmap(data_norm,
         annotation_row = my_gene_col,
         annotation_col = my_sample_col,
         cutree_cols = 6,fontsize = 8, legend_labels = 'log of normalized RNA viral transcript counts')
abund_table<-read.csv("RNA_by_cat_final.csv",row.names=1,check.names=FALSE)

##fig 5.b
setwd("~/Desktop/Project/P2_virome_multiomics/Manuscript_dir/")
library(ggplot2)
library(DESeq2)
abund_table<-read.csv("RNA_by_cat_final.csv",row.names=1,check.names=FALSE)
abund_table<-t(abund_table)
head(abund_table)
grouping_info<-data.frame(row.names=rownames(abund_table),t(as.data.frame(strsplit(rownames(abund_table),"_"))))
head(grouping_info[2])
countData = round(as(abund_table, "matrix"), digits = 0)
countData<-(t(countData+1)) 
dds <- DESeqDataSetFromMatrix(countData, grouping_info, as.formula(~ X2))
head (dds)
data_deseq_test = DESeq(dds)
res = results(data_deseq_test, cooksCutoff = FALSE)
head(res)
mcols(res, use.names=TRUE)
res_tax = cbind(as.data.frame(res), as.matrix(countData[rownames(res), ]), OTU = rownames(res))
sig = 0.05
fold = 0
plot.point.size = 2
label=T
tax.display = NULL
tax.aggregate = "OTU"
res_tax_sig = subset(res_tax, padj < sig & fold < abs(log2FoldChange))
res_tax_sig <- res_tax_sig[order(res_tax_sig$padj),]
res_tax$Significant <- ifelse(rownames(res_tax) %in% rownames(res_tax_sig) , "Yes", "No")
res_tax$Significant[is.na(res_tax$Significant)] <- "No"
p1 <- ggplot(data = res_tax, aes(x = baseMean, y = log2FoldChange, color = Significant)) +
  geom_point(size = plot.point.size) +
  scale_x_log10() +
  scale_color_manual(values=c("black", "red")) +
  labs(x = "Mean abundance", y = "abundance of Rhizobiales OTUs")+theme_bw()
if(label == T){
  if (!is.null(tax.display)){
    rlab <- data.frame(res_tax, Display = apply(res_tax[,c(tax.display, tax.aggregate)], 1, paste, collapse="; "))
  } else {
    rlab <- data.frame(res_tax, Display = res_tax[,tax.aggregate])
  }
  p1 <- p1 + geom_text(data = subset(rlab), aes(label = Display), size = 3, vjust = 1.5)
}
p1

res_tax_sig_abund = cbind(as.data.frame(countData[rownames(res_tax_sig), ]), OTU = rownames(res_tax_sig), padj = res_tax[rownames(res_tax_sig),"padj"]) 
write.table(res_tax_sig_abund, file = "RNA_virus_by_cat.txt", sep = "\t", row.names = TRUE)
data<-abund_table
data<-as.data.frame(data)
df<-NULL
for(i in res_tax[rownames(res_tax_sig),"OTU"]){
  tmp<-data.frame(data[,i],grouping_info[2],rep(paste(i," padj = ",round(res_tax[i,"padj"],5),sep=""),dim(data)[1]))
  if(is.null(df)){df<-tmp} else { df<-rbind(df,tmp)} 
}
colnames(df)<-c("Value","Type","Taxa")
df
p<-ggplot(df,aes(x=Type,y=Value,colour=Type))+ylab("abundance of Rhizobiales OTUs")+facet_wrap( ~Taxa, scales="free")
p<-p+geom_boxplot()+geom_jitter()+theme_bw()
p
p<-p+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
p

###extended data file 
####18s
setwd('/Users/wuru978/Desktop/Project/P2_virome_multiomics/Manuscript_dir/')
data<-read.delim("18s_composition_normalizedByReads_addingFungi.txt",header=T, row.names="gene")
pheatmap(data)
cal_z_score <- function(x){
  (log(x+1))
}

data_norm <- t(apply(data, 1, cal_z_score))
pheatmap(data_norm)
my_hclust_gene <- hclust(dist(data_norm), method = "complete")
library(dendextend)
as.dendrogram(my_hclust_gene) %>%
  plot(horiz = TRUE)
my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = 4)
my_gene_col <- data.frame(cluster = ifelse(test = my_gene_col == 1, yes = "cluster 1", no = "cluster 2"))
my_sample_col <- data.frame(sample = rep(c("dry", "wet"), c(6,6)))
row.names(my_sample_col) <- colnames(data_norm)
my_sample_col
pheatmap(data_norm, annotation_row = my_gene_col, annotation_col = my_sample_col)
pheatmap(data_norm,
         annotation_row = my_gene_col,
         annotation_col = my_sample_col,
         cutree_cols = 6,fontsize = 8, legend_labels = 'log of normalized 18s rRNA counts')

####16s 
setwd('/Users/wuru978/Desktop/Project/P2_virome_multiomics/Manuscript_dir/')
library(ggplot2)
library(DESeq2)
#abund_table<-read.csv("/Users/wuru978/Desktop/Project/P2_virome_multiomics/",row.names=1,check.names=FALSE)
abund_table<-read.csv("/Users/wuru978/Desktop/Project/P2_virome_multiomics/16s_dir/OTU_16s_forR.csv",row.names=1,check.names=FALSE)
abund_table<-t(abund_table)
head(abund_table)
grouping_info<-data.frame(row.names=rownames(abund_table),t(as.data.frame(strsplit(rownames(abund_table),"_"))))
head(grouping_info[2])
countData = round(as(abund_table, "matrix"), digits = 0)
head(countData)
countData<-(t(countData+1)) 
dds <- DESeqDataSetFromMatrix(countData, grouping_info, as.formula(~ X2))
data_deseq_test = DESeq(dds, test="Wald", fitType="local")
res = results(data_deseq_test, cooksCutoff = FALSE)
head(res)
mcols(res, use.names=TRUE)
res_tax = cbind(as.data.frame(res), as.matrix(countData[rownames(res), ]), OTU = rownames(res))
sig = 0.05
fold = 0
plot.point.size = 2
label=T
tax.display = NULL
tax.aggregate = "OTU"
res_tax_sig = subset(res_tax, padj < sig & fold < abs(log2FoldChange))
res_tax_sig <- res_tax_sig[order(res_tax_sig$padj),]
res_tax$Significant <- ifelse(rownames(res_tax) %in% rownames(res_tax_sig) , "Yes", "No")
res_tax$Significant[is.na(res_tax$Significant)] <- "No"
p1 <- ggplot(data = res_tax, aes(x = baseMean, y = log2FoldChange, color = Significant)) +
  geom_point(size = plot.point.size) +
  scale_x_log10() +
  scale_color_manual(values=c("black", "red")) +
  labs(x = "Mean abundance", y = "Log2 fold change (wet relative to dry)")+theme_bw()
if(label == T){
  if (!is.null(tax.display)){
    rlab <- data.frame(res_tax, Display = apply(res_tax[,c(tax.display, tax.aggregate)], 1, paste, collapse="; "))
  } else {
    rlab <- data.frame(res_tax, Display = res_tax[,tax.aggregate])
  }
  p1 <- p1 + geom_text(data = subset(rlab, Significant == "Yes"),aes(label = Display), size = 2, vjust = 1.5)
}
p1

