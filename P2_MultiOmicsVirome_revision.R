###Fig. 1a host virus pairing
setwd('/Users/wuru978/Desktop/2021_five_Manuscripts/PNNL_Manuscript_Revision/P2_multiomics_MS_CommBio_resubmission/')
###host virus pairing 
library(ggalluvial)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
d<-read.csv("fig1a_host_treatment_VC_ntw.csv",check.names=FALSE)
#d<-read.csv("fig1a_host_treatment_VC_ntw_wet_dry_noWetDry.csv",check.names=FALSE)
d
p1<-ggplot(data = d,
           aes(axis1 = Host, axis2=VC)) +
  scale_x_discrete(limits = c("Host","Viral contig"), expand = c(.2, .05)) +
  geom_alluvium(aes(fill = Host), alpha=0.9) +
  geom_stratum(width=1/12, alpha=.5, color='white') +
  #geom_text_repel(stat = "stratum", size=1, aes(label = after_stat(stratum))) +
  facet_wrap(~Treatment, scale='free')
mycolors<-c("#8DD3C7","#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69","#BC80BD")
fig.1a<-p1+scale_fill_manual(values=mycolors)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank())+
  theme(legend.position = "left", legend.title = element_blank(), legend.text = element_text(size = 7),
        legend.key.size = unit(0.4, "cm"),
        axis.text.x = element_text(size = 8, color='black'), 
        axis.text.y = element_text(size = 8, color='black'),
        panel.spacing=unit(1, "lines"))+
  guides(fill=guide_legend(ncol=1))+
  ylab('NO. of virus-host pairings')+
  theme(axis.title.y = element_text(size=9))

fig.1a

###Fig. 1b wet_dry_VC_Count_Venn
setwd('/Users/wuru978/Desktop/2021_five_Manuscripts/PNNL_Manuscript_Revision/P2_multiomics_MS_CommBio_resubmission/')
install.packages("remotes")
remotes::install_github("yanlinlin82/ggvenn")
library(ggvenn)
library(grid)
d<-read.csv("VC_wet_dryForVenn.csv",check.names=FALSE)
#d <- tibble(value   = c(1, 2, 3, 5, 6, 7, 8, 9, 10, 12, 13),
#            `Set 1` = c(T, F, T, T, F, T, F, T, F,  F,  F),
#            `Set 2` = c(T, F, F, T, F, F, F, T, F,  F,  T),
#            `Set 3` = c(T, T, F, F, F, F, T, T, F,  F,  F),
#            `Set 4` = c(F, F, F, F, T, T, F, F, T,  T,  F))
#d

p2 <- ggvenn(d, c("wet", "dry"), fill_color = c("#00BFC4", "#F8766D"), stroke_size = 0.2, set_name_size = 5) + theme_void() + coord_fixed()
p2

###Fig1.cd barchat with sig
library(ggplot2)
library(ggsignif)
library(tidyverse)
library(rstatix)
library(ggpubr)
setwd("~/Desktop/Project/P2_virome_multiomics/1_mapping_transcript_to_VC_contigs/")
d <-read.csv("Fig1c_fecet.csv", header = TRUE)
d
stat.test <- d %>%
  group_by(Variable) %>%
  t_test(Value ~ Treatment) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test
stat.test <- stat.test %>% add_xy_position(x = "Treatment")
stat.test
p.adj<- scientific(stat.test$p.adj, digits = 3)
p.adj
bxp<- ggplot(d,aes(x=Treatment,y=Value))+
  geom_boxplot(fill=NA,alpha=0.9, aes(color=Treatment))+
  geom_jitter(aes(shape=Site, color=Treatment),width=0.2, alpha=0.7)+
  facet_wrap(~Variable,scales = "free", strip.position='left')+
  #  stat_pvalue_manual(stat.test, size=3,bracket.nudge.y = -200, tip.length = 0.01)+
  scale_y_continuous(expand = expansion(mult = c(0.1,0.1)))+
  ylab(NULL)+
  theme_bw()+
  theme(axis.text = element_text(size=10),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.placement='outside',
        axis.title.x = element_blank(),
        axis.text.y=element_text(size = 8),
        panel.spacing=unit(1, "lines"))

bxp



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

library(ggforce)
library(ggrepel)
library(ggplot2)
library(ggforce)
library(grid)
library(ggforce)

# define facet_zoom2 function to use FacetZoom2 instead of FacetZoom
# (everything else is the same as facet_zoom)
facet_zoom2 <- function(x, y, xy, zoom.data, xlim = NULL, ylim = NULL, 
                        split = FALSE, horizontal = TRUE, zoom.size = 2, 
                        show.area = TRUE, shrink = TRUE) {
  x <- if (missing(x)) if (missing(xy)) NULL else lazyeval::lazy(xy) else lazyeval::lazy(x)
  y <- if (missing(y)) if (missing(xy)) NULL else lazyeval::lazy(xy) else lazyeval::lazy(y)
  zoom.data <- if (missing(zoom.data)) NULL else lazyeval::lazy(zoom.data)
  if (is.null(x) && is.null(y) && is.null(xlim) && is.null(ylim)) {
    stop("Either x- or y-zoom must be given", call. = FALSE)
  }
  if (!is.null(xlim)) x <- NULL
  if (!is.null(ylim)) y <- NULL
  ggproto(NULL, FacetZoom2,
          shrink = shrink,
          params = list(
            x = x, y = y, xlim = xlim, ylim = ylim, split = split, zoom.data = zoom.data,
            zoom.size = zoom.size, show.area = show.area,
            horizontal = horizontal
          )
  )
}

# define FacetZoom as a ggproto object that inherits from FacetZoom,
# with a modified draw_panels function. the compute_layout function references
# the version currently on GH, which is slightly different from the CRAN
# package version.
FacetZoom2 <- ggproto(
  "FacetZoom2",
  ggforce::FacetZoom,
  
  compute_layout = function(data, params) {
    layout <- rbind( # has both x & y dimension
      data.frame(name = 'orig', SCALE_X = 1L, SCALE_Y = 1L),
      data.frame(name = 'x', SCALE_X = 2L, SCALE_Y = 1L),
      data.frame(name = 'y', SCALE_X = 1L, SCALE_Y = 2L),
      data.frame(name = 'full', SCALE_X = 2L, SCALE_Y = 2L),
      data.frame(name = 'orig_true', SCALE_X = 1L, SCALE_Y = 1L),
      data.frame(name = 'zoom_true', SCALE_X = 1L, SCALE_Y = 1L)
    )
    if (is.null(params$y) && is.null(params$ylim)) { # no y dimension
      layout <- layout[c(1,2, 5:6),]
    } else if (is.null(params$x) && is.null(params$xlim)) { # no x dimension
      layout <- layout[c(1,3, 5:6),]
    }
    layout$PANEL <- seq_len(nrow(layout))
    layout
  },
  
  draw_panels = function(panels, layout, x_scales, y_scales, ranges, coord,
                         data, theme, params) {
    
    if (is.null(params$x) && is.null(params$xlim)) {
      params$horizontal <- TRUE
    } else if (is.null(params$y) && is.null(params$ylim)) {
      params$horizontal <- FALSE
    }
    if (is.null(theme[['zoom']])) {
      theme$zoom <- theme$strip.background
    }
    if (is.null(theme$zoom.x)) {
      theme$zoom.x <- theme$zoom
    }
    if (is.null(theme$zoom.y)) {
      theme$zoom.y <- theme$zoom
    }
    axes <- render_axes(ranges, ranges, coord, theme, FALSE)
    panelGrobs <- ggforce:::create_panels(panels, axes$x, axes$y)
    panelGrobs <- panelGrobs[seq_len(length(panelGrobs) - 2)]
    if ('full' %in% layout$name && !params$split) {
      panelGrobs <- panelGrobs[c(1, 4)]
    }
    
    # changed coordinates in indicator / lines to zoom from 
    # the opposite horizontal direction
    if ('y' %in% layout$name) {
      if (!inherits(theme$zoom.y, 'element_blank')) {
        zoom_prop <- scales::rescale(
          y_scales[[2]]$dimension(ggforce:::expansion(y_scales[[2]])),
          from = y_scales[[1]]$dimension(ggforce:::expansion(y_scales[[1]])))
        indicator <- polygonGrob(
          x = c(0, 0, 1, 1), # was x = c(1, 1, 0, 0), 
          y = c(zoom_prop, 1, 0), 
          gp = gpar(col = NA, fill = alpha(theme$zoom.y$fill, 0.5)))
        lines <- segmentsGrob(
          x0 = c(1, 1), x1 = c(0, 0), # was x0 = c(0, 0), x1 = c(1, 1)
          y0 = c(0, 1), y1 = zoom_prop,
          gp = gpar(col = theme$zoom.y$colour,
                    lty = theme$zoom.y$linetype,
                    lwd = theme$zoom.y$size,
                    lineend = 'round'))
        indicator_h <- grobTree(indicator, lines)
      } else {
        indicator_h <- zeroGrob()
      }
    }
    
    if ('x' %in% layout$name) {
      if (!inherits(theme$zoom.x, 'element_blank')) {
        zoom_prop <- scales::rescale(x_scales[[2]]$dimension(ggforce:::expansion(x_scales[[2]])),
                                     from = x_scales[[1]]$dimension(ggforce:::expansion(x_scales[[1]])))
        indicator <- polygonGrob(c(zoom_prop, 1, 0), c(1, 1, 0, 0), 
                                 gp = gpar(col = NA, fill = alpha(theme$zoom.x$fill, 0.5)))
        lines <- segmentsGrob(x0 = c(0, 1), y0 = c(0, 0), x1 = zoom_prop, y1 = c(1, 1), 
                              gp = gpar(col = theme$zoom.x$colour,
                                        lty = theme$zoom.x$linetype,
                                        lwd = theme$zoom.x$size,
                                        lineend = 'round'))
        indicator_v <- grobTree(indicator, lines)
      } else {
        indicator_v <- zeroGrob()
      }
    }
    
    if ('full' %in% layout$name && params$split) {
      space.x <- theme$panel.spacing.x
      if (is.null(space.x)) space.x <- theme$panel.spacing
      space.x <- unit(5 * as.numeric(convertUnit(space.x, 'cm')), 'cm')
      space.y <- theme$panel.spacing.y
      if (is.null(space.y)) space.y <- theme$panel.spacing
      space.y <- unit(5 * as.numeric(convertUnit(space.y, 'cm')), 'cm')
      
      # change horizontal order of panels from [zoom, original] to [original, zoom]
      # final <- gtable::gtable_add_cols(panelGrobs[[3]], space.x)
      # final <- cbind(final, panelGrobs[[1]], size = 'first')
      # final_tmp <- gtable::gtable_add_cols(panelGrobs[[4]], space.x)
      # final_tmp <- cbind(final_tmp, panelGrobs[[2]], size = 'first')
      final <- gtable::gtable_add_cols(panelGrobs[[1]], space.x)
      final <- cbind(final, panelGrobs[[3]], size = 'first')
      final_tmp <- gtable::gtable_add_cols(panelGrobs[[2]], space.x)
      final_tmp <- cbind(final_tmp, panelGrobs[[4]], size = 'first')
      
      final <- gtable::gtable_add_rows(final, space.y)
      final <- rbind(final, final_tmp, size = 'first')
      final <- gtable::gtable_add_grob(final, list(indicator_h, indicator_h),
                                       c(2, 6), 3, c(2, 6), 5,
                                       z = -Inf, name = "zoom-indicator")
      final <- gtable::gtable_add_grob(final, list(indicator_v, indicator_v), 
                                       3, c(2, 6), 5, 
                                       z = -Inf, name = "zoom-indicator")
      heights <- unit.c(
        unit(max_height(list(axes$x[[1]]$top, axes$x[[3]]$top)), 'cm'),
        unit(1, 'null'),
        unit(max_height(list(axes$x[[1]]$bottom, axes$x[[3]]$bottom)), 'cm'),
        space.y,
        unit(max_height(list(axes$x[[2]]$top, axes$x[[4]]$top)), 'cm'),
        unit(params$zoom.size, 'null'),
        unit(max_height(list(axes$x[[2]]$bottom, axes$x[[4]]$bottom)), 'cm')
      )
      
      # swop panel width specifications according to the new horizontal order
      widths <- unit.c(
        # unit(max_width(list(axes$y[[3]]$left, axes$y[[4]]$left)), 'cm'),
        # unit(params$zoom.size, 'null'),
        # unit(max_height(list(axes$y[[3]]$right, axes$y[[4]]$right)), 'cm'),
        # space.x,
        # unit(max_width(list(axes$y[[1]]$left, axes$y[[2]]$left)), 'cm'),
        # unit(1, 'null'),
        # unit(max_height(list(axes$y[[1]]$right, axes$y[[2]]$right)), 'cm')        
        unit(max_width(list(axes$y[[1]]$left, axes$y[[2]]$left)), 'cm'),
        unit(1, 'null'),
        unit(max_height(list(axes$y[[1]]$right, axes$y[[2]]$right)), 'cm'),
        space.x,
        unit(max_width(list(axes$y[[3]]$left, axes$y[[4]]$left)), 'cm'),
        unit(params$zoom.size, 'null'),
        unit(max_height(list(axes$y[[3]]$right, axes$y[[4]]$right)), 'cm')
        
      )
      final$heights <- heights
      final$widths <- widths
    } else {
      if (params$horizontal) {
        space <- theme$panel.spacing.x
        if (is.null(space)) space <- theme$panel.spacing
        space <- unit(5 * as.numeric(convertUnit(space, 'cm')), 'cm')
        heights <- unit.c(
          unit(max_height(list(axes$x[[1]]$top, axes$x[[2]]$top)), 'cm'),
          unit(1, 'null'),
          unit(max_height(list(axes$x[[1]]$bottom, axes$x[[2]]$bottom)), 'cm')
        )
        
        # change horizontal order of panels from [zoom, original] to [original, zoom]
        # first <- gtable::gtable_add_cols(panelGrobs[[2]], space)
        # first <- cbind(final, panelGrobs[[1]], size = 'first')
        final <- gtable::gtable_add_cols(panelGrobs[[1]], space) 
        final <- cbind(final, panelGrobs[[2]], size = "first") 
        
        final$heights <- heights
        
        # swop panel width specifications according to the new horizontal order
        # unit(c(params$zoom.size, 1), 'null')
        final$widths[panel_cols(final)$l] <- unit(c(1, params$zoom.size), 'null') 
        
        final <- gtable::gtable_add_grob(final, indicator_h, 2, 3, 2, 5, 
                                         z = -Inf, name = "zoom-indicator")
      } else {
        space <- theme$panel.spacing.y
        if (is.null(space)) space <- theme$panel.spacing
        space <- unit(5 * as.numeric(convertUnit(space, 'cm')), 'cm')
        widths <- unit.c(
          unit(max_width(list(axes$y[[1]]$left, axes$y[[2]]$left)), 'cm'),
          unit(1, 'null'),
          unit(max_height(list(axes$y[[1]]$right, axes$y[[2]]$right)), 'cm')
        )
        final <- gtable::gtable_add_rows(panelGrobs[[1]], space)
        final <- rbind(final, panelGrobs[[2]], size = 'first')
        final$widths <- widths
        final$heights[panel_rows(final)$t] <- unit(c(1, params$zoom.size), 'null')
        final <- gtable::gtable_add_grob(final, indicator_v, 3, 2, 5, 
                                         z = -Inf, name = "zoom-indicator")
      }
    }
    final
  }
)


fig2a <- ggplot(data = res_tax, aes(y = baseMean, x = log2FoldChange, color = Significant)) +
  geom_jitter(position = position_jitter(width = 0.5, height=0.5), size=1) +
  scale_color_manual(values=c("black", "red")) +
  labs(y = "Mean transcript abundance", x = "Log2 fold change (wet relative to dry)")+theme_bw()+
  theme(axis.title = element_text(size=11), axis.text=element_text(size=8), 
        legend.text = element_text(size=8), legend.title = element_text(size=10))+
  facet_zoom2(ylim = c(0, 22), zoom.size = 1)
#+scale_y_log10() 
fig2a

res_tax_sig_abund
res_tax_sig_abund = cbind(as.data.frame(countData[rownames(res_tax_sig), ]), OTU = rownames(res_tax_sig), padj = res_tax[rownames(res_tax_sig),"padj"]) 
write.table(res_tax_sig_abund, file = "KS_wet_dry_normalized_normalized_sig_DNA_V.txt", sep = "\t", row.names = TRUE)
data<-log(abund_table)
data<-as.data.frame(data)

df<-NULL
for(i in res_tax[rownames(res_tax_sig),"OTU"]){
  tmp<-data.frame(data[,i],grouping_info[2],grouping_info[1],rep(paste(i," padj = ",scientific(round(res_tax[i,"padj"],5),sep="")),dim(data)[1]))
  if(is.null(df)){df<-tmp} else { df<-rbind(df,tmp)} 
}
colnames(df)<-c("Value","Treatment",'Site',"Taxa")
df

p<-ggplot(df,aes(x=Treatment,y=Value,colour=Treatment))+ylab("Normalized transcript counts")+facet_wrap( ~Taxa, scales="free")
p<-p+geom_boxplot(alpha=0.9,outlier.shape=NA)+geom_jitter(aes(shape=Site),alpha=0.7)+theme_bw()
p<-p+theme(axis.text.x=element_blank(), axis.title.x=element_blank(), strip.text = element_text(size=9),panel.spacing=unit(1, "lines"))+
  xlab(NULL)
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
  ungroup() %>%
  group_by(Treatment, Length) %>% 
  summarise(q05 = quantile(dens, 0.025),
            q50 = quantile(dens, 0.5),
            q95 = quantile(dens, 0.975),
            mean_d=mean(dens)) 
head(d)
#densities.qtiles

ggplot(densities.qtiles, aes(Length, mean_d)) +
  facet_wrap(~Treatment, ncol=1,strip.position = c("right")) +
  theme_bw()+
  theme(strip.text = element_text(size =10, face = "bold"))+
  geom_ribbon(aes(ymin = q05, ymax = q95, fill='0.05-0.95 quantile'), alpha = 0.5) +
  scale_fill_manual('', values='grey50')+
  geom_line(size = 0.8, aes(color = Treatment)) + 
  labs(y = "Reads coverage")+
  scale_x_continuous(breaks= seq(0,8915,by=400))+
  theme(panel.background = element_blank(), panel.grid.minor = element_blank(),legend.key.size = unit(0.4, "cm"), axis.text = element_text(size = 6))+
  xlab('Position (Length: 8915 bp)')
#+facet_zoom(xlim = c(1527,2700),zoom.size = 0.5)
#  facet_zoom(xlim = c(1250,2750),zoom.size = 1)

###Fig3.ab, add noncoding regions) 
setwd('/Users/wuru978/Desktop/2021_five_Manuscripts/PNNL_Manuscript_Revision/P2_multiomics_MS_CommBio_resubmission/')
#d <-read.csv("CodingNoncoding_bar_errBar_fig3a.csv", header=TRUE)
d <-read.csv("codingNonCoding_absuluteABD.csv", header=TRUE)
library(ggplot2)
library(Rmisc) 
library(plyr)
library(ggsignif)
library(tidyverse)
library(rstatix)
library(ggpubr)

stat.test <- d %>%
  group_by(Part, Site) %>%
  t_test(Percentage ~ Treatment) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

stat.test <- stat.test %>% add_xy_position(x = "Treatment",fun = "mean_sd")
stat.test

fig3a<-ggbarplot(d, x = "Treatment", y = "Percentage", fill = 'Treatment', 
                 add = c("mean_sd", "jitter"),
                 facet = c("Part", "Site")) +
  stat_pvalue_manual(stat.test, hide.ns = TRUE)+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))+
  theme(legend.position = 'none', axis.title.x = element_blank(), 
        axis.title=element_text(size=11), axis.text = element_text(size=10), 
        strip.text = element_text(size=9))+
  ylab('Percent of transcript abundance (%)')

fig3a


setwd('/Users/wuru978/Desktop/2021_five_Manuscripts/PNNL_Manuscript_Revision/P2_multiomics_MS_CommBio_resubmission/AMG_doubleCheck_dir/')
abund_table<-read.csv("test.csv")

library(reshape2)
library(wesanderson)
library(ggplot2)
colnames(abund_table)<-c("Gene","Treatment","Value")
pal <- wes_palette("Zissou1", 50, type = "continuous")
tail(abund_table)
fig3b<-ggplot(abund_table, aes(x=Treatment,y=Gene))+ 
  geom_tile(aes(fill = log(Value))) + 
  scale_fill_gradientn(colours = pal, na.value = 'white')+
  theme_bw()+
  theme(legend.position = "left",
        axis.text.y = element_text(size=8),
        axis.ticks = element_blank(),
        axis.text.x = element_text(hjust = 1, vjust=0.5,size=8),
        strip.text = element_text(size = 8), 
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size=8))
fig3b

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
  tmp<-data.frame(data[,i],grouping_info[2],grouping_info[1], rep(paste(i," padj = ",scientific(round(res_tax[i,"padj"],5)),sep=""),dim(data)[1]))
  if(is.null(df)){df<-tmp} else { df<-rbind(df,tmp)} 
}
colnames(df)<-c("Value","Treatment","Site","Taxa")
df
p<-ggplot(df,aes(x=Treatment,y=Value, colour=Treatment))+
  ylab("Normalized reads coverage of RNA viral contig")+
  facet_wrap( ~Taxa, scales="free", ncol=1)+
  geom_boxplot(alpha=0.9,outlier.shape=NA)+geom_jitter(aes(shape=Site),alpha=0.7)+
  theme_bw()+
  theme(axis.title.x = element_blank(), axis.text.x=element_blank(),strip.text = element_text(size=9),panel.spacing=unit(1, "lines"))+
  xlab(NULL)
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
