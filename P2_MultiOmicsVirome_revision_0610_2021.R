###setwd to where the Source data you downloaded
###Fig. 1a host virus pairing
setwd('/Users/wuru978/Desktop/2021_five_Manuscripts/PNNL_Manuscript_Revision/P2_multiomics_MS_CommBio_resubmission/')
###host virus pairing 
library(ggalluvial)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
d<-read.csv("Fig1a_data.csv",check.names=FALSE)
p1<-ggplot(data = d,
           aes(axis1 = Host, axis2=VC)) +
  scale_x_discrete(limits = c("Host","Viral contig"), expand = c(.2, .05)) +
  geom_alluvium(aes(fill = Host), alpha=0.9) +
  geom_stratum(width=1/12, alpha=.5, color='white') +
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
install.packages("remotes")
remotes::install_github("yanlinlin82/ggvenn")
library(ggvenn)
library(grid)
d<-read.csv("Fig1b_data.csv",check.names=FALSE)
fig1b <- ggvenn(d, c("wet", "dry"), fill_color = c("#00BFC4", "#F8766D"), stroke_size = 0.2, set_name_size = 5) + theme_void() + coord_fixed()
fig1b

###Fig1.cd barchat with sig
library(ggplot2)
library(ggsignif)
library(tidyverse)
library(rstatix)
library(ggpubr)
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
fig1cd<- ggplot(d,aes(x=Treatment,y=Value))+
  geom_boxplot(fill=NA,alpha=0.9, aes(color=Treatment))+
  geom_jitter(aes(shape=Site, color=Treatment),width=0.2, alpha=0.7)+
  facet_wrap(~Variable,scales = "free", strip.position='left')+
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

fig1cd

####fig2.a,b
library(ggplot2)
library(DESeq2)
abund_table<-read.csv("Fig2ab_data.csv",row.names=1,check.names=FALSE)
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
      

      widths <- unit.c(
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
        final <- gtable::gtable_add_cols(panelGrobs[[1]], space) 
        final <- cbind(final, panelGrobs[[2]], size = "first") 
        
        final$heights <- heights
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

fig2b<-ggplot(df,aes(x=Treatment,y=Value,colour=Treatment))+ylab("Normalized transcript counts")+facet_wrap( ~Taxa, scales="free")
fig2b<-fig2b+geom_boxplot(alpha=0.9,outlier.shape=NA)+geom_jitter(aes(shape=Site),alpha=0.7)+theme_bw()
fig2b<-fig2b+theme(axis.text.x=element_blank(), axis.title.x=element_blank(), strip.text = element_text(size=9),panel.spacing=unit(1, "lines"))+
  xlab(NULL)
fig2b

###fig2.c
library(dplyr)
library(broom)
library(tidyr)
library(ggplot2)
d <-read.csv("Fig2c_data.csv", header = TRUE)
densities.qtiles <-
  d %>%
  ungroup() %>%
  group_by(Treatment, Length) %>% 
  summarise(q05 = quantile(dens, 0.025),
            q50 = quantile(dens, 0.5),
            q95 = quantile(dens, 0.975),
            mean_d=mean(dens)) 
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

###Fig3.ab, add noncoding regions) 
d <-read.csv("Fig3ab_data.csv", header=TRUE)
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

####Fig3b_partial
abund_table<-read.csv("Fig3b_data.csv")
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
info <- read.csv("Fig3c_data1.csv")
tree <- read.tree("Fig3c_data2.tree")
colnames(info)<-c("id1","group",'id3')
tail(info)
p<-ggtree(tree, layout='circular') %<+% info + 
  geom_tippoint(aes(color=group),size=0.5)+ 
  theme(legend.position="right", legend.text = element_text(size=10))+
  guides(color = guide_legend(override.aes = list(size = 2)))
fig3c<- msaplot(p,'Fig3c_data3.fasta', window=c(1,6), width=0.2)+ 
  theme(legend.position="right", legend.text = element_text(size=10))+
  guides(color = guide_legend(override.aes = list(size = 2)))
fig3c

###fig 4
#BiocManager::install("ggtree")
library(ggtree)
library(ggplot2)
####pss tree_1
info <- read.csv("Fig4_data1.csv")
tree <- read.tree("Fig4_data2.txt")
colnames(info)<-c("id","taxa",'Dry','Wet')
head(info)
library(RColorBrewer)
library(ggnewscale)
library(pheatmap)
library(viridis)
nb.cols <- 17
mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols)
p<-ggtree(tree) %<+% info+
  geom_cladelabel(node=271, label="Leviviridae(U)",color="blue", offset=1.5, align=TRUE,fontsize=3)+
  geom_cladelabel(node=243, label="Narnaviridae",color="blue", offset=1.5, align=TRUE,fontsize=3)+
  geom_cladelabel(node=242, label="Potyviridae",color="blue", offset=1.5, align=TRUE,fontsize=3)+
  geom_cladelabel(node=65, label="Hypoviridae",color="blue", offset=1.5, align=TRUE,fontsize=3)+
  geom_cladelabel(node=238, label="Virgaviridae",color="blue", offset=1.5, align=TRUE,fontsize=3)+
  geom_cladelabel(node=263, label="Flaviviridae",color="blue", offset=1.5, align=TRUE,fontsize=3)+
  geom_cladelabel(node=235, label="Alphaflexiviridae",color="blue", offset=1.5, align=TRUE,fontsize=3)+
  geom_cladelabel(node=59, label="Solemoviridae",color="blue", offset=1.5, align=TRUE,fontsize=3)+
  geom_cladelabel(node=58, label="Polycipiviridae",color="blue", offset=1.5, align=TRUE,fontsize=3)+
  geom_cladelabel(node=92, label="Ifaviridae",color="blue", offset=1.5, align=T,fontsize=3)+
  geom_cladelabel(node=269, label="Closteroviridae",color="blue", offset=1.5, align=T,fontsize=3)+
  geom_cladelabel(node=96, label="Bromoviridae",color="blue", offset=1.5, align=TRUE,fontsize=3)+
  geom_cladelabel(node=93, label="Picornaviridae",color="blue", offset=1.5, align=TRUE,fontsize=3)+
  geom_cladelabel(node=267, label="Tombusviridae",color="blue", offset=1.5, align=TRUE,fontsize=3)+
  geom_cladelabel(node=56, label="Botourmiaviridae",color="blue", offset=1.5, align=TRUE,fontsize=3)+
  geom_cladelabel(node=57, label="Secoviridae",color="blue", offset=1.5, align=TRUE,fontsize=3)+
  geom_cladelabel(node=174, label="Leviviridae(L)",color="blue", offset=1.5, align=TRUE,fontsize=3)+
  geom_cladelabel(node=1, label="Alphaproteobacteria",color="blue", offset=1.5, align=TRUE,fontsize=3)

heatmapData=read.csv("Fig4_data3.csv", row.names=1)

p<-rotate(p, 227)%>% rotate(271)  
p
p1<-collapse(p, node=271)+geom_point2(aes(subset=(node==271)), shape=15, size=3, fill='grey')
p1<-collapse(p1, node=174)+geom_point2(aes(subset=(node==174)), shape=15, size=3, fill='grey')+
  geom_point2(aes(subset=(node==1)), shape=16, size=3, fill='grey')+ xlim(NA, 10)
p1
plot1<-gheatmap(p1, log(heatmapData), offset = 0.1,colnames=FALSE, width=0.2,color = FALSE)+
  scale_fill_viridis_c(option="A",na.value = '#000000')+
  theme(legend.position='bottom',plot.margin = unit(c(1, 0, 0, 0), "lines"))
plot1
###optional view subclade
p2<-ggtree(tree) %<+% info+
  geom_cladelabel(node=271, label="Leviviridae(U)",color="blue", offset=0.6, align=TRUE,fontsize=3)+
  geom_cladelabel(node=174, label="Leviviridae(L)",color="blue", offset=0.6, align=TRUE,fontsize=3)+xlim(NA, 20)
p2
p21<-viewClade(p2, node=174)
p21
p22<-viewClade(p2, node=271)
p22
plot2<-gheatmap(p21, log(heatmapData), offset = -1.5,colnames=FALSE, width=0.3,color = FALSE)+xlim(NA, 15)+
  scale_fill_viridis_c(option="A",na.value = '#000000', name="log")+
  theme(legend.position='NA',plot.margin = unit(c(1, 0, 0, 0), "lines"))
plot3<-gheatmap(p22, log(heatmapData), offset = -1.5,colnames=FALSE, width=0.3,color = FALSE)+xlim(NA, 15)+
  scale_fill_viridis_c(option="A",na.value = '#000000', name="log")+
  theme(legend.position='NA',plot.margin = unit(c(1, 0, 0, 0), "lines"))

library(gridExtra)
gridExtra::grid.arrange(plot1, arrangeGrob(plot2,plot3), ncol=2)


####ds
info <- read.csv("Fig4_data4.csv")
tree <- read.tree("Fig4_data5.txt")
ggtree(tree)%<+% info+geom_text(aes(label=node), hjust=-.3, size=4)+geom_text(aes(label=taxa), size=4,nudge_x = 0.5)
plot4<-ggtree(tree) %<+% info+
  geom_cladelabel(node=10, label="Reoviridae",color="blue", offset=0.6, align=TRUE,fontsize=3)+
  geom_cladelabel(node=2, label="Partitiviridae",color="blue", offset=0.6, align=TRUE,fontsize=3)+
  geom_cladelabel(node=3, label="Partitiviridae",color="blue", offset=0.6, align=TRUE,fontsize=3)+
  geom_cladelabel(node=1, label="Alphaproteobacteria",color="blue", offset=0.6, align=TRUE,fontsize=3)+
  geom_point2(aes(subset=(node==1)), shape=16, size=3, fill='grey')+xlim(NA, 5)
plot4
heatmapData=read.csv("dsRNA_heatmap_dryWet_info.csv", row.names=1)
plot4<-gheatmap(plot4, log(heatmapData), offset = 0.1,colnames=FALSE, width=0.2, high=0.001,color = FALSE)+
  scale_fill_viridis_c(option="A",na.value = '#000000')+theme(legend.position='NA',plot.margin = unit(c(1, 1, 1, 1), "lines"))
heatmapData
###nss
info <- read.csv("Fig4_data6.csv")
tree <- read.tree("Fig4_data7.txt")
heatmapData=read.csv("Fig4_data8.csv", row.names=1)
ggtree(tree)%<+% info+geom_text(aes(label=node), hjust=-.3, size=4)+geom_text(aes(label=taxa), size=1.5,nudge_x = 0.4)
plot5<-ggtree(tree) %<+% info+
  geom_cladelabel(node=8, label="Mymonaviridae",color="blue", offset=0.8, align=TRUE,fontsize=3)+
  geom_cladelabel(node=7, label="Unclassified Mononegavirales",color="blue", offset=0.8, align=TRUE,fontsize=3)+
  geom_cladelabel(node=6, label="Paramyxoviridae",color="blue", offset=0.8, align=TRUE,fontsize=3)+
  geom_cladelabel(node=5, label="Nairoviridae",color="blue", offset=0.8, align=TRUE,fontsize=3)+
  geom_cladelabel(node=4, label="Unclassified Bunyavirales",color="blue", offset=0.8, align=TRUE,fontsize=3)+
  geom_cladelabel(node=15, label="Peribunyaviridae",color="blue", offset=0.8, align=TRUE,fontsize=3)+
  geom_cladelabel(node=21, label="Arenaviridae",color="blue", offset=0.8, align=TRUE,fontsize=3)+
  geom_cladelabel(node=9, label="Hantaviridae",color="blue", offset=0.8, align=TRUE,fontsize=3)+
  geom_cladelabel(node=1, label="Alphaproteobacteria",color="blue", offset=0.8, align=TRUE,fontsize=3)+
  geom_point2(aes(subset=(node==1)), shape=16, size=3, fill='grey')+xlim(NA, 8)
plot5

plot5<-gheatmap(plot5, log(heatmapData), offset = 0.05,colnames=FALSE, hjust=0, font.size=3,width=0.2,colnames_offset_y = 1, color = FALSE)+
  scale_fill_viridis_c(option="A",na.value = '#000000')+theme(legend.position='bottom',plot.margin = unit(c(1, 0, 0, 0), "lines"))
plot5
library(gridExtra)
gridExtra::grid.arrange(arrangeGrob(plot4,plot5),plot1, arrangeGrob(plot3,plot2), ncol=3)

####fig 5.a
library(ggplot2)
library(pheatmap)
data<-read.delim("Fig5a_data.txt",header=T, row.names="taxa")
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

##fig 5.b
library(ggplot2)
library(DESeq2)
abund_table<-read.delim("Fig5b_data.txt",row.names=1,check.names=FALSE)
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
  labs(x = "Mean abundance", y = "Log2 fold change")+theme_bw()

###only show the significant
if(label == T){
  if (!is.null(tax.display)){
    rlab <- data.frame(res_tax, Display = apply(res_tax[,c(tax.display, tax.aggregate)], 1, paste, collapse="; "))
  } else {
    rlab <- data.frame(res_tax, Display = res_tax[,tax.aggregate])
  }
  #p1 <- p1 + geom_text(data = subset(rlab, Significant == "Yes"), aes(label = Display), size = 3, vjust = 1, hjust=1)
  p1 <- p1 + geom_text(data = subset(rlab), aes(label = Display), size = 3, vjust = 1, hjust=1)
}
p1
###or show all the label
library(ggrepel)
if(label == T){
  if (!is.null(tax.display)){
    rlab <- data.frame(res_tax, Display = apply(res_tax[,c(tax.display, tax.aggregate)], 1, paste, collapse="; "))
  } else {
    rlab <- data.frame(res_tax, Display = res_tax[,tax.aggregate])
  }
  p1 <- p1 + geom_text(data = subset(rlab), aes(label = Display), size = 3, vjust = 1, hjust=1)
}

p1
p1+geom_text_repel()

res_tax_sig_abund = cbind(as.data.frame(countData[rownames(res_tax_sig), ]), OTU = rownames(res_tax_sig), padj = res_tax[rownames(res_tax_sig),"padj"]) 
write.table(res_tax_sig_abund, file = "RNA_virus_by_cat.txt", sep = "\t", row.names = TRUE)
data<-abund_table
data<-as.data.frame(data)
df<-NULL
for(i in res_tax[rownames(res_tax_sig),"OTU"]){
  #tmp<-data.frame(data[,i],grouping_info[2],grouping_info[1], rep(paste(i," padj = ",scientific(round(res_tax[i,"padj"],5)),sep=""),dim(data)[1]))
  #tmp<-data.frame(data[,i],grouping_info[2],grouping_info[1], rep(paste(i," padj = ",res_tax[i,"padj"],sep=""),dim(data)[1]))
  tmp<-data.frame(data[,i],grouping_info[2],grouping_info[1], rep(paste(i,", p=",formatC(res_tax[i,"padj"], format='e',digits=2),sep=""),dim(data)[1]))
  if(is.null(df)){df<-tmp} else { df<-rbind(df,tmp)} 
}
colnames(df)<-c("Value","Treatment","Site","Taxa")
df
fig5b<-ggplot(df,aes(x=Treatment,y=Value, colour=Treatment))+
  ylab("Normalized reads coverage of RNA viral contig")+
  facet_wrap( ~Taxa, scales="free", ncol=4)+
  geom_boxplot(alpha=0.9,outlier.shape=NA)+geom_jitter(aes(shape=Site),alpha=0.7)+
  theme_bw()+
  theme(axis.title.x = element_blank(), axis.text.x=element_blank(),strip.text = element_text(size=9),panel.spacing=unit(1, "lines"))+
  xlab(NULL)
fig5b

##################################################################################################################################################################
###extended data file 
####Supplementary data 2a: 16s 
library(ggplot2)
library(DESeq2)
abund_table<-read.csv("SupplementaryData2a.csv",row.names=1,check.names=FALSE)
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
SupplementaryData2a <- ggplot(data = res_tax, aes(x = baseMean, y = log2FoldChange, color = Significant)) +
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
  SupplementaryData2a <- SupplementaryData2a + geom_text(data = subset(rlab, Significant == "Yes"),aes(label = Display), size = 2, vjust = 1.5)
}
SupplementaryData2a


#SupplemtaryData2b:18s
data<-read.delim("SupplementaryData2b.txt",header=T, row.names="gene")
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




####SupplementaryData2cd:16s 18s correlate to DNA and RNA viruses
library(ggplot2)
abdCorr<-read.csv("SupplementaryData2cd.csv", header = TRUE)
ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(data=abdCorr, aes(y=BacABD, x=DNAVtranscripts, grouping(Treatment)),fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
    geom_point(aes(colour =Treatment, shape=Site), size = 2) +
    stat_smooth(method = "lm", se=TRUE, fullrange=FALSE, level=0.95) +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 4),
                       "Intercept =",signif(fit$coef[[1]],4 ),
                       " Slope =",signif(fit$coef[[2]], 4),
                       " P =",signif(summary(fit)$coef[2,4], 4)))
} 
fit1 <- lm(DNAVtranscripts~BacABD, data = abdCorr)
ggplotRegression(fit1)

ggplotRegression(fit1)+facet_wrap(~Treatment)

library("ggplot2")
library("datasets")
library("plyr")
df_16s<-read.csv("SupplementaryData2cd.csv", header = TRUE)
regression=function(df_16s){
  reg_fun<-lm(formula=df_16s$BacABD~df_16s$DNAVtranscripts) #regression function
  slope<-round(coef(reg_fun)[2],3)  
  intercept<-round(coef(reg_fun)[1],3) 
  R2<-round(as.numeric(summary(reg_fun)[8]),3)
  R2.Adj<-round(as.numeric(summary(reg_fun)[9]),3)
  c(slope,intercept,R2,R2.Adj)
}
regressions_data_16s<-ddply(df_16s,"Treatment",regression)
colnames(regressions_data_16s)<-c ("Treatment","slope","intercept","R2","R2.Adj")
df

ggplot(data = df_16s, aes(x = BacABD, y = DNAVtranscripts)) +
  geom_smooth(method = "lm") +
  geom_point()+ facet_grid(Treatment~., scales = "free")+
  xlab('16s rRNA gene abundance')+
  ylab('DNA viral transcript abundance')
label_16s<-paste("y=",regressions_data_16s$slope,"x+",regressions_data_16s$intercept,', ',"R^2=",regressions_data_16s$R2,", ","R^2.Adj=",regressions_data_16s$R2.Adj)
label_16s

df_18s<-read.csv("SupplementaryData2cd.csv", header = TRUE)
regression=function(df_18s){
  reg_fun<-lm(formula=df_18s$EukABD~df_18s$RNAVabd) #regression function18sABD	RNAVabd
  slope<-round(coef(reg_fun)[2],3)  
  intercept<-round(coef(reg_fun)[1],3) 
  R2<-round(as.numeric(summary(reg_fun)[8]),3)
  R2.Adj<-round(as.numeric(summary(reg_fun)[9]),3)
  c(slope,intercept,R2,R2.Adj)
}
regressions_data_18s<-ddply(df_18s,"Treatment",regression)
colnames(regressions_data_18s)<-c ("Treatment","slope","intercept","R2","R2.Adj")
df

ggplot(data = df_18s, aes(x = EukABD, y = RNAVabd)) +
  geom_smooth(method = "lm") +
  geom_point()+ facet_grid(Treatment~., scales = "free")+
  xlab('18s rRNA gene transcript abundance')+
  ylab('RNA viral abundance')
label_18s<-paste("y=",regressions_data_18s$slope,"x+",regressions_data_18s$intercept,', ',"R^2=",regressions_data_18s$R2,", ","R^2.Adj=",regressions_data_18s$R2.Adj)
label_18s

regressions_data_16s<-ddply(df_16s,'Select', regression)
colnames(regressions_data_16s)<-c ('Select',"slope","intercept","R2","R2.Adj")
ggplot(data = df_16s, aes(x = BacABD, y = DNAVtranscripts)) +
  geom_smooth(method = "lm") +
  geom_point()+
  xlab('16s rRNA gene abundance')+
  ylab('DNA viral transcript abundance')
label_16s<-paste("y=",regressions_data_16s$slope,"x+",regressions_data_16s$intercept,', ',"R^2=",regressions_data_16s$R2,", ","R^2.Adj=",regressions_data_16s$R2.Adj)
label_16s

regressions_data_18s<-ddply(df_18s,'Select', regression)
colnames(regressions_data_18s)<-c ('Select',"slope","intercept","R2","R2.Adj")
ggplot(data = df_18s, aes(x = EukABD, y = RNAVabd)) +
  geom_smooth(method = "lm") +
  geom_point()+ 
  xlab('18s rRNA gene transcript abundance')+
  ylab('RNA viral abundance')
label_18s<-paste("y=",regressions_data_18s$slope,"x+",regressions_data_18s$intercept,', ',"R^2=",regressions_data_18s$R2,", ","R^2.Adj=",regressions_data_18s$R2.Adj)
label_18s