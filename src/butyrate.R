# ##############################################################################
#
## Check the butyrate production capacity of the microbiome in CRC
#
# ##############################################################################

library("tidyverse")
library("here")
library("ggpubr")
library("ggthemes")
library("cowplot")

# ##############################################################################
## PART 1: Check for butyrate production in the Gut Metabolic Modules (GMMs),
## see for that Vieira-Silva et al. (https://doi.org/10.1038/nmicrobiol.2016.88)
## computed as part of our CRC meta-analysis 
## (http://dx.doi.org/10.1038/s41591-019-0406-6)

data.location <- 'https://www.embl.de/download/zeller/crc_meta/'

gmm.data <- read.table(paste0(data.location, 
                              'functional_profiles/GMM_profiles.tsv'), 
                       sep='\t', stringsAsFactors = FALSE, check.names = FALSE)
gmm.data <- as.matrix(gmm.data)

gmm.pvals <- read.table(paste0(data.location, 'results/GMM_p_adj.tsv'),
                        sep='\t', stringsAsFactors = FALSE, check.names = FALSE)

meta.all <- read_tsv(paste0(data.location, 'metadata/meta_all.tsv')) %>% 
  filter(Group %in% c('CTR', 'CRC')) %>% 
  filter(Sample_ID %in% colnames(gmm.data))

df.plot <- meta.all %>% 
  filter(Sample_ID %in% colnames(gmm.data)) %>% 
  select(Sample_ID, Group) %>% 
  mutate(`butyrate production I`=log10(gmm.data['MF0088', Sample_ID])+1e-09) %>% 
  mutate(`butyrate production II`=log10(gmm.data['MF0089', Sample_ID])+1e-09)


# butyrate production I
g1 <- df.plot %>% 
  mutate(Group=factor(Group, levels=c('CTR', 'CRC'))) %>% 
  ggplot(aes(x=Group, y=`butyrate production I`, fill=Group)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha=0.3) +
  theme_gdocs() +
  theme(plot.background = element_blank(),
        plot.title = element_text(size=rel(1))) + 
  scale_fill_manual(values=c("#C2C2C290", "#d9544d90"), guide=FALSE) +
  xlab('') +
  ylab('log10(relative abundance)') + 
  ggtitle('butyrate production I') +
  annotate(x=1.5, y=-3.4, geom='text', 
           label=paste0('blocked Wilcoxon, p = ', 
                        formatC(format='e', digits=2, 
                                gmm.pvals['MF0088', 'all'])))

# butyrate production II
g2 <- df.plot %>% 
  mutate(Group=factor(Group, levels=c('CTR', 'CRC'))) %>% 
  ggplot(aes(x=Group, y=`butyrate production II`, fill=Group)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha=0.3) +
  theme_gdocs() +
  theme(plot.background = element_blank(),
        plot.title = element_text(size=rel(1))) + 
  scale_fill_manual(values=c("#C2C2C290", "#d9544d90"), guide=FALSE) +
  xlab('') +
  ylab('log10(relative abundance)') + 
  ggtitle('butyrate production II') +
  annotate(x=1.5, y=-3.4, geom='text', 
           label=paste0('blocked Wilcoxon, p = ', 
                        sprintf(fmt='%.2f', 
                                gmm.pvals['MF0089', 'all'])))

# combine
g <- plot_grid(g1, g2)
ggsave(g, filename = here('figures', 'butyrate_modules.pdf'),
       width = 5, height = 4, useDingbats=FALSE)


# ##############################################################################
## PART 2: As shotgun metagenomic data also allows for a direct quantification 
## of gene and pathway enrichments (summarised across all community members 
## encoding these), we analysed butyrate production in the metagenomic data 
## from our meta-analysis and  asked whether these pathways are enriched in 
## control relative to CRC metagenomes. Pathway definitions were based on 
## Vital et al. (https://msystems.asm.org/content/2/6/e00130-17)

## Explanation: I downloaded the database of butyrate-pathway enzymes 
## published with the paper above. Created HMMs based on MUSCLE-alignments 
## of these proteins for each of the 27 different genes.
## I then used the HMMs to search for fitting genes in the IGC gene catalogue.
## The e-values were chosen for each gene individually by looking at the 
## resulting e-value distribution.
## The downloaded files are the HMM hits (butyrate_hits.tsv) and the 
## quantification of these genes from the functional profiling that was done 
## for the meta-analysis (butyrate_genes.tsv).

library("matrixStats")
library("ComplexHeatmap")
library("circlize")
library("coin")

gene.mat.but <- read.table(paste0(data.location, 
                                  'functional_profiles/butyrate_genes.tsv'), 
                           sep='\t', quote='', check.names = FALSE)
gene.mat.but <- as.matrix(gene.mat.but)
# apply abundance filtering
f.idx <- which(rowMaxs(gene.mat.but) > 1e-07)
gene.mat.but <- gene.mat.but[f.idx,]

but.hits <- read_tsv(paste0(data.location, 'files/butyrate_hits.tsv')) %>% 
  filter(gene %in% rownames(gene.mat.but))

# sum up pathway abundance
pathways <- list(
  'Acetyl-CoA' = c('thl', 'bhbd', 'cro', 'but', 'buk'),
  'Glutarate' = c('gctA', 'gctB', 'hgCoAdA', 'hgCoAdB', 'hgCoAdC', 'gcdA', 
                  'gcdB'),
  '4Aminobutyrate' = c('abfD', 'abfH', '4hbt'),
  'Lysine' =  c('kamA', 'kamD', 'kamE', 'kdd', 'kce', 'kal',  'atoA', 'atoD'),
  'Central-Complex'=c('bcd', 'etfA', 'etfB'))

df.temp <- meta.all %>% 
  select(Sample_ID, Study, Group)
df.plot <- tibble(pathway=character(0), study=character(0), gFC=double(0))

for (pathway.name in names(pathways)){
  but.hits.p <- but.hits %>% 
    filter(query %in% pathways[[pathway.name]])
  
  temp <- gene.mat.but[but.hits.p$gene,df.temp$Sample_ID]
  
  df.temp <- df.temp %>% 
    mutate(!!pathway.name:=log10(colSums(temp) + 1e-09))
  for (s in unique(df.temp$Study)){
    g.ctr <- df.temp %>% 
      filter(Study==s) %>% 
      filter(Group=='CTR') %>% 
      pull(pathway.name) %>% 
      quantile(probs=seq(0.05, 0.95, 0.05))
    g.crc <- df.temp %>% 
      filter(Study==s) %>% 
      filter(Group=='CRC') %>% 
      pull(pathway.name) %>% 
      quantile(probs=seq(0.05, 0.95, 0.05))
    df.plot <- df.plot %>% 
      add_row(pathway=pathway.name, study=str_remove(s, '-CRC'), 
              gFC=mean(g.ctr-g.crc))
  }
}

# general preparations
colors <- ggthemes_data$tableau$`color-palettes`$
  `ordered-diverging`$`Red-Blue Diverging`$value

pathway.colours <- c("Acetyl-CoA"='#e87f28',
                     "Glutarate"='#314975',
                     "4Aminobutyrate"='#286a34', 
                     "Lysine"='#995229',
                     "Central-Complex"='slategrey')

matrix.plot <- df.plot %>% 
  mutate(study=factor(study, levels=c('FR', 'AT', 'CN', 'US', 'DE'))) %>% 
  mutate(pathway=factor(pathway, levels=names(pathways))) %>% 
  spread(key=study, value=gFC) %>% 
  as.data.frame()
rownames(matrix.plot) <- matrix.plot$pathway
matrix.plot$pathway <- NULL

bt.annot <- HeatmapAnnotation(
  df=data.frame(pathway=df.plot %>% 
                  mutate(pathway=factor(pathway, levels=names(pathways))) %>% 
                  spread(key=study, value=gFC) %>% 
                  pull(pathway)), 
  col=list(pathway=pathway.colours),
  show_legend = FALSE)


h <- Heatmap(t(matrix.plot), cluster_rows = FALSE, cluster_columns = FALSE,
             col = colorRamp2(colors=colors,breaks=seq(-0.3, 0.3, 0.1)), 
             cell_fun = function(j, i, x, y, width, height, fill){
               grid.text(sprintf("%.2f", -t(matrix.plot)[i, j]), x, y, 
                         gp = gpar(fontsize = 10))},
             show_heatmap_legend = FALSE,
             bottom_annotation = bt.annot)
pdf(here('figures', 'butyrate_pathways_metaG_gFC.pdf'), 
    width = 5, height = 4, useDingbats = FALSE)
print(h)
dev.off()

# for single genes
df.temp <- meta.all %>% 
  select(Sample_ID, Study, Group)
df.plot <- tibble(gene=character(0), study=character(0), gFC=double(0))
df.p <- tibble(gene=character(0), p.val=double(0))

for (q in unique(but.hits$query)){
  but.hits.p <- but.hits %>% 
    filter(query == q)
  
  temp <- gene.mat.but[but.hits.p$gene,df.temp$Sample_ID]
  
  df.temp <- df.temp %>% 
    mutate(!!q:=log10(colSums(temp) + 1e-09))
  
  for (s in unique(df.temp$Study)){
    g.ctr <- df.temp %>% 
      filter(Study==s) %>% 
      filter(Group=='CTR') %>% 
      pull(q) %>% 
      quantile(probs=seq(0.05, 0.95, 0.05))
    g.crc <- df.temp %>% 
      filter(Study==s) %>% 
      filter(Group=='CRC') %>% 
      pull(q) %>% 
      quantile(probs=seq(0.05, 0.95, 0.05))
    df.plot <- df.plot %>% 
      add_row(gene=q, study=str_remove(s, '-CRC'), 
              gFC=mean(g.ctr-g.crc))
  }
  t <- wilcox_test(df.temp[[q]]~as.factor(df.temp$Group)|
                     as.factor(df.temp$Study))
  df.p <- df.p %>% 
    add_row(gene=q, p.val=as.double(pvalue(t)))
}

pathway.info <- str_remove(names(unlist(pathways)), '[0-9]$')

matrix.plot <- df.plot %>% 
  mutate(study=factor(study, levels=c('FR', 'AT', 'CN', 'US', 'DE'))) %>% 
  mutate(gene=factor(gene, levels=unlist(pathways))) %>% 
  na.omit() %>% 
  spread(key=study, value=gFC) %>% 
  as.data.frame()
rownames(matrix.plot) <- matrix.plot$gene
matrix.plot$gene <- NULL

bt.annot <- HeatmapAnnotation(
  df=data.frame(pathway=pathway.info),
  col=list(pathway=pathway.colours),
  show_legend = FALSE)


h <- Heatmap(t(matrix.plot), cluster_rows = FALSE, cluster_columns = FALSE,
             col = colorRamp2(colors=colors,
                              breaks=seq(-0.51, 0.51, 
                                         length.out = length(colors))), 
             cell_fun = function(j, i, x, y, width, height, fill){
               grid.text(sprintf("%.2f", -t(matrix.plot)[i, j]), x, y, 
                         gp = gpar(fontsize = 10))},
             show_heatmap_legend = FALSE,
             bottom_annotation = bt.annot)
pdf(here('figures', 'butyrate_genes_metaG_gFC.pdf'), 
    width = 10, height = 5, useDingbats = FALSE)
print(h)
dev.off()

g <- df.p %>% 
  mutate(gene=factor(gene, levels=unlist(pathways))) %>% 
  na.omit() %>% 
  ggplot(aes(x=gene, y=-log10(p.val))) + 
  geom_bar(stat='identity') + 
  theme_base() + 
  xlab('') + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(g, filename = here('figures', 'butyrate_genes_metaG_pval.pdf'), 
       width = 10, height = 3, useDingbats=FALSE)    
