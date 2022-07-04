# ##############################################################################
#
## Analyze the 16S mouse data
#
# ##############################################################################

# packages
library("tidyverse")
library("here")
library("vegan")
library("labdsv")
library("progress")
library("matrixStats")
library("pheatmap")
library("lmerTest")
library("RColorBrewer")
library("ggrepel")

# ##############################################################################
# functions
.f_alpha <- function(mat, min.depth=NULL){
  if (is.null(min.depth)){
    min.depth <- min(range(colSums(mat)))
    message("Rarefaction depth: ", min.depth)
  }
  
  df.res <- list()
  pb <- progress_bar$new(total=100)
  for (r in seq_len(100)){
    temp <- rrarefy(t(mat), sample=min.depth)
    for (ind in c('shannon', 'simpson', 'invsimpson', 'richness')){
      if (ind == 'richness'){
        div <- enframe(rowSums(temp>3),
                       name='Sample_ID', value='Diversity') %>% 
          mutate(type=ind, n.repeat=r)
      } else {
        div <- diversity(temp, index=ind) %>% 
          enframe(name='Sample_ID', value='Diversity') %>% 
          mutate(type=ind, n.repeat=r)
      }
      df.res[[(length(df.res)+1)]] <- div
    }
    pb$tick()
  }
  df.res <- bind_rows(df.res)
  df.res <- df.res %>% 
    group_by(Sample_ID, type) %>% 
    summarise(div=mean(Diversity)) %>% 
    ungroup()
  return(df.res)
}
.f_beta <- function(mat, method='euclidean'){
  dist <- vegdist(t(mat), method = method)
  pco.res <- pco(dist)
  df.plot <- as_tibble(pco.res$points, rownames='Sample_ID', 
                       .name_repair=make.names)
  colnames(df.plot) <- c('Sample_ID', 'Axis1', 'Axis2')
  axis.names <- paste0(
    c('PCo 1 [', 'PCo 2 ['),
    sprintf(fmt='%.2f', c(pco.res$eig[1:2]/sum(pco.res$eig))*100),
    c('%]', '%]'))
  return(list(df.plot, axis.names))
}
.f_summarize <- function(mat, taxa, level='Genus'){
  stopifnot(nrow(mat) == nrow(taxa))
  taxa <- taxa %>% 
    as_tibble(rownames='ID')
  taxa$lvl <- case_when(is.na(taxa[[level]])~'unclassified', 
                        TRUE~taxa[[level]])
  
  bins <- unique(taxa$lvl)
  mat.res <- matrix(NA, ncol = ncol(mat), nrow=length(bins),
                    dimnames = list(bins, colnames(mat)))
  for (g in bins){
    mat.res[g,] <- colSums(mat[taxa %>% filter(lvl==g) %>% pull(ID),,
                               drop=FALSE])
  }
  return(mat.res)
}

# ##############################################################################
# load data
feat.1 <- read.table(here('files', 'feat_exp_1.tsv'), sep='\t', 
                     stringsAsFactors = FALSE, check.names = FALSE, 
                     quote = '', comment.char = '')
feat.1 <- as.matrix(feat.1)
pheatmap(log10(prop.table(feat.1, 2) + 1e-05), show_rownames = FALSE)
# size filter
hist(nchar(rownames(feat.1)), 100)
feat.1 <- feat.1[nchar(rownames(feat.1)) > 380,]
# feat.1 <- feat.1[nchar(rownames(feat.1)) < 425,]
pheatmap(log10(prop.table(feat.1, 2) + 1e-05), show_rownames = FALSE)
feat.rel.1 <- prop.table(feat.1, 2)

# abundance filter
# maybe at least 3 samples over 1e-03
hist(log10(rowMaxs(feat.rel.1) + 1e-05), 100)
feat.rel.1 <- feat.rel.1[rowSums(feat.rel.1 > 1e-03) > 3,]
# prevalence filter
hist(rowMeans(feat.rel.1!=0), 100)
feat.rel.1 <- feat.rel.1[rowMeans(feat.rel.1 != 0) > 0.1,]
pheatmap(log10(feat.rel.1 + 1e-05), show_rownames = FALSE)

# metadata
meta.1 <- read_tsv(here('files', 'meta_exp_1.tsv'))

# taxa information
taxa.1 <- read.table(here('files', 'taxonomy_1.tsv'), sep='\t', 
                     stringsAsFactors = FALSE, check.names = FALSE,
                     quote = '', comment.char = '')
taxa.1 <- taxa.1[rownames(feat.1),]
# convert to genus-level abundances
feat.genus.1 <- .f_summarize(feat.1, taxa.1)
feat.family.1 <- .f_summarize(prop.table(feat.rel.1, 2), 
                              taxa.1[rownames(feat.rel.1),], level='Family')
feat.order.1 <- .f_summarize(feat.1, taxa.1, level='Order')
feat.class.1 <- .f_summarize(feat.1, taxa.1, level='Class')

# ##############################################################################
# alpha diversity
df.plot <- full_join(meta.1, .f_alpha(feat.1, min.depth = 35000))

temp <- df.plot %>% 
  filter(type=='shannon') %>% 
  mutate(div=exp(div)) %>% 
  mutate(type='eff_richness')

df.plot <- rbind(df.plot, temp)
pdf(here('figures', 'ASV_alpha_diversity_all.pdf'), 
    useDingbats = FALSE,  width = 4, height = 3.5)
for (x in unique(df.plot$type)){
  g <- df.plot %>% 
    filter(type==x) %>% 
    mutate(Group=case_when(Treatment=='T0'~' ', 
                           Group=='Bp'~'CC4',
                           TRUE~Group)) %>% 
    mutate(Group=factor(Group, levels = c(' ', 'C', 'Ps', 'CC4'))) %>% 
    mutate(Treatment=factor(Treatment, 
                            levels = c('T0', 'T1', 'T2', 'T3', 'TK'))) %>% 
    ggplot(aes(x=Group, y=div, fill=Group)) + 
    geom_bar(stat='summary', fun='mean') +
    geom_jitter(aes(col=Group), position = position_jitterdodge(), pch=21, 
                fill='white') +
    facet_grid(~Treatment, scales = 'free', space='free_x') + 
    xlab('') + 
    ylab("Alpha diversity") + 
    scale_fill_manual(values=c('#333333', '#0024b1', 
                               '#da4432','#4a7a31'), guide=FALSE) +
    scale_colour_manual(values=c('#333333', '#0024b1', 
                                 '#da4432','#4a7a31'), guide=FALSE) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5),
          panel.grid = element_blank(),
          strip.background = element_blank(), 
          panel.border = element_blank(),
          axis.line = element_line(size=0.4),
          axis.ticks.x = element_blank())
  
  
  temp <- df.plot %>% 
    filter(type==x)
  
  fit <- lmerTest::lmer(div~Group+(1|Treatment), data=temp)
  res <- anova(fit)
  
  g <- g + ggtitle(paste0(x, ', P=', sprintf(fmt='%.2f', res$`Pr(>F)`[1])))
  print(g)
}
dev.off()

# ##############################################################################
# beta diversity
dis <- 'euclidean'
res <- .f_beta(log10(feat.rel.1[,meta.1 %>% 
                                  filter(Treatment=='TK') %>% 
                                  pull(Sample_ID)] + 1e-05),
               method='euclidean')
mat.dist <- vegdist(t(log10(feat.rel.1[,meta.1 %>% 
                                         filter(Treatment=='TK') %>% 
                                         pull(Sample_ID)] + 1e-05)), 
                    method = 'euclidean')
adonis(mat.dist~meta.1 %>% 
         filter(Treatment=='TK') %>% 
         pull(Group))
df.plot <- full_join(meta.1, res[[1]], by='Sample_ID') %>% 
  filter(!is.na(Axis1))
group.mean <- df.plot %>% 
  group_by(Group) %>% 
  summarise(x=mean(Axis1), y=mean(Axis2))
g <- df.plot %>% 
  left_join(group.mean, by='Group') %>% 
  mutate(Group=factor(Group, levels = c('C', 'Bp', 'Ps'))) %>% 
  ggplot(aes(x=Axis1, y=Axis2, col=Group)) + 
  geom_point() + 
  geom_point(shape=17, data=group.mean, aes(x=x, y=y, col=Group)) + 
  geom_segment(aes(x=x, y=y, xend=Axis1, yend=Axis2), alpha=0.3) + 
  xlab(res[[2]][1]) + 
  ylab(res[[2]][2]) + 
  scale_colour_manual(values=c('#0024b1', '#4a7a31', '#da4432'),
                      labels=c('C', 'CC4', 'Ps')) + 
  theme_bw() + 
  NULL
ggsave(g, filename = here('figures','ASV_beta_diversity.pdf'),
       width = 4, height = 2.5, useDingbats=FALSE)
dev.off()
# ##############################################################################
# heatmap
mat.log <- log10(feat.rel.1[,meta.1 %>% 
                              filter(Treatment=='TK') %>% 
                              pull(Sample_ID)] + 1e-05)

pdf(here('figures','ASV_heatmap_TK.pdf'), 
    width = 4, height = 6, useDingbats = FALSE)
annot.col <- meta.1 %>% 
  filter(Treatment=='TK') %>% 
  mutate(Group=case_when(Group=='Bp'~"CC4",TRUE~Group)) %>% 
  select(Sample_ID, Group) %>% 
  as.data.frame()
rownames(annot.col) <- annot.col$Sample_ID
annot.col$Sample_ID <- NULL
pheatmap(mat.log, show_rownames = FALSE, 
         clustering_method = 'ward.D2',
         annotation_col = annot.col, 
         border_color = NA, 
         annotation_colors = list('Group'=c('C'='#0024b1', 
                                            'CC4'='#4a7a31', 'Ps'='#da4432')))
dev.off()

# ##############################################################################
# fold change at TK
fold.change.ps <- vapply(rownames(feat.rel.1),FUN=function(x){
  ctr <- feat.rel.1[x,meta.1 %>% 
                      filter(Treatment=='TK', Group=='C') %>% 
                      pull(Sample_ID)] %>% 
    +1e-04 %>% 
    log10 %>% sort
  case <- feat.rel.1[x,meta.1 %>% 
                       filter(Treatment=='TK', Group=='Ps') %>% 
                       pull(Sample_ID)] %>%   
    +1e-04 %>% 
    log10 %>% sort
  fc <- mean(case-ctr)
  p.val <- wilcox.test(ctr, case, exact=FALSE)
  return(c(fc, p.val$p.val))},FUN.VALUE = double(2))
fold.change.bp <- vapply(rownames(feat.rel.1),FUN=function(x){
  ctr <- feat.rel.1[x,meta.1 %>% 
                      filter(Treatment=='TK', Group=='C') %>% 
                      pull(Sample_ID)] %>% 
    +1e-04 %>% 
    log10 %>% sort
  case <- feat.rel.1[x,meta.1 %>% 
                       filter(Treatment=='TK', Group=='Bp') %>% 
                       pull(Sample_ID)] %>%   
    +1e-04 %>% 
    log10 %>% sort
  fc <- mean(case-ctr)
  p.val <- wilcox.test(ctr, case, exact=FALSE)
  return(c(fc, p.val$p.val))},FUN.VALUE = double(2))

rownames(fold.change.ps) <- c('fc_ps', 'pval_ps')
rownames(fold.change.bp) <- c('fc_bp', 'pval_bp')

fold.change.ps <- t(fold.change.ps) %>% as_tibble(rownames='ASV')
fold.change.bp <- t(fold.change.bp) %>% as_tibble(rownames='ASV')

df.plot <- full_join(fold.change.bp, fold.change.ps) %>% 
  mutate(family=taxa.1[ASV, 'Family']) %>% 
  mutate(diff=fc_bp-fc_ps) %>% 
  mutate(select=abs(diff) > 1) %>% 
  filter(!is.na(family)) %>% 
  mutate(label=paste0('ASV-', seq_len(nrow(.)))) %>% 
  mutate(label=case_when(select~label, TRUE~''))

g <- df.plot %>% 
  mutate(diff=factor(-sign(diff))) %>% 
  mutate(fc.type=case_when(diff==-1~fc_bp, diff==1~fc_ps)) %>%  
  mutate(label=case_when(fc.type > 0.1~label, TRUE~'')) %>%
  mutate(select=label!='') %>% 
  ggplot(aes(x=fc_bp, y=fc_ps, col=select, alpha=select)) + 
  geom_hline(yintercept = 0, colour='darkgrey') +
  geom_vline(xintercept = 0, colour='darkgrey') +
  geom_abline(slope = 1, intercept = 0, colour='darkgrey', lty=3) +
  geom_point(pch=19, stroke=0) + 
  theme_bw() + 
  scale_colour_manual(values=c('#707372', '#FFA300'), guide=FALSE) + 
  scale_alpha_manual(values=c(0.75, 1), guide=FALSE) +
  xlab('Generalized fold change CC4 vs Ctrl') +
  ylab('Generalized fold change Ps vs Ctrl') + 
  geom_text_repel(aes(label=label), size=1, segment.size = 0.3)
ggsave(g, filename = here('figures', 'ASV_fc_plot.pdf'), 
       width = 3, height = 3, useDingbats=FALSE)


# heatmap
g <- log10(feat.rel.1 + 1e-04) %>% 
  as_tibble(rownames = 'ASV') %>% 
  pivot_longer(-ASV, names_to = 'Sample_ID', values_to = 'ab') %>% 
  left_join(df.plot, by='ASV') %>% 
  filter(select) %>% 
  mutate(diff=factor(-sign(diff))) %>% 
  mutate(fc=case_when(diff==-1~fc_bp, diff==1~fc_ps)) %>% 
  filter(fc > 0.1) %>% 
  arrange(desc(fc)) %>% 
  mutate(label=factor(label, levels = unique(label))) %>% 
  separate(Sample_ID, into=c('Group', 'Treatment', 'ID'), sep = '-') %>% 
  mutate(Group=case_when(Treatment=='T0'~'Initial', TRUE~Group)) %>% 
  mutate(split.group=case_when(fc_ps>fc_bp~'Ps',TRUE~'Bp')) %>% 
  group_by(Group, Treatment, label) %>% 
  summarise(m=mean(ab), split.group=unique(split.group)) %>% 
  ungroup() %>% 
  mutate(Group=factor(Group, levels = c('Initial', 'C', 'Bp', 'Ps'))) %>% 
  ggplot(aes(x=Treatment, y=label, fill=m)) + 
  geom_tile() + 
  facet_grid(split.group~Group, scales = 'free', space = 'free') + 
  xlab('') + 
  ylab('') + 
  scale_fill_gradientn(
    colours=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100),
    limits=c(-4, -1), name='log10(Relative\n abundance)') + 
  theme_bw() + 
  theme(strip.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank())
ggsave(g, filename = here('figures','ASV_full_heatmap.pdf'), 
       width = 6, height = 4, useDingbats=FALSE)

df.plot.2 <- feat.rel.1[df.plot %>% 
                          filter(select) %>% 
                          pull(ASV), 
                        meta.1 %>% 
                          filter(Treatment=='TK') %>% 
                          pull(Sample_ID)] %>% 
  as_tibble(rownames = 'ASV') %>% 
  pivot_longer(-ASV, names_to = 'Sample_ID', values_to = 'ab') %>% 
  mutate(ab=log10(ab+1e-04)) %>% 
  left_join(meta.1, by='Sample_ID') %>% 
  left_join(df.plot %>% select(ASV, family, diff, fc_bp, fc_ps, label), 
            by='ASV') %>% 
  mutate(genus=taxa.1[ASV,'Genus'])

g <- df.plot.2 %>% 
  mutate(diff=factor(-sign(diff))) %>% 
  mutate(fc=case_when(diff==-1~fc_bp, diff==1~fc_ps)) %>% 
  filter(fc > 0.1) %>% 
  arrange(desc(fc)) %>% 
  mutate(label=factor(label, levels = unique(label))) %>% 
  mutate(Group=factor(Group, levels = c('C', 'Bp', 'Ps'))) %>% 
  ggplot(aes(x=label, y=ab, fill=Group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(col=Group), fill='white', pch=21,
              position = position_jitterdodge()) + 
  facet_grid(~diff+family, space = 'free', scales = 'free') + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)) +
  xlab('Amplicon sequence variants') + 
  ylab('log10(relative abundance)') + 
  scale_fill_manual(values=c('C'='#0024b1', 'Bp'='#4a7a31', 'Ps'='#da4432'),
                    labels=c('C', 'CC4', 'Ps')) + 
  scale_colour_manual(values=c('C'='#0024b1', 'Bp'='#4a7a31', 'Ps'='#da4432'),
                      labels=c('C', 'CC4', 'Ps'))
ggsave(g, filename = here('figures','ASV_family_plot.pdf'),
       width = 9, height = 4, useDingbats=FALSE)