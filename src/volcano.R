# ##############################################################################
#
## Check the mOTUs profiles from the CRC meta-analysis for CTR-enriched bacteria
#
# ##############################################################################

# packages
library("tidyverse")
library("here")
library("ggrepel")
library("coin")
library("pROC")
library('progress')

# ##############################################################################
# read in data
data.location <- 'https://www.embl.de/download/zeller/crc_meta/'
# features
feat.ab <- read.table(paste0(data.location, 
                              "mOTU_profiles/crc_meta_g2_l75_motus2.0.0.motus"),
                       sep='\t', stringsAsFactors = FALSE, check.names = FALSE, 
                       row.names = 1, header = TRUE, quote = '', 
                       comment.char = '', skip=2)
feat.rel <- prop.table(as.matrix(feat.ab), 2)

# metadata
meta.all <- read_tsv(paste0(data.location, 'metadata/meta_all.tsv')) %>% 
  filter(Sample_ID %in% colnames(feat.rel)) %>% 
  filter(Group %in% c('CTR', 'CRC'))
feat.rel <- feat.rel[,meta.all$Sample_ID]

# taxonomy
motus.taxonomy <- read_tsv(paste0(data.location, 'files/motus_2.1.tax'))

# calculate mean abundance in controls
mean.log.ab <- rowMeans(log10(feat.rel[,meta.all %>% 
                                         filter(Group=='CTR') %>% 
                                         pull(Sample_ID)] + 1e-05)) %>% 
  enframe(name='species', value='mean.log.ab')

# read in results from differential abundance testing
all.res <- read_tsv(paste0(data.location, 'results/all_measures.tsv')) %>% 
  left_join(mean.log.ab) %>% 
  mutate(mOTU_ID=str_extract(species, '(ref|meta)_mOTU_v2_[0-9]*')) %>% 
  left_join(motus.taxonomy) %>% 
  # filter(!str_detect(order, '^NA')) %>% 
  mutate(p.adj=p.adj+1e-14) %>% 
  mutate(g.col=ifelse(order=='186802 Clostridiales', 
                      'Clostridiales', 'other')) %>% 
  mutate(name=case_when(mOTU_ID=='ref_mOTU_v2_1427'~'R. intestinalis',
                        mOTU_ID=='ref_mOTU_v2_4207'~'E. hallii',
                        mOTU_ID=='ref_mOTU_v2_4211'~'F. prausnitzii',
                        mOTU_ID=='ref_mOTU_v2_1381'~'A. caccae',
                        mOTU_ID=='ref_mOTU_v2_4614'~'P. stomatis',
                        TRUE~'')) %>% 
  mutate(highlight=ifelse(name=='', FALSE, TRUE)) %>% 
  mutate(alpha.highlight=case_when(highlight~TRUE,
                                   p.adj < 0.005 & fc < -0.25~TRUE,
                                   TRUE~FALSE)) %>% 
  mutate(g.col.2=case_when(g.col=='other'~'other',
                           genus=='841 Roseburia'~'Roseburia',
                           genus=='216851 Faecalibacterium'~'Faecalibacterium',
                           genus=='1730 Eubacterium'~'Eubacterium',
                           genus=='207244 Anaerostipes'~'Anaerostipes',
                           TRUE~'other Clostridiales')) %>% 
  mutate(g.col.2=factor(g.col.2, 
                        levels = c('Anaerostipes', 'Eubacterium', 
                                   'Faecalibacterium', 'Roseburia',
                                   'other Clostridiales', 'other')))

df.sp <- all.res %>% 
  filter(highlight)

g <- all.res %>% 
  ggplot(aes(x=fc, y=-log10(p.adj), col=g.col, size=mean.log.ab)) + 
  geom_vline(xintercept = 0, colour='lightgrey', lty=3) + 
  geom_vline(xintercept = -0.25, colour='lightgrey', lty=2) + 
  geom_hline(yintercept = -log10(0.005), 
             colour='lightgrey') + 
  geom_point(alpha=0.4, stroke=0.3) + 
  xlab('Generalized fold change') + 
  ylab('-log10(Adj. P-value)') + 
  xlim(-0.65, 0.65) + 
  theme_bw() + 
  scale_colour_manual(values=c('#307FE2', 'darkgrey'), name='',
                      label=c('Clostridiales', 'other')) +
  geom_point(data=df.sp, shape=21, col='black', fill='#307FE2') +
  geom_text_repel(aes(label=name), show.legend = FALSE, 
                  colour='black', size=4) + 
  scale_size(range=c(1, 6), trans='exp')


ggsave(g, width = 6, height = 5, useDingbats=FALSE,
       filename = here('figures', "volcano.pdf"))

# ##############################################################################
# test differences on Order-level
rownames(feat.rel) <- str_extract(rownames(feat.rel), 
                                  pattern = '(ref|meta)_mOTU_v2_[0-9]*')

tax.table <- motus.taxonomy %>% 
  filter(mOTU_ID %in% rownames(feat.rel)) %>% 
  filter(mOTU_name!='-1') %>% 
  select(-mOTU_name) 

# combine on order level
o.unique <- unique(tax.table$order)
feat.order <- matrix(NA, nrow = length(o.unique), ncol=ncol(feat.rel),
                     dimnames = list(o.unique, colnames(feat.rel)))
for (o in o.unique){
  feat.order[o,] <- colSums(feat.rel[tax.table %>% 
                                       filter(order==o) %>% 
                                       pull(mOTU_ID),,drop=FALSE])
}
feat.order <- feat.order[rowMeans(feat.order!=0) > 0.05,]


# test all orders
meta.test <- meta.all %>% 
  filter(!is.na(Sampling_rel_to_colonoscopy)) %>%
  mutate(block=ifelse(Study!='CN-CRC', Study, 
                      paste0(Study, '_', Sampling_rel_to_colonoscopy)))
studies <- unique(meta.test$Study)
feat.test <- feat.order[,meta.test$Sample_ID]

p.val <- matrix(NA, nrow=nrow(feat.test), ncol=length(studies)+1, 
                dimnames=list(row.names(feat.test), c(studies, 'all')))
fc <- p.val
aucs.mat <- p.val
aucs.all  <- vector('list', nrow(feat.test))

log.n0 <- 1e-05

cat("Calculating effect size for every feature...\n")
pb <- progress_bar$new(total = nrow(feat.test))

# caluclate wilcoxon test and effect size for each feature and study
for (f in row.names(feat.test)) {
  
  # for each study
  for (s in studies) {
    
    x <- feat.test[f, meta.test %>% filter(Study==s) %>% 
                     filter(Group=='CRC') %>% pull(Sample_ID)]
    y <- feat.test[f, meta.test %>% filter(Study==s) %>% 
                     filter(Group=='CTR') %>% pull(Sample_ID)]
    
    # Wilcoxon
    p.val[f,s] <- wilcox.test(x, y, exact=FALSE)$p.value
    
    # AUC
    aucs.all[[f]][[s]] <- c(roc(controls=y, cases=x, 
                                direction='<', ci=TRUE, auc=TRUE)$ci)
    aucs.mat[f,s] <- c(roc(controls=y, cases=x, 
                           direction='<', ci=TRUE, auc=TRUE)$ci)[2]
    
    # FC
    q.p <- quantile(log10(x+log.n0), probs=seq(.1, .9, .05))
    q.n <- quantile(log10(y+log.n0), probs=seq(.1, .9, .05))
    fc[f, s] <- sum(q.p - q.n)/length(q.p)
  }
  
  # calculate effect size for all studies combined
  # Wilcoxon + blocking factor
  d <- data.frame(y=feat.test[f,], 
                  x=as.factor(meta.test$Group), block=as.factor(meta.test$block))
  p.val[f,'all'] <- pvalue(wilcox_test(y ~ x | block, data=d))
  # other metrics
  x <- feat.test[f, meta.test %>% filter(Group=='CRC') %>% pull(Sample_ID)]
  y <- feat.test[f, meta.test %>% filter(Group=='CTR') %>% pull(Sample_ID)]
  # FC
  fc[f, 'all'] <- mean(fc[f, studies])
  # AUC
  aucs.mat[f,'all'] <- c(roc(controls=y, cases=x, 
                             direction='<', ci=TRUE, auc=TRUE)$ci)[2]
  
  # progressbar
  pb$tick()
}


# multiple hypothesis correction
p.adj <- data.frame(apply(p.val, MARGIN=2, FUN=p.adjust, method='fdr'),
                    check.names = FALSE)
df.all <- tibble(species=rownames(feat.test),
                 fc=fc[,'all'],
                 auroc=aucs.mat[,'all'],
                 p.val=p.val[,'all'],
                 p.adj=p.adj[,'all'])

df.all %>% 
  mutate(p.adj=p.adj+1e-08) %>% 
  mutate(l=case_when(p.adj < 1e-06~species, TRUE~'')) %>% 
  ggplot(aes(x=fc, y=-log10(p.adj))) + 
  geom_point() + 
  geom_text_repel(aes(label=l)) + 
  xlab('Generalized fold change') + 
  ylab('-log10(Adj. p-value)')

# ##############################################################################
# boxplots of selected species

df.plot <- meta.all %>% 
  select(Sample_ID, Group, AJCC_stage) %>% 
  mutate(`R. intestinalis`=feat.rel['ref_mOTU_v2_1427',Sample_ID]) %>% 
  mutate(`E. hallii`=feat.rel['ref_mOTU_v2_4207',Sample_ID]) %>% 
  mutate(`F. prausnitzii`=feat.rel['ref_mOTU_v2_4211',Sample_ID]) %>% 
  mutate(`A. caccae`=feat.rel['ref_mOTU_v2_1381',Sample_ID]) %>% 
  mutate(Clostridiales=feat.order['186802 Clostridiales',Sample_ID]) %>% 
  mutate_if(is_double, .funs = function(x){log10(x+1e-05)})

g1 <- df.plot %>% 
  pivot_longer(cols=c(`R. intestinalis`, `E. hallii`, 
                      `F. prausnitzii`, `A. caccae`)) %>% 
  mutate(name=factor(name, 
                     levels = c('R. intestinalis', 'E. hallii', 
                                'F. prausnitzii', 'A. caccae'))) %>% 
  ggplot(aes(x=name, y=value, fill=Group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(col=Group), 
              position = position_jitterdodge(jitter.width = 0.1)) + 
  coord_flip() + 
  facet_wrap(~name, scales='free_y', ncol = 1) + 
  xlab('') + ylab('log. rel. ab.') + 
  scale_fill_manual(values=c('#E4004675', alpha('darkgrey', alpha=0.75)), 
                    guide='none') + 
  scale_colour_manual(values=c('#E40046', 'darkgrey'), guide='none') + 
  theme_bw() +
  theme(strip.text = element_blank(), 
        panel.grid.major.x = element_line(colour='lightgrey'),
        axis.ticks.y=element_blank(), 
        axis.text.y = element_text(angle = 90, hjust=1))

g2 <- df.plot %>% 
  ggplot(aes(x=1, y=Clostridiales, fill=Group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(col=Group), 
              position = position_jitterdodge(jitter.width = 0.1)) + 
  coord_flip() + 
  xlab('') + ylab('log. rel. ab.') + 
  scale_fill_manual(values=c('#E4004675', alpha('darkgrey', alpha=0.75)), 
                    guide='none') +
  scale_colour_manual(values=c('#E40046', 'darkgrey'), guide='none') +
  theme_bw() +
  theme(strip.text = element_blank(), 
        panel.grid.major.x = element_line(colour='lightgrey'),
        axis.ticks.y=element_blank(), 
        axis.text.y = element_text(angle = 90, hjust=1))


# along the stages?
g <- df.plot %>% 
  pivot_longer(cols=c(`R. intestinalis`, `E. hallii`, 
                      `F. prausnitzii`, `A. caccae`)) %>% 
  mutate(name=factor(name, 
                     levels = c('R. intestinalis', 'E. hallii', 
                                'F. prausnitzii', 'A. caccae'))) %>% 
  mutate(AJCC_stage=case_when(Group=='CTR'~'CTR', 
                              AJCC_stage %in% c('0', 'I', 'II')~'0/I/II',
                              AJCC_stage %in% c('IV', 'III')~'III/IV')) %>% 
  filter(!is.na(AJCC_stage)) %>% 
  mutate(AJCC_stage=factor(AJCC_stage, 
                           levels = rev(c("CTR", '0/I/II', 'III/IV')))) %>% 
  ggplot(aes(x=name, y=value, fill=AJCC_stage)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(col=AJCC_stage), 
              position = position_jitterdodge(jitter.width = 0.1)) + 
  coord_flip() + 
  facet_wrap(~name, scales='free_y', ncol = 1) + 
  xlab('') + ylab('log. rel. ab.') + 
  scale_fill_manual(values=alpha(rev(c('darkgrey', '#E58F9E', '#E40046')), 0.75),
                    guide='none') +
  scale_colour_manual(values=rev(c('darkgrey', '#E58F9E', '#E40046'))) + 
  theme_bw() +
  theme(strip.text = element_blank(), 
        panel.grid.major.x = element_line(colour='lightgrey'),
        axis.ticks.y=element_blank(), 
        axis.text.y = element_text(angle = 90, hjust=1))
ggsave(g, width = 85, height = 90, units = 'mm',
       filename = here('figures', 'volcano_boxes.pdf'))

g <- df.plot %>% 
  mutate(AJCC_stage=case_when(Group=='CTR'~'CTR', 
                              AJCC_stage %in% c('0', 'I', 'II')~'0/I/II',
                              AJCC_stage %in% c('IV', 'III')~'III/IV')) %>% 
  filter(!is.na(AJCC_stage)) %>% 
  mutate(AJCC_stage=factor(AJCC_stage, 
                           levels = rev(c("CTR", '0/I/II', 'III/IV')))) %>% 
  ggplot(aes(x=1, y=Clostridiales, fill=AJCC_stage)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(col=AJCC_stage), 
              position = position_jitterdodge(jitter.width = 0.1)) + 
  coord_flip() + 
  xlab('') + ylab('log. rel. ab.') +
  scale_fill_manual(values=alpha(rev(c('darkgrey', '#E58F9E', '#E40046')), 0.75),
                    guide='none') +
  scale_colour_manual(values=rev(c('darkgrey', '#E58F9E', '#E40046'))) + 
  theme_bw() +
  theme(strip.text = element_blank(), 
        panel.grid.major.x = element_line(colour='lightgrey'),
        axis.ticks.y=element_blank(), 
        axis.text.y = element_text(angle = 90, hjust=1))
ggsave(g, width = 85, height = 30, units = 'mm',
       filename = here('figures', 'volcano_boxes_Clostridiodes.pdf'))

# ##############################################################################
# Other boxplots for the supplementary figure
g <- df.plot %>% 
  pivot_longer(cols=c(`R. intestinalis`, `E. hallii`, 
                      `F. prausnitzii`, `A. caccae`)) %>% 
  mutate(name=factor(name, 
                     levels = c('R. intestinalis', 'E. hallii', 
                                'F. prausnitzii', 'A. caccae'))) %>% 
  left_join(meta.all %>% select(Sample_ID,Study), by='Sample_ID') %>% 
  mutate(Study=case_when(Study=='AT-CRC'~"Feng 2015",
                         Study=='CN-CRC'~'Yu 2017',
                         Study=='FR-CRC'~'Zeller 2014',
                         Study=='US-CRC'~'Vogtmann 2016',
                         Study=='DE-CRC'~'Wirbel 2019')) %>% 
  ggplot(aes(x=name, y=value, fill=Group)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(col=Group), 
              position = position_jitterdodge(jitter.width = 0.1)) + 
  facet_grid(Study~name, scales = 'free') +
  xlab('') + ylab('log10(Relative abundance)') + 
  scale_fill_manual(values=c('#E4004675', alpha('darkgrey', alpha=0.75)), 
                    guide='none') +
  scale_colour_manual(values=c('#E40046', 'darkgrey'), guide='none') +
  theme_bw() +
  theme(axis.ticks.x=element_blank(), 
        strip.text.x = element_blank(),
        panel.grid.major.y = element_line(colour='lightgrey', linetype=3))
ggsave(g, width = 6, height = 8, useDingbats=FALSE,
       filename = here('figures', "volcano_boxes_study.pdf"))

# ##############################################################################
# Top control-enriched species + boxplots

# which species are enriched in controls?
top10.ctr <- all.res %>% 
  filter(fc < 0) %>% 
  filter(!str_detect(species, 'meta_mOTU')) %>% 
  arrange(p.val) %>% 
  slice_head(n=10) %>% 
  mutate(species=factor(species, levels = rev(species)))

g.pval <- top10.ctr %>% 
  ggplot(aes(x=species, y=-log10(p.adj))) +
    geom_bar(stat='identity') + 
    coord_flip() + 
    theme_bw() + 
    theme(panel.grid = element_blank()) + 
    ylab('-log10(q-value)') + xlab('')
ggsave(g.pval, width = 5, height = 5, useDingbats=FALSE,
       filename = here('figures', 'top10_species_ctr_p_val.pdf'))

g.ab <- prop.table(as.matrix(feat.ab), 2)[as.character(top10.ctr$species),] %>% 
  as_tibble(rownames='species') %>% 
  pivot_longer(-species, names_to ='Sample_ID', values_to='rel.ab') %>% 
  mutate(log.rel.ab=log10(rel.ab+1e-05)) %>% 
  right_join(meta.all) %>% 
  mutate(species=factor(species, levels = levels(top10.ctr$species))) %>% 
  mutate(Study=str_remove(Study, '-CRC')) %>% 
  mutate(Study=factor(Study, levels = c('FR', 'AT', 'CN', 'US', 'DE'))) %>% 
  ggplot(aes(x=species, y=log.rel.ab, fill=Group)) + 
    geom_boxplot(outlier.shape = NA) + 
    geom_jitter(position = position_jitterdodge(jitter.width = 0.1), 
                alpha=0.3) +
    facet_grid(~Study) + 
    coord_flip() + 
    xlab('') + ylab('log10(relative abundance)') + 
    theme_bw() + 
    theme(panel.grid = element_blank(), axis.text.y = element_blank(), 
          axis.ticks.y = element_blank(), strip.background = element_blank()) + 
    scale_fill_manual(values=c('#d9544d90', "#C2C2C290"))
ggsave(g.ab, width = 12, height = 6, useDingbats=FALSE,
       filename = here('figures', 'top10_species_abundance.pdf'))  
