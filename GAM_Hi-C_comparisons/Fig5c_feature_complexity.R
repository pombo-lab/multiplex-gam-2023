library(ggplot2)
library(tidyverse)
library(ggbeeswarm)

options(stringsAsFactors=F)
options(scipen=99999)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#setwd("/data/pombo/christoph/1808_GAM_mESC_diffcontacts/scripts/publication/")


mytheme = theme_bw()+ theme(text = element_text(size=16))


featuretable=read_tsv("./data/featuretable_clean.tsv")

###################
# mean complexity for given features

flist.figgrpset=c(
  "CTCF",
  "p300",
  "Oct4",
  "RNAPII-S5p",
  "RNAPII-S7p",
  "RNAPII-S2p",
  "Enhancer",
  "Super-enhancers (SE)",
  "H4K20me3",
  "H3K9me3")

cat.colors=c("HiC specific cat"="#4a739f",
             "Strong-and-common cat"="#ff990cdd",
             "GAM specific cat"="#c83737",
             "Unspecific, many"="#5788af",
             "Unspecific, few"="#70a4ea")


complexcounts = tibble()
for(feat in flist.figgrpset){
  print(feat)
  complexcounts = complexcounts %>% bind_rows(
    featuretable  %>%  select(chrom, start, end, !!feat, complexity_bincount, "GAM specific cat", "HiC specific cat", "Strong-and-common cat") %>% 
      mutate(feature = !!feat, complexity_bincount=as.numeric(complexity_bincount)) %>%
      filter(!is.na(!!sym(feat)), !!sym(feat)==T, !complexity_bincount==".") %>%
      select(-!!feat)
  )
}


complexcounts.long = complexcounts %>% 
  mutate("genomewide" = T)%>% pivot_longer(cols=c("GAM specific cat", "HiC specific cat", "Strong-and-common cat", "genomewide"), names_to ="bin_class", values_to ="bin_member") %>% 
  filter(bin_member==T)

complexcounts.u = complexcounts.long %>% distinct(feature, chrom, start, complexity_bincount, bin_class, .keep_all=T)


#get feature order by mean
feature.ranks = complexcounts.u  %>% filter(bin_class=="genomewide") %>% 
  group_by(feature) %>%
  dplyr::summarise(fmean.compl.count = mean(complexity_bincount, na.rm = T)) %>% 
  dplyr::arrange(fmean.compl.count)

# Fig 5c
ggplot(complexcounts.u, aes(x=factor(feature, levels=feature.ranks$feature), y=complexity_bincount, color=bin_class))+
  geom_quasirandom(size=0.5)+
  geom_violin(fill=NA, draw_quantiles = c(0.25,0.5,0.75), color="black", scale="count")+
  geom_text(data= complexcounts.u %>% count(feature, bin_class), aes(x=feature, y=4500, label=n), color="black")+
  mytheme+
  facet_grid(factor(bin_class, levels=c("genomewide", "HiC specific cat",  "Strong-and-common cat", "GAM specific cat"))~.)+
  xlab("")+ylab("")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+ 
  theme(legend.position = "none")+
  scale_fill_manual(values = cat.colors)+
  scale_color_manual(values = cat.colors)+
  ylim(c(0,4700))


