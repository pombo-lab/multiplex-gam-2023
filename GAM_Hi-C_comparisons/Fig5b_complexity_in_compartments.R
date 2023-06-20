library(ggplot2)
library(tidyverse)
library(ggbeeswarm)

options(stringsAsFactors=F)
options(scipen=99999)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#setwd("/data/pombo/christoph/1808_GAM_mESC_diffcontacts/scripts/publication/")


mytheme = theme_bw()+ theme(text = element_text(size=16))

cat.colors=c("Mostly Hi-C specific"="#4a739f",
             "Most Strong-and-common"="#ff990cdd",
             "Mostly GAM specific"="#c83737",
             "Unspecific, many"="#5788af",
             "Unspecific, few"="#70a4ea")
cat.order = c("Mostly Hi-C specific", "Most Strong-and-common", "Mostly GAM specific")



featuretable=read_tsv("./data/featuretable_clean.tsv")


####################
# mean complexity at preferred regions

compl.gam = featuretable %>% dplyr::filter(`GAM specific cat`==T) %>% dplyr::select(chrom, start, `complexity_bincount`, `GAM eigenvalues`) %>% mutate(ds = "Mostly GAM specific", comp = ifelse(`GAM eigenvalues`>0, "A", "B"))
compl.hic = featuretable %>% dplyr::filter(`HiC specific cat`==T) %>% dplyr::select(chrom, start, `complexity_bincount`, `GAM eigenvalues`) %>% mutate(ds = "Mostly Hi-C specific", comp = ifelse(`GAM eigenvalues`>0, "A", "B"))
compl.sac = featuretable %>% dplyr::filter(`Strong-and-common cat`==T) %>% dplyr::select(chrom, start, `complexity_bincount`, `GAM eigenvalues`) %>% mutate(ds = "Most Strong-and-common", comp = ifelse(`GAM eigenvalues`>0, "A", "B"))
compl.data = dplyr::bind_rows(compl.hic, compl.sac, compl.gam)

# Fig 5b
ggplot(compl.data, aes(x=factor(ds, levels=cat.order), y=as.numeric(complexity_bincount), fill=ds))+
  geom_violin(scale = "count", draw_quantiles = c(0.25,.5,0.75))+
  #geom_jitter(height = 0)+
  mytheme+
  facet_grid(.~comp)+
  ylab("Complexity (Mean PiABC)")+
  xlab("Bin assignment")+
  scale_fill_manual(values = cat.colors)+ 
  theme(legend.position = "none")
