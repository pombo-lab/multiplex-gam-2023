library(tidyverse)
library(ggplot2)
library(patchwork)

options(stringsAsFactors=F)
options(scipen=99999)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


mytheme = theme_bw()+ theme(text = element_text(size=16))


cat.colors=c("HiC specific cat"="#4a739f",
             "Strong-and-common cat"="#ff990cdd",
             "GAM specific cat"="#c83737",
             "genomewide"="#dddddd")

cat.order = names(cat.colors)


### transform genelist to 40kb binned version
#genes.positions.file="./data/Ferrai2017_ESC_genes_mm9.bed"
#system(sprintf("bedtools intersect -wo -a '../40000/windows.bed' -b '%s' > ../data/Ferrai2017_ESC_genes_mm9_40000.tsv", genes.positions.file))
genoverlap.win = read_tsv("./data/genes/Ferrai2017_ESC_genes_mm9_40000.tsv", col_names = c("chrom","start","end","g.chrom","g.start","g.end", "clusterID", "day0_ID", "cov"))

#get gene activity as TPM counts
genes.activitylist=read.table("./data/genes/Ferrai2017_msb177754-sup-0004-datasetev2.tsv", header=T, sep="\t", quote="\"", stringsAsFactors = F) %>%
  dplyr::mutate(tpm_quant = ntile(log(ESC_TPM), 10)) # add quantile for TPMs here - need to do before adding duplicates by the LO-join on the binned windows

#annotate mm9 gene positions with activity data
geneslist = merge(genoverlap.win, genes.activitylist, by.x=c("clusterID", "day0_ID", "chrom"), by.y=c("clusterID", "ESC_ID", "chrom")) %>% 
  dplyr::select(chrom, start, end, ESC_TPM, tpm_quant, geneSymbol, description)

table(geneslist$tpm_quant)



featuretable=read_tsv("./data/featuretable_clean.tsv")

window.categories = featuretable %>% 
  tibble() %>%
  dplyr::select(1:3, ends_with("cat")) %>% 
  mutate(genomewide=T) %>%
  arrange(chrom, start) %>% #add group with all bins
  pivot_longer(cols = all_of(cat.order)) %>% 
  filter(!is.na(value), value==T) %>% 
  left_join(geneslist,  by=c("chrom", "start", "end"), multiple = "all")

window.categories %>% count(name)


###################
# Fig 4c - barplot with gene counts
plot.genecounts = window.categories %>% filter(name!="genomewide") %>% distinct(name, ESC_TPM, geneSymbol) %>%  group_by(name) %>% summarize(genecount = n()) %>%
  ggplot(aes(x=factor(name, levels=cat.order), y=genecount, fill=factor(name, levels=cat.order)))+
  geom_col(position=position_dodge(.7), width=0.4)+
  mytheme+
  ggtitle("Gene count")+
  xlab("")+ylab("")+
  #theme(aspect.ratio = 5/4)+
  scale_fill_manual(values = cat.colors)+ 
  theme(legend.position = "bottom")


###################
# Fig 4d - gene TPM values
plot.expression = window.categories  %>% distinct(name, ESC_TPM, geneSymbol) %>% 
  ggplot(aes(x=factor(name, levels=cat.order), y=log2(ESC_TPM), fill=factor(name, levels=cat.order)))+
  geom_violin(draw_quantiles = c(0.25,0.5,0.75))+
  mytheme+
  ggtitle("Gene count")+
  xlab("")+ylab("")+
  #theme(aspect.ratio = 5/4)+
  scale_fill_manual(values = cat.colors)+ 
  theme(legend.position = "bottom")



###################
# Fig 4e - LAD assignment
window.lads = featuretable %>% 
  group_by(chrom) %>% 
  dplyr::select(1:3, matches("*LAD*"), matches("*cat*")) %>% 
  dplyr::select(1:7) %>% 
  mutate(genomewide=T) #add group with all bins

laddata = window.lads %>% 
  pivot_longer(cols = cat.order) %>% 
  filter(!is.na(value), value==T) %>% 
  group_by(name) %>% 
  summarize(wincount = n(), ladcount =sum(`MESC-LADs`>0), ladratio = sum(`MESC-LADs`>0)/n(), 
            lbl = sprintf("%s/%s", ladcount, wincount))

plot.lads = ggplot(laddata, aes(x=factor(name, levels=cat.order),y=ladratio, fill=name, label=lbl))+
  geom_bar(stat="identity")+
  geom_text(vjust=-1)+
  mytheme+
  geom_hline(yintercept=(laddata %>% dplyr::filter(name=="genomewide"))$ladratio, linetype=3, alpha=0.5)+
  #theme(aspect.ratio = 5/4)+
  ylab("Frequency")+
  xlab("Bin assignment")+
  ylim(c(0,1))+
  scale_fill_manual(values = cat.colors)+ 
  theme(legend.position = "none")




###################
# Fig 4f - feature presence
cat.order.3g = names(cat.colors)[1:3]

flist.figgrpset=c("CTCF",
                  "p300",
                  "Oct4",
                  "RNAPII-S5p",
                  "RNAPII-S7p",
                  "RNAPII-S2p",
                  "Enhancer",
                  "Super-enhancers (SE)",
                  "H4K20me3",
                  "H3K9me3"
                  )

featuredata.groups = featuretable %>% mutate(genomewide=T)  %>%  #add group with all bins
  select(1:3, all_of(cat.order.3g), all_of(flist.figgrpset)) %>% 
  pivot_longer(cols = all_of(cat.order.3g)) %>% 
  filter(!is.na(value), value==T) %>% 
  pivot_longer(-c("chrom", "start", "end", "name", "value"), names_to="features", values_to="feature_presence") %>% 
  filter(feature_presence>0) 
  
plot.feats = ggplot(featuredata.groups %>% filter(name!="genomewide"), aes(x=factor(name, levels=cat.order), fill=factor(name, levels=cat.order)))+
  geom_bar(position=position_dodge(.7), width=0.4)+
  facet_wrap(factor(features, levels=flist.figgrpset) ~., scales = "free_y", ncol=4)+ 
  mytheme+
  xlab("")+ylab("")+
  scale_fill_manual(values = cat.colors)+ 
  theme(legend.position = "None")




###################
# Fig 4f - compartment assignment
flist.figgrpset.c=c("HiC-Acompartment",
                  "HiC-Bcompartment",
                  "GAM-Acompartment",
                  "GAM-Bcompartment")

featuredata.groups.c = featuretable  %>%  #add group with all bins
  select(1:3, all_of(cat.order.3g), all_of(flist.figgrpset.c)) %>% 
  pivot_longer(cols = all_of(cat.order.3g)) %>% 
  filter(!is.na(value), value==T) %>% 
  pivot_longer(-c("chrom", "start", "end", "name", "value"), names_to="features", values_to="feature_presence") %>% 
  filter(feature_presence>0) 

plot.comps = ggplot(featuredata.groups.c, aes(x=factor(name, levels=cat.order), fill=factor(name, levels=cat.order)))+
  geom_bar(position=position_dodge(.7), width=0.4)+
  facet_wrap(features ~., scales = "free_y", ncol=2)+ 
  mytheme+
  xlab("")+ylab("")+
  ylim(c(0,6000))+
  scale_fill_manual(values = cat.colors)+ 
  theme(legend.position = "None")



((plot.genecounts / plot.expression / plot.lads) | plot.feats | plot.comps) & theme(axis.text.x = element_text(angle = 90, hjust = 1))


