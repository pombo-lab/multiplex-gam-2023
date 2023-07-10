library(dplyr)
library(tidyr)
library(readr)
library(forcats)
library(ggplot2)
library(ggrepel)
library("RColorBrewer")
library(patchwork)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

options(scipen=9)
options(digits=2)

cGAM="#c83737"
cHiC="#4a739f"


bcolors = colorRampPalette(brewer.pal(9, "Blues"))(10)

mytheme = 
  theme_bw() + 
  theme(legend.position = "right")+
  theme(text = element_text(size=14))


catcolors.spec=c("strong"="#828282",
                 "common"="#828282",
                 "GAM"="#c83737",
                 "HiC"="#4a739f",
                 "GAM.strong"="#FF6900",
                 "HiC.strong"="#779AFF",
                 "strongcommon"="#ff990cdd"
)


###################
## 10% sets

enrichments = read_tsv("./data/ST9 - Feature pair enrichments.tsv") %>%
  dplyr::mutate(
    Enrichment_in_GAM=`% covered in GAM`/`% covered in Hi-C`,
    Enrichment_in_HiC=1/Enrichment_in_GAM
  ) 

enrichments_dtop10 = enrichments  %>% 
  mutate(dselected = rank( desc(`Mean discriminatory power (Gini)`)), 
         selected = rank(desc(`Mean discriminatory power (Gini)`))<11 | `HiC rank enrichment`<4) %>% 
  filter(selected) %>%
  mutate(FeatureF = fct_reorder(`Feature pair`, (`Mean discriminatory power (Gini)`)))


# Ext. Data Fig 8b
ggplot(enrichments, aes(x=Enrichment_in_GAM))+
  geom_segment(aes(color=`Mean discriminatory power (Gini)`, x=Enrichment_in_GAM,xend=Enrichment_in_GAM, y=`% covered in GAM`,  yend=`% covered in Hi-C`), linewidth=1) +
  geom_text_repel(data=enrichments_dtop10, aes(label=`Feature pair`, y=`% covered in GAM`), min.segment.length=0) +
  geom_point(aes(y=`% covered in GAM`), color=cHiC)+
  geom_point(aes( y=`% covered in Hi-C`), color=cGAM)+
  ylab("Feature frequency")+
  xlab("Enrichment in GAM/Hi-C for 10% contact sets")+
  geom_vline(xintercept = 1)+
  mytheme + 
  scale_color_gradient(low = "yellow", high = "red")

