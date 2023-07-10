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


# Ext. Data Fig 8c
ggplot(enrichments_dtop10) + 
  geom_point(aes(x = FeatureF, y = `% covered in GAM`, size = `covered in GAM`), alpha = 0.9, color = cGAM) +
  geom_point(aes(x = FeatureF, y = `% covered in Hi-C`, size = `covered in Hi-C`), alpha = 0.9, color = cHiC) +
  scale_size(range = c(1, 30))+  # Adjust the range of points size
  geom_text(aes(label=sprintf("%.3f", `Mean discriminatory power (Gini)`), x=FeatureF), y=0)+
  mytheme+
  xlab("Feature pair") +
  ylab("Percentage of contacts covered in 10% sets") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  ylim(c(0,0.6))+ 
  coord_flip() 

