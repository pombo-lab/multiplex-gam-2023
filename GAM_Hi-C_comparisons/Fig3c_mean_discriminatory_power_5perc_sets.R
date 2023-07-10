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
## 5% sets

enrichments = read_tsv("./data/ST8 - Feature pair enrichments.tsv") %>% 
  left_join(read_tsv("./data/feature_enrichments/sampled_005_pairfeatures.tsv")) #add frequencies from randomized background data

enrichments = enrichments  %>% 
  mutate(dselected = rank( desc(`Mean discriminatory power (Gini)`)), selected = rank(desc(`Mean discriminatory power (Gini)`))<11 | `HiC rank enrichment`<4)

enrichments_dtop10 = enrichments %>% 
  filter(selected) %>%
  mutate(FeatureF = fct_reorder(`Feature pair`, (`Mean discriminatory power (Gini)`)))

# Fig 3c
ggplot(enrichments_dtop10) + 
  geom_point(aes(x = FeatureF, y = `% covered in GAM`, size = `covered in GAM`), alpha = 0.9, color = cGAM) +
  geom_point(aes(x = FeatureF, y = `% covered in Hi-C`, size = `covered in Hi-C`), alpha = 0.9, color = cHiC) +
  scale_size(range = c(1, 30))+  # Adjust the range of points size
  geom_text(aes(label=sprintf("%.3f", `Mean discriminatory power (Gini)`), x=FeatureF), y=0)+
  mytheme+
  xlab("Feature pair") +
  ylab("Percentage of contacts covered in 5% contact sets") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  ylim(c(0,0.6))+ 
  #theme(legend.position="bottom")+
  coord_flip() 
