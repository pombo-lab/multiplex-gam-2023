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


# number of observations to adjust for set sizes
NUMRAND=250173
NUMGAM=265167
NUMHIC=231165

enrichments = enrichments %>% mutate(transp = as.factor(ntile(`Mean discriminatory power (Gini)`,10)),
                                     `Feature pair` = fct_reorder(`Feature pair`, -(GAM_random1))) 


a= ggplot(enrichments, aes(x=`Feature pair`))+
  geom_point(aes(y=`% covered in GAM`), color=cGAM)+
  geom_point(aes(y=`% covered in Hi-C`), color=cHiC)+
  
  geom_point(aes(y=GAM_random1/NUMGAM), color= "#dddddd", alpha=0.4)+
  geom_point(aes(y=GAM_random2/NUMGAM), color= "#dddddd", alpha=0.4)+
  geom_point(aes(y=GAM_random3/NUMGAM), color= "#dddddd", alpha=0.4)+
  
  geom_point(aes(y=`Hi-C_random1`/NUMHIC), color="#dddddd", alpha=0.4)+
  geom_point(aes(y=`Hi-C_random2`/NUMHIC), color="#dddddd", alpha=0.4)+
  geom_point(aes(y=`Hi-C_random3`/NUMHIC), color="#dddddd", alpha=0.4)+
  mytheme+
  theme(axis.text.x = element_blank())+
  scale_color_manual(values = bcolors)+
  ylab("Percentage of contacts covered") + xlab("") +
  geom_text_repel(data=enrichments %>% filter( selected), aes(label=`Feature pair`, y=pmax(`% covered in GAM`, `% covered in Hi-C`)), min.segment.length=0) +
  theme(
    plot.background = element_blank(), #element_rect(fill = "transparent",colour = NA),
    panel.grid=element_blank(),
    panel.background=element_blank(),
  )


b= ggplot(enrichments, aes(x=`Feature pair`, y=0, fill=transp))+
   geom_tile()+
   mytheme+
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
   scale_fill_manual(values = bcolors) +
   theme(axis.text.y = element_blank())+
   ylab("Gini") +
  xlab("Feature pair") +
  theme(
  plot.background = element_blank(),
  panel.grid=element_blank(),
  panel.background=element_blank(),
  axis.text.x=element_text(size=unit(0,"npc")))  

# Fig 3b
(a/b)+ plot_layout(heights = c(19, 1))

