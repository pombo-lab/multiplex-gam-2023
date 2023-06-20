library(tidyverse)
library(ComplexUpset)
library(RColorBrewer)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

options(scipen=9)
options(digits=2)

colors=brewer.pal(10,"RdBu")


# get fixed order for set classes based on observations in background
refclasses = read_tsv("./data/upset_data/random_005_1.upset.tsv") %>% pivot_longer(-"Window pair") %>%
  group_by(name) %>% summarise(n=sum(value)) %>% arrange(n)
refclasses = refclasses$name


getuplots = function(fn, cl, iss=15){
  contacts_all = read_tsv(fn)

  #filter empty contacts
  contacts_all[refclasses] = contacts_all[refclasses] > 0 # turn into binary table
  contacts_all = contacts_all %>% filter_at(refclasses, any_vars(.==T)) #remove empty groups
  
  up.full = upset(contacts_all, refclasses,
                  #name='GAM',
                  min_size=100,
                  n_intersections=iss,
                  keep_empty_groups=T,
                  base_annotations=list(
                    'Intersection size'=intersection_size(
                      mapping=aes(fill='bars_color')
                    ) + scale_fill_manual(values=c('bars_color'=cl), guide='none') +
                      ylim(c(0,8000))
                  ),
                              set_sizes=(
                                upset_set_size()
                                + geom_text(aes(label=..count..), hjust=1.1, stat='count')
                                + theme(axis.text.x=element_text(angle=90))
                                +expand_limits(y=120000)
                                #ylim(c(0,))
                              ),                  
                  
                  height_ratio=c(1,1),
                  sort_sets=FALSE #use given order
  )
  return(up.full)
}


ym=10000 #max value on y
ts=1.6 #font size
iss=11
u1p = getuplots("./data/upset_data/GAM_05.upset.tsv", "#c83737", iss)
u2p = getuplots("./data/upset_data/HiC_05.upset.tsv", "#4a739f", iss)
u3p = getuplots("./data/upset_data/common_005.upset.tsv", "#ff990c", iss)
u4p = getuplots("./data/upset_data/random_005_1.upset.tsv", "grey30", iss)


# Fig 3e
(u1p | u2p) / (u3p | u4p)



# define colors based on group sizes
getcolors = function(windows){
  wpercentiles = windows %>% tidyr::pivot_longer(cols=2:ncol(windows)) %>% 
    dplyr::select(name, value) %>%
    dplyr::group_by(name) %>%
    dplyr::summarize(cnt = sum(value)) %>%
    dplyr::ungroup() %>%
    bind_rows(list("name"="dummy", "cnt"=0)) %>%  #ntiles should start with zero counts
    dplyr::mutate(ptile = ntile(cnt, 10), cl = colors[ptile]) %>% 
    filter(name!="dummy") %>%
    dplyr::arrange(cnt)
  return(wpercentiles$cl)
}





# get color codes
contacts_all = read_tsv("./data/upset_data/GAM_05.upset.tsv") %>% mutate(ds="GAM") %>%
  bind_rows( read_tsv("./data/upset_data/HiC_05.upset.tsv") %>% mutate(ds="HiC")) %>%
  bind_rows( read_tsv("./data/upset_data/common_005.upset.tsv") %>% mutate(ds="common")) %>%
  bind_rows( read_tsv("./data/upset_data/random_005_1.upset.tsv") %>% mutate(ds="random"))

contacts_all %>% count(ds)

wpercentiles = contacts_all %>% tidyr::pivot_longer(cols=-c("Window pair", "ds")) %>% 
  dplyr::select(name, ds, value) %>%
  dplyr::group_by(name, ds) %>%
  dplyr::summarize(cnt = sum(value)) %>%
  dplyr::ungroup() %>%
  bind_rows(list("name"="dummy", "cnt"=0)) %>%  #ntiles hould start with zero counts
  dplyr::mutate(ptile = ntile(cnt, 10), 
                #fract = ceiling(cnt/ntotal*10), 
                lab=sprintf("%s", cnt), 
                lab2=sprintf("%s (%s)", name, cnt)) %>%
  #dplyr::mutate(ptile = ntile(cnt, 10), cl = colors[ptile]) %>% 
  filter(name!="dummy") %>%
  dplyr::arrange(cnt)

tst = ggplot(wpercentiles, aes(y=factor(name, levels=refclasses),x=ds, fill=factor(ptile, levels=1:10), label=lab )) +
  geom_tile()+
  geom_text(size=3)+
  ylab("")+
  xlab("Set Size")+
  scale_fill_brewer(palette = "RdBu", direction = -1)
  
plot(tst)

  



