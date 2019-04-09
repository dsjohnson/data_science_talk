library(sf)
library(spdep)
library(tidyverse)
library(mgcv)
library(ptolemy)
library(ggspatial)

chl <- readRDS("chl_data.rds")
ak <- ptolemy::alaska()

## list of spatial neighbors
nb <- chl %>% spdep::poly2nb() %>% `names<-`(chl$cell) %>% `class<-`("list")

chl <- chl %>% 
  mutate(
    cell = factor(cell),
    weight = ifelse(is.na(log10chl), 0, 1),
    log10chl = ifelse(is.na(log10chl), mean(log10chl, na.rm=T), log10chl)
  )

plt_chl <- ggplot() + layer_spatial(chl, aes(fill=log10(chl))) + scale_fill_viridis_c()+
  annotation_spatial(ak, fill=1, color=1) + 
  theme(legend.title = element_blank())

### Fit a GMRF model
fit <- mgcv::gam(log10chl ~ s(cell, bs="mrf",xt=list(nb=nb)),
                 data=chl, weights=weight)
### Predict missing values
chl <- bind_cols(chl, predict(fit, se.fit=T) %>% as.data.frame())

plt_chl_pred <- ggplot() + layer_spatial(chl, aes(fill=fit)) + scale_fill_viridis_c()+
annotation_spatial(ak, fill=1, color=1) + 
  theme(legend.title = element_blank())
