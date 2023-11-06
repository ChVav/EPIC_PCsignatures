
# Copyright (c) 2023 Charlotte D. Vavourakis. Licensed under the MIT license (see LICENSE.md).

library(tidyverse)
library(patchwork)
library(rcartocolor)

palette = carto_pal(8,"Safe")

## note that also Lasso penalized regression was performed, but not shown here
## note that all models were on all CpGs, we tested also removing flagged unreliable probes (low normalized intensit, not on EPIC v2.0), not shown here

## Experiment 450K blood ##----

accuracy <- readRDS("../1-src/analysis/3-training-size-experiment/8-output/accuracy_slope.Rds")
precision <- readRDS("../1-src/analysis/3-training-size-experiment/8-output/precision_ICC.Rds")

#add variable so that can have different shades of colors for lines/points
accuracy <- accuracy %>%
  mutate(colPoints = case_when(method == "ElNet" ~ palette[7], TRUE ~ palette[2]),
         colLines = case_when(method == "ElNet" ~ palette[4], TRUE ~ palette[6]))

precision <- precision %>%
  mutate(colPoints = case_when(method == "ElNet" ~ palette[7], TRUE ~ palette[2]),
         colLines = case_when(method == "ElNet" ~ palette[4], TRUE ~ palette[6]))

#### plot all results and modeled curves
plota <- accuracy %>%
  filter(method %in% c("ElNet","PC")) %>% 
  droplevels() %>%
  ggplot(aes(x=train.sample.size, shape = method, linetype = method)) +
  geom_point(aes(y=val.slope.regress, color = colPoints), size=0.8) +
  geom_line(aes(y=val.slope.regress.fit, color = colLines), size=0.8) +
  scale_color_identity() + #guide="legend" adds a legend that is not useful, modify manually in Inscape
  scale_linetype_manual(values=c("solid", "dashed"))+
  labs(color="", shape = "", linetype = "") +
  ylab("Slope regression age vs predicted age") +
  xlab('n training samples \n from "BloodFull_450K"') +
  ylim(0, 1) +
  theme_minimal() +
  ggtitle('Accuracy in "Hannum" test set')

plotb <- precision %>%
  filter(method %in% c("ElNet","PC")) %>% 
  droplevels() %>%
  ggplot(aes(x=train.sample.size, shape = method, linetype = method)) +
  geom_point(aes(y=rep.ICC, color = colPoints), size=0.8) +
  geom_line(aes(y=rep.ICC.fit, color = colLines), size=0.8) +
  scale_color_identity() +
  scale_linetype_manual(values=c("solid", "dashed"))+
  labs(color="", shape = "", linetype = "") +
  ylab("ICC") +
  xlab('n training samples \n from "BloodFull_450K"') +
  ylim(0, 1) +
  theme_minimal() +
  ggtitle('Reliability in "BloodRep_450K" test set')

## Experiment EPIC cervical smear samples  ##----
accuracy <- readRDS("../1-src/analysis/3-training-size-experiment/4-output/accuracy_slope.Rds")
precision <- readRDS("../1-src/analysis/3-training-size-experiment/4-output/precision_ICC.Rds")

#add variable so that can have different shades of colors for lines/points
accuracy <- accuracy %>%
  mutate(colPoints = case_when(method == "ElNet" ~ palette[7], TRUE ~ palette[2]),
         colLines = case_when(method == "ElNet" ~ palette[4], TRUE ~ palette[6]))
 
precision <- precision %>%
  mutate(colPoints = case_when(method == "ElNet" ~ palette[7], TRUE ~ palette[2]),
         colLines = case_when(method == "ElNet" ~ palette[4], TRUE ~ palette[6]))

#### plot all results and modeled curves
plotc <- accuracy %>%
  filter(method %in% c("ElNet","PC")) %>% 
  droplevels() %>%
  ggplot(aes(x=train.sample.size, shape = method, linetype = method)) +
  geom_point(aes(y=val.slope.regress, color = colPoints), size=0.8) +
  geom_line(aes(y=val.slope.regress.fit, color = colLines), size=0.8) +
  scale_color_identity() + #guide="legend" adds a legend that is not useful, modify manually in Inscape
  scale_linetype_manual(values=c("solid", "dashed"))+
  labs(color="", shape = "", linetype = "") +
  ylab("Slope regression age vs predicted age") +
  xlab('n training samples \n from "3CDisc"') +
  ylim(0, 1) +
  theme_minimal() +
  ggtitle('Accuracy in "3CExtVal" test set')

plotd <- precision %>%
  filter(method %in% c("ElNet","PC")) %>% 
  droplevels() %>%
  ggplot(aes(x=train.sample.size, shape = method, linetype = method)) +
  geom_point(aes(y=rep.ICC, color = colPoints), size=0.8) +
  geom_line(aes(y=rep.ICC.fit, color = colLines), size=0.8) +
  scale_color_identity() +
  scale_linetype_manual(values=c("solid", "dashed"))+
  labs(color="", shape = "", linetype = "") +
  ylab("ICC") +
  xlab('n training samples \n from "3CDisc"') +
  ylim(0, 1) +
  theme_minimal() +
  ggtitle('Reliability in "repeatability" test set')


## combined plot  ##----

des <- c("
         AB
         CD")

plot <- plota + plotb + plotc + plotd + plot_layout(guides="collect", design=des) + plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')')

ggsave(filename='Figure2.pdf', plot, width=11, height=8, dpi = 600) #legend needs to be edited manually unfortunately
