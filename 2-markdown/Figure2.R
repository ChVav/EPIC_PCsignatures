
# Copyright (c) 2023 Charlotte D. Vavourakis. Licensed under the MIT license (see LICENSE.md).

library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)
library(rcartocolor)

palette = carto_pal(8,"Safe")

## note that also Lasso penalized regression was performed
## note that all models were also trained without removing flagged probes (low normalized intensit, not on EPIC v2.0), but will only show results after removal

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
plota <- accuracy %>%
  filter(method %in% c("ElNet","PC")) %>% 
  filter(input.CpGs == "unreliable.removed") %>%
  droplevels() %>%
  ggplot(aes(x=train.sample.size, shape = method, linetype = method)) +
  geom_point(aes(y=val.slope.regress, color = colPoints), size=0.8) +
  geom_line(aes(y=val.slope.regress.fit, color = colLines), size=0.8) +
  scale_color_identity() + #guide="legend" adds a legend that is not useful, modify manually in Inscape
  scale_linetype_manual(values=c("solid", "dashed"))+
  labs(color="", shape = "", linetype = "") +
  ylab("Accuracy 3CExt test set \n [slope]") +
  xlab("n training samples") +
  theme_minimal()

plotb <- precision %>%
  filter(method %in% c("ElNet","PC")) %>% 
  filter(input.CpGs == "unreliable.removed") %>%
  droplevels() %>%
  ggplot(aes(x=train.sample.size, shape = method, linetype = method)) +
  geom_point(aes(y=rep.ICC, color = colPoints), size=0.8) +
  geom_line(aes(y=rep.ICC.fit, color = colLines), size=0.8) +
  scale_color_identity() +
  scale_linetype_manual(values=c("solid", "dashed"))+
  labs(color="", shape = "", linetype = "") +
  ylab("Precision repeatability test set \n [ICC]") +
  xlab("n training samples") +
  theme_minimal()

plot <- plota + plotb + plot_layout(guides="collect") + plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')')

ggsave(filename='Figure2.pdf', plot, width=11, height=4, dpi = 600) #legend needs to be edited manually unfortunately
