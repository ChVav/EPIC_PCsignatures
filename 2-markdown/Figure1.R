# Copyright (c) 2023 Charlotte D. Vavourakis. Licensed under the MIT license (see LICENSE.md).

library(patchwork)
library(tidyverse)
library(rcartocolor)

palette = carto_pal(8,"Safe")
cols <- c(palette[1],palette[5])

# source helper functions
source("../1-src/functions/regresAge.R") # helper function for calculating root mean squared error

### a/b) performance Hannum original and Hannum_PC, Hannum PC_proxy in 450K blood data, full population study ##----
result_BloodFull <- readRDS("../1-src/analysis/1-Hannum-clocks-Higgins/3-output/pheno_BloodFull_450K_Hannum_PC.Rds")

# Original Hannum vs chronological age
fit <- lm(round(Hannum, digits=4) ~ round(age, digits=4), data = result_BloodFull)
anno <- paste0("R2 = ", signif(summary(fit)$r.squared, 3),
               "\ny = ", round(fit$coefficients[2],2) , "x + ", round(fit$coefficients[1],2))

#CI for slope
#fit <- lm(Hannum ~ age, data = result_BloodFull)
#confint(fit,'age',level=0.95) #0.8121643 - 0.8431847

plota <- result_BloodFull %>%
  ggplot(aes(x = round(age, digits=4), y = round(Hannum, digits=4))) +
  geom_point(color=cols[1], shape=15, size=0.3) +
  theme(panel.background = element_blank()) +
  geom_smooth(method = "lm",
              formula = y ~ x,
              size = 0.5,
              se = FALSE,
              alpha = 0.3,
              colour = "gray40") +
  xlab("Age") +
  ylab("Hannum predicted age")  +
  xlim(35,85) +
  ylim(35,85) +
  annotate("text",
           x=50,
           y = 40,
           label = anno,
           hjust = 0,
           size = 3) 

# PC Age vs chronological age
fit <- lm(round(Hannum.PC2, digits=4) ~ round(age, digits=4), data = result_BloodFull)
anno <- paste0("R2 = ", signif(summary(fit)$r.squared, 3),
               "\ny = ", round(fit$coefficients[2],2) , "x + ", round(fit$coefficients[1],2))

#CI for slope
#fit <- lm(Hannum.PC2 ~ age, data = result_BloodFull)
#confint(fit,'age',level=0.95) # 0.5725942 - 0.5983141

plotb <- result_BloodFull %>%
  ggplot(aes(x = round(age, digits=4), y = round(Hannum.PC2, digits=4))) +
  geom_point(color=cols[2], shape=16, size=0.3) +
  theme(panel.background = element_blank(),plot.margin = margin(0, 0, 0, 0, "pt")) +
  geom_smooth(method = "lm",
              formula = y ~ x,
              size = 0.5,
              se = FALSE,
              alpha = 0.3,
              colour = "gray40") +
  xlab("Age") +
  ylab("Hannum_PC2 predicted age")  +
  xlim(35,85) +
  ylim(35,85) +
  annotate("text",
           x=50,
           y = 40,
           label = anno,
           hjust = 0,
           size = 3) 

## combined plot ##----

des <- "
AB"

plot <- plota +  plotb + plot_layout(design = des) + plot_annotation(tag_levels = 'a', tag_prefix = '(', tag_suffix = ')')

ggsave(filename='Figure1.pdf', plot, width=6, height=2.5, dpi = 600)
ggsave(filename= 'Figure1.jpeg', plot, width=120, height=60, units= "mm", dpi=600)

