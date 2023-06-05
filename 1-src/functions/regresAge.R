# Copyright (c) 2023 Charlotte D. Vavourakis. Licensed under the MIT license (see LICENSE.md).

## helper function for plotting and calculating fit for regression predicted age ~ chronological age
regresAge <- function(df,ypar,title,mode){
  require(dplyr)
  require(tidyverse)
  require(ggplot2)
  
  fit <- lm(get(ypar)~age, data=df)
  anno <- paste0("R2 = ", signif(summary(fit)$r.squared, 3),
                 "\ny = ", round(fit$coefficients[2],2) , "x + ", round(fit$coefficients[1],2),
                 "\nRMSE = ", round(sqrt(mean(fit$residuals^2)),3))
  plot <- df %>%
    ggplot(aes(x = age, y = get(ypar))) +
    geom_point() +
    theme(panel.background = element_blank(),aspect.ratio = 1) +
    geom_smooth(method = "lm",
                formula = y ~ x,
                size = 0.5,
                se = FALSE,
                alpha = 0.3,
                colour = "gray40") +
    xlim(20,90) +
    ylim(20,90) +
    ylab(ypar) +
    annotate("text",
             x=25,
             y = 85,
             label = anno,
             hjust = 0,
             size = 3) +
    ggtitle(title)
  
  ifelse(mode=="plotting",return(plot),return(fit))
}