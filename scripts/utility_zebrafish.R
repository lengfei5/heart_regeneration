# global functions ####
combine_plots <- function(...,
                          title.text = NULL,
                          title.color = "black",
                          title.size = 16,
                          title.vjust = 0.5,
                          title.hjust = 0.5,
                          title.fontface = "bold",
                          caption.text = NULL,
                          caption.color = "black",
                          caption.size = 10,
                          caption.vjust = 0.5,
                          caption.hjust = 0.5,
                          caption.fontface = "plain",
                          sub.text = NULL,
                          sub.color = "black",
                          sub.size = 12,
                          sub.vjust = 0.5,
                          sub.hjust = 0.5,
                          sub.fontface = "plain",
                          sub.x = 0.5,
                          sub.y = 0.5,
                          sub.vpadding = ggplot2::unit(1, "lines"),
                          sub.angle = 0,
                          sub.lineheight = 0.9,
                          title.rel.heights = c(0.1, 1.2),
                          caption.rel.heights = c(1.2, 0.1),
                          title.caption.rel.heights = c(0.1, 1.2, 0.1)) {
  # preparing the basic plot
  plot <- cowplot::plot_grid(...)
  
  # preparing the title
  if (!is.null(title.text)) {
    title <-
      cowplot::ggdraw() +
      cowplot::draw_label(
        label = title.text,
        colour = title.color,
        size = title.size,
        vjust = title.vjust,
        hjust = title.hjust,
        fontface = title.fontface
      )
  }
  # preparing the caption
  if (!is.null(caption.text)) {
    caption <-
      cowplot::ggdraw() +
      cowplot::draw_label(
        label = caption.text,
        colour = caption.color,
        size = caption.size,
        vjust = caption.vjust,
        hjust = caption.hjust,
        fontface = caption.fontface
      )
  }
  
  # combining the basic plot with the either title or caption or title and
  # caption
  if (!is.null(title.text)) {
    if (!is.null(caption.text)) {
      # if both title and caption are needed
      plot <-
        cowplot::plot_grid(title,
                           plot,
                           caption,
                           ncol = 1,
                           rel_heights = title.caption.rel.heights
        )
    } else {
      # if only title is needed
      plot <-
        cowplot::plot_grid(title,
                           plot,
                           ncol = 1,
                           rel_heights = title.rel.heights
        )
    }
  } else if (!is.null(caption.text)) {
    # if only caption is needed
    plot <-
      cowplot::plot_grid(plot,
                         caption,
                         ncol = 1,
                         rel_heights = caption.rel.heights
      )
  }
  
  # finally adding sub if it's needed
  if (!is.null(sub.text)) {
    plot <-
      cowplot::ggdraw(
        cowplot::add_sub(
          plot = plot,
          label = sub.text,
          x = sub.x,
          y = sub.y,
          vpadding = sub.vpadding,
          colour = sub.color,
          size = sub.size,
          vjust = sub.vjust,
          hjust = sub.hjust,
          fontface = sub.fontface,
          angle = sub.angle,
          lineheight = sub.lineheight
        )
      )
  }
  
  # return the final, combined plot
  return(plot)
}

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
  
}
