library("ggplot2")

cleanPlot2 <- function() {theme_bw()+
    theme(
      panel.spacing = unit(0.2, "lines"),
      strip.background = element_blank(),
      strip.text.y = element_text(angle = 0, size = 16),
      strip.text.x = element_text(size = 16),
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 20),
      legend.text = element_text(size = 16),
      legend.title = element_text(size = 16),
      #legend.position = 'top',
      axis.line = element_blank(),
      #axis.ticks = element_blank(),
      panel.border = element_rect(fill = NA,colour = "grey20"),
      panel.grid = element_blank(),
      panel.grid.major = element_line(color = 'grey85', linetype = 'dashed'),
      #panel.border = element_blank(),
      #panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    )
}

guide_textcolourguide <- function(...) {
  guide <- guide_legend(...)
  class(guide) <- c("guide", "textcolourguide", "legend")
  guide
}

guide_gengrob.textcolourguide <- function(guide, theme) {
  legend <- NextMethod()
  
  # Figure out what are keys and labels
  keys <- grep("^key(?!.*bg)", legend$layout$name, perl = TRUE)
  labels <- grep("^label", legend$layout$name)
  
  # Recolour the labels based on keys, assumes parallel ordering
  newlabels <- mapply(function(key, lab) {
    colour <- legend$grobs[[key]]$gp$col
    lab <- legend$grobs[[lab]]
    lab$children[[1]]$children[[1]]$gp$col <- colour
    return(lab)
  }, key = keys, lab = labels, SIMPLIFY = FALSE)
  
  # Replace labels
  legend$grobs[labels] <- newlabels
  
  # Purge keys
  gtable::gtable_filter(legend, "key", invert = TRUE)
}