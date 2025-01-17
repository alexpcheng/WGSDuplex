theme_bw2 <- function() {
  theme_bw()+
    theme(axis.text.x = element_text(family = "Helvetica", size=6, color = "black"),
          axis.text.y = element_text(family = "Helvetica", size=6, color = "black"),
          axis.title = element_text(family = "Helvetica", size=8, color = "black"),
          strip.text = element_text(family = "Helvetica", size=8, color = "black"),
          strip.background = element_rect(fill="White", colour = "black", size = 0.25),
          panel.border = element_rect(color = "black", size = 0.25),
          axis.ticks = element_line(color = "black", size = 0.25),
          panel.grid.major = element_line(size=0.25),
          panel.grid.minor = element_blank())
}
