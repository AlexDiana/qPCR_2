library(ggplot2)

fun1 <- function(x){
  
  (exp(x) - 1) * ifelse(x > 0, 1, 0)
  
}

ggplot() + 
  xlim(c(- 1, 2)) +
  stat_function(fun = fun1) + 
  xlab("") + ylab("") + 
  theme(
    axis.text = element_text(size = 12),
    panel.background = element_rect(fill = "white", 
                                        colour = NA), panel.border = element_rect(fill = NA, 
                                                                                  colour = "grey20"), panel.grid = element_line(colour = "grey92"), 
        panel.grid.minor = element_line(size = rel(0.5)), 
        strip.background = element_rect(fill = "grey85", 
                                        colour = "grey20"), legend.key = element_rect(fill = "white", 
                                                                                      colour = NA), complete = TRUE)
