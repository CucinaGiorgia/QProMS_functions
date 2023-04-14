theme_cuc <- function(){
  font <- "NimbusRom"
  theme_light() %+replace%
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, vjust = 2), 
      legend.position = 'right',
      legend.title = element_blank(),
      axis.title.x = element_text(size=15),
      axis.title.y = element_text(size = 15, angle=90, vjust = 2),
      axis.text.x = element_text(size=10),
      axis.text.y = element_text(size = 10)
  
  )
    
}
