# 
plot_adjust <- function(plots){
  library(gridExtra)
  # Your code, but using p1 and p2, not the plots with adjusted margins
  gl <- lapply(plots, ggplotGrob)
  widths <- do.call(grid::unit.pmax, lapply(gl, function (x) x$widths)) #"[[", "widths"))
  heights <- do.call(grid::unit.pmax, lapply(gl, function (x) x$heights))
  lg <- lapply(gl, function(g) {g$widths <- widths; g$heights <- heights; g})
  
  # # New code
  # library(gtable)
  # gt = cbind(lg[[1]], lg[[2]][, -(1:3)], size = "first")
  # 
  # gt$widths[5] = unit(2, "lines")
  # 
  # # Draw the plot
  # grid::grid.newpage()
  # grid::grid.draw(gt)
  fig <- do.call("grid.arrange", c(lg, nrow = 2))
  
  
}
