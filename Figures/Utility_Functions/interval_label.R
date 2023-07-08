# add bar labels to indicate time periods
interval_label <- function(fig,
                           x_interval_list = list(c(-2,-0.5), c(-0.5, 1)),
                           x_interval_text = NULL, #list(c(-2), c(-0.5)),
                           text_list = list("Baseline Period", "Reward Period"),
                           scl = 1.01, # percent for lines above original ylim values
                           x_scl = 0.01, # percent for gap in lines
                           txt_adjust = 0.03, # percent text is above line
                           txt_size = 3.75,
                           col_seq = c("#ca0020", "#0868ac", "#E69F00", "#525252"),
                           ylim = NULL,
                           y_val = NULL,
                           alpha = 1
){
  
  if(is.null(ylim))     ylim = layer_scales(fig)$y$range$range
  xlim = layer_scales(fig)$x$range$range
  y_range <- diff(ylim)
  x_range <- diff(xlim)
  x_gap <- c(x_range * x_scl, -x_range * x_scl) # add a small bit to lower end and subtract from upper end
    
  x_interval_list <- lapply(x_interval_list, function(xx) xx + x_gap )
  
  # text_list = list("Baseline Period", "Reward Period")
  # x_interval_list = list(c(-2,-0.5), c(-0.5, 1))
  
  if(is.null(y_val)){
    y_val <- ylim[2] - y_range * (1-scl)
  }
  
  y_val_txt <- ylim[2] - y_range * (1-scl - txt_adjust)
  
  len <- length(x_interval_list)
  
  # use mean if not specified
  if(is.null(x_interval_text))   x_interval_text  <- lapply(x_interval_list, mean)
  
  # col_seq <- seq(30, 300, by = 20)
  
  # fig = plot.f4[[prdtr]]
  for(i in 1:len){
    col <- col_seq[i] #paste("grey", col_seq[i])
    fig <- fig + 
      annotate(geom="text", x=x_interval_text[[i]], size = txt_size,
               y=y_val_txt, label=text_list[[i]], color = col, alpha = alpha) +
      geom_segment(aes_string(x=x_interval_list[[i]][1],xend=x_interval_list[[i]][2],y=y_val,yend=y_val), 
                   inherit.aes = TRUE, color = col, size = 1, alpha = alpha)
  }
  
  return(fig)
  
}
