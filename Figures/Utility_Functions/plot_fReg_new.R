# plot freg models

plot.freg <- function(fit_dat,
                      r, 
                      fig = NULL, 
                      align, 
                      Hz, 
                      var_name = NULL, 
                      title = NULL,
                      y_val_lim = 1.1,
                      ylim = NULL,
                      y_scal_orig = 0.05){

  library(gridExtra)
  name = NULL
  if(is.null(var_name))    var_name <- paste0("Variable", r)

  decimal <- c(2,2,2,2,3)
  if(class(fit_dat) == "fosr"){
    beta.hat.plt <- data.frame(s = seq(1, nrow(fit_dat$est.func), length.out = nrow(fit_dat$est.func)), 
                               beta = fit_dat$est.func[,r],
                               lower.joint = fit_dat$est.func[,r] - 1.96 * fit_dat$se.func[,r],
                               upper.joint = fit_dat$est.func[,r] + 1.96 * fit_dat$se.func[,r])
  }else{
    # bayes fosr
    beta.hat.plt <- data.frame(s = seq(1, ncol(fit_dat$beta.hat), length.out = ncol(fit_dat$beta.hat)), 
                               beta = fit_dat$beta.hat[r,],
                               lower.joint = fit_dat$beta.LB[r,],
                               upper.joint = fit_dat$beta.UB[r,] )
  }
  
  p.CI <- ggplot() +
    #theme_bw() +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    geom_ribbon(aes(x = s / Hz - align/Hz - 1/Hz, ymax = upper.joint, ymin = lower.joint), 
                data = beta.hat.plt, fill = "gray10", alpha = 0.2) +
    geom_line(aes(x = s / Hz - align/Hz - 1/Hz, y = beta, color = "Estimate"), # , color = "Estimate"
              data = beta.hat.plt, alpha = 1, size = 1) + # , lty = 5
    scale_colour_manual(name="", values=c("Estimate"="black")) + # "blue3"
    #scale_colour_manual(name="", values=c("black")) + # "blue3"
    scale_y_continuous(labels=function(x) sprintf(paste0("%.", decimal[r], "f"), x)) #+
  
  p.CI <- p.CI + 
    labs(x = "Time (s)", y = bquote(paste(beta[.(r-1)], "(s)")), 
         title = var_name[r]) +
    theme(legend.position = "none") #+
  
  # make x and y intercepts
  if(!is.null(ylim)){
    p.CI <- p.CI + coord_cartesian(ylim = ylim)
  }else{  
    ylim <- c(min(beta.hat.plt$lower.joint), max(beta.hat.plt$upper.joint)) #layer_scales(p.CI)$y$range$range
    y_adjust <- y_scal_orig * (max(beta.hat.plt$upper.joint) - min(beta.hat.plt$lower.joint)) #layer_scales(p.CI)$y$range$range
    ylim[1] <- ylim[1] - y_adjust # just scale bottom because top is scaled below
  }     
  
  xlim = layer_scales(p.CI)$x$range$range
  
  x_range <- diff(xlim) * 0.1
  y_range <- diff(ylim) * 0.1
  y_range_up <- diff(ylim) * 0.02
  
  # extend upper limit
  y_val_lim <- c(1, y_val_lim)
  y_top <- (0.975) * diff(ylim*y_val_lim) + ylim[1]*y_val_lim[1]
  
  p.CI <- p.CI + 
    coord_cartesian(ylim = ylim*y_val_lim,
                    xlim = xlim) + 
    geom_segment(aes(x=xlim[1] - x_range, xend=xlim[2] + x_range,
                     y=0,yend=0),
                 color = "black", lwd=0.5, alpha = 0.75, linetype = "dashed") + # x - intercept
    geom_segment(aes(y=ylim[1] - y_range, yend=y_top,#ylim[2]+y_range_up,
                     x=0,xend=0), # don't extend up
                 color = "black", lwd=0.5, alpha = 0.75, linetype = "dashed") 
  
  return(p.CI)
}
