# plot fLME models

## Figure 4
plot.FUI <- function(r, 
                     fig = NULL, 
                     align, 
                     Hz, 
                     var_name = NULL, 
                     title = NULL,
                     y_val_lim = 1.1,
                     ylim = NULL,
                     xlim = NULL,
                     y_scal_orig = 0.05){
  library(gridExtra)
  name = NULL
  if(!is.null(fig))     fit_dat <- fig
  
  # var_name <- c("Intercept", "Case", "Scan Date (yr)", "Sex", "Age at Baseline (yr)")
  if(is.null(var_name))    var_name <- rownames(fit_dat$betaHat)
  if(nrow(fit_dat$betaHat) != length(var_name) )  var_name <- rownames(fit_dat$betaHat)
  
  decimal <- c(2,2,2,2,3)
  beta.hat.plt <- data.frame(s = seq(1, ncol(fit_dat$betaHat), length.out = ncol(fit_dat$betaHat)), 
                             beta = fit_dat$betaHat[r,],
                             lower = fit_dat$betaHat[r,] - 2*sqrt(diag(fit_dat$betaHat.var[,,r])),
                             upper = fit_dat$betaHat[r,] + 2*sqrt(diag(fit_dat$betaHat.var[,,r])),
                             lower.joint = fit_dat$betaHat[r,] - fit_dat$qn[r]*sqrt(diag(fit_dat$betaHat.var[,,r])),
                             upper.joint = fit_dat$betaHat[r,] + fit_dat$qn[r]*sqrt(diag(fit_dat$betaHat.var[,,r])))
  p.CI <- ggplot() +
    #theme_bw() +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    geom_ribbon(aes(x = s / Hz - align/Hz - 1/Hz, ymax = upper.joint, ymin = lower.joint), data = beta.hat.plt, fill = "gray20", alpha = 0.2) +
    geom_ribbon(aes(x = s / Hz - align/Hz - 1/Hz, ymax = upper, ymin = lower), data = beta.hat.plt, fill = "gray10", alpha = 0.4) +
    geom_line(aes(x = s / Hz - align/Hz - 1/Hz, y = beta, color = "Estimate"), # , color = "Estimate"
              data = beta.hat.plt, alpha = 1, size = 1) + # , lty = 5
    scale_colour_manual(name="", values=c("Estimate"="black")) + # "blue3"
    #scale_colour_manual(name="", values=c("black")) + # "blue3"
    scale_y_continuous(labels=function(x) sprintf(paste0("%.", decimal[r], "f"), x)) #+

  p.CI <- p.CI + 
    labs(x = "Time (s)", y = bquote(paste(beta[.(r-1)], "(s)")), 
         title = var_name[r]) +
    theme(legend.position = "none") #+
  
  # if(r == 1){
  #   p.CI <- p.CI + 
  #     labs(x = "Time (s)", y = expression(beta[0](s)), 
  #          title = var_name[r]) +
  #     theme(legend.position = "none")
  # }else{
  #   p.CI <- p.CI + 
  #     labs(x = "Time (s)", y = bquote(paste(beta[.(r-1)], "(s)")), 
  #                       title = var_name[r]) +
  #     theme(legend.position = "none") #+
  # }
  
  # make x and y intercepts
  if(!is.null(ylim)){
    p.CI <- p.CI + coord_cartesian(ylim = ylim)
  }else{  
    ylim <- c(min(beta.hat.plt$lower.joint), max(beta.hat.plt$upper.joint)) #layer_scales(p.CI)$y$range$range
    y_adjust <- y_scal_orig * (max(beta.hat.plt$upper.joint) - min(beta.hat.plt$lower.joint)) #layer_scales(p.CI)$y$range$range
    ylim[1] <- ylim[1] - y_adjust # just scale bottom because top is scaled below
    }     

  if(is.null(xlim)){
    xlim = layer_scales(p.CI)$x$range$range
  }
  
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
  
  #ggplot2::layer_scales(p.CI)$y$range$range <- ylim*y_val_lim
  
  return(p.CI)
}
