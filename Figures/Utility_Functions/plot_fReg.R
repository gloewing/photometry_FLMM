# plot fLME models

## Figure 4
plot.freg <- function(fit_dat, r, align, Hz, var_name = NULL, title = NULL){
  library(gridExtra)
  name = NULL
  # var_name <- c("Intercept", "Case", "Scan Date (yr)", "Sex", "Age at Baseline (yr)")
  if(is.null(var_name))    var_name <- paste0("Variable", r)
  # if(nrow(fit_dat$beta.hat) != length(var_name) )  var_name <- rownames(fit_dat$beta.hat)
  
  decimal <- c(2,2,2,2,3)
  if(class(fit_dat) == "fosr"){
    beta.hat.plt <- data.frame(s = seq(1, nrow(fit_dat$est.func), length.out = nrow(fit_dat$est.func)), 
                               beta = fit_dat$est.func[,r],
                               lower.joint = fit_dat$est.func[,r] + 1.96 * fit_dat$se.func[,r],
                               upper.joint = fit_dat$est.func[,r] - 1.96 * fit_dat$se.func[,r])
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
    geom_ribbon(aes(x = s / Hz - align/Hz - 1/Hz, ymax = upper.joint, ymin = lower.joint), data = beta.hat.plt, fill = "gray30", alpha = 0.2) +
    geom_line(aes(x = s / Hz - align/Hz - 1/Hz, y = beta, color = "Estimate"), data = beta.hat.plt, alpha = 1, lty = 5) +
    scale_colour_manual(name="", values=c("Estimate"="blue3")) +
    scale_y_continuous(labels=function(x) sprintf(paste0("%.", decimal[r], "f"), x)) #+
  #scale_x_continuous(breaks = c(1, 24, 47, 70, 93))
  
  if(r == 1){
    p.CI <- p.CI + labs(x = "Time (s)", y = expression(beta[0](s)), title = var_name[r]) +
      geom_hline(yintercept=0, color = "black", lwd=0.5, linetype = "dotted") +
      theme(legend.title=element_blank(),
            legend.position = c(0.15, 0.99),
            legend.justification = c("left", "top"),
            legend.box.just = "right",
            legend.margin = margin(0, 0, 0, 0),
            legend.background = element_rect(fill=alpha('white', 0)))
  }else{
    p.CI <- p.CI + labs(x = "Time (s)", y = bquote(paste(beta[.(r-1)], "(s)")), 
                        title = var_name[r]) +
      theme(legend.position = "none") +
      geom_hline(yintercept=0, color = "black", lwd=0.5, linetype = "dotted") +
      geom_vline(xintercept=0, color = "black", lwd=0.5, linetype = "dashed", alpha = 0.5)
    
  }
  
  # if(!is.null(title)){
  #   p.CI <- p.CI + ggtitle(title)
  # }
  return(p.CI)
}
