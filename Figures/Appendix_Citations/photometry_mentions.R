# https://scholar.google.com/scholar?q=%22fiber+photometry%22&hl=en&as_sdt=0%2C21&as_ylo=2014&as_yhi=2014

years <- 2014:2022
t <- 1:length(years)
refs <- c(6, 27, 71, 137, 240, 372, 532, 614, 879) # 549

add_nls <- nls(refs ~ a*exp(r*t), 
               start = list(a = 20, r = 0.5)) # estimates are very stable for reasonable initializations
b_hat <- coef(add_nls)
y_hat <- b_hat[1] * exp(b_hat[2] * t) # fitted values
#plot(refs, resid(add_nls))

pdf(file="~/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/citations/citations_per_year.pdf", 
    width=6, height=6)
barplot(refs, names = years, 
        ylab = "Google Scholar ''Fiber Photometry'' Mentions",
        xlab = "Year")
lines(y_hat, lwd = 3, col = "blue")
dev.off()
rm(add_nls, b_hat, y_hat)
# total results ~ 3,600


# web of science citation report
# https://www.webofscience.com/wos/woscc/summary/9f8ebd37-0a5a-45a5-88b1-be1267152ab6-933d613d/relevance/1(overlay:export/exc)
refs <- c(1, 3, 11, 19, 33, 45, 75, 70, 123) # 50

add_nls <- nls(refs ~ a*exp(r*t), 
               start = list(a = 20, r = 0.5)) # estimates are very stable for reasonable initializations
b_hat <- coef(add_nls)
y_hat <- b_hat[1] * exp(b_hat[2] * t) # fitted values

pdf(file="~/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/citations/citations_per_year_wos.pdf", 
    width=6, height=6)
barplot(refs, names = years, 
        ylab = "Web of Science ''Fiber Photometry'' Citations",
        xlab = "Year")
lines(y_hat, lwd = 3, col = "blue")
dev.off()
