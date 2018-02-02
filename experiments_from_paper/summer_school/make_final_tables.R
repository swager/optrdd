rm(list = ls())
library(optrdd)
library(RColorBrewer)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

data = read.csv("output/result_summaries.csv")
data = data[data$cate.at.pt == data$center,]

idx = c(1, 3)

pretty = data.frame(stringsAsFactors = FALSE,
	CI=paste0("$", round(data$tau.hat, 3), ' \\pm ', round(data$tau.plusminus, 3), "$"),
	max.bias=as.character(round(data$max.bias, 3)),
	samp.se=as.character(round(data$sampling.se, 3)))
	
for(offs in c(0, 2, 4, 6)) {
	cat(Reduce(function(a, b) paste(a, b, sep = " & "), unlist(c(pretty[offs+1,], pretty[offs+2,]))))
	print("")
}
