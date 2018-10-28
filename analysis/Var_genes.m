			means = apply(filtered.normalized,1,function(v) log(mean(exp(v) - 1) + 1))
			frac.cells = rowSums(filtered.normalized>0)
			vars = apply(filtered.normalized,1,function(v) log(var(exp(v) - 1) + 1))
			dispersion = apply(filtered.normalized,1,function(v) log(var(exp(v) - 1) / mean(exp(v) - 1)))
			dispersion[is.na(x = dispersion)] = 0
			means[is.na(x = means)] = 0

			library(ggplot2,quietly=T)
			pdf(file.path(environment$baseline.work.path,paste(dataset,'VariableGenes.pdf',sep='.')))
			plot.data = data.frame(gene = names(means),means,dispersion)#head(plot.data)
			print(ggplot(plot.data, aes(x=means, y=dispersion, label = gene)) + geom_text(check_overlap = TRUE,size=2))
			smoothScatter(means,dispersion)
			dev.off()

			num.bin = 20
			bins = cut(x = means, breaks = num.bin)
			names(x = bins) = names(x = means)
			mean_y = tapply(dispersion,bins,mean)
			sd_y = tapply(dispersion,bins,sd)
			dispersion.scaled = (dispersion - mean_y[as.numeric(x = bins)]) / sd_y[as.numeric(x = bins)]
			dispersion.scaled[is.na(x = dispersion.scaled)] = 0
			names(x = dispersion.scaled) = names(x = means)
			
			criterion = means >= min.mean & frac.cells >= min.frac.cells & dispersion.scaled >= min.dispersion.scaled
			HVG = names(means)[criterion]