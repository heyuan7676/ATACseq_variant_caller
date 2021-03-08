peak = read.table('data/Test/Peak22759.txt', header = T, stringsAsFactors = F, sep=' ')
real_gt = read.table('data/Test/real_GT.txt', header = T, stringsAsFactors = F, sep=' ')
called_gt = read.table('data/Test/called_GT.txt', header = T, stringsAsFactors = F, sep=' ')
imputed_gt = read.table('data/Test/imputed_GT.txt', header = T, stringsAsFactors = F, sep=' ')
integrated_gt = read.table('data/Test/integrated_GT.txt', header = T, stringsAsFactors = F, sep=' ')

samples = colnames(integrated_gt)
samples = samples[grepl('HG', samples)]
samples = intersect(samples, colnames(peak))
samples = intersect(colnames(real_gt), samples)

real_gt[, samples] = sapply(real_gt[,samples], function(x) sum(c(as.numeric(strsplit(x, '|')[[1]][1]), as.numeric(strsplit(x, '|')[[1]][3]))))
dat = rbind(peak[,samples], real_gt[, samples], called_gt[,samples], imputed_gt[,samples], integrated_gt[,samples])
dat = as.data.frame(dat)
rownames(dat) = c("Peak", "Real","Called", "Imputed", "Integrated")

dat = as.data.frame(t(dat))
