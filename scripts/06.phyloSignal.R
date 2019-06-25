## phylogenetic signal
library(picante)
library(geiger)
tail = 2
roundTo = 3
perms = 4999

dat.2014.allTaxa.physignal <-
  multiPhylosignal(dat.2014.bySp[c('kd', 'ln.kd', paste('bio', c(1,12,4,14), sep = ''), paste('eco.MDS', 1:4, sep = ''))], tr, reps = perms)
dat.2014.angiosperms.physignal <-
  multiPhylosignal(dat.2014.bySp.angio[c('kd', 'ln.kd', paste('bio', c(1,12,4,14), sep = ''), paste('eco.MDS', 1:4, sep = ''))], tr.angio, reps = perms)
dat.2014.allTaxa.physignal$lambda <-
  fitContinuous(tr, dat.2014.bySp[c('kd', 'ln.kd', paste('bio', c(1,12,4,14), sep = ''), paste('eco.MDS', 1:4, sep = ''))], model = 'lambda') %>%
  sapply(., function(x) x$opt$lambda)
dat.2014.angiosperms.physignal$lambda <-
  fitContinuous(tr.angio, dat.2014.bySp.angio[c('kd', 'ln.kd', paste('bio', c(1,12,4,14), sep = ''), paste('eco.MDS', 1:4, sep = ''))], model = 'lambda') %>%
  sapply(., function(x) x$opt$lambda)


dat.2014.physignal.formatted <- cbind(
  Angiosperms = apply(round(dat.2014.angiosperms.physignal, roundTo), 1, function(x) paste('K = ', x[1], ', p = ', tail*x[4], ', lambda = ', x[6], sep = '')),
  allTaxa = apply(round(dat.2014.allTaxa.physignal, roundTo), 1, function(x) paste('K = ', x[1], ', p = ', tail*x[4], ', lambda = ', x[6], sep = ''))
)

dimnames(dat.2014.physignal.formatted)[[2]] <- c('Angiosperms tree', 'All taxa tree')
write.csv(dat.2014.physignal.formatted, '../out/TABLE.01.phylogenetic.signal.csv')
