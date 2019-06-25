## redoing auteur analyses
library(geiger)
do.rbm = F
do.jump.bm = F
do.kd = T

if(do.rbm){
  dat.2014.rjmcmc.rates.rbm.noSE <- rjmcmc.bm(tr,
    structure(dat.2014.bySp$ln.kd, names = row.names(dat.2014.bySp)),
    ngen = 1E06,
    filebase = 'ln.kd.rbm.noSE')

  dat.2014.rjmcmc.rates.rbm.noSE.results <- load.rjmcmc('jump-relaxedBM.ln.kd.rbm.noSE')

  pdf('SUPPLEMENT.rjmcmc.rbm.rates.PDF', 8.5, 11)
  plot(dat.2014.rjmcmc.rates.rbm.noSE.results, burnin = 0.20, par = 'jumps',
  show.tip = F, edge.width = 2)
  dev.off()
}

if(do.jump.bm) {
  dat.2014.rjmcmc.rates.jump.bm.noSE <- rjmcmc.bm(tr,
    structure(dat.2014.bySp$ln.kd, names = row.names(dat.2014.bySp)),
    ngen = 1E06,
    type = 'jump-bm',
    filebase = 'ln.kd.noSE')

  dat.2014.rjmcmc.rates.jumpbm.noSE.results <- load.rjmcmc('jump-BM.ln.kd.noSE')

  pdf('SUPPLEMENT.rjmcmc.jump-bm.rates.PDF', 8.5, 11)
  plot(dat.2014.rjmcmc.rates.jumpbm.noSE.results, burnin = 0.20, par = 'jumps',
  show.tip = F, edge.width = 2)
  dev.off()
}

if(do.kd) {
  dat.2014.rjmcmc.rates.jump.bm.kd.noSE <- rjmcmc.bm(tr,
    structure(dat.2014.bySp$kd, names = row.names(dat.2014.bySp)),
    ngen = 1E06,
    type = 'jump-bm',
    filebase = 'kd.noSE')

  dat.2014.rjmcmc.rates.jump.bm.kd.noSE.results <- load.rjmcmc('jump-BM.kd.noSE')

  pdf('../out/SUPPLEMENT.rjmcmc.jump-bm.rates.PDF', 8.5, 11)
  plot(dat.2014.rjmcmc.rates.jump.bm.kd.noSE.results, burnin = 0.20, par = 'jumps',
  show.tip = F, edge.width = 2)
  dev.off()

}
