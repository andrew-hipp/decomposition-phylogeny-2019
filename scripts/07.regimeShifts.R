## shifts if lnKd using OU models

library(l1ou)
library(PhylogeneticEM)
library(magrittr)

do.analyses <- TRUE

if(do.analyses) {
  dat.2014.phyEm <- PhyloEM(tr.edit,
                            dat.2014.bySp[tr$tip.label, 'ln.kd', drop = F] %>%
                              t %>%
                              as.data.frame,
                            process = 'OU', parallel_alpha = TRUE, Ncores = 14)

  dat.2014.l1ou <- adjust_data(tr.edit, dat.2014.bySp[, 'ln.kd', drop = F])
  dat.2014.l1ou.est <- estimate_shift_configuration(dat.2014.l1ou$tree, dat.2014.l1ou$Y)
  save(dat.2014.phyEm, file = ('PhyloEM.tree.Rdata'))
  save(dat.2014.l1ou.est, file = ('l1ou.tree.Rdata'))
  }

pdf('../out/FigS2a.PhyloEM.tre.pdf', 8.5, 11)
plot(dat.2014.phyEm, cex = 0.3)
dev.off()

pdf('../out/FigS2b.l1ou.tre.pdf', 8.5, 11)
plot(dat.2014.l1ou.est, cex = 0.3)
dev.off()
