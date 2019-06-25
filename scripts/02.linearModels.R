require(vegan)
require(ape)

## changed 2018-04-27 to use PCOA for the phylogenetic eigenvectors

dat.2014.lm <- list(
  allTaxa = list(
    bySp = 0,
    bySpSite = 0
  ),
  angiosperms = list(
    bySp = 0,
    bySpSite = 0
  )
)

phylo.pcoa <- pcoa(cophenetic.phylo(tr))
eco.mds <- metaMDS(scale(dat.2014.bySp[, grep('bio', names(dat.2014.bySp), value = T)]), distance = 'euclidean', k = 4, verbose = F)
eco.mds.stress <- lapply(1:10, function(x) {
  metaMDS(scale(dat.2014.bySp[, grep('bio', names(dat.2014.bySp), value = T)]), distance = 'euclidean', k = x, verbose = F)
  })
plot(1:10, sapply(eco.mds.stress, "[[", 'stress'))

dimnames(phylo.pcoa$vectors)[[2]] <-  paste('phylo', dimnames(phylo.pcoa$vectors)[[2]], sep = '.')
dimnames(eco.mds$points)[[2]] <- paste('eco', dimnames(eco.mds$points)[[2]], sep = '.')
dat.2014.bySp <-
  cbind(dat.2014.bySp,
    eco.mds$points[row.names(dat.2014.bySp), ],
    phylo.pcoa$vectors[row.names(dat.2014.bySp), ]) %>%
  as.data.frame

#added 2018-05-08 to rescale
dat.2014.bySp.scaled <- scale(dat.2014.bySp[, grep('bio|kd|phylo|MDS', colnames(dat.2014.bySp))]) %>%
  as.data.frame

dat.2014.lm$allTaxa$bySp <- list(
      bioclim.phyloPCOA = lm(ln.kd ~ bio1 + bio12 + bio4 + bio14 + phylo.Axis.1 + phylo.Axis.2 + phylo.Axis.3 + phylo.Axis.4, dat.2014.bySp.scaled),
      bioclimMDS.phyloPCOA = lm(ln.kd ~ eco.MDS1 + eco.MDS2 + eco.MDS3 + eco.MDS4 + phylo.Axis.1 + phylo.Axis.2 + phylo.Axis.3 + phylo.Axis.4, dat.2014.bySp.scaled),
      phyloPCOA = lm(ln.kd ~ phylo.Axis.1 + phylo.Axis.2 + phylo.Axis.3 + phylo.Axis.4, dat.2014.bySp.scaled),
      bioclim = lm(ln.kd ~ bio1 + bio12 + bio4 + bio14, dat.2014.bySp.scaled),
      bioclimMDS = lm(ln.kd ~ eco.MDS1 + eco.MDS2 + eco.MDS3 + eco.MDS4, dat.2014.bySp.scaled)
    )

### by species, with phylo predictors -- angiosperms only

phylo.pcoa.angio <- pcoa(cophenetic.phylo(tr.angio))
dimnames(phylo.pcoa.angio$vectors)[[2]] <- paste('phylo', dimnames(phylo.pcoa.angio$vectors)[[2]], sep = '.')

dat.2014.bySp.angio <- dat.2014.bySp.scaled[taxa.angio,
                                     -grep('phylo.Axis', names(dat.2014.bySp))]
dat.2014.bySp.angio <- cbind(dat.2014.bySp.angio, phylo.pcoa.angio$vectors)
dat.2014.bySp.angio.scaled <- scale(dat.2014.bySp.angio) %>%
  as.data.frame

dat.2014.lm$angiosperms$bySp <- list(
  bioclim.phyloPCOA = lm(ln.kd ~ bio1 + bio12 + bio4 + bio14 + phylo.Axis.1 + phylo.Axis.2 + phylo.Axis.4 + phylo.Axis.5, dat.2014.bySp.angio.scaled),
  bioclimMDS.phyloPCOA = lm(ln.kd ~ eco.MDS1 + eco.MDS2 + eco.MDS3 + eco.MDS4 + phylo.Axis.1 + phylo.Axis.2 + phylo.Axis.4 + phylo.Axis.5, dat.2014.bySp.angio.scaled),
  phyloPCOA = lm(ln.kd ~ phylo.Axis.1 + phylo.Axis.2 + phylo.Axis.4 + phylo.Axis.5, dat.2014.bySp.angio.scaled),
  bioclim = lm(ln.kd ~ bio1 + bio12 + bio4 + bio14, dat.2014.bySp.angio.scaled),
  bioclimMDS = lm(ln.kd ~ eco.MDS1 + eco.MDS2 + eco.MDS3 + eco.MDS4, dat.2014.bySp.angio.scaled)
  )

dat.2014.bySpSite <-
  cbind(dat.2014.bySpSite,
    phylo.pcoa$vectors[dat.2014.bySpSite$tidiedName, ],
    eco.mds$points[dat.2014.bySpSite$tidiedName, ]) %>%
  as.data.frame
temp.spp <- dat.2014.bySpSite$tidiedName

## not changed 2018-05-08: already scaled here.
dat.2014.bySpSite <-
  dat.2014.bySpSite[, grep('bio|eco|phylo|kd', names(dat.2014.bySpSite))] %>%
  scale %>%
  as.data.frame
attr(dat.2014.bySpSite, 'spp')  <- temp.spp
dat.2014.bySpSite$site <- as.factor(gsub(paste(c('[',letters,']'), collapse = ''), '',row.names(dat.2014.bySpSite)))
dat.2014.bySpSite$tipLabel <- gsub('[-.0123456789]', '', row.names(dat.2014.bySpSite))
dat.2014.bySpSite$tidiedName <- temp.spp

dat.2014.lm$allTaxa$bySpSite <- list(
  bioclim.phyloPCOA = lm(ln.kd ~ bio1 + bio12 + bio4 + bio14 + phylo.Axis.1 + phylo.Axis.2 + phylo.Axis.3 + phylo.Axis.4, dat.2014.bySpSite),
  bioclimMDS.phyloPCOA = lm(ln.kd ~ eco.MDS1 + eco.MDS2 + eco.MDS3 + eco.MDS4 + phylo.Axis.1 + phylo.Axis.2 + phylo.Axis.3 + phylo.Axis.4, dat.2014.bySpSite),
  phyloPCOA = lm(ln.kd ~ phylo.Axis.1 + phylo.Axis.2 + phylo.Axis.3 + phylo.Axis.4, dat.2014.bySpSite),
  bioclim = lm(ln.kd ~ bio1 + bio12 + bio4 + bio14, dat.2014.bySpSite),
  bioclimMDS = lm(ln.kd ~ eco.MDS1 + eco.MDS2 + eco.MDS3 + eco.MDS4, dat.2014.bySpSite)
  )

### by species, with phylo predictors -- angiosperms only
cols.angio <- which(attr(dat.2014.bySpSite, 'spp') %in% taxa.angio)
dat.2014.bySpSite.angio <-
  dat.2014.bySpSite[cols.angio, -grep('phylo.Axis', names(dat.2014.bySpSite))]

dat.2014.bySpSite.angio.text <- dat.2014.bySpSite.angio[c("site","tipLabel","tidiedName")]
dat.2014.bySpSite.angio[c("site","tipLabel","tidiedName")] <- NULL

dat.2014.bySpSite.angio <-
  cbind(dat.2014.bySpSite.angio,
        phylo.pcoa.angio$vectors[attr(dat.2014.bySpSite, 'spp')[cols.angio], ]) %>%
  scale %>%
  as.data.frame
attr(dat.2014.bySpSite.angio, 'spp') <-
  attr(dat.2014.bySpSite, 'spp')[which(attr(dat.2014.bySpSite, 'spp') %in% taxa.angio)]
dat.2014.bySpSite.angio <- cbind(dat.2014.bySpSite.angio, dat.2014.bySpSite.angio.text)

dat.2014.lm$angiosperms$bySpSite <- list(
  bioclim.phyloPCOA = lm(ln.kd ~ bio1 + bio12 + bio4 + bio14 + phylo.Axis.1 + phylo.Axis.2 + phylo.Axis.3 + phylo.Axis.4, dat.2014.bySpSite.angio),
  bioclimMDS.phyloPCOA = lm(ln.kd ~ eco.MDS1 + eco.MDS2 + eco.MDS3 + eco.MDS4 + phylo.Axis.1 + phylo.Axis.2 + phylo.Axis.3 + phylo.Axis.4, dat.2014.bySpSite.angio),
  phyloPCOA = lm(ln.kd ~ phylo.Axis.1 + phylo.Axis.2 + phylo.Axis.3 + phylo.Axis.4, dat.2014.bySpSite.angio),
  bioclim = lm(ln.kd ~ bio1 + bio12 + bio4 + bio14, dat.2014.bySpSite.angio),
  bioclimMDS = lm(ln.kd ~ eco.MDS1 + eco.MDS2 + eco.MDS3 + eco.MDS4, dat.2014.bySpSite.angio)
  )
