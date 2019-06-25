# added 2019-06-25
require(vegan)
require(ape)

d2pcComp <- list(
  allTaxa = list(
    bySp = 0,
    bySpSite = 0
  ),
  angiosperms = list(
    bySp = 0,
    bySpSite = 0
  )
)

eco.pca <- prcomp(scale(dat.2014.bySp[, grep('bio', names(dat.2014.bySp), value = T)]))
dimnames(eco.pca$x)[[2]] <- paste('eco', dimnames(eco.pca$x)[[2]], sep = '.')

d2p2 <-
  cbind(dat.2014.bySp,
    eco.pca$x[row.names(dat.2014.bySp), ]) %>%
  as.data.frame

#added 2018-05-08 to rescale
d2p2.scaled <- scale(d2p2[, grep('bio|kd|phylo|MDS|PC', colnames(d2p2))]) %>%
  as.data.frame

d2pcComp$allTaxa$bySp <- list(
      bioclim.phyloPCOA = lm(ln.kd ~ bio1 + bio12 + bio4 + bio14 + phylo.Axis.1 + phylo.Axis.2 + phylo.Axis.3 + phylo.Axis.4, d2p2.scaled),
      bioclimMDS.phyloPCOA = lm(ln.kd ~ eco.MDS1 + eco.MDS2 + eco.MDS3 + eco.MDS4 + phylo.Axis.1 + phylo.Axis.2 + phylo.Axis.3 + phylo.Axis.4, d2p2.scaled),
      bioclimPCA.phyloPCOA = lm(ln.kd ~ eco.PC1 + eco.PC2 + eco.PC3 + eco.PC4 + phylo.Axis.1 + phylo.Axis.2 + phylo.Axis.3 + phylo.Axis.4, d2p2.scaled),
      phyloPCOA = lm(ln.kd ~ phylo.Axis.1 + phylo.Axis.2 + phylo.Axis.3 + phylo.Axis.4, d2p2.scaled),
      bioclim = lm(ln.kd ~ bio1 + bio12 + bio4 + bio14, d2p2.scaled),
      bioclimMDS = lm(ln.kd ~ eco.MDS1 + eco.MDS2 + eco.MDS3 + eco.MDS4, d2p2.scaled),
      bioclimPCA = lm(ln.kd ~ eco.PC1 + eco.PC2 + eco.PC3 + eco.PC4, d2p2.scaled)
    )

### by species, with phylo predictors -- angiosperms only

d2p2.angio <- d2p2.scaled[taxa.angio, -grep('phylo.Axis', names(d2p2.scaled))]
d2p2.angio <- cbind(d2p2.angio, phylo.pcoa.angio$vectors)
d2p2.angio.scaled <- scale(d2p2.angio) %>%
  as.data.frame

d2pcComp$angiosperms$bySp <- list(
  bioclim.phyloPCOA = lm(ln.kd ~ bio1 + bio12 + bio4 + bio14 + phylo.Axis.1 + phylo.Axis.2 + phylo.Axis.4 + phylo.Axis.5, d2p2.angio.scaled),
  bioclimMDS.phyloPCOA = lm(ln.kd ~ eco.MDS1 + eco.MDS2 + eco.MDS3 + eco.MDS4 + phylo.Axis.1 + phylo.Axis.2 + phylo.Axis.4 + phylo.Axis.5, d2p2.angio.scaled),
  bioclimPCA.phyloPCOA = lm(ln.kd ~ eco.PC1 + eco.PC2 + eco.PC3 + eco.PC4 + phylo.Axis.1 + phylo.Axis.2 + phylo.Axis.4 + phylo.Axis.5, d2p2.angio.scaled),
  phyloPCOA = lm(ln.kd ~ phylo.Axis.1 + phylo.Axis.2 + phylo.Axis.4 + phylo.Axis.5, d2p2.angio.scaled),
  bioclim = lm(ln.kd ~ bio1 + bio12 + bio4 + bio14, d2p2.angio.scaled),
  bioclimMDS = lm(ln.kd ~ eco.MDS1 + eco.MDS2 + eco.MDS3 + eco.MDS4, d2p2.angio.scaled),
  bioclimPCA = lm(ln.kd ~ eco.PC1 + eco.PC2 + eco.PC3 + eco.PC4, d2p2.angio.scaled)
  )

temp.spp <- dat.2014.bySpSite$tidiedName
d2p2.spSite <-
  cbind(dat.2014.bySpSite,
    eco.pca$x[temp.spp, ]) %>%
  as.data.frame

## not changed 2018-05-08: already scaled here.
d2p2.spSite <-
  d2p2.spSite[, grep('bio|eco|phylo|kd', names(d2p2.spSite))] %>%
  scale %>%
  as.data.frame
attr(d2p2.spSite, 'spp')  <- temp.spp
d2p2.spSite$site <- as.factor(gsub(paste(c('[',letters,']'), collapse = ''), '',row.names(d2p2.spSite)))
d2p2.spSite$tipLabel <- gsub('[-.0123456789]', '', row.names(d2p2.spSite))
d2p2.spSite$tidiedName <- temp.spp

d2pcComp$allTaxa$bySpSite <- list(
  bioclim.phyloPCOA = lm(ln.kd ~ bio1 + bio12 + bio4 + bio14 + phylo.Axis.1 + phylo.Axis.2 + phylo.Axis.3 + phylo.Axis.4, d2p2.spSite),
  bioclimMDS.phyloPCOA = lm(ln.kd ~ eco.MDS1 + eco.MDS2 + eco.MDS3 + eco.MDS4 + phylo.Axis.1 + phylo.Axis.2 + phylo.Axis.3 + phylo.Axis.4, d2p2.spSite),
  bioclimPCA.phyloPCOA = lm(ln.kd ~ eco.PC1 + eco.PC2 + eco.PC3 + eco.PC4 + phylo.Axis.1 + phylo.Axis.2 + phylo.Axis.3 + phylo.Axis.4, d2p2.spSite),
  phyloPCOA = lm(ln.kd ~ phylo.Axis.1 + phylo.Axis.2 + phylo.Axis.3 + phylo.Axis.4, d2p2.spSite),
  bioclim = lm(ln.kd ~ bio1 + bio12 + bio4 + bio14, d2p2.spSite),
  bioclimMDS = lm(ln.kd ~ eco.MDS1 + eco.MDS2 + eco.MDS3 + eco.MDS4, d2p2.spSite),
  bioclimPCA = lm(ln.kd ~ eco.PC1 + eco.PC2 + eco.PC3 + eco.PC4, d2p2.spSite)
  )

### by species, with phylo predictors -- angiosperms only
cols.angio <- which(attr(d2p2.spSite, 'spp') %in% taxa.angio)
d2p2.spSite.angio <-
  d2p2.spSite[cols.angio, -grep('phylo.Axis', names(d2p2.spSite))]

d2p2.spSite.angio.text <- d2p2.spSite.angio[c("site","tipLabel","tidiedName")]
d2p2.spSite.angio[c("site","tipLabel","tidiedName")] <- NULL

d2p2.spSite.angio <-
  cbind(d2p2.spSite.angio,
        phylo.pcoa.angio$vectors[attr(d2p2.spSite, 'spp')[cols.angio], ]) %>%
  scale %>%
  as.data.frame
attr(d2p2.spSite.angio, 'spp') <-
  attr(d2p2.spSite, 'spp')[which(attr(d2p2.spSite, 'spp') %in% taxa.angio)]
d2p2.spSite.angio <- cbind(d2p2.spSite.angio, d2p2.spSite.angio.text)

d2pcComp$angiosperms$bySpSite <- list(
  bioclim.phyloPCOA = lm(ln.kd ~ bio1 + bio12 + bio4 + bio14 + phylo.Axis.1 + phylo.Axis.2 + phylo.Axis.3 + phylo.Axis.4, d2p2.spSite.angio),
  bioclimMDS.phyloPCOA = lm(ln.kd ~ eco.MDS1 + eco.MDS2 + eco.MDS3 + eco.MDS4 + phylo.Axis.1 + phylo.Axis.2 + phylo.Axis.3 + phylo.Axis.4, d2p2.spSite.angio),
  bioclimPCA.phyloPCOA = lm(ln.kd ~ eco.PC1 + eco.PC2 + eco.PC3 + eco.PC4 + phylo.Axis.1 + phylo.Axis.2 + phylo.Axis.3 + phylo.Axis.4, d2p2.spSite.angio),
  phyloPCOA = lm(ln.kd ~ phylo.Axis.1 + phylo.Axis.2 + phylo.Axis.3 + phylo.Axis.4, d2p2.spSite.angio),
  bioclim = lm(ln.kd ~ bio1 + bio12 + bio4 + bio14, d2p2.spSite.angio),
  bioclimMDS = lm(ln.kd ~ eco.MDS1 + eco.MDS2 + eco.MDS3 + eco.MDS4, d2p2.spSite.angio),
  bioclimPCA = lm(ln.kd ~ eco.PC1 + eco.PC2 + eco.PC3 + eco.PC4, d2p2.spSite.angio)
  )
