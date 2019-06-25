## phylogenetic community structure of sites
library(picante)
library(magrittr)

u.sites <- unique(dat.2014.bySpSite$site)
u.taxa <- unique(dat.2014.bySpSite$tipLabel)

dat.bySiteDiversity <- table(dat.2014.bySpSite$site) %>% table

dat.2014.commMat <- matrix(0, length(u.sites), length(u.taxa),
                           dimnames = list(u.sites, u.taxa))
for(i in 1:dim(dat.2014.bySpSite)[1]) {
  dat.2014.commMat[dat.2014.bySpSite[i, 'site'], dat.2014.bySpSite[i, 'tipLabel']] <- 1
}

dat.2014.commMat <- dat.2014.commMat[rowSums(dat.2014.commMat) > 1, ]
dat.2014.commMat <- dat.2014.commMat[, colSums(dat.2014.commMat) > 0]
tr.comm.temp <- tr.edit
tr.comm.temp$tip.label <- tidyName(tr.comm.temp$tip.label)
dat.2014.mntd <- ses.mntd(dat.2014.commMat, cophenetic(tr.comm.temp), runs = 9999)
dat.2014.mpd <- ses.mpd(dat.2014.commMat, cophenetic(tr.comm.temp), runs = 9999)
message(paste('MPD sign, alpha = 0.01:',
              sum(dat.2014.mpd$mpd.obs.p <= 0.005 | dat.2014.mpd$mpd.obs.p >= 0.995)))
message(paste('MNTD sign, alpha = 0.01:',
              sum(dat.2014.mntd$mntd.obs.p <= 0.005 | dat.2014.mntd$mntd.obs.p >= 0.995)))

pdf('../out/Supplement.mpd.mntd.pdf', 8.5, 11)
layout(matrix(1:6, 3))
hist(dat.2014.mpd$mpd.obs, 20, main = "MPD observed",
     xlab = 'MPD observed (x 10^6 yrs)')
hist(dat.2014.mpd$mpd.obs.z, 20, main = "MPD, std effect size",
     xlab = 'MPD s.e.s., observed Z')
hist(dat.2014.mpd$mpd.obs.p, 100, main = "MPD, std effect size (p)",
     xlab = 'MPD s.e.s., p-value')
abline(v = c(0.01, 0.99), col = 'red', lty = 'dashed')
hist(dat.2014.mntd$mntd.obs, 20, main = "MNTD observed",
     xlab = 'MNTD observed (x 10^6 yrs)')
hist(dat.2014.mntd$mntd.obs.z, 20, main = "MNTD, std effect size",
     xlab = 'MNTD s.e.s., observed Z')
hist(dat.2014.mntd$mntd.obs.p, 100, main = "MNTD, std effect size (p)",
     xlab = 'MNTD s.e.s., p-value')
abline(v = c(0.01, 0.99), col = 'red', lty = 'dashed')
dev.off()

## clades per site
site.temp <- paste(dat.2014$lat, dat.2014$long)
site.clades <- lapply(unique(site.temp), function(x) unique(dat.2014$Clade[site.temp == x]))
site.clades <- site.clades[which(sapply(site.clades, length) > 1)]
sapply(site.clades, '%in%', x='Eudicot') %>% print
