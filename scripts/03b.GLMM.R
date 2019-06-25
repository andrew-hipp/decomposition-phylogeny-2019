library(MCMCglmm)
library(phangorn)
library(ape)
library(nlme)

#source('https://raw.githubusercontent.com/andrew-hipp/morton/master/R/gls.r.squared.R')
#load('../out/')
dat.2014.glmm <- list(
  allTaxa = list(
    bySp = 0,
    bySpSite = 0
  ),
  angiosperms = list(
    bySp = 0,
    bySpSite = 0
  )
)

dat.2014.bySp.angio$tipLabel <- row.names(dat.2014.bySp.angio)

dat.2014.bySp$tipLabel <- row.names(dat.2014.bySp)

for(dataSet in c('allTaxa', 'angiosperms')) {
  for(analysisLevel in c('bySp', 'bySpSite')) {
    message(paste('... MCMC for', dataSet, 'with grouping', analysisLevel))
    if(dataSet == 'allTaxa' & analysisLevel == 'bySp') {
      datTemp <- dat.2014.bySp
      trTemp <- tr.edit}
    if(dataSet == 'angiosperms' & analysisLevel == 'bySp') {
      datTemp <- dat.2014.bySp.angio
      trTemp <- tr.angio.edit}
    if(dataSet == 'allTaxa' & analysisLevel == 'bySpSite') {
      datTemp <- datList.spSite$dat
      trTemp <- datList.spSite$tr}
    if(dataSet == 'angiosperms' & analysisLevel == 'bySpSite') {
      datTemp <- datList.angio.spSite$dat
      trTemp <- datList.angio.spSite$tr}
    if(analysisLevel == 'bySp') {
      dat.2014.glmm[[dataSet]][[analysisLevel]] <-
      list(
          noPred = MCMCglmm(ln.kd ~ 1, data = datTemp,
            nitt = 50000, burnin = 10000, verbose = F),
          eco = MCMCglmm(ln.kd ~ bio1 + bio4 + bio12 + bio14, data = datTemp,
            nitt = 50000, burnin = 10000, verbose = F),
          phy = MCMCglmm(ln.kd ~ 1, data = datTemp, random = ~tipLabel,
            ginverse = list(tipLabel = inverseA(trTemp)$Ainv),
            nitt = 50000, burnin = 10000, verbose = F),
          eco.phy = MCMCglmm(ln.kd ~ bio1 + bio4 + bio12 + bio14, data = datTemp, random = ~tipLabel,
            ginverse = list(tipLabel = inverseA(trTemp)$Ainv),
            nitt = 50000, burnin = 10000, verbose = F)
          )
        }
      else {
        dat.2014.glmm[[dataSet]][[analysisLevel]] <-
        list(
          noPred = MCMCglmm(ln.kd ~ 1, data = datTemp,
            nitt = 50000, burnin = 10000, verbose = F),
          phy = MCMCglmm(ln.kd ~ 1, data = datTemp, random = ~tipLabel,
                         ginverse = list(tipLabel = inverseA(trTemp)$Ainv),
                         nitt = 50000, burnin = 10000, verbose = F),
          eco = MCMCglmm(ln.kd ~ bio1 + bio4 + bio12 + bio14, data = datTemp,
                          nitt = 50000, burnin = 10000, verbose = F),
          site = MCMCglmm(ln.kd ~ 1, data = datTemp, random = ~site,
                          nitt = 50000, burnin = 10000, verbose = F),
          eco.site = MCMCglmm(ln.kd ~ bio1 + bio4 + bio12 + bio14, data = datTemp,
                              random = ~site,
                              nitt = 50000, burnin = 10000, verbose = F),
          phy.site = MCMCglmm(ln.kd ~ 1, data = datTemp,
                                  random = ~tipLabel + site,
                                  ginverse = list(tipLabel = inverseA(trTemp)$Ainv),
                                  nitt = 50000, burnin = 10000, verbose = F),
          phy.eco = MCMCglmm(ln.kd ~ bio1 + bio4 + bio12 + bio14, data = datTemp, random = ~tipLabel,
                              ginverse = list(tipLabel = inverseA(trTemp)$Ainv),
                              nitt = 50000, burnin = 10000, verbose = F),
          phy.eco.site = MCMCglmm(ln.kd ~ bio1 + bio4 + bio12 + bio14, data = datTemp,
                                  random = ~tipLabel + site,
                                  ginverse = list(tipLabel = inverseA(trTemp)$Ainv),
                                  nitt = 50000, burnin = 10000, verbose = F)
          )
      }
    } # close for analysisLevel
  }# close for dataset
