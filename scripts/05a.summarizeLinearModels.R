## make a table summarizing lm results

roundTo = 3

counter.coef <- c(lapply(dat.2014.lm[[1]][[1]], coef),
                  lapply(dat.2014.lm[[2]][[1]], coef)) %>%
              lapply(., names) %>%
              unlist %>%
              unique %>%
              as.character %>%
              sort %>%
              .[-grep('Intercept', .)]

counter.rows <- matrix(NA, 0, 3, dimnames = list(NULL, c('Taxon set', 'Analysis level', 'Model')))

for(taxa in c('allTaxa', 'angiosperms')) {
  for(analysis in c('bySp', 'bySpSite')) {
    for(model in names(dat.2014.lm[[1]][[1]])) {
      counter.rows <- rbind(counter.rows, c(taxa, analysis, model))
    }
  }
}

counter.cols.extra <- c('r2', 'AIC', 'delta AIC', 'Sum AICw')

dat.2014.lmSummaries <- cbind(counter.rows,
                              replicate(length(counter.coef) + length(counter.cols.extra), rep('-', dim(counter.rows)[1]))
                            )

dimnames(dat.2014.lmSummaries)[[2]] <- c(dimnames(counter.rows)[[2]], counter.cols.extra, counter.coef)

for(i in 1:dim(dat.2014.lmSummaries)[1]){
    taxa <- dat.2014.lmSummaries[i, 'Taxon set']
    analysis <- dat.2014.lmSummaries[i, 'Analysis level']
    model <- dat.2014.lmSummaries[i, 'Model']
    dat.2014.lmSummaries[i, counter.coef] <- paste(round(as.data.frame(summary(dat.2014.lm[[taxa]][[analysis]][[model]])$coefficients)[counter.coef, 1], roundTo),
                                                 ' (p = ', round(as.data.frame(summary(dat.2014.lm[[taxa]][[analysis]][[model]])$coefficients)[counter.coef, 4], roundTo),
                                                 ')', sep = '')
    dat.2014.lmSummaries[i, 'r2'] <- round(summary(dat.2014.lm[[taxa]][[analysis]][[model]])$r.squared, roundTo)
    dat.2014.lmSummaries[i, 'AIC'] <- round(AIC(dat.2014.lm[[taxa]][[analysis]][[model]]), roundTo)
    dat.2014.lmSummaries[i, grep('NA ', dat.2014.lmSummaries[i, ])] <- '-'
}

dat.2014.lmSummaries <- dat.2014.lmSummaries[order(dat.2014.lmSummaries[, 'Taxon set'],
dat.2014.lmSummaries[, 'Analysis level'],
as.numeric(dat.2014.lmSummaries[,'AIC'])), ]

dataset = apply(dat.2014.lmSummaries[, 1:2], 1, paste, collapse = '')
for(i in unique(dataset)) {
  dat.2014.lmSummaries[which(dataset == i), 'delta AIC'] <-
    as.numeric(dat.2014.lmSummaries[which(dataset == i), 'AIC']) - min(as.numeric(dat.2014.lmSummaries[which(dataset == i), 'AIC'])) %>%
    round(roundTo)
  dat.2014.lmSummaries[which(dataset == i), 'Sum AICw'] <-
    as.numeric(dat.2014.lmSummaries[which(dataset == i), 'AIC']) %>%
    aic.w %>%
    cumsum %>%
    round(roundTo)
  }

dat.sp.angio.write <- dat.2014.bySp.angio
row.names(dat.sp.angio.write) <- spp.all[row.names(dat.sp.angio.write)]
dat.sp.write <- dat.2014.bySp
row.names(dat.sp.write) <- spp.all[row.names(dat.sp.write)]
write.csv(dat.2014.lmSummaries, '../out/TABLE.02.linear.model.summaries.csv')

write.csv(dat.sp.write, '../out/dat.2014.bySp.allTaxa.csv')
write.csv(cbind(datList.spSite$dat, fullName = spp.all[datList.spSite$dat$tidiedName]),
          '../out/dat.2014.bySpSite.allTaxa.csv')
write.csv(dat.sp.angio.write, '../out/dat.2014.bySp.angiosperms.csv')
write.csv(cbind(datList.angio.spSite$dat, fullName = spp.all[datList.angio.spSite$dat$tidiedName]),
          '../out/dat.2014.bySpSite.angiosperms.csv')
