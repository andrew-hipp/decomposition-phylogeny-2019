## make a table summarizing lm results

roundTo = 3

counter.coef <- c(lapply(d2pcComp[[1]][[1]], coef),
                  lapply(d2pcComp[[2]][[1]], coef)) %>%
              lapply(., names) %>%
              unlist %>%
              unique %>%
              as.character %>%
              sort %>%
              .[-grep('Intercept', .)]

counter.rows <- matrix(NA, 0, 3, dimnames = list(NULL, c('Taxon set', 'Analysis level', 'Model')))

for(taxa in c('allTaxa', 'angiosperms')) {
  for(analysis in c('bySp', 'bySpSite')) {
    for(model in names(d2pcComp[[1]][[1]])) {
      counter.rows <- rbind(counter.rows, c(taxa, analysis, model))
    }
  }
}

counter.cols.extra <- c('r2', 'AIC', 'delta AIC', 'Sum AICw')

d2pcComp.lmSummaries <- cbind(counter.rows,
                              replicate(length(counter.coef) + length(counter.cols.extra), rep('-', dim(counter.rows)[1]))
                            )

dimnames(d2pcComp.lmSummaries)[[2]] <- c(dimnames(counter.rows)[[2]], counter.cols.extra, counter.coef)

for(i in 1:dim(d2pcComp.lmSummaries)[1]){
    message(paste('doing', i, 'of', dim(d2pcComp.lmSummaries)[1]))
    taxa <- d2pcComp.lmSummaries[i, 'Taxon set']
    analysis <- d2pcComp.lmSummaries[i, 'Analysis level']
    model <- d2pcComp.lmSummaries[i, 'Model']
    message(paste(taxa,analysis,model, sep = ' -- '))
    d2pcComp.lmSummaries[i, counter.coef] <- paste(round(as.data.frame(summary(d2pcComp[[taxa]][[analysis]][[model]])$coefficients)[counter.coef, 1], roundTo),
                                                 ' (p = ', round(as.data.frame(summary(d2pcComp[[taxa]][[analysis]][[model]])$coefficients)[counter.coef, 4], roundTo),
                                                 ')', sep = '')
    d2pcComp.lmSummaries[i, 'r2'] <- round(summary(d2pcComp[[taxa]][[analysis]][[model]])$r.squared, roundTo)
    d2pcComp.lmSummaries[i, 'AIC'] <- round(AIC(d2pcComp[[taxa]][[analysis]][[model]]), roundTo)
    d2pcComp.lmSummaries[i, grep('NA ', d2pcComp.lmSummaries[i, ])] <- '-'
}

d2pcComp.lmSummaries <- d2pcComp.lmSummaries[order(d2pcComp.lmSummaries[, 'Taxon set'],
d2pcComp.lmSummaries[, 'Analysis level'],
as.numeric(d2pcComp.lmSummaries[,'AIC'])), ]

dataset = apply(d2pcComp.lmSummaries[, 1:2], 1, paste, collapse = '')
for(i in unique(dataset)) {
  d2pcComp.lmSummaries[which(dataset == i), 'delta AIC'] <-
    as.numeric(d2pcComp.lmSummaries[which(dataset == i), 'AIC']) - min(as.numeric(d2pcComp.lmSummaries[which(dataset == i), 'AIC'])) %>%
    round(roundTo)
  d2pcComp.lmSummaries[which(dataset == i), 'Sum AICw'] <-
    as.numeric(d2pcComp.lmSummaries[which(dataset == i), 'AIC']) %>%
    aic.w %>%
    cumsum %>%
    round(roundTo)
  }

d2p2.angio.write <- d2p2.angio
row.names(d2p2.angio.write) <- spp.all[row.names(d2p2.angio.write)]
d2p2.write <- d2p2
row.names(d2p2.write) <- spp.all[row.names(d2p2.write)]
write.csv(d2pcComp.lmSummaries, '../out/TABLE.02-supp-pcaCompare.linear.model.summaries.csv')
