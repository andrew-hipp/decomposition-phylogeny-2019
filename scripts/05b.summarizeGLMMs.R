## make a table summarizing GLMM results

roundTo = 3

preds <- c('bio1', 'bio12', 'bio4', 'bio14')
gcovs <- c('site', 'tipLabel')
headers <- c('Taxon set', 'Analysis level', 'Model',
             'DIC', 'delta DIC', 'Sum DICw',
             'Residuals var',
             'Site var', 'Phylo var',
             preds)

glmmSummary <- matrix(NA, 0, length(headers), dimnames = list(NULL, headers))

for(dataSet in c('allTaxa', 'angiosperms')) {
  for(analysisLevel in c('bySp', 'bySpSite')) {
    for(model in names(dat.2014.glmm[[dataSet]][[analysisLevel]])) {
      mTemp <- summary(dat.2014.glmm[[dataSet]][[analysisLevel]][[model]])
      if(is.null(mTemp$Gcovariances)) mTemp$Gcovariances <- data.frame(x=1, y=2, z=3) # create dummy data.frame if no random effects specified
      glmmSummary <- rbind(glmmSummary,
                            c(dataSet, analysisLevel, model,
                              DIC = round(mTemp$DIC, roundTo),
                              Delta.DIC = NA,
                              Sum.DICw = NA,
                              Residual.cov = mTemp$Rcovariances[1:3] %>%
                                round(roundTo) %>%
                                format(digits = roundTo) %>%
                                (function(x) paste(x[1], ' (', x[2], ',',
                                      x[3], ')', sep = '')),
                              as.data.frame(mTemp$Gcovariances)[gcovs, 1:3] %>%
                                round(roundTo) %>%
                                format(digits = roundTo) %>%
                                apply(., 1, function(x) paste(x[1],
                                      ' (', x[2], ',',
                                      x[3],')', sep = '')),
                              as.data.frame(mTemp$solutions)[preds, c(1:3, 5)] %>%
                                round(roundTo) %>%
                                format(digits = roundTo) %>%
                                apply(., 1, function(x) paste(x[1],
                                                              ' (',
                                                              x[2],
                                                              ',',
                                                              x[3],
                                                              ', p = ',
                                                              x[4],
                                                              ')',
                                                              sep = ''))
                              ) # close c
                            ) # close rbind
    } # close for model
  } # close for analysisLevel
} # close for dataSet

glmmSummary[grep('NA ', glmmSummary)] <- '-'
glmmSummary <- gsub('p = 0)', paste('p < ', 10 ^ (-roundTo), ')', sep = ''), glmmSummary)
glmmSummary <- glmmSummary[order(glmmSummary[, 'Taxon set'],
glmmSummary[, 'Analysis level'],
as.numeric(glmmSummary[,'DIC'])), ]

dataset = apply(glmmSummary[, 1:2], 1, paste, collapse = '')
for(i in unique(dataset)) {
  glmmSummary[which(dataset == i), 'delta DIC'] <-
    as.numeric(glmmSummary[which(dataset == i), 'DIC']) - min(as.numeric(glmmSummary[which(dataset == i), 'DIC'])) %>%
    round(roundTo)
  glmmSummary[which(dataset == i), 'Sum DICw'] <-
    as.numeric(glmmSummary[which(dataset == i), 'DIC']) %>%
    aic.w %>%
    cumsum %>%
    round(roundTo)
  }

write.csv(glmmSummary, '../out/TABLE.03.glmmSummary.csv')
