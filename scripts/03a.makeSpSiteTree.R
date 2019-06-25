## make a phylogeny for the sp*site GLMMs
## creates tree with one tip per row, making spp tips a designated maximum sp age
## AH 2018-09-08; currently doesn't work under Linux, though I'm not sure why
##   -- I am currently running this script in windows and then loading data
##      in the next script

require(phytools)

# arguments to make.spSite:
#tr : phylo object
#dat : data frame with
#  (1) column tipLabel equivalent to tip labels in tr
#  (2) row.names equivalent to the new labels you want for your fancy new tree


make.spSite <- function(tr, dat,
                        maxSpAge = 1.5,
                        brScalar = 0.999,
                        makeDi = TRUE,
                        verbose = FALSE,
                        do.tidy.tree = TRUE) {
  ARGS <- as.list(match.call())
  tr.spSite <- tr
  if(do.tidy.tree) tr.spSite$tip.label <- tidyName(tr.spSite$tip.label)
  for(i in seq(dim(dat)[1])) {
    node.temp <- which(tr.spSite$tip.label == dat$tipLabel[i])
    new.label <- row.names(dat)[i]
    pos.temp <- min(tr.spSite$edge.length[which(tr.spSite$edge[, 2] == node.temp)] * brScalar,
                    maxSpAge)
    if(verbose) message(paste('adding tip', new.label, 'with age position', round(pos.temp, 1), 'Ma'))
    tr.spSite <- bind.tip(tr.spSite, new.label,
                                where = node.temp,
                                position = pos.temp)
  }
  tr.spSite <- drop.tip(tr.spSite, unique(dat$tipLabel))
  if(makeDi) tr.spSite <- multi2di(tr.spSite)
  dat$tipLabel <- row.names(dat)
  write.tree(tr.spSite, paste('../out/', ARGS$tr, '.spSite.tre', sep = ''))
  message(paste('\nOn tree', ARGS$tr, 'and dataset', ARGS$dat, '---'))
  message(paste("Tree length =", length(tr.spSite$tip.label)))
  message(paste("Dataset rows =", dim(dat)[1]))
  message(paste("Data overlap = ",
            length(intersect(tr.spSite$tip.label, row.names(dat)))))
  out <- list(tr = tr.spSite, dat = dat)
  return(out)
}

datList.angio.spSite <- make.spSite(tr = tr.angio.edit,
                              dat = dat.2014.bySpSite.angio,
                              maxSpAge = 1.5,
                              brScalar = 0.999,
                              makeDi = TRUE,
                              verbose = TRUE,
                              do.tidy.tree = TRUE)

datList.spSite <- make.spSite(tr = tr.edit,
                              dat = dat.2014.bySpSite,
                              maxSpAge = 1.5,
                              brScalar = 0.999,
                              makeDi = TRUE,
                            verbose = TRUE,
                            do.tidy.tree = TRUE)

save(datList.angio.spSite, datList.spSite, file='../out/datLists.spSite.Rdata')
