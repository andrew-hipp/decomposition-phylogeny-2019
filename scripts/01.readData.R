## command lines for Nov 2014 decomp phylogeny
## final version, 2018-09-20:
###  - add three remaining taxa lats longs (Carri to do)
###  - clean up lat longs and duplicates (Carri to do)
###  - clean up site names
###  - family averages are first from species averages

library(ape)
library(nlme)
library(magrittr)
library(sp)
library(phangorn)
library(phytools)

taxa.outsies <- c("blechnumnudum", "dacrycarpusdacrydioides", "larixdecidua",
                  "piceaabies", "pinuspinaster", "pinusponderosa", "pinusstrobus",
                  "pinussylvestris", "ptisanasalicina", "tsugacanadensis", "tsugaheterophylla",
                  'cryptomeriajaponica'
                  )

source('../scripts/decomp.fcts.R')

dat.2014 <- read.csv('../data/k-values_with_select_env_variables_CJL_ALH_9-26-18.csv', as.is = T)
tr <- tr.orig <- read.tree('../data/tree.newages.tidy.final.v2.tre')

spp.matrix <- read.csv('../data/spp.matrix.csv', as.is = T) # this data matrix includes output from TNRS


## clean up accepted names
dat.2014.meta <- read.csv('../data/k-values_with_select_env_variables.colNames.csv',
                            as.is = T,
                            row.names = 1)
names(dat.2014) <- dat.2014.meta[names(dat.2014), 'name.new']
dat.2014$tidyNameOrig <- tidyName(nameAssign(dat.2014$Scientific.Name))

spp.matrix$Accepted_name <- gsub('subsp.', 'ssp.', spp.matrix$Accepted_name)
spp.matrix$Accepted_name <- gsub(' $', ' sp.', spp.matrix$Accepted_name)

dat.2014$newName <- dat.2014$tidiedName <- nameAssign(dat.2014$Scientific.Name)

tr$tip.label <- dat.2014$newName[match(tr$tip.label, dat.2014$tidyNameOrig)]
taxa.outsies <- dat.2014$newName[match(taxa.outsies, dat.2014$tidyNameOrig)]

### add on worldclim data
if(!exists('worldclim')) {
  library(raster)
  library(sp)
  worldclim <- raster:::getData("worldclim",var="bio",res=10)
}
dat.2014 <- cbind(dat.2014, extract(worldclim, SpatialPoints(dat.2014[ , c('long', 'lat')])))

discard.noPhylo <- !dat.2014$tidiedName %in% tr$tip.label
discard.duplicate <- duplicated(dat.2014$Data.point.code)
discard.noBio <- is.na(dat.2014$bio1)
discard.all <- discard.noPhylo | discard.duplicate | discard.noBio
dat.2014.discarded <- cbind(noPhylo = discard.noPhylo,
                            duplicate_Data.point.code = discard.duplicate,
                            noBioclim = discard.noBio,
                          dat.2014)
dat.2014.discarded <- dat.2014.discarded[discard.all, ]
dat.2014 <- dat.2014[!discard.all, ]
write.csv(dat.2014.discarded, 'dat.2014.discardedRows.csv')
write.csv(dat.2014, 'dat.2014.retainedRows.csv')

dat.2014.means <- as.data.frame(t(sapply(split(dat.2014[dat.2014.meta$numeric], dat.2014$tidiedName), function(x) apply(x, 2, function(y) mean(as.numeric(y), na.rm = T)))))
dat.2014.means$ln.kd <- log(dat.2014.means$kd)

indsAll <- which(dat.2014$kd != 0 &
             !is.na(dat.2014$kd) &
               !is.na(dat.2014$bio1))

## added 2018-09-20

siteSpVec <- paste(tidyName(dat.2014$tidiedName),
                     dat.2014$lat, dat.2014$long,
                     sep = '')

### dataset by species and site
dat.2014.bySpSite <-
  cbind(dat.2014,
        siteSp = siteSpVec) %>%
  split(f = .$siteSp) %>%
  sapply(., function(x) apply(x, 2, function(y) mean(as.numeric(y), na.rm = T))) %>%
  t

dat.2014.bySpSite <- cbind(dat.2014.bySpSite,
                            dat.2014[match(row.names(dat.2014.bySpSite), siteSpVec),
                                      c('Scientific.Name', 'tidiedName', 'Genus', 'Family', 'Order')
                                      ]
                            ) %>% as.data.frame

temp.spp <- dat.2014.bySpSite$tidiedName
attr(dat.2014.bySpSite, 'spp')  <- temp.spp

dat.2014.bySpSite[which(is.na(dat.2014.bySpSite[1,]))] <- NULL
dat.2014.bySpSite <- dat.2014.bySpSite[which(dat.2014.bySpSite$kd > 0 & !is.na(dat.2014.bySpSite$bio1)), ]
dat.2014.bySpSite$ln.kd <- log(dat.2014.bySpSite$kd)
dat.2014.bySpSite$site <- as.factor(gsub(paste(c('[',letters,']'), collapse = ''), '',row.names(dat.2014.bySpSite)))
dat.2014.bySpSite$tipLabel <- gsub('[-.0123456789]', '', row.names(dat.2014.bySpSite))

## dataset by species
dat.2014.bySp <-
  split(dat.2014, dat.2014$tidiedName) %>%
  sapply(., function(x) apply(x, 2, function(y) mean(as.numeric(y), na.rm = T))) %>%
  t %>%
  as.data.frame
dat.2014.bySp <- dat.2014.bySp[-which(is.na(dat.2014.bySp[1, ]))]
if(any(is.na(dat.2014.bySp[ , 'bio1']))) {
  dat.2014.bySp <- dat.2014.bySp[-which(is.na(dat.2014.bySp[ , 'bio1'])), ]
}
dat.2014.bySp$ln.kd <- log(dat.2014.bySp$kd)

indsMult <- which(dat.2014.bySpSite$tidiedName %in% names(which(table(dat.2014.bySpSite$tidiedName)>2)))

tr <- drop.tip(tr, which(!tr$tip.label %in% row.names(dat.2014.bySp)))
taxa.angio <- setdiff(intersect(tr$tip.label, row.names(dat.2014.bySp)),
                                taxa.outsies)
tr.angio <- drop.tip(tr, which(!tr$tip.label %in% taxa.angio))
tr.edit <- nnls.tree(cophenetic(tr), tr, rooted = T)
tr.edit$node.label <- NULL
tr.angio.edit <- nnls.tree(cophenetic(tr.angio), tr.angio, rooted = T)
tr.angio.edit$node.label <- NULL

write.tree(tr.edit, '../out/SUPPL.finalSpTree.allTaxa.tre')
write.tree(tr.angio.edit, '../out/SUPPL.finalSpTree.angiosperms.tre')

## added 2018-09-09
rndDig = 4
dat.2014.splBySp <- split(dat.2014, dat.2014$tidiedName)
dat.2014.matBySp <- data.frame(
    kd.mean = sapply(dat.2014.splBySp, function(x) round(mean(x$kd), rndDig)),
    N.samples = sapply(dat.2014.splBySp, function(x) dim(x)[1]),
    kd.sd = sapply(dat.2014.splBySp, function(x) {
      ifelse(dim(x)[1] > 1, round(sd(x$kd), rndDig), '--')
      }),
    kd.SEM = sapply(dat.2014.splBySp, function(x) {
      ifelse(dim(x)[1] > 1, round(sd(x$kd) / sqrt(dim(x)[1]), rndDig), '--')
      }),
    Family = sapply(dat.2014.splBySp, function(x) x$Family[1]),
    Order = sapply(dat.2014.splBySp, function(x) x$Order[1])
  )
dat.2014.splByFam <- split(dat.2014.matBySp, dat.2014.matBySp$Family)
dat.2014.matByFam <- data.frame(
    kd.mean = sapply(dat.2014.splByFam, function(x) round(mean(x$kd.mean), rndDig)),
    N.species = sapply(dat.2014.splByFam, function(x) dim(x)[1]),
    kd.sd = sapply(dat.2014.splByFam, function(x) round(sd(x$kd.mean), rndDig)),
    kd.SEM = sapply(dat.2014.splByFam, function(x) {
      ifelse(dim(x)[1] > 1, round(sd(x$kd.mean) / sqrt(dim(x)[1]), rndDig), '--')
      }),
    Order = sapply(dat.2014.splByFam, function(x) x$Order[1]),
    Clade = factor(
      dat.2014$Clade[match(names(dat.2014.splByFam), dat.2014$Family)],
      levels = c("Fern/fern allies", "Gymnosperm", "Monocot", "Magnoliid", "Eudicot")
    )
  )
dat.2014.matByFam <- dat.2014.matByFam[with(dat.2014.matByFam,
                                            order(Clade, Order)),
                                      ]

write.csv(dat.2014.matBySp, '../out/TABLE.kd.meansBySpecies.csv')
write.csv(dat.2014.matByFam, '../out/TABLE.kd.meansByFamily.csv')

fam.table <- aggregate(kd ~ Family, data = dat.2014, FUN = function(x) c(round(mean(x), rndDig),
                                                                         round(sd(x), rndDig),
                                                                         paste(round(range(x), rndDig), collapse = '--')
                                                                       ),
                       simplify = T)
fam.table <- data.frame(fam.table[[2]], row.names = fam.table$Family)
names(fam.table) <- c('mean.by.sample', 'sd.by.sample', 'range.by.sample')
fam.table$'Num. taxa' <- sapply(row.names(fam.table), function(x) {
  length(unique(dat.2014$tidiedName[which(dat.2014$Family == x)]))
})
fam.table$'Num. samples' <- sapply(row.names(fam.table), function(x) {
  sum(dat.2014$Family == x)
})
fam.table$Order <- dat.2014$Order[match(row.names(fam.table), dat.2014$Family)]
fam.table$Clade <- factor(
  dat.2014$Clade[match(row.names(fam.table), dat.2014$Family)],
  levels = c("Fern/fern allies", "Gymnosperm", "Monocot", "Magnoliid", "Eudicot")
)
fam.table$sd <- as.character(fam.table$sd)
fam.table$sd[is.na(fam.table$sd)] <- '-'

fam.table <- fam.table[c('Clade', 'Order', 'Num. taxa', 'Num. samples', 'mean.by.sample', 'sd.by.sample', 'range.by.sample')]
fam.table <- fam.table[with(fam.table, order(Clade, Order)), ]

write.csv(fam.table, '../out/SUPPL.ONLY.FOR.COMPARISON.kdByFamily.BYSAMPLENOTSP.csv')

## create a silly thing that is required by later scripts b/c of how I wrote them originally
spp.all <- unique(dat.2014$tidiedName) %>% sort
names(spp.all) <- spp.all

## MAKE STATS FOR METHODS
dat.stats <- c(
  "TOTAL SITES",
  paste(dat.2014$lat, dat.2014$long) %>% unique %>% length,
  "TOTAL RECORDS",
  dim(dat.2014)[1],
  "TOTAL TAXA",
  unique(dat.2014$tidiedName) %>% length,
  "TOTAL ORDERS",
  unique(dat.2014$Order) %>% length,
  "TOTAL GENERA",
  unique(dat.2014$Genus) %>% length,
  "TOTAL FAMILIES",
  unique(dat.2014$Family) %>% length
)

writeLines(dat.stats, '../out/dat.stats.txt')

## summary of phylogenetic diversity by site
dat.bySiteDiversity <- table(dat.2014.bySpSite$site) %>% table
