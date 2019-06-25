library(phytools)
library(openxlsx)
library(magrittr)
source('../scripts/decomp.fcts.R')

### NEED TO ADD EDGE LENGTHS!!!

tr <- read.tree('../data/tree.newages.tidy.final.tre')
tr.names <- read.xlsx('../data/addTaxa.xlsx', 1)
for(i in seq(dim(tr.names)[1])) {
  if(tr.names$rule[i] == 'ignore') next
  tax.add <- tr.names$taxonToAdd[i]
  if(tr.names$rule[i] == 'sister') {
    message(paste('adding a sister named', tax.add))
    node <-
    tr <- bind.tip(tr,
                   tidyName(tax.add),
                   where = which(tr$tip.label == tr.names$lookup[i]),
                   position = tr.names$position[i],
                   edge.length = tr.names$position[i]
                 )
  }
  if(tr.names$rule[i] == 'mrca') {
    message(paste('adding a taxon named', tax.add, 'to an mrca'))
    nodeTax <- strsplit(tr.names$lookup[i], '|', fixed = T)[[1]]
    tr <- bind.tip(tr, tidyName(tax.add),
          where = fastMRCA(tr, nodeTax[1], nodeTax[2]),
          position = tr.names$position[i],
          edge.length = max(as.numeric(nodeHeights(tr))) -
                            fastHeight(tr, nodeTax[1], nodeTax[2]) +
                            tr.names$position[i]
                  ) # close bind.tip
  }
}

write.tree(tr, '../data/tree.newages.tidy.final.v2.tre')
