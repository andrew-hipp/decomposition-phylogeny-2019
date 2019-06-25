# make decomposition phylogeny

require(ape)
require(ggtree)

paint.nodes <- mrca(tr) %>%
  {c(.['Tamarix sp.', 'Inga spectabilis']
  #  .['phragmitesaustralis', 'phragmitesaustralis'],
  #  .['rhododendron', 'rhododendronponticum']
    )}
paint.colors <- c('blue', 'red'
  #               'orange', 'yellow'
                  )
names(paint.colors) <- as.character(0:1)
paint.labels <- c('Gymnosperms, basal dicots, monocots',
                  'Core eudicots'
                  #'Rhododendron',
                  #'Phragmites australis'
                  )
tr.plotVersion <- groupClade(tr, paint.nodes)
levels(attr(tr.plotVersion, 'group')) <- as.character(0:1)
#attr(tr.plotVersion, 'group')[which(tr$tip.label == 'phragmitesaustralis')] <- '3'

dat.2014.bySp$blank <- -0.05
dat.tr <- dat.2014.bySp[intersect(tr$tip.label, row.names(dat.2014.bySp)),
                         c('ln.kd', 'blank',
                         'bio1', 'bio12', 'bio4', 'bio14', 'blank',
                         'phylo.Axis.1', 'phylo.Axis.2', 'phylo.Axis.3', 'phylo.Axis.4')]
dat.tr$bio1[dat.tr$bio1 < 0] <- 0 ## chopping off lowest temperatures
dat.tr <- apply(dat.tr, 2, function(x) {(x - min(x, na.rm = T))/max(x-min(x, na.rm = T), na.rm = T)})
dat.tr[is.na(dat.tr)] <- -0.05
dimnames(dat.tr)[[2]] <- c('ln(Kd)',
                           ' ',
                           'Mean temperature', 'Mean precipitation',
                           'Temperature seasonality', 'Precipitation seasonality',
                           '',
                           'Phylo PCOA 1', 'Phylo PCOA 2', 'Phylo PCOA 3', 'Phylo PCOA 4')
tr.plotVersion$node.label[tr.plotVersion$node.label == ""] <- NA
tr.plotVersion$node.label[grep('ales', tr.plotVersion$node.label, invert = T)] <- NA
p <- ggtree(tr.plotVersion,
            layout='fan', open.angle = 15,
            aes(color = group))
            #colour = 'gray')
p <- p + geom_label(aes(x=branch),
                    label = c(rep(NA, length(tr.plotVersion$tip.label)),
                              tr.plotVersion$node.label),
                    size = 1.5,
                    fill= 'gray',
                    label.padding = unit(0.18, "lines"),
                    label.r = unit(0.1, "lines"),
                    color = 'black')
p <- gheatmap(p, dat.tr, width = 0.25,
              colnames_position = 'bottom',
              colnames_angle = 270,
              font.size = 1.5, hjust = 0)
p <- p + scale_fill_continuous(name = 'Scaled trait values',
                              low = 'white', high = 'black')
#                                low = 'blue', high = 'red')

p <- p + scale_colour_manual(name = 'Clade-level kd regimes',
                              values = paint.colors,
                             labels = paint.labels)
p <- p + theme(legend.title = element_text(size=7),
               legend.text = element_text(size=5),
              legend.title.align = 0,
              #legend.position = 'right',
              #legend.position = 'none',
              legend.position = c(0.45, 0.45),
              legend.key.size = unit(0.15, 'cm'),
              plot.margin = margin(t = 0.25, r = 0, b = 0, l = 0, unit = "inches")
              )

version <- max(sapply(strsplit(dir('../out', patt = 'decompTree'), '[.v]'), function(x) as.numeric(x[3]) + 1))
pdf(paste('../out/FIG01-decompTree.v', version, format(Sys.time(), ".%Y-%m-%d.pdf"), sep = ''), 7,9)
print(p)
dev.off()
