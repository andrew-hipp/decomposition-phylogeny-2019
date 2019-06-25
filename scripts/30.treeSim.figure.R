library(geiger)
library(ggtree)
library(gridExtra)
library(picante)

tr.sim <- sim.bdtree(n = 50)
a = data.frame(Brownian = sim.char(tr.sim, par = 1, nsim = 1, model = 'BM')[,,1],
          Random = sim.char(lambdaTree(tr.sim, 0.1), par = 1, nsim = 1, model = 'BM')[,,1])
a.fit <- fitContinuous(tr.sim, a, model = 'lambda')

tr <- ggtree(tr.sim, layout = 'circular')
tr1 <- gheatmap(tr, a[,'Brownian', drop = FALSE],
                width = 0.1, colnames = F)
tr1 <- tr1 + scale_fill_continuous(low = 'white', high = 'black')
tr1 <- tr1 + ggtitle('High phylogenetic signal',
                      subtitle = paste("Pagel's lambda = ", round(a.fit$Brownian$opt$lambda, 4),
                      ",\nBlomberg's K = ", Kcalc(x=a[tr.sim$tip.label, 'Brownian'], phy=tr.sim) %>% as.numeric %>% round(4),
                    sep = '')
                      )
tr1 <- tr1 + theme(legend.position = 'none')
  #                  plot.title = element_text(hjust = 0.5))

tr2 <- gheatmap(tr, a[,'Random', drop = FALSE],
                width = 0.1, colnames = F)
tr2 <- tr2 + scale_fill_continuous(low = 'white', high = 'black')
tr2 <- tr2 + labs(title = 'Low phylogenetic signal',
                  subtitle = paste("Pagel's lambda = ", round(a.fit$Random$opt$lambda, 4),
                  ",\nBlomberg's K = ", Kcalc(x=a[tr.sim$tip.label, 'Random'], phy=tr.sim) %>% as.numeric %>% round(4),
                sep = '')
                  )

tr2 <- tr2 + theme(legend.position = 'none')

pdf('../out/Fig.samplePhyloSignal-v2.pdf')
grid.arrange(tr1, tr2, nrow = 1)
dev.off()
