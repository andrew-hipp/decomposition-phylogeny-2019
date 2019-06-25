## making mop-up figures

# MDS supplement
library(vegan)
library(maps)
library(magrittr)

eco.env.12 <- envfit(eco.mds, dat.2014.bySp[, grep("bio", names(dat.2014.bySp), value = T)],
                  permutations = 4999)
eco.env.34 <- envfit(eco.mds, dat.2014.bySp[, grep("bio", names(dat.2014.bySp), value = T)],
                  permutations = 4999,
                  choices = 3:4)

pdf('../out/SUPPLEMENT.eco.mds.pdf', 8.5, 11)
  layout(matrix(1:6, 3, byrow = TRUE))
  plot(eco.mds, xlim = range(eco.mds$points[, c(1,3)]),
                ylim = range(eco.mds$points[, c(2,4)]),
                main = "NMDS ordination\nBIOCLIM variables")
  plot(eco.env.12, p.max = 0.01)

  plot(eco.mds, xlim = range(eco.mds$points[, c(1,3)]),
                ylim = range(eco.mds$points[, c(2,4)]),
                choices = 3:4,
                main = "NMDS ordination\nBIOCLIM variables4")
  plot(eco.env.34, p.max = 0.01)

  ordisurf(eco.mds ~ bio1, dat.2014.bySp, bubble = 3, choices = 1:2,
           main = "Ordination surface\nannual mean temp (BIO1)")
  ordisurf(eco.mds ~ bio1, dat.2014.bySp, bubble = 3, choices = 3:4,
           main = "Ordination surface\nannual mean temp (BIO1)")
  ordisurf(eco.mds ~ bio14, dat.2014.bySp, bubble = 3, choices = 1:2,
           main = "Ordination surface\nprecip of driest month (BIO14)")
  ordisurf(eco.mds ~ bio14, dat.2014.bySp, bubble = 3, choices = 3:4,
           main = "Ordination surface\nprecip of driest month (BIO14)")
dev.off()

dat.2014.bySp$b1corrected <- dat.2014.bySp$bio1 / 10
pdf('../out/FIG.03.eco.lnKd.pdf')
  plot(eco.mds$points, cex = dat.2014.bySp$kd * 40, xlab = 'NMDS 1', ylab = 'NMDS 2')
  ordisurf(eco.mds ~ b1corrected, dat.2014.bySp, add  = T, col = 'gray')
dev.off()

## figure maps
pdf('../out/FIG.01.map.newCols-v2.pdf', 8.5, 11)
map(col = 'gray60', lwd = 0.15)
a=paste(dat.2014$long, dat.2014$lat, sep="|") %>% table
abins <- list(a=1:2, b=3:5, c=6:10, d=11:20, e=21:50, f=51:max(a))
abins.labels = sapply(abins, function(x) paste(range(x), collapse = '-'))
for(i in names(abins)) a[which(a %in% abins[[i]])] <- i
a[1:length(a)] <- match(a, letters)
a <- sapply(a, as.numeric)
a.col=rainbow(max(a)+1, alpha=0.6, start = 1/6 , end = 4.5/6)
a.rim=rainbow(max(a)+1, alpha = 1, start = 1/6 , end = 4.5/6)
points(sapply(strsplit(names(a), "|", fixed = T), function(x) x[1]),
       sapply(strsplit(names(a), "|", fixed = T), function(x) x[2]),
       pch = 21, cex = 0.75,lwd = 0.01,
       bg = a.col[as.integer(a)],
	     #col = a.rim[as.integer(a)+1]
       col = "gray50"
)

LegLoc <- cbind(x = c(-170, -163, -163, -170),
               y = c(-30, -30, -70 , -70)
			   )
legend.gradient2(pnts=LegLoc, cols = a.col[1:max(a)], title = "Samples per site",
				limits = range(a, na.rm = T),
				title.cex = 0.8,
				labels.cex = 0.5,
        labelSet = abins.labels
				)
dev.off()
