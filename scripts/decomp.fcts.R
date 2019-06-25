tidyName <- function(x) tolower(gsub("['_ .-]", "", x))

nameAssign <- function(ecoNames = as.character(putNameVectorHere),
                       oldNames = as.character(spp.matrix$Name_submitted),
                       newNames = as.character(spp.matrix$Accepted_name)) {
## assigns the updated name for the old names in ecoMat, puts them in a new column entitled "tree label"
  matches <- newNames[match(tidyName(ecoNames), tidyName(oldNames))]
  return(matches)
  }

aic.w <- function(x) {
    e.aic <- exp(-0.5 * (x - min(x)))
    aic.w <- e.aic / sum(e.aic)
    return(aic.w)
  }

legend.gradient2 <- function (pnts, cols = heat.colors(100), limits = c(0, 1), title = "Legend",
    title.cex = 1, labels.cex = 0.5, labels.inset = 0.10, labelSet = NA,
    ...)
	## modified from SDMTools, ah 2018-09-12
{
    pnts = try(as.matrix(pnts), silent = T)
    if (!is.matrix(pnts))
        stop("you must have a 4x2 matrix")
    if (dim(pnts)[1] != 4 || dim(pnts)[2] != 2)
        stop("Matrix must have dimensions of 4 rows and 2 columms")
    if (length(cols) < 2)
        stop("You must have 2 or more colors")
    yvals = seq(min(pnts[, 2]), max(pnts[, 2]), length = length(cols) +
        1)
    for (i in 1:length(cols)) {
        polygon(x = pnts[, 1], y = c(yvals[i], yvals[i], yvals[i +
            1], yvals[i + 1]), col = cols[i], border = F)
        if(!is.na(labelSet[1])) text(max(pnts[, 1]),
                                     mean(yvals[i:(i+1)]),
                                     labelSet[i], cex = labels.cex,
                                     pos = 4)
    }
    if(is.na(labelSet[1])){
      text(max(pnts[, 1]), min(pnts[, 2]) + labels.inset * (max(pnts[, 2]) - min(pnts[, 2])),
	     labels = limits[1],
        pos = 4, cex=labels.cex, ...)
      text(max(pnts[, 1]), max(pnts[, 2]) - labels.inset * (max(pnts[, 2]) - min(pnts[, 2])),
	     labels = limits[2] ,
        pos = 4, cex=labels.cex, ...)
      }
    text(min(pnts[, 1]), max(pnts[, 2]), labels = title, adj = c(0,
        -1), cex=title.cex, ...)
}
