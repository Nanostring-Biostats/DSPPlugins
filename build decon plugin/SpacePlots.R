#### functions for drawing points in space and polygons behind them

# contents:
# - spaceplot: a fn to draw expression values in space. Calls the below two functoins
# - getBoundary: fn to infer a polygonal shape around a set of points
# - makeTissueGrid: fn to arrange all the xy points from different tissues in a study in a non-overlapping way



# dev note: what's needed here:
# - 0 and tiny values get an empty pch = 1


#' Plot a gene or score in space
#'
#' Wrapper function for other spatial plotting functions.
#' Draws plots the x-y coords of tissues, with point size giving variable value.
#' @param x X coords
#' @param y Y coords
#' @param z.mat A non-negative matrix of expression levels, pathway/cell scores, etc., with genes as rownames
#' @param point_annot data frame with one column or a clinical variable to be used to color points.
#'      Rownames should match sample names in z.mat and column name should be defined
#' @param tissue A vector of tissue IDs
#' @param tissue.order Optional, vector of tissue names giving the order in which to plot them.
#'  If NULL, then alphabetical order will be used
#' @param tissuecols Named vector of colors to use for each tissue's polygon
#' @param rescale Logical for whether to rescale to give z a mean of 1.
#' @param cex A scalar, controls point size. (Points are further scaled by the values of z.)
#' @param col either a single color, or a named vector of colors assigned to clinical annotation provided in point_annot
#' @return Draws a plot in xy space. Returns a list:
#'  x: new x coords
#'  y: new y coords,
#'  outlines: list of polygonal tissue boundaries' x-y coords
#' @example
#' # sim data from 5 tissues:
#' set.seed(0)
#' x = rnorm(50)
#' y = rnorm(50)
#' z = runif(50)
#' tissue = rep(letters[1:5], each = 10)
#' spaceplot(x, y, z, tissue, point_annot = NULL, tissuecols = NULL, rescale = FALSE,
#'           cex = 1, col = "#00008B80", expansion = 1.2)
spaceplot <- function(x, y, z.mat, point_annot = NULL, tissue, tissue.order = NULL, tissuecols = NULL, tissuecols.alpha = 0.2, rescale = FALSE,
                      cex = 2, col = "#00008B80", expansion = 1.2) {

  # first, get new xy coords:
  newxy.custom <- makeTissueGrid(
    x = x, y = y,
    tissue = tissue, tissue.order = tissue.order,
    expansion = expansion
  )
  x <- newxy.custom$x
  y <- newxy.custom$y

  # also get polygon tissue boundaries:
  boundaries <- list()
  for (tiss in unique(tissue)) {
    tempuse <- tissue == tiss
    boundaries[[tiss]] <- getBoundary(x = x[tempuse], y = y[tempuse], marg = 0.2)
  }

  # transform z:
  z.mat <- z.mat / mean(z.mat, na.rm = TRUE)
  if (rescale) {
    for (i in 1:nrow(z.mat)) {
      z.mat[i, ] <- z.mat[i, ] / mean(z.mat[i, ])
    }
  }

  # tissue boundaries:
  if (length(tissuecols) == 0) {
    tissuecols <- rep("grey50", length(unique(tissue)))
    names(tissuecols) <- unique(tissue)
  }

  # plotting:
  plot_list <- list()

  # plot points
  for (i in 1:nrow(z.mat)) {
    z <- z.mat[i, ]
    gene.name <- rownames(z.mat)[i]
    df <- as.data.frame(cbind(x, y, z))

    # plot points with colors based on point_annot if specified
    if (is.null(point_annot)) {
      df <- as.data.frame(cbind(x, y, z))

      p <- ggplot(data = df, mapping = aes(x, y)) +
        geom_point(aes(color = col), size = z * cex, shape = 16, show.legend = FALSE)
    } else {
      df <- as.data.frame(cbind(x, y, z, point_annot))
      p <- ggplot(data = df, mapping = aes(x, y)) +
        geom_point(aes(color = df[, 4]), size = z * cex, shape = 16, show.legend = TRUE)

      if (length(col) > 1) {
        p <- p + scale_color_manual(name = colnames(df)[4], values = col)
      }
    }

    p <- p + ylab(gene.name)

    # add polygon boundaries
    for (tiss in names(boundaries)) {
      p <- p + geom_polygon(
        data = as.data.frame(boundaries[[tiss]]),
        mapping = aes(x, y), color = NA,
        fill = tissuecols[[tiss]], alpha = tissuecols.alpha
      )

      # get mean of each polygon to use for centering tissue names
      boundary.means <- lapply(boundaries, lapply, mean)
      x.label.points <- as.data.frame(lapply(boundary.means, `[[`, 1))
      colnames(x.label.points) <- names(boundary.means)
    }

    # add tissue labels to the top of the first row of plots
    if (i == 1) {
      p <- p + scale_x_continuous(breaks = x.label.points, labels = colnames(x.label.points), position = "top")
    } else {
      p <- p + theme(axis.text.x = element_blank())
    }

    plot_list[[i]] <- p + theme(
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      panel.background = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_blank()
    )
  }

  # return coordinates:
  out <- list(x = x, y = y, boundaries = boundaries)
  plot_grid <- ggarrange(plotlist = plot_list, ncol = 1, common.legend = TRUE, legend = "right")
  return(plot_grid)
  # return(out)
}






#' Function to define a polygon boundary around a set of points in xy space
#'
#' Fits a convex hull polygon around the points, leaving a bit of margin beyond the points
#' @param x vector of x coords
#' @param y vector of y coords
#' @param marg Amount of extra margin to draw. Default 10%
#' @return a list: x: coords of the polygon. y: y coords on the polygon
#' @example
#'  x = rnorm(30)
#'  y = rnorm(30)
#'  plot(x, y)
#'  bound = getBoundary(x, y)
#'  polygon(bound, col = rgb(0,0,1,0.5))
getBoundary <- function(x, y, marg = 0.2) {
  # get center of shape:
  mx <- mean(x)
  my <- mean(y)
  # expand all points away from center
  x2 <- mx + (x - mx) * (1 + marg)
  y2 <- my + (y - my) * (1 + marg)

  # get convex hull around expanded points:
  ch <- chull(x2, y2)

  out <- list(
    x = x2[ch],
    y = y2[ch]
  )
  return(out)
}


#' Define non-overlapping x-y coords for plotting multiple tissues in the same frame
#'
#' Given xy coords for segments from several tissues, shifts each tissue so they don't overlap
#' @param x Vector of x coords
#' @param y Vector of y coords
#' @param tissue Vector of tissue IDs corresponding to x and y
#' @param tissue.order Optional, vector of tissue names giving the order in which to plot them.
#'  If NULL, then alphabetical order will be used
#' @param expansion A constant scaling factor, for how much to expand the margins between tissues
#' @return A list: x = new x coords. y = new y coords.
#' @examples
#' # sim data from 5 tissues:
#' set.seed(0)
#' x <- rnorm(50)
#' y <- rnorm(50)
#' tissue <- rep(letters[1:5], each = 10)
#'
#' # arrange in default layout:
#' newxy <- makeTissueGrid(
#'   x = x, y = y, tissue = tissue,
#'   tissue.order = NULL, expansion = 1.2
#' )
#' plot(newxy, pch = 16, col = 0)
#' text(newxy$x, newxy$y, tissue, col = as.numeric(as.factor(tissue)))
#'
#' # specify a layout:
#' newxy.custom <- makeTissueGrid(x = x, y = y, tissue = tissue, expansion = 1.2)
#' plot(newxy.custom, pch = 16, col = 0)
#' text(newxy.custom$x, newxy.custom$y, tissue, col = as.numeric(as.factor(tissue)))
makeTissueGrid <- function(x, y, tissue,
                           tissue.order = NULL, expansion = 1.2) {

  # define the tissue ids and their order:
  if (length(tissue.order) > 0) {
    tissues <- tissue.order
  }
  else {
    tissues <- sort(unique(tissue))
  }

  # get tissue spans and centers:
  xspans <- yspans <- xcenters <- ycenters <- c()
  for (tiss in tissue) {
    rx <- range(x[tissue == tiss])
    ry <- range(y[tissue == tiss])
    xspans[tiss] <- diff(rx)
    yspans[tiss] <- diff(ry)
    xcenters[tiss] <- median(rx)
    ycenters[tiss] <- median(ry)
  }

  # define how far apart tissues should be in x and y space:
  xmarg <- max(xspans) * expansion
  ymarg <- max(yspans) * expansion

  # now get xy offsets for each tissue
  xnew <- ynew <- replace(x, TRUE, NA)
  for (tiss in tissues) {
    tempxoffset <- which(tissues == tiss) * xmarg - xcenters[tiss]
    tempyoffset <- ymarg - ycenters[tiss] # tissues will always be in one row
    xnew[tissue == tiss] <- x[tissue == tiss] + tempxoffset
    ynew[tissue == tiss] <- y[tissue == tiss] + tempyoffset
  }

  out <- list(x = xnew, y = ynew)
  return(out)
}
