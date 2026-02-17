library(dplyr)
library(plotly)

# ---- Function: draw_landscape_plotly ----
# Interactive 3D surface using plotly
draw_landscape_plotly <- function(x, y, f, zlim = NULL, showscale = TRUE) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' required for interactive plot. install.packages('plotly')")
  }
  stopifnot(is.numeric(x), is.numeric(y))
  nx <- length(x); ny <- length(y)
  if (is.function(f)) {
    vf <- Vectorize(f, vectorize.args = c("x", "y"))
    Z <- outer(x, y, vf)
  } else if (is.matrix(f)) {
    if (!all(dim(f) == c(nx, ny))) stop("Matrix f must have dimensions length(x) x length(y)")
    Z <- f
  } else stop("f must be either a function or a matrix")
  
  if (is.null(zlim)) zlim <- range(Z, na.rm = TRUE)
  plotly::plot_ly(x = ~x, y = ~y, z = ~Z, type = "surface", showscale = showscale, alpha = .75) %>%
    plotly::layout(scene = list(xaxis = list(title = "x"),
                                yaxis = list(title = "y"),
                                zaxis = list(title = "z", range = zlim)))
}



# example vectors
x <- y <- seq(-2, 2, length.out = 1e2)

# example function (sinc-like)
f <- function(x, y) {
  10 * cos(pi * x) - 7 * y ^ 2
}

draw_landscape_plotly(x, y, f)