\name{plot.DBEST}
\alias{plot.DBEST}
\title{
Methods for objects of class "DBEST".
}
\description{
Plot methods for objects of class "DBEST".
}
\usage{
\method{plot}{DBEST}(x, figure = 1, ...)
}
\arguments{
  \item{x}{
DBEST object.
}
  \item{figure}{
If omitted, the plot function displays 'figure 1' and 'figure 2'. If this option is set to 1, the plot function will display the input data and the estimated trend, plus the trend local change (figure 1). If this option is set to 2, the plot function will display a graph with the decomposition of the time-series, including the actual data, the trend, the seasonal component and the remainder.
}
  \item{...}{further arguments passed to the \code{\link[graphics]{plot}} function.}
}
\value{
No return value, called for side effects. The function draws the plots
described under \code{figure} on the current graphics device, and the
graphical parameters of the user are restored when it exits.
}
\author{
Sadegh Jamali, Hristo Tomov
}
\seealso{
See also as \code{\link[DBEST]{DBEST}}
}
\examples{
# Create a "DBEST" object (NDVI data for Site 1, used in Fig. 5, Jamali et al. 2015)
data(NDVI.Site1)
NDVI.Site1 <- ts(NDVI.Site1, start=1982, frequency=12)

DBEST.obj <- DBEST(data=NDVI.Site1, data.type="cyclical",
                   seasonality=12, algorithm="change detection",
                   breakpoints.no=3, first.level.shift=0.1,
                   second.level.shift=0.2, duration=24,
                   distance.threshold="default", alpha=0.05)

## Figure 1: the input data and the estimated trend, plus the trend local change.
plot(DBEST.obj, figure = 1)

## Figure 2: the decomposition of the time-series, including the actual data,
## the trend, the seasonal component and the remainder.
plot(DBEST.obj, figure = 2)
}
\keyword{ts}
