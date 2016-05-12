\name{DBEST}
\alias{DBEST}
\encoding{UTF-8}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Detecting Breakpoints and Estimating Segments in Trend (DBEST)
}
\description{
An algorithm for analysing vegetation time series with two main application domains:
generalising vegetation trends and characterizing vegetation trend changes. The algorithm
simplifies the trend into linear segments using one of three parameters: a generalisation threshold
parameter, the largest changes, or a threshold for the magnitude of changes. The outputs are the
simplified trend, the change type (abrupt or non-abrupt), and the time and magnitude of the change.
}
\usage{
DBEST(data, data.type, seasonality = -1, algorithm, breakpoints.no = -1, 
generalization.percent = -1, change.magnitude = -1, first.level.shift, 
second.level.shift, duration, distance.threshold, alpha, plot = -1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
univariate time-series to be analysed. This should be an object of class 'ts'/'zoo' with a
frequency greater than one without NA's or a vector without NA's. If the input data is a vector the
algorithm will automatically assign a start year.
}
  \item{data.type}{
the seasonal component. There are two options: "cyclical" or "non-cyclical".
}
  \item{seasonality}{
the seasonality period as a number. If the input data type is non-cyclical this variable should
be omitted or set to 'none'/'null'. This parameter will overwrite the frequency value of an input
object of class 'ts'/'zoo'. However, if the input data is an object of class 'ts'/'zoo' and the seasonality
period is omitted, then the algorithm will use as a seasonality period the frequency value from the
input data.
}
  \item{algorithm}{
the algorithm mode. There are two options: "change detection" and "generalisation".
}
  \item{breakpoints.no}{
the number of changes or breakpoints of interest for detection. This parameter should be
omitted if 'generalization.percent' or 'change.magnitude' is in use.
}
  \item{generalization.percent}{
the generalization threshold value in percent as a number. This parameter should be omitted
if 'breakpoints.no' or 'change.magnitude' is in use.
}
  \item{change.magnitude}{
the change magnitude threshold value as a number. This parameter should be omitted if
'breakpoints.no' or 'change.magnitude' is in use.
}
  \item{first.level.shift}{
the first level-shift-threshold value. This parameter corresponds to the lowest absolute
difference in the time-series data between the level-shift point and the next data point.
}
  \item{second.level.shift}{
the second level-shift-threshold value. This parameter corresponds to the lowest absolute
difference in the means of the data calculated over the duration period before and after the level-
shift point.
}
  \item{duration}{
the duration threshold value. This parameter corresponds to the the lowest time period (time
steps) within which the shift in the mean of the data level, before and after the level-shift
point, persists; and, the lowest spacing (time steps) between successive level-shift points.
}
  \item{distance.threshold}{
the distance threshold value. This correspond to the the lowest perpendicular distance from
farthest data point to the straight line passing through every pair of successive peak and valley
points. The algorithm will estimate a distance threshold if this parameter is set to 'default'.
}
  \item{alpha}{
the statistical significance level value used for testing the significance of detected changes.
}
  \item{plot}{
display figures. This parameter could be omitted or set to: "on", "fig1", "fig2" and "off".
The 'on' option will display 'figure 1' and 'figure 2'. The "fig1" option will display the input data and
the estimated trend, plus the trend local change. The "fig2" option will display a graph with the
decomposition of the time-series, including the actual data, the trend, the seasonal component and
the reminder.
}
}
\details{
An object of the class "DBEST" is a list with elements depending on whether the algorithm
is used for change detection or generalisation.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
\item{BreakpointNo}{the number of break points estimated by the algorithm.}
\item{SegmentNo}{the number of segments estimated by the algorithm.}
\item{Start}{a list with numbers representing the starting points of the changes as time-steps.}
\item{Duration}{a list with numbers representing the durations of the changes as time-steps.}
\item{End}{a list with numbers representing the ending points of the changes as time-steps.}
\item{Change}{a list with the values of the changes.}
\item{ChangeType}{a list with the types of the changes as numbers which could be 0 or 1. The number 0
correspond to a non-abrupt change and 1 to an abrupt change.}
\item{Significance}{a list with the statistical significances of the changes as numbers which could be 0 or
1. The number 0 correspond to a statistically insignificant change and 1 to a significant
change.}
\item{RMSE}{the calculated Root Mean Squared Error of the fit.}
\item{Change}{the calculated Mean Absolute Deviation of the fit.}

}
\references{
Jamali, S, Eklundh, L, Ardö, J, Seaquist, J, & Jönsson, P 2015, 'Detecting changes in vegetation trends using
time series segmentation', Remote Sensing Of Environment, 156, p. 182-195, Scopus®, EBSCOhost,
viewed 27 March 2015.

Tomov, H. 2016, 'Automated temporal NDVI analysis over the Middle East for the period 1982 - 2010', http://lup.lub.lu.se/student-papers/record/8871893

}
\author{
Sadegh Jamali, Hristo Tomov
}

\seealso{
\code{\link[DBEST]{plot.DBEST}} for plotting of DBEST() results.
}

\examples{

\dontrun{
data(NDVI.Example1)
NDVI.Example1 <- ts(NDVI.Example1, start=1982, frequency=12)

data(NDVI.Example2)
NDVI.Example2 <- ts(NDVI.Example2, start=1982, frequency=12)

# Examples of a trend generalisation (Site 1 & Site 2, Jamali et al. 2015)
#(a, f) no-breakpoint
DBEST.a <- DBEST(data=NDVI.Example1, data.type="cyclical", 
                   seasonality=12, algorithm="generalisation", 
                   breakpoints.no=0, first.level.shift=0.1, 
                   second.level.shift=0.2, duration=24, 
                   distance.threshold="default", alpha=0.05, plot="fig1")

DBEST.f <- DBEST(data=NDVI.Example2, data.type="cyclical", 
                   seasonality=12, algorithm="generalisation", 
                   breakpoints.no=0, first.level.shift=0.1, 
                   second.level.shift=0.2, duration=24, 
                   distance.threshold="default", alpha=0.05, plot="fig1")

#(b, g) major breakpoint
DBEST.b <- DBEST(data=NDVI.Example1, data.type="cyclical", 
                 seasonality=12, algorithm="generalisation", 
                 breakpoints.no=1, first.level.shift=0.1, 
                 second.level.shift=0.2, duration=24, 
                 distance.threshold="default", alpha=0.05, plot="fig1")

DBEST.g <- DBEST(data=NDVI.Example2, data.type="cyclical", 
                 seasonality=12, algorithm="generalisation", 
                 breakpoints.no=1, first.level.shift=0.1, 
                 second.level.shift=0.2, duration=24, 
                 distance.threshold="default", alpha=0.05, plot="fig1")

#(c, h) change magnitude 0.1
DBEST.c <- DBEST(data=NDVI.Example1, data.type="cyclical", 
                 seasonality=12, algorithm="generalisation", 
                 change.magnitude=0.1, first.level.shift=0.1, 
                 second.level.shift=0.2, duration=24, 
                 distance.threshold="default", alpha=0.05, plot="fig1")

DBEST.h <- DBEST(data=NDVI.Example2, data.type="cyclical", 
                 seasonality=12, algorithm="generalisation", 
                 change.magnitude=0.1, first.level.shift=0.1, 
                 second.level.shift=0.2, duration=24, 
                 distance.threshold="default", alpha=0.05, plot="fig1")

#(d, i) all-inclusive breakpoints
DBEST.d <- DBEST(data=NDVI.Example1, data.type="cyclical", 
                 seasonality=12, algorithm="generalisation", 
                 generalization.percent=0, first.level.shift=0.1, 
                 second.level.shift=0.2, duration=24, 
                 distance.threshold="default", alpha=0.05, plot="fig1")

DBEST.i <- DBEST(data=NDVI.Example2, data.type="cyclical", 
                 seasonality=12, algorithm="generalisation", 
                 generalization.percent=0, first.level.shift=0.1, 
                 second.level.shift=0.2, duration=24, 
                 distance.threshold="default", alpha=0.05, plot="fig1")

# Examples of a change detection (Site 1 & Site 2, Jamali et al. 2015)

#(a, b) detection the three major changes
DBEST.a <- DBEST(data=NDVI.Example1, data.type="cyclical", 
                 seasonality=12, algorithm="change detection", 
                 breakpoints.no=3, first.level.shift=0.1, 
                 second.level.shift=0.2, duration=24, 
                 distance.threshold="default", alpha=0.05, plot="on")

DBEST.b <- DBEST(data=NDVI.Example2, data.type="cyclical", 
                 seasonality=12, algorithm="change detection", 
                 breakpoints.no=3, first.level.shift=0.1, 
                 second.level.shift=0.2, duration=24, 
                 distance.threshold="default", alpha=0.05, plot="on")

}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ts}% use one of  RShowDoc("KEYWORDS")
