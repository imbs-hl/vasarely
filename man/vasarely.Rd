% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vasarely.R
\name{vasarely}
\alias{vasarely}
\title{vasarely chart}
\usage{
vasarely(data, color = grey.colors(256, start = 1, end = 0),
  name_xaxis = "allele 2", name_yaxis = "allele 1", lower_color_value = 0,
  upper_color_value = 1)
}
\arguments{
\item{data}{The input data to create the chart,
e.g. genetic snp data. The input data must have
exactly to columns.}

\item{color}{The optional colors which should be used
for the chart, e.g. blues9. Color must be a character
vector with at least two colors. If color is NULL grey
values will be taken.}

\item{name_xaxis}{The optional title for the x-axis of
the plot. Otherwise the x-axis is called "allele 2"
as a genetic input is expected.}

\item{name_yaxis}{The optional title for the y-axis of
the plot. Otherwise the y-axis is called "allele 1"
as a genetic input is expected.}

\item{lower_color_value}{The lower limit for
spreading the color over the relative frequencies values.
It must be a number between 0 and 1.}

\item{upper_color_value}{The upper limit for
spreading the color over the relative frequencies values.
It must be a number between 0 and 1.}
}
\value{
returns a list of the calculated relative
frequencies, the number of observations, the vasarely
chart and a list of the p-value and statistic of
a chi-squared test.
}
\description{
vasarely produces a so called vasarely chart.
E.g. for some (genetic) input data expected and observed
relative frequencies of allele combinations can be calculated
and plotted. For all possible combinations
the expected relative frequency is described by a geom_raster
plot. Each cell of the plot has a color which depends
on the value of the expected relative frequcency. The real
relative frequency is described by a dot in the corresponding
cell which also has a color that depends on the
calculated relative frequency. The chart and its colors for
the two relative frequencies help to check how expected
and observed relative frequencies differ from each other.
}
\examples{
# create data
a1 <- c(rep("A", each = 25), rep("B", each = 75))
a2 <- c(rep("A", each = 50), rep("B", each = 50))
data <- data.frame(a1, a2)

# use function
vasarely(data = data, color = c("yellow", "red"),
name_xaxis = "a1", name_yaxis = "a2")

}
