% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/refs_functions.R
\name{plot_refs}
\alias{plot_refs}
\title{Plot reference timelines}
\usage{
plot_refs(datapkg, ...)
}
\arguments{
\item{datapkg}{the name of the data-package to plot references from}

\item{...}{extra arguments passed on to the datapkg's plot function}
}
\description{
Plots a timeline chart for references of the input data-package.
This corresponds to a \code{.plot_refs} function in the data-package.
}
\examples{
\donttest{
plot_refs("wormRef")
}

}
