
\name{mbo_recommended}
\alias{mbo_recommended}
\title{Fit recommended regularized tram based on model based optimization output}
\usage{
mbo_recommended(mbo_obj, m0, x, ...)
}
\arguments{
\item{mbo_obj}{object return by \code{mbo_tramnet}}

\item{m0}{null model of class \code{"tram"}}

\item{x}{matrix of covariables}

\item{...}{additional arguments to \code{tramnet()}}
}
\description{
Extracts the "optimal" tuning parameters from an object of class
\code{"MBOSingleObjResult"} and fits the corresponding \code{tramnet} model
}
\value{Object of class \code{"tramnet"}}
