
\name{NEWS}
\title{NEWS file for the tram package}

\section{Changes in Version 0.6-1 (2021-03-10)}{
  \itemize{
    \item{Allow diagonals of the inverse Cholesky factor to be modelled
          as a function of covariates in \code{mmlt}.}
  }
}


\section{Changes in Version 0.6-0 (2021-03-08)}{
  \itemize{
    \item{Add \code{simulate} method.}
    \item{Enable nonparametric likelihood estimation by \code{Polr}.}
    \item{New argument \code{sparse_nlevels} to \code{tram}, requesting
          a model \code{Matrix} for ordered factors with at least
          \code{sparse_nlevels} levels.}
  }
}


\section{Changes in Version 0.5-2 (2020-10-23)}{
  \itemize{
    \item{Add experimental \code{score_test} and \code{perm_test} for
          \code{survival::coxph}.}
    \item{Update simulation code.}
  }
}


\section{Changes in Version 0.5-1 (2020-09-03)}{
  \itemize{
    \item{Robustify computation of score intervals.}
    \item{Add Spearman correlation to \code{mmlt} example.}
  }
}

\section{Changes in Version 0.5-0 (2020-08-27)}{
  \itemize{
    \item{\code{Lehmann()} parameter is now labelled as proportional reverse time hazard ratio.}
    \item{Add \code{mmlt} for multivariate conditional transformation models
          (continuous responses only for now).}
    \item{Reproducibility material for <arxiv:1906.03151> added to package.}
  }
}


\section{Changes in Version 0.4-0 (2020-05-11)}{
  \itemize{
    \item{Force support to start at time zero and constrain cumulative
          hazard to zero for time zero in \code{Aareg}.}
  }
}


\section{Changes in Version 0.4-0 (2020-05-11)}{
  \itemize{
    \item{Add \code{Aareg} for Aalen additive hazard models.}
  }
}


\section{Changes in Version 0.3-3 (2020-01-16)}{
  \itemize{
    \item{\code{score_test} generic and method added for computing 
          score tests and corresponding confidence intervals.}
    \item{\code{perm_test} generic and method added for computing 
          permutation score tests.}
    \item{\code{PI} converts log-odds ratios to probabilistic indices and 
          vice versa.}
  }
}


\section{Changes in Version 0.3-2 (2020-01-08)}{
  \itemize{
    \item{\code{profile} method added for computing 
          profile likelihood confidence intervals.}
    \item{Allow additional linear inequality constraints for regression
          coefficients.}
  }
}

\section{Changes in Version 0.3-1 (2019-11-05)}{
  \subsection{Bugfixes}{
    \itemize{
      \item{Fix problem in C header file.}
      \item{Fix problems with responses of class \code{R}, spotted by Balint Tamasi.}
    }
  }
}

\section{Changes in Version 0.3-2 (2019-12-02)}{
  \itemize{
    \item{Add new demo.}
  }
  \subsection{Bugfixes}{
    \itemize{
      \item{Deal with \code{model_only} in \code{Survreg}, spotted by Balint Tamasi.}
    }
  }
}


\section{Changes in Version 0.3-0 (2019-10-09)}{
  \subsection{Bugfixes}{
    \itemize{
      \item{\code{y | s1:s2 ~ x} lead to rank-deficit models, the intercept is
            now removed.}
      \item{Not really a bug, but time-consuming examples were removed
            (or replaced by static output) in package vignettes to make
            CRAN happy.}
      \item{\code{as.mlt} returned an object of class \code{mlt} but not
            \code{mlt_fit}.}
    }
  }
  \subsection{New Features}{
    \itemize{
      \item{New marginally interpretable models for clustered data
            in \code{mtram}, see \code{mtram} vignette. This is still
            experimental and very rough, so please be careful.}
      \item{Compute \code{vcov.tram} via Schur complement when possible.}
  }
  }
}

\section{Changes in Version 0.2-6 (2019-06-20)}{
  \itemize{
    \item{Add regression test for count data and \code{log_first = TRUE}.}
  }
}


\section{Changes in Version 0.2-5 (2019-02-25)}{
  \itemize{
    \item{\code{tram} objects are now in line with \code{multcomp::glht},
          allowing post-hoc multiple comparisons.}
  }
}

\section{Changes in Version 0.2-4 (2019-01-22)}{
  \itemize{
    \item{\code{log_first = TRUE} is now accepted by \code{tram} and all 
          model functions calling \code{tram}. This is still experimental.}
    \item{Silence regression tests on Windows.}
  }
}


\section{Changes in Version 0.2-3 (2018-09-17)}{
  \itemize{
    \item{New \code{extrapolate = TRUE} option for all models constraints
          the second derivate of the Bernstein transformation function to be
          zero at \code{support}. This helps to get rid of possible
          discontinuous densities at \code{support}.}
    \item{Add \code{Lehmann} for fitting shift transformation models those
          linear predictor can be interpreted as log-Lehmann-alternative.}
    \item{Add \code{method = "loglog"} to \code{Polr}.}
    \item{New \code{update} method.}
  }
}

\section{Changes in Version 0.2-2 (2018-06-22)}{
  \itemize{
    \item{Use default \code{bounds = c(0, Inf)} for \code{Surv} objects (for example
          in \code{Coxph}.}
    \item{Update to \pkg{mlt} 1.0-2.}
  }
}


\section{Changes in Version 0.2-1 (2018-04-20)}{
  \itemize{
    \item{Likelihood ratio tests with fixed parameters gave an error.}
    \item{Remove strata from linear predictor (experimental).}
    \item{Try to implement a better guess for the \code{add} argument
          later defining the axes corresponding to the response
          in plots.}
    \item{All models accept \code{add} and \code{bounds} arguments.}
  }
}

\section{Changes in Version 0.2-0 (2018-02-25)}{
  \itemize{
    \item{\pkg{tram} first published on CRAN.}
  }
}
