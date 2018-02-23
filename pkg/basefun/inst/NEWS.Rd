
\name{NEWS}
\title{NEWS file for the basefun package}

\section{Changes in Version 0.0-40 (2018-02-23)}{
  \itemize{
    \item{\code{nparm} sometimes generated data with sanitised variable
          names, we don't want this to happen.}
    \item{Allow positive increasing Bernstein polynomials.}
  }
}

\section{Changes in Version 0.0-39 (2017-12-18)}{
  \itemize{
    \item{Add \code{intercept_basis}, ie an intercept only model matrix.}
  }
}


\section{Changes in Version 0.0-38 (2017-05-15)}{
  \itemize{
    \item{\code{model.matrix.cbind_bases} can now deal with 
          lists as \code{newdata} argument (given appropriate 
          dimensions in \code{dim}).}
  }
}



\section{Changes in Version 0.0-37 (2017-04-04)}{
  \itemize{
    \item{Deal with integer variables when \code{sumconstr = TRUE}.}
    \item{Remove unnecessary warning when \code{sumconstr = TRUE} and 
          \code{negative = TRUE}.}
  }
}


\section{Changes in Version 0.0-36 (2017-03-01)}{
  \itemize{
    \item{Register C routines.}
  }
}


\section{Changes in Version 0.0-35 (2017-02-14)}{
  \itemize{
    \item{Fix bug in \code{predict} when \code{newdata} was a data frame.}
  }
}


\section{Changes in Version 0.0-34 (2016-12-14)}{
  \itemize{
    \item{Fix bug in \code{predict} when \code{newdata} was a list. The
          ordering of the resulting array may have been incorrect under
          rare circumstances.}
  }
}

\section{Changes in Version 0.0-33 (2016-09-12)}{
  \itemize{
    \item{Approximate Bernstein polynomial f outside its \code{support} (l, u)
          by a linear function f(u + x) = f(u) + f'(u) * (u - x), now also
          for \code{deriv > 0} (\code{maxderiv} argument in 
          \code{model.matrix.Bernstein_basis}. The default (linear extrapolation
          for \code{deriv = 0} only) did not change.}
  }
}

\section{Changes in Version 0.0-32 (2016-05-31)}{
  \itemize{
    \item{Contraints in \code{log_basis} were potentially wrong.}
  }
}


\section{Changes in Version 0.0-31 (2016-04-25)}{
  \itemize{
    \item{\code{log_basis} has better checks for positivity of
          the variable to be log-transformed.}
  }
}

\section{Changes in Version 0.0-30 (2016-02-19)}{
  \itemize{
    \item{The \pkg{basefun} package was first 
          published on CRAN.}
  }
}
