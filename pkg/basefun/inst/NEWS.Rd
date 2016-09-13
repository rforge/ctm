
\name{NEWS}
\title{NEWS file for the basefun package}

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
