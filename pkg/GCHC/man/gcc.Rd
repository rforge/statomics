% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gchc.R
\name{gcc}
\alias{gcc}
\title{Genetic Correlation heterogeneity test for Causality (GCC)}
\usage{
gcc(rg1.log = NULL, rg2.log = NULL, exposure.sumstats = NULL,
  outcome.sumstats1 = NULL, outcome.sumstats2 = NULL,
  exposure = "Exposure", outcome = "Outcome", heading = TRUE)
}
\arguments{
\item{rg1.log}{LDSC log file for the 1st genetic correlation estimation between the exposure and outcome phenotypes.}

\item{rg2.log}{LDSC log file for the 2nd genetic correlation estimation between the exposure and outcome phenotypes.}

\item{exposure.sumstats}{The gzipped file of exposure phenotype GWAS summary statistics, munged by LDSC.}

\item{outcome.sumstats1}{The gzipped file of outcome phenotype GWAS (1st population) summary statistics, munged by LDSC.}

\item{outcome.sumstats2}{The gzipped file of outcome phenotype GWAS (2nd population) summary statistics, munged by LDSC.}

\item{exposure}{Name of the exposure phenotype.}

\item{outcome}{Name of the outcome phenotype.}

\item{heading}{logical value that specified whether the software intro heading is printed.}
}
\value{
A list is returned with:
\itemize{
\item{rg.diff }{The estimated difference between two genetic correlation estimates.}
\item{se }{The standard error of \code{rg.diff}.}
\item{p.value }{The p-value testing the null hypothesis of \code{rg.diff} = 0.}
\item{r.rg }{The estimated correlation between two genetic correlation estimates.}
}
}
\description{
This function tests the difference between genetic correlation estimates 
reported by LD Score regression (LDSC), to infer causality.
}
\note{
GWAS of the exposure phenotype is required to be done in only one population, to keep its heritability fixed. 
Two different populations/sources are required for GWAS of the outcome phenotype.
}
\examples{
\dontrun{
gchc(rg1.log = 'EDU_BMI1.log', 
                rg2.log = 'EDU_BMI2.log', 
                exposure.sumstats = 'EDU.sumstats.gz', 
                outcome.sumstats1 = 'BMI1.sumstats.gz',
                outcome.sumstats2 = 'BMI2.sumstats.gz',
                exposure = 'EA', 
                outcome = 'BMI')
}
}
\references{
Shen X, Ning Z, Joshi PK, Lee Y, Wilson JF, Pawitan Y (2017). Heterogeneity in genetic correlation detects causal factors
for complex traits. \emph{Submitted}.
}
\seealso{
GCC homepage: http://gcc.shen.se
}
\author{
Xia Shen
}
\keyword{causal}
\keyword{correlation}
\keyword{genetic}
\keyword{inference,}
