\name{gchc}
\alias{gchc}
\title{
Genetic Correlation Heterogeneity for Causality (GCHC)
}
\description{
This function tests the difference between genetic correlation estimates reported by LD Score regression (LDSC), to infer causality.
}
\usage{
gchc(rg1.log = NULL, rg2.log = NULL, 
     exposure.sumstats = NULL, 
     outcome.sumstats1 = NULL, outcome.sumstats2 = NULL, 
     exposure = "Exposure", outcome = "Outcome")
}
\arguments{
  \item{rg1.log}{
LDSC log file for the 1st genetic correlation estimation between the exposure and outcome phenotypes.
}
  \item{rg2.log}{
LDSC log file for the 2nd genetic correlation estimation between the exposure and outcome phenotypes.
}
  \item{exposure.sumstats}{
The gzipped file of exposure phenotype GWAS summary statistics, munged by LDSC.
}
  \item{outcome.sumstats1}{
The gzipped file of outcome phenotype GWAS (1st population) summary statistics, munged by LDSC.
}
  \item{outcome.sumstats2}{
The gzipped file of outcome phenotype GWAS (2nd population) summary statistics, munged by LDSC.
}
  \item{exposure}{
Name of the exposure phenotype.
}
  \item{outcome}{
Name of the outcome phenotype.
}
}
\details{
GWAS of the exposure phenotype is required to be done in only one population, to keep its heritability fixed. 
Two different populations/sources are required for GWAS of the outcome phenotype.
}
\value{
A list is returned with:
  \item{rg.diff }{The estimated difference between two genetic correlation estimates.}
  \item{se }{The standard error of \code{rg.diff}.}
  \item{p.value }{The p-value testing the null hypothesis of \code{rg.diff} = 0.}
  \item{r.rg }{The estimated correlation between two genetic correlation estimates..}
}
\references{
Shen et al. (2017)
}
\author{
Xia Shen
}
\note{
(c) 2017
}

\examples{
\dontrun{
gchc(rg1.log = 'EDU_BMI.log', 
                rg2.log = 'EDU_st02_off_bmi_res_rn.log', 
                exposure.sumstats = 'EDU.sumstats.gz', 
                outcome.sumstats1 = 'BMI.sumstats.gz',
                outcome.sumstats2 = 'st02_off_bmi_res_rn.sumstats.gz',
                exposure = 'EA', 
                outcome = 'BMI')
}
}