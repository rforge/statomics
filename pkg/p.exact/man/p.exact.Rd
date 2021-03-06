% Generated by roxygen2 (4.0.2): do not edit by hand
\docType{package}
\name{p.exact}
\alias{p.exact}
\alias{p.exact-package}
\title{Exact p-values}
\description{
p.exact: Exact P-values for Genome-Wide Association Analyses in Inbred Populations
}
\details{
The p.exact package implements for calculating exact dataset-specific p-values for
genome-wide association (GWA) analyses in inbred populations. The module is
compatible with existing GenABEL/gwaa.data data formats.

For converting data from other formats, see

\code{\link{convert.snp.illumina}} (Illumina/Affymetrix-like format).
\code{\link{convert.snp.text}} (conversion from human-readable GenABEL format),
\code{\link{convert.snp.ped}} (Linkage, Merlin, Mach, and similar files),
\code{\link{convert.snp.mach}} (Mach-format),
\code{\link{convert.snp.tped}} (from PLINK TPED format),
\code{\link{convert.snp.affymetrix}} (BRML-style files).

For converting of GenABEL's data to other formats, see
\code{\link{export.merlin}} (MERLIN and MACH formats),
\code{\link{export.impute}} (IMPUTE, SNPTEST and CHIAMO formats),
\code{\link{export.plink}} (PLINK format, also exports phenotypic data).

To load the data, see \code{\link{load.gwaa.data}}.

For data managment and manipulations see
\code{\link{merge.gwaa.data}},
\code{\link{merge.snp.data}},
\code{\link{gwaa.data-class}},
\code{\link{snp.data-class}},
\code{\link{snp.names}},
\code{\link{snp.subset}}.
}
\author{
Xia Shen
}
\references{
If you use p.exact package in your analysis, please cite the following work:

Xia Shen (2015). Flaw or discovery? Calculating exact p-values for
genome-wide association studies in inbred populations. \emph{Submitted}.
}
\seealso{
\code{p.exact.gaussian}, \code{p.exact.binary}
}
\keyword{package}

