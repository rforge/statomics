#' @name sum.stat.validation
#' @title UK Biobank GWAS summary statistics on height around rs11090631 on chromosome 22
#' @description This data set is the GWAS summary statistics on height computed using UK Biobank Wave 1 data. 
#' 582 SNPs within a 1-Mb centred at rs11090631 are included in the data.
#' @docType data
#' @usage data(sum.stat.validation)
#' @format 1 row per SNP. The variables are as follows:
#' \itemize{
#'    \item SNP, SNP ID
#'    \item A1, effect allele
#'    \item A2, reference allele
#'    \item Freq1, the allele frequency of Allele1
#'    \item b, estimate of marginal effect in GWAS
#'    \item se, standard error of the estimates of marginal effects in GWAS
#'    \item N, sample size
#' }
#' @source Height GWAS results for 120,286 genetically British individuals in UK Biobank Wave 1 data. 
#' The phenotypes were age and sex adjusted and then inverse-gaussian transformed.  
#' @author Zheng Ning, 2018-05-29
NULL