pkgname <- "p.exact"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "p.exact-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('p.exact')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("p.exact.binary")
### * p.exact.binary

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: p.exact.binary
### Title: Exact p-values for genome-wide association analysis of
###   case-control data in inbred populations
### Aliases: p.exact p.exact.binary
### Keywords: discovery exact false genome-wide p-value, p.exact, rate

### ** Examples

## Not run: 
##D ## loading example gwaa.data of data from Atwell et al. (2010) Nature
##D data(arab)
##D 
##D ## running a regular GWA analysis for AvrRPM1
##D cc1 <- ccfast('X33_.i.avrRpm1..i.', data = arab)
##D 
##D ## check the top finding using the exact p-value
##D top <- which.min(cc1[,'P1df'])
##D or <- cc1[top,'effB']
##D f <- summary(arab[,top])$Q.2
##D maf <- min(f, 1 - f)
##D exact <- p.exact.binary(pheno = 'X33_.i.avrRpm1..i.', gwaa.object = arab, or = or, or.maf = maf, con.table = 'tab.AvrRPM1.RData')
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("p.exact.binary", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("p.exact.gaussian")
### * p.exact.gaussian

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: p.exact.gaussian
### Title: Exact p-values for genome-wide association analysis in inbred
###   populations
### Aliases: p.exact p.exact.gaussian
### Keywords: discovery exact false genome-wide p-value, p.exact, rate

### ** Examples

## Not run: 
##D ## loading example gwaa.data in GenABEL
##D data(ge03d2ex.clean)
##D 
##D ## running a regular GWA analysis for height
##D qt1 <- qtscore(height, data = ge03d2ex.clean)
##D 
##D ## running the multivariate GWAS again
##D exact <- p.exact.gaussian(ge03d2ex.clean, beta = qt1[,'effB'])
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("p.exact.gaussian", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
