\name{drawPedigree}
\alias{drawPedigree}
\title{Produce a graphical representation of a pedigree}
\description{
  Plots a pedigree, with options specific to considerations for pedigrees used to for quantitative genetic inferencein natural populations.  Pedigrees containing only those individuals that are informative with respect to (genetic) variation in an arbitrary trait can be plotted, potentially overlain on a complete pedigree.  Functions also exist to plot various types of pedgiree links associated with focal individuals.
}
\usage{
drawPedigree(Ped, cohorts = NULL, sex = NULL, dat = NULL, dots = "n", 
              plotfull = "y", writeCohortLabels = "n", links = "all", 
              sexInd = c(0, 1), dotSize = 0.001, dataDots = "n", 
              dataDots.cex = 2, cohortLabs.cex = 1, retain="informative", 
              focal=NULL, sexColours=c('red','blue'), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Ped}{An ordered pedigree with 3 columns: id, dam, sire}
  \item{cohorts}{An optional numeric vector of the same length as the pedigree designating, for example cohort afinities or birth years}
  \item{sex}{An optional numeric vector of the same length as the pedigree containing the sexes (may be unknown) of all indivduals with entries in the pedigree.  Defaults (modifiable with \code{sexInd}) are 0=male and 1=female}
  \item{dat}{An optional vector or data frame containg indicators of data availablility.  If dat contains only ones and zeros, then any individual with any entry of one will be considered as having data records.  If data contains values other than ones and zeros, individuals in the pedigree with rows in data that contain at least one available record, i.e., one data record is not NA, will be treated as having data.}
  \item{dots}{If 'y', then a dot will be printed representing each individual in the pedigree.  If sexes are available, dots will be colour coded by sex.}
  \item{plotfull}{To be used when dat is supplied.  If 'y' (the default), individuals in the pedigree that are uninformative with respect to the available data have their pedigree links plotted in gray.}
  \item{writeCohortLabels}{To be used when cohorts is used.  Will plot the cohort values on the left hand side of the pedigree image.}
  \item{links}{Default is 'all', other values are 'mums' to print only maternal pedigree links and 'dads' to print only paternal pedigree links.}
  \item{sexInd}{To be used with if sex is supplied and if the vector of sex specifiers differs from the default.}
  \item{dotSize}{Set the dot size bigger or smaller}
  \item{dataDots}{Will print dots over the dots denoting indivdiuals, but denoting individuals with available data as indicated by dat.}
  \item{dataDots.cex}{controls the size of dataDots relative to dots.}
  \item{cohortLabs.cex}{controls the size of cohort lables.}
  \item{retain}{When those pedigree links only informative relative to phenotypic data availability are to be plotted, this controls whether or not a pruned pedigree based on phenotypic data is plotted (if set to "pruned"), or whether strictly only those informative pedigree links are plotted (the default)}
  \item{focal}{An optional list containing the id of an individual and the kinds of relatives of the focal individual to which to plot pedigree links.  Available types are 'offspring','descendents','parents',,ancestors', and 'kin'.}
  \item{sexColours}{The colours that will be used to draw points and or lines associated with males and females.}
  \item{...}{Additional graphical parameters.}
}


\references{Morrissey, M.B, and A.J. Wilson, 2009. pedantics, an R package for pedigree-based genetic simulation, and pedigree manipulation, characterisation, and viewing. Molecular Ecology Resources.}
\author{Michael Morrissey \email{michael.morrissey@st-andrews.ac.uk}}

\seealso{ \code{\link{fixPedigree}} to prepare pedigrees that may not explicitly contain records for all individuals (i.e., where founding individuals may only appear in the dam or sire column).)}
\examples{
data(gryphons)
pedigree<-fixPedigree(gryphons[,1:3])

## draw the gryphon pedigree by pedigree depth
drawPedigree(pedigree)

## draw the gryphon pedigree by cohort
# drawPedigree(pedigree,cohorts=gryphons$cohort,writeCohortLabels='y',
#                                                   cohortLabs.cex=1)

\dontrun{

## draw the gryphon pedigree by cohort with only maternal links
drawPedigree(pedigree,cohorts=gryphons$cohort,links='mums')

## draw the gryphon pedigree by cohort with colour only for those
## indiduals that are informative relative to the quantitative
## genetics of a hypothetical trait only measured for individuals
## in the last two cohorts, emphasize the phenotyped individuals
## with large black dots, and all other individuals with dots 
## colour coded by sex:
dataAvailability<-(gryphons$cohort>=(max(gryphons$cohort)-1))+0

# not run
# drawPedigree(pedigree,cohorts=gryphons$cohort,sex=gryphons$sex,
#   dots='y',dat=dataAvailability,writeCohortLabels='y',dataDots='y')

}

}
\keyword{hplot}
