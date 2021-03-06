\name{phensim}
\alias{phensim}
\title{A function to simulated phenotypic data}
\description{
 Simulates phenotypic data across arbitrary pedigrees.  \code{phensim} simulate direct, maternal and paternal genetica and environmental effects for an arbitrary number of traits with arbitrary patterns of missing data.
}
\usage{
phensim(pedigree, traits = 1, randomA = NULL, randomE = NULL, 
        parentalA = NULL, parentalE = NULL, sampled = NULL, 
              records = NULL, returnAllEffects = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pedigree}{A pedigree}
  \item{traits}{The number of traits for which data should be simulated.}
  \item{randomA}{An additive genetic covariance matrix, with dimensions a multiple of traits - see details}
  \item{randomE}{An additive environmental covariance matrix, with dimensions a multiple of traits - see details}
  \item{parentalA}{A vector indicating which effects in \code{randomA} (if any) to treat as parental effects}
  \item{parentalE}{A vector indicating which effects in \code{randomE} (if any) to treat as parental effects}
  \item{sampled}{A vector indicating which individuals are sampled}
  \item{records}{A single value, array of matrix specifying data record availability - see details}
  \item{returnAllEffects}{If \code{TRUE} then all individual breeding values and environmental effects are returned}
}
\details{
\code{randomA} and \code{randomE} are square matrices with dimension equal to the sum of the number direct and indirect effects.  This must be a multiple of the number of traits, i.e. if an indirect effect is to be simulated for only one of multiple traits, those traits with no indirect effect should be included with (co)variances of zero.

\code{parentalA} and \code{parentalE} are optional vectors of characters indicating which trait positions in \code{randomA} and \code{randomE} are to be treated as indirect effects, and which effects to treat as maternal or paternal.  Valid values are 'd', 'm', and 'p', for direct, maternal indirect and paternal indirect effects, respectively.

\code{records} can be specified either (1) as asingle value to be applied to all individuals and traits, (2) as a vector the same length as the number of traits, representing trait-specific rates to be applied uniformly to all individuals, or (3) as data frames with rows for each individual and columns for each trait.  In the third option, observed patterns of data availability can be simulated by supplying 0s and 1s for missing and available individual genotypes, respectively.
}
\value{
  \item{phenotypes}{A dataframe containing phenotypes for all indivduals specified to have records.}
  \item{allEffects}{(optional) A dataframe with all direct and indirect genetic and environmental effects.}
}

\references{Morrissey et al. 2007. Journal of Evolutionary Biology 20:2309-2321., Morrissey, M.B, and A.J. Wilson, 2009. pedantics, an R package for pedigree-based genetic simulation, and pedigree manipulation, characterisation, and viewing. Molecular Ecology Resources.}
\author{Michael Morrissey \email{michael.morrissey@st-andrews.ac.uk}}

\examples{
## make up a pedigree
id<-   c("a1","a2","a3","a4","a5","a6","a7","a8","a9")
dam<-  c(NA,NA,NA,"a1","a1","a1","a4","a4","a4")
sire<- c(NA,NA,NA,"a2","a2","a2","a5","a6","a6")
pedigree<-as.data.frame(cbind(id,sire,dam))

traits<-2
## no correlations
randomA<-diag(4)
randomE<-diag(4)
parentalA<-c("d","d","m","m")
parentalE<-c("d","d","m","m")

## generate phenoypic data based on this architecture
phensim(pedigree=pedigree,traits=2,randomA=randomA,randomE=randomE,
          parentalA=parentalA,parentalE=parentalE)

## let's do it again but see how the phenotypes were composed
phensim(pedigree=pedigree,traits=2,randomA=randomA,randomE=randomE,
          parentalA=parentalA,parentalE=parentalE,returnAllEffects=TRUE)
}
\keyword{datagen}
