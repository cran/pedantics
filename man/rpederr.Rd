\name{rpederr}
\alias{rpederr}
\title{Permutes a pedigree to create a plausible complete pedigree}
\description{
  Given estimates of indivdiual life histories and rates and patterns of errors in pedigree links, \code{rpederr} probabilistically assigns "true" parents given an incomplete and potentially eroneous pedigree.
}
\usage{
rpederr(assumedPedigree, founders = NULL, sex = NULL, samp = NULL, 
         sireE = NULL, damE = NULL, sireS = NULL, damS = NULL, 
         cohort = NULL, first = NULL, last = NULL, monoecey = 0, 
         modifyAssumedPedigree = 0)
}
\arguments{
  \item{assumedPedigree}{A pedigree}
  \item{founders}{A vector of indicator variables denoting founder status (1=founder, 0=non-founder)}
  \item{sex}{A vector of indicator variables denoting sex (0=male,1=female,anything else=unknown)}
  \item{samp}{A vector denoting whether or not individuals are sampled (1), or dummy indivdiuals (0) added to the pedigree for the purpose of simulating potential "true" pedigree links that go outside the sampled population}
  \item{sireE}{Sire assignment error rates, see details}
  \item{damE}{Dam assignment error rates, see details}
  \item{sireS}{Proportion of "true" simulated sires that are to be taken from the unsampled portion of the pedgiree.}
  \item{damS}{Proportion of "true" simulated dams that are to be taken from the unsampled portion of the pedgiree.}
  \item{cohort}{A numeric vector the same length as the pedigree containing cohorts}
  \item{first}{A numeric vector the same length as the pedigree indicating the first cohort for which an individual is to be considered a potential parent}
  \item{last}{A numeric vector the same length as the pedigree indicating the last cohort for which an individual is to be considered a potential parent}
  \item{monoecey}{An indicator specifying whether or not bisexuality is allowed (0=no (default), 1=yes)}
  \item{modifyAssumedPedigree}{An indicator variable specifying whether or not an assumed pedigree with errors but no missing links should be supplied.}
}
\value{
  \item{assumedPedigree}{echos the supplied pedigree}
  \item{truePedigree}{A plausible pedigree with no errors and no missing links}
  \item{supplementalPedigree}{A plausible pedigree with errors but no missing links}
}
\references{Morrissey et al. 2007. Journal of Evolutionary Biology 20:2309-2321., Morrissey, M.B, and A.J. Wilson, 2009. pedantics, an R package for pedigree-based genetic simulation, and pedigree manipulation, characterisation, and viewing. Molecular Ecology Resources.}
\author{Michael Morrissey \email{michael.morrissey@st-andrews.ac.uk}}

\seealso{\code{\link{rpederrBird}}}
\examples{
id<-   c("a1","a2","a3","a4","a5","a6","a7","a8","a9")
dam<-  c(NA,NA,NA,"a1","a1","a1","a4","a4","a4")
sire<- c(NA,NA,NA,NA,NA,NA,"a5","a5","a5")
found<-c(1,1,1,0,0,0,0,0,0)
samp<- c(1,1,1,1,1,1,1,1,1)
sex<-  c(1,0,0,1,0,0,1,0,0)
dade<- rep(0,9)
mume<- rep(0,9)
dads<- rep(1,9)
mums<- rep(1,9)
cohort<-c(1,1,1,2,2,2,3,3,3)
first<-c(2,2,2,3,3,3,4,4,4)
last<-c(2,2,2,3,3,3,4,4,4)
pedigree<-as.data.frame(cbind(id,sire,dam))

### don't simulate any errors, just fill in the missing sires
rpederr(assumedPedigree=pedigree,founders=found,sex=sex,
                  samp=samp,cohort=cohort,first=first,last=last)

## fill in the missing sires, and additionally simulate a problem 
## with the second maternal sibship note that it is probabilistic, 
## so this example may need to be run a couple of times before the
## error comes up, given the very small example dataset
fatherErrors<-c(0,0,0,0,0,0,0.8,0.8,0.8)
rpederr(assumedPedigree=pedigree,founders=found,sex=sex,samp=samp,
           sireE=fatherErrors,cohort=cohort,first=first,last=last)
}
\keyword{datagen}
