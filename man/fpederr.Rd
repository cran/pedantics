\name{fpederr}
\alias{fpederr}
\title{Simulates a pedigree with errors and missing data from a complete pedigree.}
\description{
 Implements the 'forward' approach to producing pairs of pedigrees for power and sensitivity analyses.
}
\usage{
fpederr(truePedigree, founders = NULL, sex = NULL, samp = NULL, 
          sireE = NULL, damE = NULL, sireA = NULL, damA = NULL, 
          cohort = NULL, first = NULL, last = NULL, monoecey = 0, 
          modifyAssumedPedigree = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{truePedigree}{ A complete pedigree with records for all individuals and parental ID's for all non-founders }
  \item{founders}{A vector the same length as the pedigree containing indicator variables 1 = founder, 0 = non-founder}
  \item{sex}{A vector the same length as the pedigree indicating sex, 0=male, 1=female, any other value = unknown sex}
  \item{samp}{A vector the same length as the pedigree indicating whether or not each individual is sampled (1), or an unsampled dummy individual (0).}
  \item{sireE}{Value(s) indicating the paternal error rate.  If it is a single number (between 0 and 1), it is applied to the entire pedigree; if it is a vector the length of the pedigree, then probabilities can vary among individuals.}
  \item{damE}{Value(s) indicating the maternal error rate.  If it is a single number (between 0 and 1), it is applied to the entire pedigree; if it is a vector the length of the pedigree, then probabilities can vary among individuals.}
  \item{sireA}{Value(s) indicating the paternal pedigree link assignment rate.  If it is a single number (between 0 and 1), it is applied to the entire pedigree; if it is a vector the length of the pedigree, then probabilities can vary among individuals.}
  \item{damA}{Value(s) indicating the maternal pedigree link assignment rate.  If it is a single number (between 0 and 1), it is applied to the entire pedigree; if it is a vector the length of the pedigree, then probabilities can vary among individuals.}
  \item{cohort}{A numeric vector the same length as the pedigree containing cohorts}
  \item{first}{A numeric vector the same length as the pedigree indicating the first cohort for which an individual is to be considered a potential parent}
  \item{last}{A numeric vector the same length as the pedigree indicating the last cohort for which an individual is to be considered a potential parent}
  \item{monoecey}{An indicator specifying whether or not bisexuality is allowed (0=no (default), 1=yes)}
  \item{modifyAssumedPedigree}{An indicator variable specifying whether or not an assumed pedigree with errors but no missing links should be supplied.}
}
\value{
  \item{assumedPedigree}{A pedigree differing from the supplied pedigree
so as to mimic patterns of pedigree errors and missing data that 
might occur in a real study.}
  \item{truePedigree}{Echos the pedigree supplied.}
  \item{supplementalPedigree}{(optional) a 'assumed' pedigree containing
 errorsbut no missing links.}
}
\references{Morrissey et al. 2007. Journal of Evolutionary Biology 20:2309-2321., Morrissey, M.B, and A.J. Wilson, 2009. pedantics, an R package for pedigree-based genetic simulation, and pedigree manipulation, characterisation, and viewing. Molecular Ecology Resources.}
\author{Michael Morrissey \email{michael.morrissey@ed.ac.uk}}
\seealso{\code{\link{rpederr}},\code{\link{fpederr}}}
\examples{
testData<-as.data.frame(matrix(c(
1,	NA,	NA,	1,	1,	1,	2,	2,
2,	NA,	NA,	1,	1,	1,	2,	2,
3,	NA,	NA,	1,	1,	1,	2,	2,
4,	NA,	NA,	1,	0,	1,	2,	2,
5,	NA,	NA,	1,	0,	1,	2,	2,
6,	1,	4,	0,	-1,	2,	3,	3,
7,	1,	4,	0,	-1,	2,	3,	3,
8,	1,	4,	0,	-1,	2,	3,	3,
9,	1,	4,	0,	-1,	2,	3,	3,
10,	2,	5,	0,	-1,	2,	3,	3,
11,	2,	5,	0,	-1,	2,	3,	3,
12,	2,	5,	0,	-1,	2,	3,	3,
13,	2,	5,	0,	-1,	2,	3,	3,
14,	3,	5,	0,	-1,	2,	3,	3,
15,	3,	5,	0,	-1,	2,	3,	3,
16,	3,	5,	0,	-1,	2,	3,	3,
17,	3,	5,	0,	-1,	2,	3,	3),
17,8,byrow=TRUE))

names(testData)<-c("id","dam","sire","founder","sex",
                                  "cohort","first","last")
pedigree<-as.data.frame(cbind(testData$id,testData$dam,
                                           testData$sire))
for(x in 1:3) pedigree[,x]<-as.factor(pedigree[,x])
names(pedigree)<-c("id","dam","sire")
pedigree

## some missing sire links:
fpederr(truePedigree=pedigree,founders=testData$founder,
        sex=testData$sex,sireA=0.5,cohort=testData$cohort,
        first=testData$first,last=testData$last)$assumedPedigree

## some erroneous sire links:
fpederr(truePedigree=pedigree,founders=testData$founder,
        sex=testData$sex,sireE=0.5,cohort=testData$cohort,
        first=testData$first,last=testData$last)$assumedPedigree
}
\keyword{datagen}
