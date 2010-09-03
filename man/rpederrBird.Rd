\name{rpederrBird}
\alias{rpederrBird}
\title{Permutes a bird pedigree to create a plausible complete pedigree}
\description{
  Permutes a pedigree to create a plausible complete pedigree, with errors accounted for in terms of extra-pair paternity.
}
\usage{
rpederrBird(assumedPedigree, founders = NULL, sex = NULL, samp = NULL, 
             sireS = NULL, damS = NULL, cohort = NULL, first = NULL, 
             last = NULL, broods = NULL, monoecey = 0, broodData = NULL, 
             EPPsireData = NULL, propEPPbroods = 1, 
             propEPPchicksGivenEPPbrood = 0.5, propEPPbroodsTwoFathers = 0, 
             modifyAssumedPedigree = 0,EPPlambda=0.001,
             malePhenotypes=NULL,EPPbeta=0,EPPgamma=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{assumedPedigree}{A pedigree}
  \item{founders}{A vector of indicator variables denoting founder status (1=founder, 0=non-founder)}
  \item{sex}{A vector of indicator variables denoting sex (0=male,1=female,anything else=unknown)}
  \item{samp}{A vector denoting whether or not individuals are sampled (1), or dummy indivdiuals (0) added to the pedigree for the purpose of simulating potential "true" pedigree links that go outside the sampled population}
  \item{sireS}{Proportion of "true" simulated sires that are to be taken from the unsampled portion of the pedgiree.}
  \item{damS}{Proportion of "true" simulated dams that are to be taken from the unsampled portion of the pedgiree.}
  \item{cohort}{A numeric vector the same length as the pedigree containing cohorts}
  \item{first}{A numeric vector the same length as the pedigree indicating the first cohort for which an individual is to be considered a potential parent}
  \item{last}{A numeric vector the same length as the pedigree indicating the last cohort for which an individual is to be considered a potential parent}
  \item{broods}{A numeric vector the same length as the pedigree indicating the brood to which each individual belongs.  Can be any integer value for indivduals that do not belong to broods as spedified in broodData (see below).}
  \item{monoecey}{An indicator specifying whether or not bisexuality is allowed (0=no (default), 1=yes)}
  \item{broodData}{A five column data frame containing (1) brood IDs corresponding to \code{broods}, (2) years, (3) brood latitudes, (4) brood longitudes, and (5) an indicator (0=no or 1=yes) specifying whether or not the particular brood potentially contains extra-pair chicks.}
  \item{EPPsireData}{A four column data frame with (1) ids of all sires that may be extra-pair fathers, (2) years (cohorts) in which an indivdual might have been an extra pair father (one separate line for each year), and  (3) latitude and (4) longitude in each year.}
  \item{propEPPbroods}{Proportion of broods with extra-pair chicks.}
  \item{propEPPchicksGivenEPPbrood}{Mean proportion of chicks that have extra-pair fathers in extra-pair broods.}
  \item{propEPPbroodsTwoFathers}{Proportion of extra-pair broods that have two extra-pair fathers}
  \item{modifyAssumedPedigree}{An indicator variable specifying whether or not an assumed pedigree with errors but no missing links should be supplied.}
  \item{EPPlambda}{The rate of decline of the probability of extra pair paternity with distance between males and nests.}
  \item{malePhenotypes}{(optional) if specified, zp will take these values in the calculation of male relative fertilities, see details}
  \item{EPPbeta}{(optional) coefficient of linear selection acting through undetected extra-pair breeding success, see details}
  \item{EPPgamma}{(optional) coefficient of quadratic selection acting through undetected extra-pair breeding success, see details}
}

\value{
  \item{assumedPedigree}{echos the supplied pedigree}
  \item{truePedigree}{A plausible pedigree with no errors and no missing links}
  \item{supplementalPedigree}{A plausible pedigree with errors but no missing links}
}
\references{Morrissey, M.B, and A.J. Wilson, 2009. pedantics, an R package for pedigree-based genetic simulation, and pedigree manipulation, characterisation, and viewing. Molecular Ecology Resources.}
\author{Michael Morrissey \email{michael.morrissey@ed.ac.uk}}

\seealso{ \code{\link{rpederr}} for a simpler module that may be adequate when pedigree errors are (assumably) independent between offspring.}
\examples{
## test data for a single brood
id<-   c("a1","a2","a3","a4","a5","a6","a7","a8","a9")
dam<-  c(NA,NA,NA,NA,"a1","a1","a1","a1","a1")
sire<-  c(NA,NA,NA,NA,"a2","a2","a2","a2","a2")
broods<- c(NA,NA,NA,NA,"r1","r1","r1","r1","r1")
founders<-c(1,1,1,1,0,0,0,0,0)
samp<- c(1,1,1,1,1,1,1,1,1)
sex<-  c(1,0,0,0,0,0,1,0,0)
sireS<- rep(1,9)
damS<- rep(1,9)
cohort<-c(1,1,1,1,2,2,2,2,2)
first<-c(2,2,2,3,3,3,3,3,3)
last<-c(2,2,2,3,3,3,3,3,3)

broodIDs<-c("r1")
broodYears<-2
broodN<-c(10)
broodE<-c(10)
broodEPPpotential<-c(1)
EPPsireIDs<-c("a2","a3","a4")
EPPsireYears<-c(2,2,2)
sireN<-c(10,11,5)
sireE<-c(10,11,5)
EPPlambda<- 0.4

## put the test data in the required forms
assumedPedigree<-as.data.frame(cbind(id,sire,dam))
broodData<-as.data.frame(cbind(broodIDs,broodYears,broodN,
                                            broodE,broodEPPpotential))
EPPsireData<-as.data.frame(cbind(EPPsireIDs,EPPsireYears,sireN,sireE))

## simulate extra-pair paternity (try a few times and see how often
## the paterity goes to the nearer candidate EPP father):
rpederrBird(assumedPedigree=assumedPedigree,founders=founders,sex=sex,
        samp=samp,sireS=sireS,damS=damS,cohort=cohort,first=first,
        last=last,broods=broods,broodData=broodData,EPPsireData=EPPsireData,
        propEPPbroods=1,propEPPchicksGivenEPPbrood=0.5,
        propEPPbroodsTwoFathers=0,modifyAssumedPedigree=0)

}
\keyword{datagen}
