\name{rnormMix}
\alias{rnormMix}

\title{rnormMix}
\description{
Simulate a K-states HMM with mixture of multidimensional Gaussian as emission distribution.
}
\usage{
rnormMix(n,Mat.transition,init,mean,varcov,eta)
}


\arguments{
  \item{n}{An integer corresponding to the number of simulated data}  
  \item{Mat.transition}{A matrix containing the transition probabilities between the states of the Markov chain.}
\item{init}{A numeric vector specifying the initial distribution of the HMM. The sum must be equal to 1.. }
\item{mean}{List of vector of means. Each element of the list is a matrix whose number of column corresponds to the number of components for the state}
\item{varcov}{List of variance-covariance matrices.  Each element of the list is a list of matrices.}
\item{eta}{List of vector of proportions. Each element of the list is a numeric vector whose size corresponds to the number of components for the state}
}

\value{
\code{rnormMix} returns a list with items:

  \item{data}{A numeric matrix of simulated data}
   \item{class}{A vector containing the theoretical class of each observation}
   \item{Tau}{A numeric matrix containing the theoretical posterior probabilities.}
}



\examples{
## Number of simulated data
n=800

## List of the means
mean=list()
mean[[1]] =  matrix(c(1,5,1,5),2,2)
mean[[2]] = matrix(c(8,0,8,0),2,2) 

## List of proportions
eta=list()
eta[[1]] = c(0.5,0.5)
eta[[2]] = c(0.5,0.5)

## Transition matrix
Mat.transition = matrix(c(0.8,0.2,0.2,0.8),2,2,byrow=TRUE)

## List of covariance matrices
varcov=list()
varcov[[1]] = list(matrix(c(0.1,0,0,1),2,2,byrow=TRUE),matrix(c(1,0,0,0.1),2,2,byrow=TRUE))
varcov[[2]] = list(matrix(c(0.1,0,0,1),2,2,byrow=TRUE),matrix(c(1,0,0,0.1),2,2,byrow=TRUE))

## Initial distribution
init = c(0.5,0.5)

## Simulate the data
simul = rnormMix(n,Mat.transition,init,mean,varcov,eta)

## Plot data according to their class
data = simul$data
class = simul$class
plot(data,col=class)

}

\keyword{internal}



