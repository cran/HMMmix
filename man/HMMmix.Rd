\name{HMMmix}
\alias{HMMmix}

\title{HMMmix}
\description{
Compute a hierarchical algorithm which combine the K components of a HMM. It leads to a D-states HMM with mixture of Gaussian as emission distribution. 
}
\usage{
HMMmix(data,K=5, D=4, IterationMax=1000, Threshold=1e-3,pruning.th=1,
reestimateEM=1) 
}


\arguments{
  \item{data}{A numeric matrix of observations}
  \item{K}{An integer corresponding to the total number of components. The default is 5.}  
  \item{D}{An integer corresponding to the number of cluster. The default is 4.}  
  \item{IterationMax}{An integer that corresponds to the maximum number of iterations allowed for the EM algorithm. The default is 1000.}
  \item{Threshold}{The threshold of convergence for the EM algorithm. The default is 1e-3.}
  \item{pruning.th}{A real between 0 and 1. Approximation of the merging criteria to speed up the calculations. DEFAULT VALUE is equal to 1 corresponding to not using the pruning.}
  \item{reestimateEM}{An integer that corresponds to the number of EM iterations for reestimating the parameters after combination. The default is 1.}

}

\value{

\code{HMMmix} returns a list with items:

\item{DstatesHMM}{A list containing the posterior probabilities (Tau), the transition matrix (Mat.transition) and the class of each observation (class) for the D-states HMM}
\item{MaxICLclusters}{A list containing the posterior probabilities (Tau), the transition matrix (Mat.transition) and the class of each observation (class) for the HMM where the number of states is selected by ICLs.}
\item{loglik}{A numeric vector containing the likelihood calculated after each combination of clusters.}
\item{estimation}{A list of estimated parameters for the K-states HMM, i.e. the initial HMM.}
\item{crit.mat}{A numeric matrix composed by 2 columns, the first corresponds to the BIC and the other to the ICL.}
\item{data}{A matrix n*d containing the initial data}
\item{class}{A matrix n*K containing the class of each observation for the HMM. The column of the matrix represents the number of states of the HMM. The matrix has NA values for the colums 1 to (D-1).}
\item{TAU}{A list with K elements. The i-th element of the list contains the posterior probabilities for the i-states HMM. NULL for the elements 1 to (D-1).}
\item{MAT}{A list with K elements. The i-th element of the list contains the transition matrix for the i-states HMM. NULL for the elements 1 to (D-1).}

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
data = simul$data

## Run HMMmix
result = HMMmix(data,4,2)
plotHMMmix(result)
}




