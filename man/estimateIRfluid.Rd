\name{estimateIRfluid}
\alias{estimateIRfluid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate parameters in Inversion Recovery MRI experiments model for CSF voxel 
}
\description{
The Inversion Recovery MRI signal in voxel containing only CSF follows is modeled as  $S_{InvTime} = par[1] * abs( 1 - 2 * exp(-InvTime*par[2]) )$ dependings on two parameters. These parameters are assumed to be tissue (and scanner) dependent. 
}
\usage{
estimateIRfluid(IRdataobj, TEScale = 100, dataScale = 1000,
method = c("NLR", "QL"), varest = c("RSS", "data"),
verbose = TRUE, lower = c(0, 0), upper = c(2, 2))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{IRdataobj}{ Object of class \code{"IRdata"} as generated by function \code{\link{readIRData}}.
}
  \item{TEScale}{
Internal scale factor for Echo Times. This influences parameter scales in numerical calculations.
}
  \item{dataScale}{
Internal scale factor for MR signals. This influences parameter scales in numerical calculations.
}
  \item{method}{
  Either \code{"NLS"} for nonlinear least squares (ignores Rician bias) or \code{"QL"} for Quasi-Likelihood. The second option is more accurate but requires additional information and is computationally more expensive.  
}
  \item{varest}{
  Method to, in case of \code{method="QR"}, estimate \code{sigma}if not provided. Either from residual sums of squares (\code{"RSS"}) or MR signals (\code{"data"}) using function \code{varest} from package {\code{aws}}. Only to be used in case that no image registration was needed as preprocessing.
}
  \item{verbose}{
Logical. Provide some runtime diagnostics.
}
  \item{lower}{
Lower bounds for parameter values.
}
  \item{upper}{
Upper bounds for parameter values.
}
}
\details{
The Inversion Recovery MRI signal in voxel containing only CSF follows is modeled as  $S_{InvTime} = par[1] * abs( 1 - 2 * exp(-InvTime*par[2]) )$ dependings on two parameters. These parameters are assumed to be tissue (and scanner) dependent.
}
\value{
List of class \code{IRfluid} with components
\item{IRdata }{4D array containing the IRMRI data, first dimension refers to inversion times}
\item{InvTimes }{vector of inversion times}
\item{segm}{segmentation codes, 1 for CSF, 2 for GM, 3 for WM, 0 for out of brain}
\item{sigma}{noise standard deviation, if not specified estimated fron CSF areas in image with largest inversion time}
\item{L}{effective number of coils}
\item{Sf}{Global estimate of maximal fluid signal}
\item{Rf}{Global estimate of fluid relaxation rate}
\item{Sx}{Array of maximal signals}
\item{Rx}{Array of relaxation rates}
\item{sigma}{Array of provided or estimated noise standard deviations }
\item{Convx}{Array of convergence indicators}
\item{method}{\code{"NLS"} for nonlinear regression or \code{"QL"} for quasi likelihood.}
\item{varest}{Method used for variance estimation}
The arrays only contain entries for fluid voxel.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{
Karsten Tabelow \email{tabelow@wias-berlin.de}\cr
J\"org Polzehl \email{polzehl@wias-berlin.de}
%%  ~~who you are~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{estimateIR}}, \code{\link{estimateIRsolid}}, \code{\link{estimateIRsolidfixed}},\code{\link{smoothIRSolid}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (IRdata, InvTimes, segments, TEScale = 100, dataScale = 1000, 
    method = c("NLR", "QL"), sigma = NULL, L = 1, maxR2star = 50, 
    varest = c("RSS", "data"), verbose = TRUE, lower = c(0, 0), 
    upper = c(2, 2)) 
{
    mask <- segments == 1
    nvoxel <- sum(mask)
    ntimes <- length(InvTimes)
    itmin <- order(InvTimes)[1]
    itmax <- order(InvTimes)[ntimes]
    InvTimes[InvTimes == Inf] <- 50 * max(InvTimes[InvTimes != 
        Inf])
    dimdata <- dim(IRdata)
    if (dimdata[1] != ntimes) 
        stop("estimateIRfluid: incompatible length of InvTimes")
    if (any(dimdata[-1] != dim(mask))) 
        stop("estimateIRfluid: incompatible dimension of segments")
    InvTimesScaled <- InvTimes/TEScale
    npar <- 2
    Rx <- Sx <- Conv <- array(0, dim(mask))
    isConv <- array(0, nvoxel)
    isThresh <- array(FALSE, nvoxel)
    modelCoeff <- array(0, c(npar, nvoxel))
    if (varest[1] == "data") {
        if (verbose) 
            cat("estimating variance maps from data\n")
        ind <- (InvTimes == max(InvTimes))[1]
        ddata <- IRdata[ind, , , ]
        shat <- aws::awsLocalSigma(ddata, steps = 16, mask = (segments == 
            1), ncoils = 1, hsig = 2.5, lambda = 6, family = "Gauss")$sigma
        dim(shat) <- dimdata[-1]
        shat <- shat[segments == 1]
        shat[shat == 0] <- quantile(shat, 0.8)
        shat <- shat
        if (is.null(sigma)) 
            sigma <- median(shat)
        else shat <- NULL
    }
    if (method[1] == "QL") {
        if (is.null(sigma)) {
            method <- "NLR"
            warning("estimateIRfluid: method QL needs sigma estimated or supplied")
        }
        sig <- sigma/dataScale
        CL <- sig * sqrt(pi/2) * gamma(L + 0.5)/gamma(L)/gamma(1.5)
    }
    dim(IRdata) <- c(dimdata[1], prod(dim(segments)))
    IRdataFluid <- IRdata[, segments == 1]
    thetas <- matrix(0, 2, nvoxel)
    thetas[1, ] <- IRdataFluid[itmax, ]/dataScale
    thetas[2, ] <- -log((IRdataFluid[itmin]/dataScale + thetas[1, 
        ])/2)/InvTimes[itmin] * TEScale
    if (verbose) {
        cat("Start estimation in", nvoxel, "voxel at", format(Sys.time()), 
            "\n")
        pb <- txtProgressBar(0, nvoxel, style = 3)
    }
    for (xyz in 1:nvoxel) {
        ivec <- IRdataFluid[, xyz]/dataScale
        th <- thetas[, xyz]
        res <- if (method[1] == "NLR") 
            try(nls(ivec ~ IRhomogen(par, InvTimesScaled), data = list(InvTimesScaled), 
                start = list(par = th), control = list(maxiter = 200, 
                  warnOnly = TRUE)), silent = TRUE)
        else try(nls(ivec ~ IRhomogenQL(par, InvTimesScaled, 
            CL, sig, L), data = list(InvTimesScaled, CL = CL, 
            sig = sig, L = L), start = list(par = th), control = list(maxiter = 200, 
            warnOnly = TRUE)), silent = TRUE)
        if (class(res) != "try-error") {
            thhat <- coef(res)
            outofrange <- any(thhat != pmin(upper, pmax(lower, 
                thhat)))
        }
        if (class(res) == "try-error" || outofrange) {
            th <- pmin(upper, pmax(lower, th))
            res <- if (method[1] == "NLR") 
                try(nls(ivec ~ IRhomogen(par, InvTimesScaled), 
                  data = list(InvTimes = InvTimesScaled), start = list(par = th), 
                  algorithm = "port", control = list(maxiter = 200, 
                    warnOnly = TRUE), lower = lower, upper = upper), 
                  silent = TRUE)
            else try(nls(ivec ~ IRhomogenQL(par, InvTimesScaled, 
                CL, sig, L), data = list(InvTimesScaled = InvTimesScaled, 
                CL = CL, sig = sig, L = L), start = list(par = th), 
                algorithm = "port", control = list(maxiter = 200, 
                  warnOnly = TRUE), lower = lower, upper = upper), 
                silent = TRUE)
        }
        if (class(res) != "try-error") {
            sres <- if (varest[1] == "RSS") 
                getnlspars(res)
            else getnlspars2(res, shat[, xyz], sind)
            isConv[xyz] <- as.integer(res$convInfo$isConv)
            modelCoeff[, xyz] <- sres$coefficients
        }
        if (verbose) 
            if (xyz\%/\%1000 * 1000 == xyz) 
                setTxtProgressBar(pb, xyz)
    }
    Rx[mask] <- modelCoeff[2, ]
    Sx[mask] <- modelCoeff[1, ]
    Conv[mask] <- isConv
    Sf <- median(modelCoeff[1, ], na.rm = TRUE)
    Rf <- median(modelCoeff[2, ], na.rm = TRUE)
    if (verbose) {
        close(pb)
        cat("Finished estimation", format(Sys.time()), "\n", 
            "Sf", Sf, "Rf", Rf, "\n")
    }
    list(Sf = Sf, Rf = Rf, Sx = Sx, Rx = Rx, sigma = sigma, Conv = Conv)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
 \keyword{ models }
 \keyword{ regression }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
 \concept{ IRMRI }
% \concept{ ~cpt2 }
% Use only one concept per line.
