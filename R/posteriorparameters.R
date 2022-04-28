#' Augment a DAG with parameters
#'
#' \code{DAGparameters} takes a DAG and augments it with parameters.
#' For binary data these are the parameters of the posterior beta
#' distributions and its mean. For continuous data, these are parameters 
#' of the posterior distributions of the edge coefficients from arXiv:2010.00684.
#' There is support for user-defined augmentation, with the caveat that it
#' must match the output format of either the binary or continuous cases. 
#'
#' @param incidence a single adjacency matrix with entry [i,j] equal to 1 
#' when a directed edge exists from node i to node j
#' @param dataParams the data and parameters used to learn the DAGs derived from the
#' \code{\link[BiDAG]{scoreparameters}} function of the BiDAG package
#' @param unrollDBN logical indicating whether to unroll a DBN to a full DAG over
#' all time slices (TRUE, default) or to use the compact representation (FALSE)
#'
#' @return the DAG and a list of parameters for each node given
#' its parents
#'
#' @examples
#'
#' scoreParam <- BiDAG::scoreparameters("bde", BiDAG::Asia)
#' AsiaParam <- DAGparameters(BiDAG::Asiamat, scoreParam)
#'
#' @export

DAGparameters <- function(incidence, dataParams, unrollDBN = TRUE) {
  # add parameters to an incidence matrix

  if (dataParams$DBN && !dataParams$type %in% c("bde", "bge")) {
    stop("The implementation for DBNs is currently only for the BDe and BGe score.")
  }

  bg_flag <- FALSE
  if (dataParams$type == "usr"){
    localtype <- dataParams$pctesttype
    if (!is.null(dataParams$bgremove)){
      if (dataParams$bgremove && dataParams$bgn > 0) {
        bg_flag <- TRUE
      }
    }
  } else {
    localtype <- dataParams$type
  }
    
  if (!localtype %in% c("bde", "bge")) {
    stop("The implementation is currently only for the BDe and BGe score, or user-defined scores with the same structure.")
  }

  if (dataParams$DBN) {

    n <- dataParams$n # number of nodes including background
    nsmall <- dataParams$nsmall # number of nodes without background
    bgn <- dataParams$bgn # number of background nodes
    slices <- dataParams$slices # number of time slices

### First slice parameters    
# remember background nodes were placed at the end of first slice, which we need to undo!
    
    dataParams_first <- dataParams$firstslice
    if(bgn > 0){ # place backgound nodes at the start again
      reorder <- c(1:bgn + nsmall, 1:nsmall)
      if (dataParams$type == "bde") { # bde score
        dataParams_first$data <- dataParams_first$data[, reorder]
        dataParams_first$d1 <- dataParams_first$d1[, reorder]
        dataParams_first$d0 <- dataParams_first$d0[, reorder]
      } else { # bge score
        dataParams_first$TN <- dataParams_first$TN[reorder, reorder]
      }

    }
    params_first <- DAGparameters(incidence[1:n, 1:n], dataParams_first)

## Other slices parameters    
# remember the later times were placed first, which we need to undo!     

    dataParams_other <- dataParams$otherslices

    reorder <- c(1:n + nsmall, 1:nsmall) # put them back in order
    if (dataParams$type == "bde") { # bde score
      dataParams_other$data <- dataParams_other$data[, reorder]
      dataParams_other$d1 <- dataParams_other$d1[, reorder]
      dataParams_other$d0 <- dataParams_other$d0[, reorder]
    } else { # bge score
      dataParams_other$TN <- dataParams_other$TN[reorder, reorder]
    }

    params <- DAGparameters(incidence, dataParams_other)

### Combine parameters from the slices
    
    if (dataParams$type == "bde") { # bde score
      allalphas <- params$alphas
      allalphas[1:n] <- params_first$alphas # update the first slice
      allbetas <- params$betas
      allbetas[1:n] <- params_first$betas # update the first slice
      allpmeans <- params$pmeans
      allpmeans[1:n] <- params_first$pmeans # update the first slice
    } else { # bge score
      allmus <- params$mus
      allmus[1:n] <- params_first$mus # update the first slice
      allsigmas <- params$sigmas
      allsigmas[1:n] <- params_first$sigmas # update the first slice
      alldfs <- params$dfs
      alldfs[1:n] <- params_first$dfs # update the first slice
    }
    
### Need to unroll the DBN, if more than 2 slices and we choose to unroll
    
    if (slices > 2 && unrollDBN) {
      nbig <- n + nsmall*(slices - 1)
      
      incidence_unroll <- matrix(0, nbig, nbig)
      incidence_unroll[1:nrow(incidence), 1:ncol(incidence)] <- incidence
      
      inc_names_unroll <-  paste(rep(colnames(incidence)[bgn+1:nsmall], (slices - 2)), rep(3:slices, each=nsmall), sep="_")
      colnames(incidence_unroll) <- c(colnames(incidence), inc_names_unroll)
      rownames(incidence_unroll) <- c(colnames(incidence), inc_names_unroll)
      
      if (dataParams$type == "bde") { # bde score
        allalphas_unroll <- vector("list", nbig)
        allbetas_unroll <- vector("list", nbig)
        allpmeans_unroll <- vector("list", nbig)
        
        allalphas_unroll[1:ncol(incidence)] <- allalphas
        allbetas_unroll[1:ncol(incidence)] <- allbetas
        allpmeans_unroll[1:ncol(incidence)] <- allpmeans
      } else { # bge score
        allmus_unroll <- vector("list", nbig)
        allsigmas_unroll <- vector("list", nbig)
        alldfs_unroll <- vector("list", nbig)
        
        allmus_unroll[1:ncol(incidence)] <- allmus
        allsigmas_unroll[1:ncol(incidence)] <- allsigmas
        alldfs_unroll[1:ncol(incidence)] <- alldfs
      }
      
      for (ii in 1:(slices - 2)) {
        block_rows <- n - nsmall + 1:(2*nsmall)
        block_cols <- n + 1:nsmall

        incidence_unroll[block_rows + nsmall*ii, block_cols + nsmall*ii] <- incidence[block_rows, block_cols]

        if (dataParams$type == "bde") { # bde score
          allalphas_unroll[block_cols + nsmall*ii] <- allalphas[block_cols]
          allbetas_unroll[block_cols + nsmall*ii] <- allbetas[block_cols]
          allpmeans_unroll[block_cols + nsmall*ii] <- allpmeans[block_cols]
        } else { # bge score
          allmus_unroll[block_cols + nsmall*ii] <- allmus[block_cols]
          allsigmas_unroll[block_cols + nsmall*ii] <- allsigmas[block_cols]
          alldfs_unroll[block_cols + nsmall*ii] <- alldfs[block_cols]
        }

        if (bgn > 0) { # if there are background nodes, repeat across slices
          block_rows <- 1:bgn
          incidence_unroll[block_rows, block_cols + nsmall*ii] <- incidence[block_rows, block_cols]
        }
      }
      
      incidence <- incidence_unroll
      
      if (dataParams$type == "bde") { # bde score
        allalphas <- allalphas_unroll
        allbetas <- allbetas_unroll
        allpmeans <- allpmeans_unroll
      } else { # bge score
        allmus <- allmus_unroll
        allsigmas <- allsigmas_unroll
        alldfs <- alldfs_unroll
      }      

    }
 
  } else {
  
  n <- nrow(incidence) # number of nodes in DAG

  if (localtype == "bde") { # bde score
    allalphas <- vector("list", n)
    allbetas <- vector("list", n)
    allpmeans <- vector("list", n)
  } else { # bge score
    allmus <- vector("list", n)
    allsigmas <- vector("list", n)
    alldfs <- vector("list", n)
  }

  for (j in 1:n) {
    parentNodes <- which(incidence[, j]==1)
    if (dataParams$type == "usr") {
      tempResult <- usrDAGparametersCore(j, parentNodes, dataParams)
    } else {
      tempResult <- DAGparametersCore(j, parentNodes, dataParams)
    }
    if (localtype == "bde") { # bde score
      allalphas[[j]] <- tempResult$alphas
      allbetas[[j]] <- tempResult$betas
      allpmeans[[j]] <- tempResult$pmeans
    } else { # bge score
      allmus[[j]] <- tempResult$mus
      allsigmas[[j]] <- tempResult$sigmas
      alldfs[[j]] <- tempResult$dfs
    }
  }
  
  if (bg_flag) { # remove all background nodes
    to_keep <- setdiff(1:n, dataParams$bgnodes)
    incidence <- incidence[to_keep, to_keep]
    if (localtype == "bde") { # bde score
      allalphas <- allalphas[to_keep]
      allbetas <- allbetas[to_keep]
      allpmeans <- allpmeans[to_keep]
    } else { # bge score
      allmus <- allmus[to_keep]
      allsigmas <- allsigmas[to_keep]
      alldfs <- alldfs[to_keep]
    }
  }

  }
  
  posteriorParams <- list()
  posteriorParams$DAG <- incidence
  if (localtype == "bde") { # bde score
    posteriorParams$alphas <- allalphas
    posteriorParams$betas <- allbetas
    posteriorParams$pmeans <- allpmeans
  } else { # bge score
    posteriorParams$mus <- allmus
    posteriorParams$sigmas <- allsigmas
    posteriorParams$dfs <- alldfs
  }

  return(posteriorParams)
}


usrDAGparametersCore <- function(j, parentNodes, param) {
  # this is a template function for computing the parameters 
  # and their posterior distribution. It requires the output
  # to be in the corresponding BDe or BGe format
  param$type <- param$pctesttype
  # for the template we just use the BDe or BGe
  DAGparametersCore(j, parentNodes, param)
}


DAGparametersCore <- function(j, parentNodes, param) {
  # this function computes the parameters and their posterior distribution
  # for a given node with parent set and for the data

  switch(param$type,
    "bge" = {
      coreParams <- list()
      lp <- length(parentNodes) # number of parents
      if (lp > 0) {# otherwise no regression coefficients
        df <- param$awpN - param$n + lp + 1
        R11 <- param$TN[parentNodes, parentNodes]
        R12 <- param$TN[parentNodes, j]
      
        R11inv <- solve(R11) # could avoid inversions, but here for simplicity
        mb <- R11inv %*% R12 # mean part
        divisor <- param$TN[j, j] - R12 %*% mb

        coreParams$mus <- as.vector(mb)
        coreParams$sigmas <- as.numeric(divisor/df) * R11inv
        coreParams$dfs <- df
      } else {
        coreParams$mus <- NA
        coreParams$sigmas <- NA
        coreParams$dfs <- NA
      }
      return(coreParams)
    },
    "bde" = {
      lp <- length(parentNodes) # number of parents
      noParams <- 2^lp # number of binary states of the parents
      chi <- param$chi

      alphas <- rep(NA, noParams)
      betas <- rep(NA, noParams)

      switch(as.character(lp),
        "0"={ # no parents
          N1 <- sum(param$d1[, j])
          N0 <- sum(param$d0[, j])
          NT <- N0 + N1
          alphas <- (N1 + chi/(2*noParams))
          betas <- (N0 + chi/(2*noParams))
        },
        "1"={ # one parent
          summys <- param$data[, parentNodes]
          for (i in 1:noParams-1) {
            totest <- which(summys==i)
            N1 <- sum(param$d1[totest, j])
            N0 <- sum(param$d0[totest, j])
            NT <- N0 + N1
            alphas[i+1] <- (N1 + chi/(2*noParams))
            betas[i+1] <- (N0 + chi/(2*noParams))
           }
         },
         { # more parents
           summys <- colSums(2^(c(0:(lp-1)))*t(param$data[, parentNodes]))
           N1s <- collectC(summys, param$d1[, j], noParams)
           N0s <- collectC(summys, param$d0[, j], noParams)
           NTs <- N1s + N0s
           alphas <- (N1s + chi/(2*noParams))
           betas <- (N0s + chi/(2*noParams))
         }
      )

      coreParams <- list()
      coreParams$alphas <- alphas
      coreParams$betas <- betas
      coreParams$pmeans <- alphas/(alphas + betas)

      return(coreParams)
    }
  )
}



SampleParameters <- function(DAGparams, type = "bde") {
  # this function resamples the probability parameters from the posterior 
  # beta distributions or from the posterior edge coefficient distribution
  # for an unrolled DBN they are sampled at each slice
  # rather than sampled once and copied over the slices

  if (type == "bde") {
    sampledps <- DAGparams$pmeans
  } else { # bge version
    sampledps <- DAGparams$mus
  }

  n <- length(sampledps)

  for(jj in 1:n){
    if (type == "bde") {
      as <- DAGparams$alphas[[jj]]
      bs <- DAGparams$betas[[jj]]
      ps <- rep(0,length(as))
      for(ii in 1:length(as)){
        ps[ii] <- stats::rbeta(1, as[ii], bs[ii])
      }
    } else { # bge version
      if (!is.na(DAGparams$dfs[[jj]])) {
        ps <- mvtnorm::rmvt(1, sigma = as.matrix(DAGparams$sigmas[[jj]]), 
                                         df = DAGparams$dfs[[jj]], delta = DAGparams$mus[[jj]])
      } else {
        ps <- NA
      }
    }
    sampledps[[jj]] <- ps
  }

  return(sampledps)
}



BinaryScoreAgainstDAG <- function(DAGparams, dataToScore) {
  # score a set of binary vectors against a DAG with parameters

  n <- nrow(DAGparams$DAG) # number of nodes

  logscoresagainstDAG<-matrix(NA,nrow(dataToScore),n)
  for (j in 1:n)  {
    parentNodes <- which(DAGparams$DAG[, j]==1)
    logscoresagainstDAG[, j] <- BinaryScoreAgainstDAGcore(j, parentNodes, DAGparams, dataToScore)
  }

  return(logscoresagainstDAG)
}



BinaryScoreAgainstDAGcore <- function(j, parentNodes, DAGparams, dataToScore) {
  # score of a single node of binary vectors against a DAG

  sampleNodeScores <- rep(NA, nrow(dataToScore)) # store the log scores

  lp <- length(parentNodes) # number of parents
  noParams <- 2^lp # number of binary states of the parents

  switch(as.character(lp),
    "0"={ # no parents
      theta <- DAGparams$pmeans[[j]] # the probability of each state
      sampleNodeScores[which(dataToScore[, j]==1)] <- log(theta) # log scores of 1s
      sampleNodeScores[which(dataToScore[, j]==0)] <- log(1-theta) # log scores of 0s
    },
    "1"={ # one parent
      summysfull<-dataToScore[,parentNodes]

      for (i in 1:noParams-1) {
        theta <- DAGparams$pmeans[[j]][i+1] # the probability of each state
        toScore <- which(summysfull==i)
        sampleNodeScores[toScore[which(dataToScore[toScore,j]==1)]] <- log(theta) # log scores of 1s
        sampleNodeScores[toScore[which(dataToScore[toScore,j]==0)]] <- log(1-theta) # log scores of 0s
      }
    },
    { # more parents
      summysfull <- colSums(2^(c(0:(lp-1)))*t(dataToScore[, parentNodes]))

      # find the entries where the child is 1
      toScore<-which(dataToScore[, j]==1)

      #Ns <- collectC(summysfull[toScore], rep(1, length(toScore)), noParams) # works like the table command

      Ns <- tabulate(summysfull[toScore]+1,noParams) # can use tabulate instead of collectC, but we need to add one
      tempScoreVec <- rep(log(DAGparams$pmeans[[j]]), Ns) # make relevant number of copies of each log score
      sampleNodeScores[toScore] <- tempScoreVec[rank(summysfull[toScore], ties.method="first")] # use the rank function to map scores to entries

      # find the entries where the child is 0
      toScore<-which(dataToScore[,j]==0)

      Ns<-tabulate(summysfull[toScore]+1,noParams) # again we need to add one
      tempScoreVec<-rep(log(1-DAGparams$pmeans[[j]]),Ns) # make relevant number of copies of each log score
      sampleNodeScores[toScore]<-tempScoreVec[rank(summysfull[toScore],ties.method="first")] # use the rank function to map scores to entries
    }
  )

  return(sampleNodeScores)
}


