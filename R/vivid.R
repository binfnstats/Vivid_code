
#' Title
#'
#' @param x input matrix, rows represent observations and columns represent
#' features. column names and row names are recommended.
#' @param y response variable, restricted to two classes.
#' @param bootstraps integer, the number of bootstraps to run.
#' @param cores integer, number of cores to use for parallel processing. VIVID
#' makes uses of the "future" and "furrr" packages and will set the plan to
#' "multiprocess". The cores parameter is directly fed to workers arguement of
#' future::plan.
#' @param seed integer, passed to set.seed function, a full L'Ecuyer-CMRG RNG seed (vector of 7 integers).
#' @param minSize minimum cluster size to consider when clustering.
#' @param lambda Used in cv.glmnet.
#' @param compareMethod AIC,  AICc, BIC or EBIC.
#' @param gamma gamma value for EBIC.
#' @importFrom furrr future_options
#' @importFrom magrittr %>%
#'
#' @return
#' @export
#'
#' @examples

vivid =
  function(x,
           y,
           bootstraps = 100,
           cores = parallel::detectCores() - 1,
           seed = 1234567,
           minSize = 2,
           lambda = 'lambda.1se',
           compareMethod = 'EBIC',
           gamma = 1) {
    # string to factor function
    if (base::is.vector(y)) {
      y = base::factor(y)
    }

    nClass = length(levels(y))

    if (nClass != 2) {
      stop(
        base::paste(
          "The response vector contains",
          nClass,
          "classes and should only contain 2"
        )
      )
    }

    # params
    nObs = base::nrow(x)

    base::RNGkind("L'Ecuyer-CMRG")
    base::set.seed(seed)

    # Exponential weights of Bootstrap
    mWeights = base::matrix(
      stats::rexp(nObs * bootstraps),
      nrow = bootstraps,
      ncol = nObs,
      byrow = T
    )
    mWeights = mWeights / base::rowSums(mWeights) * nObs

    lWeights = base::split(x = mWeights,
                           f = base::seq.int(bootstraps))

    # processing with future/furrr
    future::plan(strategy = future::multiprocess,
                 workers = cores)

    bootstraps = furrr::future_map(lWeights,
                                   ~ vivid_reg(.x, x, y, crossfold = 10, lambda = lambda),
                                   .options = future_options(seed = TRUE))

    bootVal = vivid_df(bootstraps)

    # Change the values to ranks
    bootRank = bootVal %>%
      dplyr::select(-1) %>%
      base::abs() %>%
      dplyr::transmute_all(base::list(rank)) %>%
      as.matrix()

    base::rownames(bootRank) = base::colnames(x)

    # Creates variance rank difference matrix using coefficient input
    bootVar = matrixStats::rowVars(bootRank)
    rankVarDiff = base::outer(bootVar, bootVar, '+') - 2 * stats::cov(x = t(bootRank))
    base::colnames(rankVarDiff) = base::rownames(bootRank)
    base::rownames(rankVarDiff) = base::rownames(bootRank)

    #Clustering for rank difference matrix
    varDist = stats::dist(x = rankVarDiff)
    varClust = stats::hclust(d = varDist)

    # Search across all possible heights of the clustering and compares values inside difference matrix
    search = cluster_search(hClust = varClust,
                                    varMat = rankVarDiff,
                                    sizeMin = minSize)
    # Reduce the dataframe to unique rows
    if (base::is.vector(search$features)) {
      selectionPath = search$features
      compareValues = 0
      optModel = search$features
    }

    if (base::is.vector(search$features) == FALSE) {
      selectionPath = search$features %>%
        base::data.frame() %>%
        dplyr::distinct()

      # Find the best set of nested models using a selection criteria

      compareValues = apply(selectionPath == 1,
                            1,
                            inf_criterion,
                            x = x,
                            y = y,
                            family = "binomial",
                            lambda = lambda,
                            metric = compareMethod)

      optModel = selectionPath[which.min(compareValues),]

    }

    optFeatures = names(optModel)[unlist(optModel) == 1]

    # Output the following:
    # a) The coefficient values from each bootstrap
    # b) The variance of rank differences matrix
    # c) The hierachical clustering
    # d) The set of nested models
    # e) The optimum model selected
    output = base::list(
      coefficients = bootVal,
      varMat = rankVarDiff,
      varClust = varClust,
      selection = selectionPath,
      sizes = unique(search$size),
      compareMethod = compareMethod,
      compareValues = compareValues,
      optModel = optModel,
      optFeatures = optFeatures,
      vividSplit = FALSE
    )

    output

  }
