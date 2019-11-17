
#' Search of heirachical clustering
#'
#' @param hClust Heirachical clustering.
#' @param varMat Variance matrix
#' @param sizeMin Minimum size of cluster
#'
#' @return
#' @export
#'
#' @examples

cluster_search = function(hClust, varMat, sizeMin) {
  p = length(x = hClust$order + 1) # Dimension of the distance matrix

  # Initializing variables
  hCutMat = base::matrix(data = NA_real_,
                         nrow = p - 2,
                         ncol = p)
  hCutVal = base::list()
  whichGroup = base::rep(x = NA_real_,
                         times = p - 2)
  whichFeatures = base::matrix(data = NA_real_,
                               nrow = p - 2,
                               ncol = p)
  minVal = base::rep(x = NA_real_,
                     times = p - 2)
  size = base::rep(x = NA_real_,
                   times = p - 2)
  checkNotNA = base::c()
  currentFeatures = base::rep(x = 1,
                              times = p)

  for (k in 2:(p - 1)) {
    # Cut the tree at every single height
    hCut = stats::cutree(tree = hClust,
                         k = k)
    # Matrix where each row contains the grouping at each tree height
    hCutMat[k - 1, ] = hCut

    #Calculate the average values
    avgValues = rep(x = NA,
                    times = k)

    for (i in 1:k) {
      # Identify which features belong to group i in cut k
      group = hCut * currentFeatures == i
      # Identify number of features in the
      count = sum(group)
      # Set defult value to 0
      avgVal = 0
      # Calculate the average value if the minimum number of features is met
      if (count > sizeMin) {
        varMatGroup = varMat[group, group]
        # Calculates the average value in the matrix excluding the diagonal
        avgVal = sum(varMatGroup) / (count * (count - 1))
      }

      avgValues[i] = avgVal
    }

    hCutVal[[k - 1]] = avgValues
    if (base::sum(avgValues) > 0) {
      # Find the smallest non-zero value
      minVal[k - 1] = base::min(avgValues[avgValues != 0])
      # Pick the group of features with the smallest average value
      whichGroup[k - 1] = base::which(avgValues == minVal[k - 1])
      whichFeatures[k - 1, ] = hCut == whichGroup[k - 1]
      currentFeatures = whichFeatures[k - 1, ] * 1
      size[k - 1] = base::sum(x = whichFeatures[k - 1, ])
      checkNotNA = base::c(checkNotNA, k - 1)
    }

  }

  whichFeatures = whichFeatures[checkNotNA, ]
  whichFeatures = matrix(whichFeatures,
                         ncol = p)
  base::colnames(whichFeatures) = hClust$labels
  minVal = minVal[checkNotNA]
  size = size[checkNotNA]

  # Output the list of features selected at each cut, the corresponding average value and the size of that group
  output = list(features = whichFeatures,
                value = minVal,
                size = size)

  output
}
