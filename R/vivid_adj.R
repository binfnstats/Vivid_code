
#' Adjust variables in a VIVID object
#'
#' @param vividObj An object passed from the vivid() funciton.
#' @param minFinalFeatures Integer value for the minimum number of final features selected.
#'
#' @return
#' @export
#'
#' @examples

vivid_adj = function(vividObj, minFinalFeatures) {
  features = vividObj$selection
  sizes = vividObj$sizes
  compareValues = vividObj$compareValues

  possModels = base::which(sizes >= minFinalFeatures)

  if (length(possModels) == 0) {
    stop("Size on minimum features exceeds that of possible feature groups found.")
  }

  optModel = features[base::which.min(compareValues[possModels]), ]

  optFeatures = base::names(optModel)[base::unlist(optModel)]

  optFeatures
}
