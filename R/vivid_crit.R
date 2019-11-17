
#' Change the information criteria used in VIVID
#'
#' @param vividObj An object passed from the vivid() funciton.
#' @param x input matrix, rows represent observations and columns represent
#' features. column names and row names are recommended.
#' @param y response variable, restricted to two classes.
#' @param metric AIC,  AICc, BIC or EBIC.
#' @param lambda Used in cv.glmnet
#'
#' @return
#' @export
#'
#' @examples

vivid_crit = function(vividObj, x, y, metric, lambda = 'lambda.1se') {
  lambda = "lambda.min"
  selectionPath = vividObj$selection

  compareValues = base::apply(
    X = selectionPath,
    MARGIN = 1,
    FUN = inf_criterion,
    x = x,
    y = y,
    lambda = lambda,
    metric = metric
  )
  optModel = selectionPath[which.min(compareValues), ]

  optFeatures = base::names(x = optModel)[base::unlist(x = optModel)]

  vividObj$compareMethod = metric
  vividObj$compareValues = compareValues
  vividObj$optModel = optModel
  vividObj$optFeatures = optFeatures

  vividObj
}
