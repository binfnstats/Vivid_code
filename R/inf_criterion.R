
#' Information Criterion
#'
#' @param model Regression model.
#' @param x input matrix, rows represent observations and columns represent
#' features. column names and row names are recommended.
#' @param y response variable, restricted to two classes.
#' @param family GLM family.
#' @param lambda Used in cv.glmnet.
#' @param metric AIC,  AICc, BIC or EBIC.
#'
#' @return
#' @export
#'
#' @examples

inf_criterion = function(model,
                         x,
                         y,
                         family = "binomial",
                         lambda = 'lambda.1se',
                         metric = "EBIC") {
  n = length(y)
  P = NCOL(x)

  df = sum(model)

  xNew = x[, model]

  CVGlmFit = glmnet::cv.glmnet(
    x = xNew,
    y = y,
    standardize = TRUE,
    alpha = 0,
    family = family
  )

  GlmFit = glmnet::glmnet(
    x = xNew,
    y = y,
    standardize = TRUE,
    alpha = 0,
    family = family,
    lambda = CVGlmFit[[lambda]]
  )

  deviance = (1L - GlmFit$"dev.ratio") * GlmFit$"nulldev"

  base::switch(
    metric,
    "AIC" = deviance + 2 * df,
    "AIC" = deviance + 2 * df + 2 * (df ^ 2 + df) / (n - df - 1),
    "BIC" = deviance + log(n) * df,
    "EBIC" = deviance + log(n) * df + 2 * gamma * lchoose(n = P, k = df),
    NA
  )

}
