
#' Matrix to Data Frame
#'
#' @param data A matrix
#'
#' @return
#' @export
#'
#' @examples

vivid_df = function(data) {
  # Changes data to dataframe structure
  df = base::as.data.frame(data)
  # Sets coloumn names
  base::colnames(df) = base::paste0("bootstrap_",
                                    base::seq.int(ncol(df)))
  df = base::cbind(id = rownames(df),
                   df)
  base::rownames(df) = NULL
  df

}
