
#' A funciton that creates graphics for the VIVID method
#'
#' @param vividObj An object passed from the vivid() funciton.
#' @param log A TRUE or FALSE variable that defines whether the values plotted are log-transformed.
#'
#' @return
#' @export
#'
#' @examples

vivid_plot = function(vividObj, log = TRUE) {
  vividSplit  = vividObj$vividSplit
  if (vividSplit == TRUE) {
    vividObj = vividObj[[length(vividObj)]]
  }

  vividObj$varClust = dendsort::dendsort(d = vividObj$varClust)

  vividObj$varClust$labels = stringr::str_trunc(
    string = colnames(vividObj$varClust$labels),
    width = 12,
    side = "right"
  )

  ddRow = stats::as.dendrogram(vividObj$varClust)
  rowOrd = stats::order.dendrogram(ddRow)
  ddCol = ddRow
  colOrd = rowOrd

  colnames(vividObj$varMat) = stringr::str_trunc(
    string = colnames(vividObj$varMat),
    width = 12,
    side = "right"
  )
  rownames(vividObj$varMat) = stringr::str_trunc(
    string = rownames(vividObj$varMat),
    width = 12,
    side = "right"
  )

  cc = grDevices::colorRampPalette(c("green", "beige", "red"))
  lattice::trellis.par.set(regions = list(col = cc(100)))


  base::switch(
    log,
    "TRUE" = lattice::levelplot(
      log(vividObj$varMat[rowOrd, colOrd]),
      col.regions = cc,
      scales = list(cex = 0.5, x = list(rot =
                                          90)),
      aspect = "fill",
      colorkey = list(space = "left"),
      legend = list(right = list(
        fun = latticeExtra::dendrogramGrob,
        args = list(
          x = ddCol,
          ord = colOrd,
          side = "right",
          size = 6
        )
      ))
    ),
    "FALSE" = lattice::levelplot(
      vividObj$varMat[rowOrd, colOrd],
      col.regions = cc,
      scales = list(cex = 0.5, x = list(rot =
                                          90)),
      aspect = "fill",
      colorkey = list(space = "left"),
      legend = list(right = list(
        fun = latticeExtra::dendrogramGrob,
        args = list(
          x = ddCol,
          ord = colOrd,
          side = "right",
          size = 6
        )
      ))
    )
  )

}
