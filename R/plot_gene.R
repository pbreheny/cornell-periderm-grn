#' Plot a gene (beeswarm) or two genes (scatter)
#'
#' @param Data     Expression data to base plot on
#' @param x,y      Name of genes to plot (character); x is mandatory, y is optional
#'
#' @examples
#' Data <- import_expr()
#' plot_gene(Data, 'irf6')
#' plot_gene(Data, 'irf6', 'gata3')
#' Data <- import_expr('farrell')
#' subset(Data, Tissue == 'EVL') |> plot_gene('irf6', 'gata3')
#' @export

plot_gene <- function(Data, x, y) {
  if (missing(y)) {
    xx <- Data$Exp[Gene == x][,-1] |> unlist
    DT <- data.table(xx=xx, Sample=names(xx))
    setnames(DT, 'xx', x)
    if (!is.null(Data$Tissue)) {
      Tissue <- Data$Tissue
      Tissue[, Sample := tolower(Sample)]
      Tissue[, Time := factor(Time)]
      DT <- merge(DT, Data$Tissue, by='Sample')
      p <- ggplot(DT, aes_string('Time', x, color='Time')) + geom_point() + ggbeeswarm::geom_beeswarm()
    } else {
      DT[, a := 1]
      p <- ggplot(DT, aes_string('a', x)) + geom_point() + ggbeeswarm::geom_beeswarm() + xlab('') +
        theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
    }
  } else {
    xx <- Data$Exp[Gene == x][,-1] |> unlist
    yy <- Data$Exp[Gene == y][,-1] |> unlist
    DT <- data.table(xx=xx, yy=yy, Sample=names(xx))
    setnames(DT, 'xx', x)
    setnames(DT, 'yy', y)
    if (!is.null(Data$Tissue)) {
      Tissue <- Data$Tissue
      Tissue[, Sample := tolower(Sample)]
      Tissue[, Time := factor(Time)]
      DT <- merge(DT, Data$Tissue, by='Sample')
      p <- ggplot(DT, aes_string(x, y, color='Time')) + geom_point()
    } else {
      p <- ggplot(DT, aes_string(x, y)) + geom_point()
    }
  }
  p
}
