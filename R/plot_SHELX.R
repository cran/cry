#
# This file is part of the cry package
#

# Functions connected to reflections data.
# Useful S3-type functions for handy use outside the main S4 framework


#' Plot SHELXC log files
#'
#' @param filename A data frame in output from \code{{read_SHELX_log}}.
#' @param filename_e A data frame with the inverted hand from shelxe
#' @param var the variable to be plotted vs the resolution
#' @param type indicate the type of file, possible value are "shelxc",
#'  "shelxd" and "shelxe".
#' @param title_chart title of the chart.
#' @return A graphical object from ggplot2 class that contains
#' the solution founded by SHELX log file.
#' @examples
#' datadir <- system.file("extdata",package="cry")
#' ## SHELXC
#' shelxc_log <- file.path(datadir,"shelxc.log")
#' shelxc <- read_SHELX_log(shelxc_log)
#' plot_shelxc <- plot_SHELX(filename = shelxc, var = shelxc$I_sig,
#' type = "shelxc", title_chart = "SHELXC")
#' plot_shelxc
#' ## SHELXD
#' shelxd_log <- file.path(datadir,"shelxd.log")
#' shelxd <- read_SHELX_log(shelxd_log)
#' plot_shelxd <- plot_SHELX(filename = shelxd, type = "shelxd",
#' title_chart = "SHELXD")
#' plot_shelxd
#' ## SHELXE
#' filename_i <- file.path(datadir,"shelxe_i.log")
#' shelxe_i <- read_SHELX_log(filename_i)
#' filename_o <- file.path(datadir,"shelxe_o.log")
#' shelxe_o <- read_SHELX_log(filename_o)
#' plot_shelxe <- plot_SHELX(filename = shelxe_i,
#' filename_e = shelxe_o, type = "shelxe", title_chart = "SHELXE")
#' plot_shelxe
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_y_continuous
#' @importFrom ggplot2 labs
#' @export

plot_SHELX <-function(filename, filename_e, var, type, title_chart)
{
  ## set to NULL global variables
  Res <- CCall <- CCweak <- ncycle <- Contrast <- dataset <- NULL
  ## Plot SHELXC data
  if(type == "shelxc") {
    gg <- ggplot(filename, aes(1/(Res)^2, var)) +
      geom_point() + geom_line() +
      xlab(expression(h^2 * (ring(A)^-2))) +
      labs(title = title_chart) +
      scale_x_continuous(expand = c(0, 0),
                         limits = c(0,
                                    max(1/(filename$Res)^2) + min(1/(filename$Res)^2))) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, max(var) + min(var)))
      #theme_bw()
    ifelse(var == filename$d_sig,
           gp <- gg +
             ylab(expression(Delta*F/sig(Delta*F))),
           ifelse(var == filename$Chi_s,
                  gp <- gg +
                    ylab(expression(Chi^2)),
                  ifelse(var == filename$I_sig,
                         gp <- gg +
                           ylab(expression(I/sig(I))),
                         ifelse(var == filename$Complete,
                                gp <- gg +
                                  ylab(expression(Completeness)),
                                gp <- gg +
                                  ylab(expression(CC_(1/2)))))))
    gp
  } else if(type == "shelxd") {
    gg <- ggplot(filename, aes(CCall, CCweak)) +
      geom_point() +
      scale_x_continuous(expand = c(0, 0),
                         limits = c(0,
                                    max(filename$CCall) + min(filename$CCall))) +
      scale_y_continuous(expand = c(0, 0),
                         limits = c(0, max(filename$CCweak) + min(filename$CCweak))) +
      #theme_bw() +
      xlab('CCall') +
      ylab('CCweak')
    gg
  } else {
    # create a new column which will be the total number of cycles
    lcycle_i <- length(filename$CYCLE[,1])
    lcycle_o <- length(filename_e$CYCLE[,1])
    ncycle_i <- seq(1,lcycle_i)
    ncycle_o <- seq(1,lcycle_o)
    CYCLE_i <- cbind(filename$CYCLE, ncycle_i)
    # Change the name of the new variable
    names(CYCLE_i)[names(CYCLE_i) == 'ncycle_i'] <- 'ncycle'
    CYCLE_o <- cbind(filename_e$CYCLE, ncycle_o)
    names(CYCLE_o)[names(CYCLE_o) == 'ncycle_o'] <- 'ncycle'

    # Since the DF can have different length, join the inverted and original DF
    # to create a new data frame. All the column must have the same name.
    df <- rbind(CYCLE_o, CYCLE_i)
    df$dataset <- c(rep("Original", nrow(CYCLE_o)),
                    rep("Inverted", nrow(CYCLE_i)))


    # Plot Contrast vs. Cycle
    gg <- ggplot(df, aes(ncycle, Contrast, color  =  dataset,
                         group = dataset)) +
      geom_line() + geom_point() +
      #theme_bw() +
      ggplot2::scale_color_manual(values = c("#4DAC26", "#2C7BB6"),
                         name = "") +
      xlab("Cycle")
    return(gg)
  }
}
