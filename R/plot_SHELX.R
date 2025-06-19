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
plot_SHELX <- function(filename, filename_e, var, type, title_chart) {
  ## Patch version for ggplot2 release (2025/06/) (J Foadi)

  ## Avoid NOTE in R CMD check
  Res <- CCall <- CCweak <- ncycle <- Contrast <- dataset <- NULL

  ## SHELXC plot
  if (type == "shelxc") {
    gg <- ggplot(filename, aes(1 / (Res)^2, var)) +
      geom_point() +
      geom_line() +
      xlab(expression(h^2 * ring(A)^-2)) +
      labs(title = title_chart) +
      scale_x_continuous(
        expand = c(0, 0),
        limits = c(0, max(1 / (filename$Res)^2) + min(1 / (filename$Res)^2))
      ) +
      scale_y_continuous(
        expand = c(0, 0),
        limits = c(0, max(var) + min(var))
      )

    # Better comparison logic
    if (identical(var, filename$d_sig)) {
      gp <- gg + ylab(expression(Delta * F / sig(Delta * F)))
    } else if (identical(var, filename$Chi_s)) {
      gp <- gg + ylab(expression(Chi^2))
    } else if (identical(var, filename$I_sig)) {
      gp <- gg + ylab(expression(I / sig(I)))
    } else if (identical(var, filename$Complete)) {
      gp <- gg + ylab(expression(Completeness))
    } else {
      gp <- gg + ylab(expression(CC[1/2]))
    }

    return(gp)

    ## SHELXD plot
  } else if (type == "shelxd") {
    gg <- ggplot(filename, aes(CCall, CCweak)) +
      geom_point() +
      scale_x_continuous(
        expand = c(0, 0),
        limits = c(0, max(filename$CCall) + min(filename$CCall))
      ) +
      scale_y_continuous(
        expand = c(0, 0),
        limits = c(0, max(filename$CCweak) + min(filename$CCweak))
      ) +
      xlab("CCall") +
      ylab("CCweak")

    return(gg)

    ## SHELXE plot
  } else if (type == "shelxe") {
    # Add cycle count
    lcycle_i <- nrow(filename$CYCLE)
    lcycle_o <- nrow(filename_e$CYCLE)

    filename$CYCLE$ncycle <- seq_len(lcycle_i)
    filename_e$CYCLE$ncycle <- seq_len(lcycle_o)

    df <- rbind(
      transform(filename_e$CYCLE, dataset = "Original"),
      transform(filename$CYCLE, dataset = "Inverted")
    )

    gg <- ggplot(df, aes(ncycle, Contrast, color = dataset, group = dataset)) +
      geom_line() +
      geom_point() +
      ggplot2::scale_color_manual(values = c("Original" = "#4DAC26", "Inverted" = "#2C7BB6"), name = "") +
      xlab("Cycle")

    return(gg)
  }

  stop("Unknown type: must be 'shelxc', 'shelxd' or 'shelxe'")
}
