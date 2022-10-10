#' cry theme for ggplot2
#'
#' @name theme_cry
#' @examples
#' \dontrun{
#'  plot_SHELX(obj_shelxc, var = obj_shelxc$Chi_sq, type = "shelxc",
#'  title_chart = "Chis ^2") + theme_cry
#' }
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 margin
#' @importFrom ggplot2 element_blank
#' @importFrom ggplot2 element_line
#' @importFrom ggplot2 element_rect
#' @export


theme_cry <- function() {
  theme(axis.title.x = element_text(size = 20),
        axis.text.x  = element_text(size = 16),
        axis.title.y = element_text(size = 20),
        axis.text.y  = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12),
        plot.title = element_text(size = 24,
                                  face = "bold",
                                  margin = margin(r = 10),
                                  hjust = 0.5),
        text = element_text(color = "navy"),
        panel.background = element_blank(),
        panel.grid.major.x = element_line(color = "grey"),
        panel.grid.minor.x = element_line(color = "grey",
                                          linetype = "dashed"),
        panel.grid.major.y = element_line(color = "grey"),
        panel.grid.minor.y = element_line(color = "grey",
                                          linetype = "dashed"),
        #axis.ticks = element_blank(),
        strip.background = element_rect(fill = "white", color ="grey"),
        legend.position ="bottom",
        legend.direction ="horizontal")
}
