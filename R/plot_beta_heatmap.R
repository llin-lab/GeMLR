#' Plot the heatmap of the beta of the best cluster number
#'
#' @param beta_result the regression beta coefficient of each cluster
#' @param output_file If you just want to have a look, let it be NULL.
#' @param width default
#' @param height default
#' @param res default
#'
#' @return a heatmap
#' @export
#'
#' @examples plot_beta_heatmap(beta_result), plot_beta_heatmap(beta_result, output_file = "beta_heatmap.png")
plot_beta_heatmap <- function(beta_result,
                              output_file = NULL,
                              width = 5.5,
                              height = 5,
                              res = 300) {

  # Check and install required packages
  if (!require(ComplexHeatmap, quietly = TRUE)) install.packages("ComplexHeatmap")
  if (!require(circlize, quietly = TRUE)) install.packages("circlize")
  library(ComplexHeatmap)
  library(circlize)

  # Create color mapping function
  max_val <- max(abs(beta_result))
  col_fun <- colorRamp2(
    c(-max_val, 0, max_val),
    c("blue", "white", "red")
  )

  # Initialize graphics device if output file specified
  if (!is.null(output_file)) {
    png(output_file, width = width, height = height, units = "in", res = res)
    on.exit(dev.off())  # Ensure device is closed
  }

  # Create main heatmap with modified display settings
  ht <- Heatmap(
    beta_result,
    name = "Coefficient",  # Legend title
    col = col_fun,
    cluster_columns = FALSE,
    cluster_rows = FALSE,
    border = TRUE,

    # Column names settings (horizontal display)
    column_names_side = "top",
    column_names_centered = TRUE,
    column_names_rot = 0,  # No rotation
    column_names_gp = gpar(fontsize = 10),

    # Row names settings (only on right side)
    row_names_side = "right",
    row_names_gp = gpar(fontsize = 10),
    show_row_names = TRUE,

    # Titles
    row_title = "Variables",
    column_title = "Clusters",

    # Legend settings
    heatmap_legend_param = list(
      title_position = "leftcenter-rot",
      legend_height = unit(4, "cm")
    )
  )

  # Draw the heatmap
  draw(ht)

  # Return the heatmap object for further manipulation
  return(invisible(ht))
}
