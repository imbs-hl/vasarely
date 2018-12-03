plot_vasarely <- function(prob, prob_ex, color, upper_color_value, lower_color_value,
                          name_xaxis, name_yaxis, chi_statistic, chi_p_value){
  library("ggplot2")
  library("forcats")

  # create plot and turn arount y-axis, so the origin
  # of the plot is top left
  # changed real_probability to prob_real
  p <- ggplot(prob, aes(x = prob_ex$allel2,
                        y = forcats::fct_rev(prob_ex$allel1))) +

    # create plot with squares for the expected relative frequency,
    geom_raster(aes(fill = prob$expected_prob),
                hjust = 0.5,
                vjust = 0.5) +

    # create plot with dots for the observed relative frequency,
    # binwidth: regulate size of dots,
    # binaxis: direction to group dots
    # stackdir: dots in the center of squares#
    # color: regulate the lines around the dots
    geom_dotplot(aes(fill = prob$real_prob),
                 binwidth = 0.90,
                 binaxis = "y",
                 stackdir = 'center',
                 color = 0.001) +
    # take color and spread it from lower to upper color value
    scale_fill_gradientn(colours = color,
                         limits = c(lower_color_value,
                                    upper_color_value)) +

    # title of legend, title of y-axis
    labs(fill = "probability",
         x = name_xaxis,
         y = name_yaxis) +

    # put x-axis to the top of the plot
    scale_x_discrete(position = "top",
                     expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    # x-axis and y-axis should always have the same length
    # so the fields are squares
    coord_fixed(ratio = 1/1) +

    # fix size of legend text and size of legend title
    # legend title in bold
    theme(legend.text = element_text(size = 8),
          legend.title = element_text(size = 10,
                                      face = "bold")) +
    # put title of the legend to the top, so it is readable
    guides(fill = guide_legend(title.position = "top"))

    # add results of chi-squared-test
    if (chi_p_value > 0.00005){
      p <- p + labs(caption = paste0( "Chi-squared test: p-value ",
                                      sprintf("%.4f",
                                              x = chi_p_value),
                                      ", statistic ",
                                      sprintf("%.2f",
                                              x = chi_statistic)))
    } else {
      p <- p + labs(caption = paste0( "Chi-squared test: p-value ",
                                      format(x = chi_p_value,
                                            scientific = T,
                                            digits=4),
                                      ", statistic ",
                                      sprintf("%.2f",
                                              x = chi_statistic)))}
  return(p)
}
