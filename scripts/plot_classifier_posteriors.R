plot_classifier_posteriors <- function(result, tree) {
  
  node_posteriors <- result$node_posteriors %>%
    dplyr::transmute(
      node = node_id,
      query_posterior = posterior
    )
  
  ggtree::ggtree(tree) %<+% node_posteriors +
    ggtree::geom_tree(
      ggplot2::aes(linewidth = query_posterior),
      lineend = "round"
    ) +
    ggtree::geom_tiplab(size = 3) +
    ggplot2::scale_linewidth_continuous(
      range = c(0.2, 5),
      limits = c(0, 1),
      name = "Posterior"
    ) +
    ggtree::theme_tree2()
}
