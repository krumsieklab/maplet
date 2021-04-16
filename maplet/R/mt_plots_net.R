#' Creates a network object
#'
#' Creates a network visualization object using the visNetwork package.
#'
#' @param D \code{SummarizedExperiment} input.
#' @param stat_name Name of the test to take correlations from.
#' @param cor_filter Filter for correlation values to plot. Default: p.value < 0.05.
#' @param node_coloring Name of the test to use for node coloring.
#' @param html_outfile Name of visnetwork html file. If empty, no html saved.
#' @param height OPTIONAL. Size in pixel of the plotting window size. Default: 500.
#'
#' @return $result$output: network ggplot
#' @return $result$output2: visnetwork plot
#'
#' @examples
#' \dontrun{#' # in the context of a SE pipeline
#' ... %>% mt_plots_net(stat_name = "xxx") %>% ...    # standard call
#' ... %>% mt_plots_net(stat_name = "xxx", cor_filter = p.adj < 0.5, node_coloring="Li's", html_outfile="Network.html", height=800) %>% ...    # filters only significant correlations and colors the nodes according to the results in the indicated test, saves visnetwork to file
#'}
#'
#' @author EB
#'
#' @importFrom ggnetwork theme_blank geom_nodetext geom_nodes geom_edges
#' @import ggplot2
#' @import network
#'
#' @export
mt_plots_net <- function(D,
                        stat_name,
                        cor_filter = p.value < 0.05,
                        node_coloring,
                        html_outfile,
                        height = 500) {

  ## check input
  stopifnot("SummarizedExperiment" %in% class(D))
  if(missing(stat_name))
    stop("stat_name must be given to plot the network")

  ## rowData
  rd1 <- subset(rowData(D), select=which(names(rowData(D))=="name")) %>%
        as.data.frame() %>%
        dplyr::mutate(var1 = rownames(D))
  colnames(rd1)[colnames(rd1)=="name"] <- "name1"
  rd1$name1 %<>% make.names()
  rd2 <- subset(rowData(D), select=which(names(rowData(D))=="name")) %>%
    as.data.frame() %>%
    dplyr::mutate(var2 = rownames(D))
  colnames(rd2)[colnames(rd2)=="name"] <- "name2"
  rd2$name2 %<>% make.names()

  ## stat
  data_plot <- mtm_get_stat_by_name(D, stat_name) %>%
    dplyr::inner_join(rd1, by = "var1") %>%
    dplyr::inner_join(rd2, by = "var2")

  # define node attributes
  nodes <- as.data.frame(unique(rbind(cbind(ids=data_plot$var1,label=data_plot$name1),cbind(ids=data_plot$var2,label=data_plot$name2))))

  if(!(missing(node_coloring))) {
    test <- mtm_get_stat_by_name(D, node_coloring);
    test <- test[match(nodes$ids, test$var),];
    nodes$map <- sign(test$statistic)*log10(test$p.value);
    nodes$node_color <- grDevices::colorRampPalette(c('blue', 'white', 'red'))(length(nodes$map))[rank(nodes$map)]
  } else nodes$node_color <- rep("lightblue",times=length(nodes$ids))

  ## apply filter on correlations
  if(!missing(cor_filter)){
    mti_logstatus("filter correlations")
    cor_filter_q <- dplyr::enquo(cor_filter)
    data_plot <- data_plot %>%
      dplyr::filter(!!cor_filter_q)
  }

  ## define edge attributes
  # rescale correlation to [0 10]
  cor.scaled <- (abs(data_plot[,which(colnames(data_plot)=="statistic")])-min(abs(data_plot[,which(colnames(data_plot)=="statistic")])))/(max(abs(data_plot[,which(colnames(data_plot)=="statistic")]))-min(abs(data_plot[,which(colnames(data_plot)=="statistic")])))*10

  edges <- data.frame(from = data_plot[,which(colnames(data_plot)=="name1")], to = data_plot[,which(colnames(data_plot)=="name2")], value = cor.scaled)
  edges$from <- as.character(edges$from)
  edges$to <- as.character(edges$to)
  edges$color <- "black"
  edges$color[edges$value<0] <- "red"

  ## plot
  e <- edges
  n <- data.frame(id=nodes$label, label= nodes$label, color=nodes$node_color)
  p_vis <- visNetwork::visNetwork(n,e, height = height, width = "100%") %>%
    # disable physics
    visNetwork::visIgraphLayout(layout = "layout_nicely", physics = F, smooth = F)

  # create ggnetwork object
  df <- list()
  df$edges <- data.frame(from = data_plot[,which(colnames(data_plot)=="var1")], to = data_plot[,which(colnames(data_plot)=="var2")], weight = cor.scaled/10)
  df$vertices <- nodes$label

  adj <- as.matrix(dils::AdjacencyFromEdgelist(df$edges, check.full = TRUE))
  colnames(adj[[1]])<- as.character(nodes$label[match(adj[[2]],nodes$ids)])
  rownames(adj[[1]])<- as.character(nodes$label[match(adj[[2]],nodes$ids)])
  mm.net <- network::network(adj[[1]], layout = "kamadakawai", directed = FALSE)

  if(!missing(node_coloring)){
    test <- mtm_get_stat_by_name(D, node_coloring)
    test <- test[match(adj[[2]], test$var),]
    map <- sign(test$statistic)*log10(test$p.value)
    mm.net %v% "strength" <- map
  } else {
    mm.net %v% "strength" <- 1
  }

  mm.col <- c("positive" = "#000000", "negative" = "#0000FF")
  x <- data_plot$statistic
  x[x>=0] <- "positive"
  x[x<0] <- "negative"
  mm.net %e% "pcor" <- abs(data_plot[,which(colnames(data_plot)=="statistic")])
  mm.net %e% "pos" <- x

  # add edge color (positive/negative)
  # add black circle around nodes
  p <- ggplot(mm.net, aes(x, y, xend = xend, yend = yend)) +
    geom_edges(color="black",aes(size=pcor)) +
    geom_nodes(aes(color = strength), size = 7) +
    scale_color_gradient2(low = "#0000FF", mid="white",high = "#FF0000") +
    geom_nodetext(color="grey50",aes(label = vertex.names),
                  size = 3, vjust = -0.6) +
    theme_blank() +
    theme(legend.position = "bottom")


  # if html_outfile given, save visnetwork to html
  if (!missing(html_outfile)) {
    # due to odd visSave path handling behavior, we need to export to a tmp file first and the move to final location
    tmpfile = sprintf("tmp_%s.html", uuid::UUIDgenerate())
    visNetwork::visSave(graph = p_vis, file = tmpfile)
    file.rename(tmpfile, html_outfile)
  }

  ## add status information & plot
  funargs <- mti_funargs()
  D %<>% 
    mti_generate_result(
      funargs = funargs,
      logtxt = sprintf("Correlation Network, aes: %s", stat_name),
      output = list(p),
      output2 = p_vis
    )
  ## return
  D
}
