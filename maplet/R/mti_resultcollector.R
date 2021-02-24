#' Class for combining pipelines
#'
#' @docType class
#'
#'
#' @field lst Combined list of all unique (as identified by uuid) metadata(D)$results entries
#' @field graph An igraph graph capturing the directional relationships between result entries
#'
#' @noRd

MTResultCollector <- R6::R6Class("MTResultCollector", list(
  lst = list(),
  graph = igraph::make_empty_graph(),

  # add single SE dataset
  add = function(D) {

    # add to graph, loop over results
    r <- metadata(D)$results
    for (i in 1:length(r)) {
      # add to list, if it doesn't exist yet
      if (!(r[[i]]$uuid %in% names(self$lst))) {
        self$lst[[r[[i]]$uuid]] = r[[i]]
      }
      # add node if it doesn't exist yet
      na <- self$graph %>% igraph::vertex_attr("name")
      if (!(r[[i]]$uuid %in% na)) {
        self$graph <- self$graph + igraph::vertex(r[[i]]$uuid, label=paste0(r[[i]]$fun,collapse="_"))
        print(sprintf('adding node %s',r[[i]]$uuid))
      }
      # add edge
      if (i>1) {
        self$graph <- self$graph + igraph::edge(r[[i-1]]$uuid,r[[i]]$uuid)
        print(sprintf('adding edge between %s to %s',r[[i-1]]$uuid,r[[i]]$uuid))
      }
    }
    #
    invisible(self)
  },

  # add multiple SE datasets
  addMultiple = function(...) {
    sapply(unlist(list(...)),self$add)
    invisible(self)
  },

  ### graph traversing helper functions

  # root node(s)
  graph_roots = function() {
    i <- which(sapply(sapply(igraph::V(self$graph), function(x) igraph::neighbors(self$graph,x, mode="in")), length) == 0)
    (self$graph %>% igraph::vertex_attr("name"))[i]
  },

  # get next node for a given vertex ID
  graph_next = function(id) {
    stopifnot(length(id)==1)
    ne <- igraph::neighbors(self$graph, id, 'out')
    unique(attr(ne,"names"))
  }

))

