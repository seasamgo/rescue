
#' Nearest Network
#'
#' Construct a nearest neighbour network based on previously computed PCs
#'
#' @param reduced_object PC reduction matrix
#' @param k_neighbors Number of k neighbors to use
#' @param minimum_shared Minimum shared neighbors
#' @param top_shared Keep at ...
#' @param verbose Be verbose
#' @param \dots Additional parameters
#'
#' @importFrom data.table := .N
#'
#' @return NN network as igraph object
#'

constructNN <- function(
  reduced_object,
  k_neighbors = 30,
  minimum_shared = 5,
  top_shared = 3,
  verbose = F,
  ...
  ) {


  # data.table variables
  from_cell_ID = from = to_cell_ID = to = shared = weight = distance = NULL

  # vector for cell_ID
  cell_names = rownames(reduced_object)
  names(cell_names) = 1:nrow(reduced_object)

  ## run nearest-neighbour algorithm ##
  if(k_neighbors >= nrow(reduced_object)) {
    k_neighbors = (nrow(reduced_object)-1)
    if(verbose == TRUE) cat('\n k is too high, adjusted to nrow(matrix)-1 \n')
  }

  nn_network = dbscan::kNN(x = reduced_object, k = k_neighbors, sort = TRUE, ...)
  nn_network_dt = data.table::data.table(
    from = rep(1:nrow(nn_network$id), k_neighbors),
    to = as.vector(nn_network$id),
    weight = 1/(1 + as.vector(nn_network$dist)),
    distance = as.vector(nn_network$dist)
    )

  nn_network_dt[, from_cell_ID := cell_names[from]]
  nn_network_dt[, to_cell_ID := cell_names[to]]

  snn_network = dbscan::sNN(x = nn_network, k = k_neighbors, kt = NULL, ...)
  snn_network_dt = data.table::data.table(
    from = rep(1:nrow(snn_network$id), k_neighbors),
    to = as.vector(snn_network$id),
    weight = 1/(1 + as.vector(snn_network$dist)),
    distance = as.vector(snn_network$dist),
    shared = as.vector(snn_network$shared)
    )

  snn_network_dt = snn_network_dt[stats::complete.cases(snn_network_dt)]
  snn_network_dt[, from_cell_ID := cell_names[from]]
  snn_network_dt[, to_cell_ID := cell_names[to]]

  # rank snn
  data.table::setorder(snn_network_dt, from, -shared)
  snn_network_dt[, rank := 1:.N, by = from]

  # filter snn
  snn_network_dt = snn_network_dt[rank <= top_shared | shared >= minimum_shared]

  ## convert to igraph object
  all_index = unique(x = c(nn_network_dt$from_cell_ID, nn_network_dt$to_cell_ID))

  missing_indices = all_index[!all_index %in% unique(snn_network_dt$from)]
  nn_network_igraph = igraph::graph_from_data_frame(snn_network_dt[, list(from_cell_ID, to_cell_ID, weight, distance, shared, rank)], directed = TRUE, vertices = all_index)

  return(nn_network_igraph)

}


#' Cluster Cells via Louvain Algorithm
#'
#' Cluster cells using a NN-network and the Louvain algorithm from the community module in Python
#'
#' @param nn_network Constructed nearest neighbor network to use
#' @param python_path Specify specific path to python if required
#' @param resolution Resolution
#' @param weight_col Weight column
#' @param louv_random Random
#' @param set_seed Set seed
#' @param seed_number Number for seed
#' @param ... Additional parameters
#'
#' @importFrom data.table :=
#'
#' @return A character vector of cluster labels
#'

clusterLouvain <- function(
  nn_network,
  python_path = NULL,
  resolution = 1,
  weight_col = NULL,
  louv_random = F,
  set_seed = T,
  seed_number = 0,
  ...
  ) {


  ## check or make paths
  # python path
  if(is.null(python_path)) {

    # this will initialize python for reticulate and suggest to install miniconda if not present
    reticulate::py_config()

    if(.Platform[['OS.type']] == 'unix') {
      python_path = try(system('which python', intern = T))
    } else if(.Platform[['OS.type']] == 'windows') {
      python_path = try(system('where python', intern = T))
      if(class(python_path) == "try-error") {
        cat('\n no python path found, set it manually when needed \n')
        python_path = '/set/your/python/path/manually/please/'
      }
    }
  }
  python_path = as.character(python_path)

  # prepare python path and louvain script
  reticulate::use_python(required = T, python = python_path)
  python_louvain_function = system.file("python", "python_louvain.py", package = 'rescue')
  reticulate::source_python(file = python_louvain_function)

  # start seed
  if(set_seed == TRUE) {
    seed_number = as.integer(seed_number)
  } else {
    seed_number = as.integer(sample(x = 1:10000, size = 1))
  }

  network_edge_dt = data.table::as.data.table(igraph::as_data_frame(x = nn_network, what = 'edges'))

  # data.table variables
  weight = NULL

  if(!is.null(weight_col)) {

    if(!weight_col %in% colnames(network_edge_dt)) {
      stop('\n weight column is not an igraph attribute \n')
    } else {
      # weight is defined by attribute of igraph object
      network_edge_dt = network_edge_dt[,c('from', 'to', weight_col), with = F]
      data.table::setnames(network_edge_dt, weight_col, 'weight')
    }
  } else {
    # weight is the same
    network_edge_dt = network_edge_dt[,c('from', 'to'), with = F]
    network_edge_dt[, weight := 1]
  }

  # do python louvain clustering
  if(louv_random == FALSE) {
    reticulate::py_set_seed(seed = seed_number, disable_hash_randomization = TRUE)
    pyth_louv_result = python_louvain(df = network_edge_dt, resolution = resolution, randomize = F)
  } else {
    reticulate::py_set_seed(seed = seed_number, disable_hash_randomization = TRUE)
    pyth_louv_result = python_louvain(df = network_edge_dt, resolution = resolution, random_state = seed_number)
  }
  ident_clusters_DT = data.table::data.table(cell_ID = rownames(pyth_louv_result), 'name' = pyth_louv_result[[1]])
  data.table::setnames(ident_clusters_DT, 'name', 'cluster')

  # return clustering result
  return(ident_clusters_DT)

}

