# Find the parent of a node.
# 
# Example:
# edges1 <- data.frame(Parent = c(1,1,1,3,3), Identity = 2:6)
# move_up(edges1, 3)
move_up <- function(edges, identity) {
  if(!(identity %in% edges$Identity) & !(identity %in% edges$Parent)) stop("Invalid identity.")
  parent <- edges[which(edges$Identity == identity), "Parent"]
  if(length(parent) == 0) return(identity) # if identity is the root then don't move
  if(is.factor(parent)) parent <- levels(parent)[parent]
  return(parent)
}

# Find the root node of a tree.
# 
# Example:
# edges1 <- data.frame(Parent = c(1,1,1,3,3), Identity = 2:6)
# find_root(edges1)
find_root <- function(edges) {
  start <- edges$Parent[1] # reasonable guess
  if(is.factor(start)) start <- levels(start)[start]
  repeat {
    if(move_up(edges, start) == start) break
    start <- move_up(edges, start)
  }
  return(start)
}

# Get a list of all subtree sizes via depth-first search.
# Optionally provide root node i and adjacency list Adj if known.
# 
# Example:
# tree1 <- data.frame(Parent = c(1,1,1,1,2,3,4), 
#                     Identity = 1:7, 
#                     Population = c(1, rep(5, 6)))
# get_subtree_sizes(tree1)
get_subtree_sizes <- function(tree,i=NULL,Adj=NULL,Col=NULL,Cumul=NULL,is_leaf=NULL){
  tree <- add_root_row(tree)
  n<-length(tree$Identity)
  has_pops <- FALSE
  if("Population" %in% colnames(tree)) has_pops <- TRUE
  if(is.null(Adj)) Adj <- get_Adj(tree)
  if(is.null(i)) i <- which(tree$Identity == find_root(tree[,1:2]))
  if(is.null(Col)) {
    Col <- rep("w",n)
    names(Col) <- unique(tree$Identity)
  }
  if(is.null(Cumul)) {
    Cumul <- rep(NA,n)
    names(Cumul) <- unique(tree$Identity)
  }
  if(is.null(is_leaf)) {
    is_leaf <- rep(FALSE, n)
    names(is_leaf) <- unique(tree$Identity)
  }
  if(is.null(Adj[[i]])) is_leaf[i] <- TRUE
  for (j in Adj[[i]]){
    if (Col[j] == "w"){
      L <- get_subtree_sizes(tree,j,Adj,Col,Cumul,is_leaf)
      Col<- L$colour
      Cumul <- L$cumulative
      is_leaf <- L$is_leaf
    }
  }
  Col[i] <- "b"
  if(has_pops) {
    Cumul[i] <- tree$Population[i] + sum(Cumul[Adj[[i]]])
  } else {
    Cumul[i] <- ifelse(is_leaf[i] == TRUE, 1, 0) + sum(Cumul[Adj[[i]]])
  }
  return(list("colour"=Col,"cumulative"=Cumul,"is_leaf"=is_leaf))
}

# Get adjacency list of a tree.
# 
# Example:
# tree1 <- data.frame(Parent = c(1,1,1,1,2,3,4), 
#                     Identity = 1:7, 
#                     Population = c(1, rep(5, 6)))
# get_Adj(tree1)
get_Adj <- function(tree) {
  n<-length(tree$Identity)
  Adj <- vector(mode = "list", length = n)
  for (i in 1:n) if(tree$Parent[i] != tree$Identity[i]) {
    p <- which(tree$Identity == tree$Parent[i])
    Adj[[p]] <- append(Adj[[p]], i)
  }
  return(Adj)
}

# Add a row to the edges list to represent the root node (if not already present)
add_root_row <- function(tree) {
  start <- setdiff(tree$Parent, tree$Identity)
  if(length(start) > 0) { # add row for root node
    if("Population" %in% colnames(tree)) {
      root_row <- data.frame(Parent = start, Identity = start, Population = 0)
      message("Assigning Population = 0 to the root node")
    }
    else root_row <- data.frame(Parent = start, Identity = start)
    tree <- rbind(root_row, tree)
  }
  return(tree)
}

# Calculate tree balance index J^1 (when nonrootdomfactor = FALSE) or
# J^{1c} (when nonrootdomfactor = TRUE).
# If population sizes are missing then the function assigns
# size 0 to internal nodes, and size 1 to leaves.
# 
# Examples using phylo object as input:
# require(ape)
# phylo_tree <- read.tree(text="((a:0.1)A:0.5,(b1:0.2,b2:0.1)B:0.2);")
# J1_index(phylo_tree)
# phylo_tree2 <- read.tree(text='((A, B), ((C, D), (E, F)));')
# J1_index(phylo_tree2)
# 
# Examples using edges lists as input:
# tree1 <- data.frame(Parent = c(1,1,1,1,2,3,4),
#                     Identity = 1:7,
#                     Population = c(1, rep(5, 6)))
# J1_index(tree1)
# tree2 <- data.frame(Parent = c(1,1,1,1,2,3,4),
#                     Identity = 1:7,
#                     Population = c(rep(0, 4), rep(1, 3)))
# J1_index(tree2)
# tree3 <- data.frame(Parent = c(1,1,1,1,2,3,4),
#                     Identity = 1:7,
#                     Population = c(0, rep(1, 3), rep(0, 3)))
# J1_index(tree3)
# cat_tree <- data.frame(Parent = c(1, 1:14, 1:15, 15),
#                        Identity = 1:31,
#                        Population = c(rep(0, 15), rep(1, 16)))
# J1_index(cat_tree)
#
# If population sizes are omitted then internal nodes are assigned population size zero
# and leaves are assigned population size one:
# sym_tree1 <- data.frame(Parent = c(1, rep(1:15, each = 2)),
#                        Identity = 1:31,
#                        Population = c(rep(0, 15), rep(1, 16)))
# sym_tree2 <- data.frame(Parent = c(1, rep(1:15, each = 2)),
#                        Identity = 1:31)
# all.equal(J1_index(sym_tree1), J1_index(sym_tree1))
J1_index <- function(tree, q = 1, nonrootdomfactor = FALSE) {
  
  if(!is.na(tree)[1]) {
    if("phylo" %in% class(tree)) { # convert from phylo object to appropriate data frame
      tree <- tree$edge
      tree <- as.data.frame(tree)
      colnames(tree) <- c("Parent", "Identity")
    }
    tree <- na.omit(tree) # remove any rows containing NA
    if(is.factor(tree$Parent)) tree$Parent <- levels(tree$Parent)[tree$Parent]
    if(is.factor(tree$Identity)) tree$Identity <- levels(tree$Identity)[tree$Identity]
    tree <- add_root_row(tree)
  }
  
  n<-length(tree$Identity)
  if (n<=1) return(0)
  Adj <- get_Adj(tree) # adjacency list
  subtree_sizes <- get_subtree_sizes(tree, Adj = Adj) # get the list of all subtree sizes
  Cumul <- subtree_sizes$cumulative # subtree sizes, including the root
  eff_int_nodes <- which(!subtree_sizes$is_leaf) # vector of internal nodes
  leaves <- which(subtree_sizes$is_leaf) # vector of leaves
  # if population sizes are missing then assign size 0 to internal nodes, and size 1 to leaves:
  if(!("Population" %in% colnames(tree))) {
    tree$Population <- rep(0, n)
    tree$Population[leaves] <- 1
    message("Assigning Population = 0 to internal nodes and Population = 1 to leaves")
  }
  if(sum(tree$Population) <= 0) stop("At least one node must have Population > 0")
  J <- 0
  Star <- Cumul - tree$Population # subtree sizes, excluding the root
  for (i in 1:n){ # loop over all nodes
    if (Star[i] > 0){ # if node has at least one child with non-zero size
      K <- 0
      if(length(Adj[[i]])>1){ # otherwise i has only one child and its balance score is 0
        eff_children <- 0 # number of children with non-zero size
        for (j in Adj[[i]]){
          if (Cumul[j]>0){ # otherwise child j has a 0-sized subtree and does not count
            eff_children <- eff_children+1
            # p is the ratio of the child subtree size including the root (root = the child) 
            # to the parent subtree size excluding the root
            p <- Cumul[j]/Star[i]
            # K is the sum of the node balance scores
            if(q == 1) {
              K <- K + -p*log(p)
            } else {
              K <- K + p^q
            }
          }
        }
        # non-root dominance factor:
        if(nonrootdomfactor) {
          h_factor <- Star[i] / Cumul[i]
        } else {
          h_factor <- 1
        }
        # normalize the sum of balance scores, adjust for non-root dominance, 
        # and then add the result to the index
        if(eff_children > 1) { # exclude nodes that have only one child with size greater than zero
          if(q == 1) {
            J <- J + h_factor * Star[i] * K / log(eff_children)
          } else {
            J <- J + h_factor * Star[i] * (1 - K) * eff_children^(q - 1) / (eff_children^(q - 1) - 1)
          }
        }
      }
    }
  }
  # normalize the index by dividing by the sum of all subtree sizes:
  if (length(eff_int_nodes)>0) J <- J/sum(Star[eff_int_nodes])
  return(as.numeric(J))
}
