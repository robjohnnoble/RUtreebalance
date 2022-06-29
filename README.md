# RUtreebalance

Robust, universal tree balance indices.

Publication: https://doi.org/10.1093/sysbio/syac027

The tree must be provided either as a `phylo` object or as a dataframe with column names Parent, Identity and (optionally) Population. The latter is similar to `tree$edge`, where `tree` is a `phylo` object; the differences are in class (dataframe versus matrix) and column names.

The dataframe may (but isn't required to) include a row for the root node.

If population sizes are omitted then internal nodes will be assigned population size zero and leaves will be assigned population size one.

For example, here are four ways of inputting the same tree with internal nodes labelled 4, 5, 6 and leaves 1, 2, 3:

``` r
# phylo object read from Newick format:
require(ape)
phylo_tree <- read.tree(text="((1)5,(2,3)6);")
J1_index(phylo_tree)

# dataframe omitting population sizes:
edges_tree <- data.frame(Parent = c(4,5,4,6,6), Identity = c(5,1,6,2,3))
J1_index(edges_tree)

# dataframe omitting population sizes and including a row for the root:
edges_tree_with_root <- data.frame(Parent = c(4,4,5,4,6,6), Identity = c(4,5,1,6,2,3))
J1_index(edges_tree_with_root)

# dataframe including population sizes:
edges_tree_with_pops <- data.frame(Parent = c(4,5,4,6,6), Identity = c(5,1,6,2,3), Population = c(0,1,0,1,1))
J1_index(edges_tree_with_pops)
```
