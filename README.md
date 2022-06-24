# RUtreebalance

Robust, universal tree balance indices.

Publication: https://doi.org/10.1093/sysbio/syac027

The tree must be provided either as a `phylo` object or as a dataframe with column names Parent, Identity and (optionally) Population. The latter is similar to `tree$edge`, where `tree` is a `phylo` object; the differences are in class (dataframe versus matrix) and column names.

For example, here are four ways of inputting the same tree:

``` r
require(ape)
phylo_tree <- read.tree(text="((a:0.1)A:0.5,(b1:0.2,b2:0.1)B:0.2);")
J1_index(phylo_tree)

edges_tree <- data.frame(Parent = c(4,3,4,6,6), Identity = c(5,1,6,2,3))
J1_index(edges_tree)

edges_tree_with_root <- data.frame(Parent = c(4,4,3,4,6,6), Identity = c(4,5,1,6,2,3))
J1_index(edges_tree_with_root)

edges_tree_with_pops <- data.frame(Parent = c(4,3,4,6,6), Identity = c(5,1,6,2,3), Population = c(1,1,0,1,0))
J1_index(edges_tree_with_pops)
```
