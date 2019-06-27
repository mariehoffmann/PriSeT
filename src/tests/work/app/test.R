library(d3treeR)
library(DT)
library(shiny)
library(treemap)


# 1  - 2 - 21 - [211, 212]
#        - 22 - [221, 222, 223, 224]
#        - 23 - [231, 232]
#        - 24 - 241
#   \- 3 - 31 - [311, 312, 313]
#        - 32 - 321
#        - 33 - [331, 332]
taxid = c(2,3,21,22,23,24,31,32,33,211,212,221,222,223,224,231,232,241,311,312,313,321,331,332)
p_taxid = c(rep(1,2),rep(2,4),rep(3,3),rep(21,2),rep(22,4),rep(23,2), 24,rep(31,3),32,33,33)
tax = data.frame(taxid, p_taxid)

root = tax$p_taxid[1]
desc = tax$taxid[tax$p_taxid == root]

# augment data frame by subtree size (number of taxid in clade including root = tax$taxid)
ctrs <- table(p_taxid)
group_ctrs <- rep(1, length(taxid))
for (i in 1:length(taxid)){
    ctr = ctrs[names(ctrs) == taxid[i]]
    if (length(ctr) > 0){
            group_ctrs[i] <- ctr
    }
}

tax$clade_size = group_ctrs

#subset = data.frame(tax[tax$p_taxid == desc])
clade <- subset(tax, tax$p_taxid %in% desc, select=c(p_taxid, taxid, clade_size))

# build treemap
tree = treemap(clade, index = c("p_taxid", "taxid"), vSize = "clade_size", type = "index")

# make it interactive ("rootname" becomes the title of the plot):
tree_inter = d3tree3(tree,  rootname = paste("root", root))
