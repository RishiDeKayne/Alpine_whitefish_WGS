library(ape)

tree <- read.tree("../99_reanalysis/reviewer_responses/Tree_renamed_for_Dsuite2.nwk")

write.tree(tree)

species<-c("Outgroup1", "Outgroup2", "Outgroup3", "Outgroup4", "Outgroup5", "Outgroup6", "Outgroup7")

pruned.tree<-drop.tip(tree,tree$tip.label[match(species, tree$tip.label)])

write.tree(pruned.tree)

#############
tree <- read.tree("../99_reanalysis/RAxML/RAxML_bipartitions.99_500kb.GTRGAMMA.raxmlout")

write.tree(tree)

species<-c("0904", "0906", "0911", "0908", "0907", "0909", "0910",
           "032", "033", "064", "063", "061",
           "053", "040", "038", "0134", "039", "055",
           "093", "019", "056", "016", "058", "0999", "017", "059", "0136",
           "026", "029",
           "013", "012",
           "0123", "0131", "0121",
           "0103", "0104",
           "0110", "0111",
           "0114", "0119",
           "087", "090", "078", "076", "079", "095", "094", "096",
           "083", "081", "0101",
           "0100", "097", "0102",
           "041", "043", "02", "05", "03",
           "06", "074", "071", "010")

pruned.tree<-drop.tip(tree,tree$tip.label[match(species, tree$tip.label)])

write.tree(pruned.tree)

pruned.tree2<-pruned.tree
pruned.tree2$edge.length<-NULL
pruned.tree2$node.label <- NULL
write.tree(pruned.tree2)

