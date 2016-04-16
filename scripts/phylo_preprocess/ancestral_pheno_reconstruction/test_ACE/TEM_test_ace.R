library(ape)
library(geiger)
# Input tree (make sure it will be binary)
blactree <- read.tree("run1.26560000.nwk_int.tree")

blactree <- multi2di(blactree)
# TEM taxa resistance phenotype 
blacdata <- read.table("TEM_pheno.txt", header=T)
blacdata
# Vector of phenotypes assigned to TEM
pheno <- blacdata$pheno
names(pheno) <- blacdata$node

# remove tips with unknown phenotypes from tree
overlap <- name.check(blactree, pheno)
overlap
GAblactree <- drop.tip(blactree, overlap$tree_not_data)
GAblactree <- makeNodeLabel(GAblactree)

# remove phenotype entries not in the phylogeny
pheno <- pheno[! names(pheno) %in% overlap$data_not_tree]
phenolabel <- character(length(GAblactree$tip.label))
phenolabel[pheno==1] <- "red"
phenolabel[pheno==2] <- "black"
phenolabel[pheno==3] <- "purple"
phenolabel[pheno==4] <- "blue"

names(phenolabel) <- names(pheno)
# ML reconstruction of ancestral phenotypes: Symmetric Rates
ace.sym <- ace(pheno, GAblactree, type="discrete", model="SYM")
ace.states.sym <- apply(ace.sym$lik.anc, 1, function(x) as.numeric(c(names(x)[which.max(x)])))
names(ace.states.sym) <- GAblactree$node.label

# ML reconstruction of ancestral phenotypes: All Rates Different
ace.ard <- ace(pheno, GAblactree, type="discrete", model="ARD")
ace.states.ard <- apply(ace.ard$lik.anc, 1, function(x) as.numeric(c(names(x)[which.max(x)])))
names(ace.states.ard) <- GAblactree$node.label

# Check which model is better
likDiff <- 2*(ace.ard$loglik - ace.sym$loglik)
pchisq(q=likDiff,df=6,lower.tail=F)
# ARD is significantly better

# Plot reconstruction
plot(GAblactree,label.offset=0.05, cex=0.4)
points(rep(0.25, length(GAblactree$tip.label)), 1:length(GAblactree$tip.label), pch=21, bg=phenolabel[GAblactree$tip.label])
nodelabels(pie=ace.ard$lik.anc, piecol=c('red','black','purple','blue'), cex=0.6)
