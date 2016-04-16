library(ape)
library(geiger)


tree.ancestral <- function(treefile, phenofile){

# Input tree (make sure it will be binary)
blactree <- read.tree(treefile)
blactree <- multi2di(blactree)
blactree

# TEM taxa resistance phenotype 
blacdata <- read.table(phenofile, header=T)
# Vector of phenotypes assigned to TEM
pheno <- blacdata$pheno
names(pheno) <- blacdata$node

# remove tips with unknown phenotypes from tree
overlap <- name.check(blactree, pheno)
GAblactree <- drop.tip(blactree, overlap$tree_not_data)
# remove phenotype entries not in the phylogeny
pheno <- pheno[! names(pheno) %in% overlap$data_not_tree]
summary(pheno)

# ML reconstruction of ancestral phenotypes: All Rates Different
ace.ard <- ace(pheno, GAblactree, type="discrete", model="ARD")
ace.states.ard <- apply(ace.ard$lik.anc, 1, function(x) as.numeric(c(names(x)[which.max(x)])))
names(ace.states.ard) <- GAblactree$node.label

ace.hr.states <- ace.states.ard
ace.hr.states[ace.states.ard==1] <- "2br"
ace.hr.states[ace.states.ard==2] <- "2b"
ace.hr.states[ace.states.ard==3] <- "2ber"
ace.hr.states[ace.states.ard==4] <- "2be"
write.table(ace.hr.states,file=gsub("nwk_int.tree","internalStates.txt",treefile))
}

args <- commandArgs(trailingOnly = TRUE)
directory <- args[1]
treefiles <- list.files(path=directory, pattern=".nwk_int.tree")
"treefiles"
treefiles
phenofile <- args[2]
"phenofile"
phenofile

for (treefile in treefiles){
	"treefile"
	treefile
	tree.ancestral(treefile, phenofile)
}