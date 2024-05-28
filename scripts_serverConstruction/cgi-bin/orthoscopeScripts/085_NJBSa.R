library("ape")
#library("ape", lib.loc="/usr/share/R/library")

args            <- commandArgs()
infile          <- args[6]
outGroup        <- args[7]
selectedModel   <- args[8]
outfile.nwk     <- args[9]

#print ('  ## 085_NJBSa.R starts.')
#print(.libPaths())
#q()


#print('infile')
#print(infile)
#print("")
#print('outGroup')
#print(outGroup)
#print("")
#print('selectedModel')
#print(selectedModel)
#print("")
#print('outfile.nwk')
#print(outfile.nwk)
#print("")

#  print (paste('  infile', infile))

#outfile.nwk <- paste(outfile, '.nwk', sep = '')

## Open the 'infile' file
if(file.access(infile) != 0){
  print (paste(infile,' does not exist.'))
} else {
  infile.phy  <- read.dna(infile)
}

## Estimate distance
dist.selectedModel <- dist.dna(infile.phy, model = selectedModel, pairwise.deletion=TRUE, gamma=5)

## #Estimate NJ tree
nj.selectedModel <- njs(dist.selectedModel)

## Reroot by outGroup
nj.selectedModel <- root(nj.selectedModel,outGroup,r=T)

## ladderize from bottom
nj.selectedModel <- ladderize(nj.selectedModel,TRUE)

## BS analysis with partitioning data (1st+2nd)
#nj.boot.nj.selectedModel <- boot.phylo(nj.selectedModel, infile.phy,     function(xx) root(nj(dist.dna(xx, model = "selectedModel", pairwise.deletion = TRUE, gamma = 5)),outGroup,r=T), 100, 2)
#nj.boot.nj.selectedModel <- boot.phylo(nj.selectedModel, infile.phy,     function(xx) root(njs(dist.dna(xx, model = "selectedModel", pairwise.deletion = TRUE)), outGroup,r=T), 100, 2)
nj.boot.nj.selectedModel <- boot.phylo(nj.selectedModel, infile.phy,     function(xx) root(njs(dist.dna(xx, model = selectedModel, pairwise.deletion = TRUE)), outGroup,r=T), 100)
nj.selectedModel$node.label <- nj.boot.nj.selectedModel

## Change bs calues: NA -> 0
#for(p in 1:length(nj.selectedModel$node.label)){
#  print(p)
#  print(nj.selectedModel$node.label[p])
#}
#cat('\n\n')

nj.selectedModel$node.label[is.na(nj.selectedModel$node.label)]<-0

#for(p in 1:length(nj.selectedModel$node.label)){
#  print(p)
#  print(nj.selectedModel$node.label[p])
#}
#cat('\n')

## Write tree
write.tree(nj.selectedModel, file = outfile.nwk)

