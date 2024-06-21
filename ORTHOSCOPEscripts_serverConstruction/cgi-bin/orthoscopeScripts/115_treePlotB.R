library(ape)
args              <- commandArgs()
infileName        <- args[6]      # 100_1stAnalysisSummary.txt
analysisGroupName <- args[7]
outfileName       <- args[8]      # 115_1st


################# Node name
nodeNameLabel_change_swich <- "on"
#nodeNameLabel_change_swich <- "off"
#########################

greenPrefixes <- c(); purplePrefixes <- c(); orangePrefixes <- c(); magentaPrefixes <-c(); bluePrefixes <- c()
redPrefixes <- c();


#######
line.picker <- function(keyWord)
{
  container <- c()
  frag <- 0
  lineStock <- c()
  for(line in infile$V1) {
    #print(line)

    ### Collect lines 
    if (regexpr('^>', line) < 0) {
      if(frag == 1) {
        lineStock <- c(lineStock, line)
      }
    }
  
   if (regexpr('^>', line) > 0)
   {
      if(frag == 1)
      {
        container <- c(container, lineStock)
        lineStock <- c()
        break
      }

      #keyWord1 <- paste('>',　keyWord,　sep='')
      if (regexpr(keyWord, line) > 0)
      {
        container <- c(container, line)
        frag <- 1
      }

    }
  }
  container <- c(container, lineStock)    
  container <- container[-1]
  return(container)
}

preab.sub <- function (keyWordA)
{
  #print(keyWordA)
  keyWordA <- paste(keyWordA, "$", sep = "")
  keyWordA <- paste(">", keyWordA, sep = "")
  #print(keyWordA)
  #print("### infile$V1 START ###")
  #print(infile$V1)
  #print("### infile$V1 END ###")
  if(any(i <- grep(keyWordA, infile$V1)))
  {
    #print("Found keyWordA")
    #print(keyWordA)
    #print("")
    containerA <- line.picker(keyWordA)
  } else {
    #print("Not Found keyWordA")
    #print(keyWordA)
    #print("")
    # print (paste('  ',keyWordA,' does not exist.', sep=''))
    containerA  <- NULL
  }
  return(containerA)
}


fontNumChange <- function (tr)
{
  tipFontNums              <- rep(1, length(tr$tip.label))
  #fontNum[queryIDNum]     <- 4
  return (tipFontNums)
}


queryNameInversion <- function (tr, queryNames)
{
  queryTipNums <- c()
  for(i in 1:length(tr$tip.label)){
    for (queryName in queryNames){
      if(regexpr(queryName, tr$tip.label[i]) > 0){
        #print(queryName)
        #print(tr$tip.label[i])
        #print("")
        queryTipNums <- c(queryTipNums, i)
      }
    }
  }
  return(queryTipNums)
}


tipColorChange <- function (tr)
{
  orangeNum <- c(); blueNum <- c(); redNum <- NULL; greenNum <- c(); magentaNum <- c(); purpleNum <- c(); humanNum <- c()

  for(i in 1:length(tr$tip.label)){

    for (redPrefix in redPrefixes){
      if(regexpr(redPrefix, tr$tip.label[i]) > 0){
        redNum <- c(redNum,i)
      }
    }

    for (greenPrefix in greenPrefixes){
      if(regexpr(greenPrefix, tr$tip.label[i]) > 0){
        greenNum <- c(greenNum,i)
      }
    }

    for (purplePrefix in purplePrefixes){
      if(regexpr(purplePrefix, tr$tip.label[i]) > 0){
        purpleNum <- c(purpleNum,i)
      }
    }
  
    for (orangePrefix in orangePrefixes){
      if(regexpr(orangePrefix, tr$tip.label[i]) > 0){
        orangeNum <- c(orangeNum,i)
      }
    }
  
    for (magentaPrefix in magentaPrefixes){
      if(regexpr(magentaPrefix, tr$tip.label[i]) > 0){
        magentaNum <- c(magentaNum,i)
      }
    }  
  
    for (bluePrefix in bluePrefixes){
      #print(bluePrefix)
      if(regexpr(bluePrefix, tr$tip.label[i]) > 0){
        blueNum <- c(blueNum,i)
      }
    }
  
  }

  tipColorNums             <- rep("black",length(tr$tip.label))
  tipColorNums[redNum]     <- "red"
  tipColorNums[greenNum]   <- "darkgreen"
  tipColorNums[purpleNum]  <- "purple"
  tipColorNums[orangeNum]  <- "darkorange1"
  tipColorNums[magentaNum] <- "hotpink2"
  tipColorNums[blueNum]    <- "blue"

  return(tipColorNums)
}

make_colorPrefixes <- function (taxonSampling_colorTMP, colorFN)
{
  colorPrefixes <- c()
  for (line in taxonSampling_colorTMP)
  {
    #spPrefixFN <- sub(' +.*$', "", line)   
    spPrefixFN <- sub('_.*$', "", line)
    if(regexpr(colorFN, line) > 0){
        colorPrefixes <- c(colorPrefixes, spPrefixFN)
    }
  }
  return(colorPrefixes)
}


nodeNameLabel_change <- function (tr)
{
  for(p in 1:length(tr$node.label)){
    #print(tr$node.label[p])
    if(regexpr('D=N', tr$node.label[p])> 0){
      tr$node.label[p] <- sub('_.*$', "", tr$node.label[p])
    } else {
      tr$node.label[p] <- sub('_.*$', "D", tr$node.label[p])
    }
    #print(tr$node.label[p])
  }
  return(tr)
}


BScolorChange <- function (tr)
{
  #print("Rearrangement_BS_value_threshold")
  #print(Rearrangement_BS_value_threshold)

  BSvalueColors <- NULL
  for(p in 1:length(tr$node.label)){
    if (regexpr('r', tr$node.label[p]) > 0){
      BSvalueColors <- c(BSvalueColors, 2)
    #} else if (as.numeric(tr$node.label[p]) < as.numeric(Rearrangement_BS_value_threshold)){
    #  BSvalueColors <- c(BSvalueColors, 2)
    } else {
      BSvalueColors <- c(BSvalueColors, 1)    
    }
  }

#  if (is.null(Rearrangement_BS_value_threshold)){
#  } else {
#    for(p in 1:length(tr$node.label)){
#      if (tr$node.label[p] == "r"){
#      } else {   
#        #print(tr$node.label[p])
#        #print(Rearrangement_BS_value_threshold)
#        #print("\n")
#        #print("tr$node.label[p]")
#        #print(tr$node.label[p])
#        if (as.numeric(tr$node.label[p]) < as.integer(Rearrangement_BS_value_threshold)){
#          BSvalueColors <- c(BSvalueColors, 2) 
#        }
#      }
#    }
#  }


  return(BSvalueColors)
}


edgeWidthChange <- function (tr, OrthoLeaves)
{
  orthoBranchNums <- c()
  for(i in 1:length(tr$tip.label)){
    for(orthoLeaf in OrthoLeaves){
      if (tr$tip.label[i] == orthoLeaf){
        orthoBranchNums <- c(orthoBranchNums, i)
      }
    }
  }

  edgeWidth <- NULL
  if(is.null(orthoBranchNums)){
    edgeWidth          <- rep(1.5, dim(tr$edge)[1])
  } else {
    edgeWidth          <- rep(1.5, dim(tr$edge)[1])
    whOrtho            <- which.edge(tr, orthoBranchNums)
    edgeWidth[whOrtho] <- 4
  }
  return(edgeWidth)
}

PNG_treeDrawing <- function (tr, prefix)
{
  png.file <- paste(outfileName, prefix, sep = "")
  pngWidth <- NULL

  pngHeight <- NULL
  if (length(tr$tip.label) > 200) {
    pngWidth <- 1500
    pngHeight = 2700
  } else if (length(tr$tip.label) > 100) {
    pngWidth <- 1200
    pngHeight = 1800
  } else if (length(tr$tip.label) > 50) {
    pngWidth <- 1000
    pngHeight = 1200
  } else if (length(tr$tip.label) > 10) {
    pngWidth <- 1000
    pngHeight = 900
  } else {
    pngWidth <- 800
    pngHeight = 600
  }

  png(png.file, width = pngWidth, height = pngHeight)
  plot(tr, no.margin=TRUE, underscore = TRUE, use.edge.length=TRUE, cex = 0.9, font = tipFontNums, tip.col = tipColorNums, edge.width = edgeWidth)
  add.scale.bar()
  #####nodelabels(tr$node.label, adj = c(1.2,-0.5), frame = "n", font = nodeLabelFontNums, cex=nodeLabelFontSizeNums, col = nodeLabelFontColorNums)
  nodelabels(tr$node.label, adj = c(1.2,-0.5), frame = "n", col = nodeLabelFontColorNums)
  if(!is.null(Num_allQueries)){
    tiplabels (tr$tip.label[Num_allQueries], Num_allQueries, cex=1.0, adj = 0, bg = "gray40", col="white")
  }
  if(!is.null(Num_1stQuery)){
    tiplabels (tr$tip.label[Num_1stQuery],   Num_1stQuery,   cex=1.0, adj = 0, bg = "navyblue", col="white")
  }
  dev.off()
}

PDF_treeDrawing <- function (tr, prefix)
{
  pdf.file <- paste(outfileName, prefix, sep = "")
  pdfWidth  <- NULL
  pdfHeight <- NULL
  if (length(tr$tip.label) > 200) {
    pdfWidth  = 38
    pdfHeight = 28
  } else if (length(tr$tip.label) > 100) {
    pdfWidth  = 31
    pdfHeight = 21
  } else if (length(tr$tip.label) > 50) {
    pdfWidth  = 24
    pdfHeight = 14
  } else if (length(tr$tip.label) > 10) {
    pdfWidth  = 25
    pdfHeight = 10
  } else {
    pdfWidth  = 15
    pdfHeight = 7
  }

  pdf(pdf.file, width = pdfWidth, height = pdfHeight)
  plot      (tr, no.margin=TRUE, underscore = TRUE, use.edge.length=TRUE, cex = 0.9, font = tipFontNums, tip.col = tipColorNums, edge.width = edgeWidth)
  add.scale.bar()
  nodelabels(tr$node.label, adj = c(1.2,-0.5), frame = "n", col = nodeLabelFontColorNums)
  if(!is.null(Num_allQueries)){
    tiplabels (tr$tip.label[Num_allQueries], Num_allQueries, cex=1.0, adj = 0, bg = "gray40", col="white")
  }
  if(!is.null(Num_1stQuery)){
    tiplabels (tr$tip.label[Num_1stQuery],   Num_1stQuery, cex=1.0,   adj = 0, bg = "navyblue", col="white")
  }

  dev.off()
}

##################################################################

infile    <- read.table(infileName, na.strings = FALSE, sep = '\t')
Querys_used_in_the_analysis <- preab.sub("Queries_used_in_the_analysis")

taxonSampling_color_color <- preab.sub("taxonSampling_color")
greenPrefixes   = make_colorPrefixes(taxonSampling_color_color, "Green")
purplePrefixes  = make_colorPrefixes(taxonSampling_color_color, "Purple")
orangePrefixes  = make_colorPrefixes(taxonSampling_color_color, "Orange")
magentaPrefixes = make_colorPrefixes(taxonSampling_color_color, "Magenta")
bluePrefixes    = make_colorPrefixes(taxonSampling_color_color, "Blue")
redPrefixes     = make_colorPrefixes(taxonSampling_color_color, "Red")

queryNames <- c()
for (line in Querys_used_in_the_analysis)
{
  line <- sub(' +.*$', "", line)   
  queryNames <- c(queryNames, line)
}
Orthogroup             <- preab.sub("Orthogroup")

Rearrangement_BS_value_threshold <- c()
Rearrangement_BS_value_threshold <- preab.sub("Rearrangement_BS_value_threshold")


##################################
Gene_tree <- preab.sub("Gene_tree_newick")
Gene_tree <- read.tree(text = Gene_tree)
Gene_tree <- ladderize(Gene_tree, TRUE)
Gene_tree$edge.length[Gene_tree$edge.length<0]<-0   ### nagative branch length, replace with 0

edgeWidth                    <- edgeWidthChange(Gene_tree, Orthogroup)

tipFontNums                  <- fontNumChange(Gene_tree)
tipColorNums                 <- tipColorChange(Gene_tree)
#nodeLabelFontNums           <- rep(1,length(Gene_tree$tip.label))
#nodeLabelFontSizeNums       <- rep(0.9, length(Gene_tree$tip.label))
nodeLabelFontColorNums       <- BScolorChange(Gene_tree)
Num_allQueries               <- queryNameInversion(Gene_tree, queryNames)
Num_1stQuery                 <- queryNameInversion(Gene_tree, queryNames[1])

PNG_treeDrawing(Gene_tree, prefix="GeneTree.png")
PDF_treeDrawing(Gene_tree, prefix="GeneTree.pdf")


##################################
Rearranged_gene_tree <- preab.sub("Rearranged_gene_tree_newick")
if(is.null(Rearranged_gene_tree)) {
  q()
}

Rearranged_gene_tree <- read.tree(text = Rearranged_gene_tree)
Rearranged_gene_tree <- ladderize(Rearranged_gene_tree, TRUE)


edgeWidth              <- edgeWidthChange(Rearranged_gene_tree, Orthogroup)

tipFontNums            <- fontNumChange(Rearranged_gene_tree)
tipColorNums           <- tipColorChange(Rearranged_gene_tree)

nodeLabelFontColorNums <- c()
if (nodeNameLabel_change_swich == "on")
{
  Rearranged_gene_tree     <- nodeNameLabel_change(Rearranged_gene_tree)
  nodeLabelFontColorNums   <- BScolorChange(Rearranged_gene_tree)
}

Num_allQueries             <- queryNameInversion(Rearranged_gene_tree, queryNames)
Num_1stQuery               <- queryNameInversion(Rearranged_gene_tree, queryNames[1])

txt.file <- paste(outfileName, "Rearranged_geneTree.txt", sep = "")
write.tree(Rearranged_gene_tree, file=txt.file)

PNG_treeDrawing(Rearranged_gene_tree, prefix="Rearranged_geneTree.png")
PDF_treeDrawing(Rearranged_gene_tree, prefix="Rearranged_geneTree.pdf")
