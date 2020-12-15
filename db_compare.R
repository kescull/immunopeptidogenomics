library(VennDiagram)
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(reshape2)
library(optparse)

sessionInfo()
options(warn=1,hrbrthemes.loadfonts=TRUE)

filterThreshold <- function(csv, threshold) {
  newCsv <- csv %>% filter(X.10lgP >= threshold) %>% droplevels() 
}
filterLength <- function(csv) {
  newCsv <- csv %>% filter(Length < 16, Length > 6) %>% droplevels()
}
getOneRowPerPeptide <- function(data) {
  dat <- data[ ,c("Peptide","Length")] %>% group_by(Peptide) %>%
            filter(row_number()==1) %>% 
            droplevels()
}
makeCategoryVector <- function(data) {
  x <- names(vennPartitions)
  for (i in 1:length(x)) {
    j <- data$Peptide %in% vennPartitions[[i]]
    data[ ,x[i]] <- j
  }
  data <- melt(data, id.vars = c("Peptide","Length"))
  data <- data[data$value==T, ] %>% select(-value)
}

getPeptides <- function(dat) {
  dat[ ,"Peptide"]
}
writePartitions <- function(x) {
  if (x[1] == TRUE && x[2] == TRUE) {
   textFile <- paste0(c(prefix,names(x[1]),names(x[2]),"overlap.txt"), collapse="_")
  }
  else if (x[1] == TRUE) {
    textFile <- paste0(c(prefix,names(x[1]),"only.txt"),collapse="_")
  }
  else {
    textFile <- paste0(c(prefix,names(x[2]),"only.txt"),collapse="_")
  }
  writeLines(x[[4]],textFile,sep = "\n")
  x[[4]]
}
makePartitionNames <- function(x) {
  if (x[1] == TRUE && x[2] == TRUE) {
     part <- "overlap"
  }
  else if (x[1] == TRUE) {
    part <- paste0(c(names(x[1]),"only"),collapse="_")
  }
  else {
    part <- paste0(c(names(x[2]),"only"),collapse="_")
  }
  part
}
assessScanSet <- function(dat, y) {
  x <- nrow(dat)
  #make new dataframe, including Scan factor vector - this needs to be exactly the same in newDat anyway
  newDat <- data.frame(Score1=numeric(x), Score2=numeric(x), N1=integer(x), N2=integer(x),Same = logical(x), Multiples = logical(x))
  j <- 0
  if (x == 1) {
    j <- j + 1
    newDat$Multiples[j] <- F 
    if (!any(is.na(dat[1,]))) {
      newDat$Same[j] <- T
    } else {
      newDat$Same[j] <- F
    } 
    newDat$Score1[j] <- ifelse(!is.na(dat[[scoreNames[1]]]), dat[[scoreNames[1]]], 0)
    newDat$Score2[j] <- ifelse(!is.na(dat[[scoreNames[2]]]), dat[[scoreNames[2]]], 0)
    newDat$N1[j] <- ifelse(!is.na(dat[[scoreNames[1]]]), 1, 0)
    newDat$N2[j] <- ifelse(!is.na(dat[[scoreNames[2]]]), 1, 0)
  } else { #if more than one row for a scan
    #make key to label pairs so can order datatable appropriately
    keyedDat <- dat[ ,c(scoreNames[1],scoreNames[2])]
    keys <- integer(x)
    n1 <- 0L
    n2 <- 0L
    cnt1 <- 2L
    cnt2 <- 2L
    for (i in 1:x) { # want keys for pairs to match. 
      if (!is.na(keyedDat[[scoreNames[1]]][i]) && !is.na(keyedDat[[scoreNames[2]]][i])) {
        keys[i] <- 1
        n1 <- n1 + 1
        n2 <- n2 + 1
      } else if (!is.na(keyedDat[[scoreNames[1]]][i])) {
        keys[i] <- cnt1
        cnt1 <- cnt1 + 1
        n1 <- n1 + 1
      } else {
        keys[i] <- cnt2
        cnt2 <- cnt2 + 1
        n2 <- n2 + 1
      }
    }
    #add keys to datatable and sort by them
    keyedDat$Key <- keys
    sorted <- keyedDat[order(keyedDat$Key),]
    i <- 1
    while (i < x + 1) {
      j <- j + 1
      newDat$Multiples[j] <- T
      newDat$N1[j] <- n1
      newDat$N2[j] <- n2
      if (sorted$Key[i] == 1) {
        newDat$Same[j] <- T
        newDat$Score1[j] <- sorted[[scoreNames[1]]][i]
        newDat$Score2[j] <- sorted[[scoreNames[2]]][i]
      } else {
        newDat$Same[j] <- F
        if (i < x && sorted$Key[i] == sorted$Key[i+1]) {
          if (!is.na(sorted[[scoreNames[1]]][i])) {
            newDat$Score1[j] <- sorted[[scoreNames[1]]][i]
            newDat$Score2[j] <- sorted[[scoreNames[2]]][i + 1]
          } else {
            newDat$Score1[j] <- sorted[[scoreNames[1]]][i + 1]
            newDat$Score2[j] <- sorted[[scoreNames[2]]][i]
          }
          i <- i + 1
        } else {
          newDat$Score1[j] <- ifelse(!is.na(sorted[[scoreNames[1]]][i]), sorted[[scoreNames[1]]][i], 0)
          newDat$Score2[j] <- ifelse(!is.na(sorted[[scoreNames[2]]][i]), sorted[[scoreNames[2]]][i], 0)
        }
      }
      i <- i + 1
    }
  }
  newDat[1:j,]
}
getMaxRow <- function(dat, y) {
  dat[which.max(dat$X.10lgP), ]
}
getMaxScore <- function(dat) {
  dat <- dat %>% group_by(Peptide) %>%
    group_modify(getMaxRow) %>% droplevels()
}
getPeptidesAndScores <- function(dat,ids) {
  dat[which(dat$Peptide %in% ids),c(1,2)] %>% droplevels()
}
getOnlyScores <- function(dat, datName) {
  index <- grep(datName, names(vennPartitions))
  newVector <- dat[which(dat$Peptide %in% vennPartitions[[index]]),2]
}
getOnlyLens <- function(dat, datName) {
  dat <- dat[ ,"Length"]
}
renameScores <- function(dat, n) {
  newDat <- dat
  names(newDat) <- c(colnames(dat)[1],n)
  newDat
}
get_inset2 <- function(df){
  p <- ggplot(df) +
    geom_col(aes(x=as.factor(Var1), y=Freq, fill=Var3), position = position_dodge2(preserve = "single"), show.legend = FALSE) +
    theme_ipsum(base_size = 12, axis_title_size = 12) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill="white"), 
          panel.border = element_rect(fill = NA, colour="gray"))
  return(p)
}
annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) {
  layer(data = data, stat = StatIdentity, position = PositionIdentity, 
        geom = ggplot2:::GeomCustomAnn,
        inherit.aes = TRUE, params = list(grob = grob, 
                                          xmin = xmin, xmax = xmax, 
                                          ymin = ymin, ymax = ymax))
}
find.iRT.scans <- function(dat) {
  dat <- dat[grep("*iRT-pep*", dat$Accession), "Scan"] %>%
    droplevels()
}
get.iRT.rows <- function(dat) {
  dat <- subset(dat, (dat$Scan %in% all.iRT.scans)) %>% 
    droplevels()
}
get.noniRT.rows <- function(dat) {
  dat <- subset(dat, !(dat$Scan %in% all.iRT.scans)) %>% 
    droplevels()
}
# # # start 'main' script here
#get data

option_list <- list(
  make_option(c("-c", "--cryptic"), type="character", default=NULL, 
              help="cryptic PEAKS results filename", metavar="character"),
  make_option(c("-n", "--normal"), type="character", default=NULL, 
              help="normal (e.g. uniprot) PEAKS results filename", metavar="character"),
  make_option(c("-j", "--junction"), type="character", default=NULL, 
              help="Optional: artificial junction peptide list for discard (txt)", metavar="character"),
  make_option(c("-p", "--prefix"), type="character", default="out",
              help="prefix for all output files (e.g. myCellsRep1)", metavar="character"),
  make_option(c("-d", "--cryptic_threshold"), type="double", default=0, 
              help="threshold score for cryptic search results", metavar="number"),
  make_option(c("-m", "--norm_threshold"), type="double", default=0,
              help="threshold score for normal (uniprot) search results", metavar="number")
); 

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

if (is.null(opt$cryptic) || is.null(opt$normal) || opt$cryptic_threshold == 0 || opt$norm_threshold == 0){
  print_help(opt_parser)
  stop("Options -c, -d, -n, -m are required (input files and thresholds)", call.=FALSE)
}

filenames <- c(opt$cryptic,opt$normal)
thresholds <- c(opt$cryptic_threshold,opt$norm_threshold)
print("Files:")
print(filenames)
print("Thresholds:")
print(thresholds)
prefix <- opt$prefix
print(paste("prefix:", prefix))

discard <- NULL
if (!is.null(opt$junction)) {
  print("discarding peptides in file:")
  print(opt$junction)
  prefix <- paste0(c(prefix,"DISCARD"), collapse = "_")
  print(paste("prefix now:", prefix))
  discard <- as.factor(scan(opt$junction, what = character()))
  str(discard)
}
easyNames <- (c("cryptic","standard"))
allDat <- lapply(filenames, read.csv)
names(allDat) <- easyNames

#Next, get rid of unwanted columns - 
#Mass, charge, ppm, m.z, rt, area, fraction, id, from.chimera, Source.File, PTM, AScore, Found.By.
#i.e. keep Peptide, length, X.10lgP, Scan, Accession
paredDown <- lapply(allDat, function(dat) dat[ ,c("Peptide","X.10lgP","Length","Scan","Accession")])

#if discard option triggered, filter cryptic results to discard junction peptides
if (!is.null(discard)) {
  paredDown$cryptic <- subset(paredDown$cryptic, !(paredDown$cryptic$Peptide %in% discard)) %>% droplevels()
}

print("Dealing with any iRT peptides...")
#get scan numbers for iRT peptides and write to file, then delete rows from both sets
iRT.scans <- lapply(paredDown, find.iRT.scans)

all.iRT.scans <- as.factor(c(levels(iRT.scans$cryptic),levels(iRT.scans$standard)))
iRT.rows <- lapply(paredDown, get.iRT.rows)
merged.iRT <- merge(iRT.rows$cryptic,iRT.rows$standard, by="Scan",all = T,suffixes=easyNames)
write.csv(merged.iRT, file = paste0(c(prefix,"iRT.csv"), collapse = "_"), na="")

no.iRT <- lapply(paredDown, get.noniRT.rows)

#filter for those over threshold
filteredThresh <- mapply(filterThreshold, no.iRT, thresholds, SIMPLIFY = FALSE)

print("Preparing length comparison plots...")
#prepare to plot all lengths of good peptides
lendfs <- lapply(filteredThresh, getOneRowPerPeptide)
onlyLens <- mapply(getOnlyLens, lendfs)
names(onlyLens) <- easyNames
stackedLens <- stack(onlyLens)

#plot all lengths of good peptides
lenNumFile <- paste0(c(prefix,"idperlen.png"), collapse = "_")
png(lenNumFile, width = 760)
lenNumPlot <- ggplot(stackedLens) +
  geom_bar(aes(x=as.factor(values), fill=ind), position = position_dodge2(preserve = "single")) +
  theme_ipsum(base_size = 12, axis_title_size = 12) +
  labs(title="Total peptides per length") +
  theme(legend.title = element_blank()) +
  xlab("Ids per length") +
  ylab("Peptide count")
print(lenNumPlot)
invisible(dev.off)

#filter for length:
filtered <- lapply(filteredThresh, filterLength)

#prepare to plot lengths of good 7-15mers
lendfs2 <- lapply(filtered, getOneRowPerPeptide)
onlyLens2 <- mapply(getOnlyLens, lendfs2)
names(onlyLens2) <- easyNames
stackedLens2 <- stack(onlyLens2)

#plot lengths of good 7-15mers
lenNumFile2 <- paste0(c(prefix,"idperlen7-15.png"), collapse = "_")
png(lenNumFile2, width = 760)
lenNumPlot2 <- ggplot(stackedLens2) +
  geom_bar(aes(x=as.factor(values), fill=ind), position = position_dodge2(preserve = "single")) +
  theme_ipsum(base_size = 12, axis_title_size = 12) +
  labs(title="Total peptides per length") +
  theme(legend.title = element_blank()) +
  xlab("Ids per length") +
  ylab("Peptide count")
print(lenNumPlot2)
invisible(dev.off)

print("Preparing Venn diagram and partition lists...")
#venn diagram
peptideLists <- lapply(filtered, getPeptides)
vennPeptideLists <- lapply(peptideLists, levels)
vennFile <- paste0(c(prefix,"Venn.png"),collapse="_")

#remove the following command to get log file for Venn
invisible(futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger"))

invisible(venn.diagram(x = vennPeptideLists, 
	category.names = names(vennPeptideLists),
	filename = vennFile,
	imagetype="png",
	cat.prompts=TRUE,
  margin=0.1,
	cat.dist=rep(0.05,length(vennPeptideLists)),
	output = TRUE))
partitions <- get.venn.partitions(vennPeptideLists)

#for each row in 'partitions', write a txt file named appropriately and containing the list of peptides found in that row
#also get named character vector of each 
vennPartitions <- apply(partitions, 1, writePartitions)
names(vennPartitions) <- apply(partitions,1,makePartitionNames)

print("Comparing length comparisons for the partitions (csv output)")
#length comparison for cryptic, standard and overlap partitions
#First merge all lists
mergedLists <- merge(lendfs2[[1]], lendfs2[[2]], by=c("Peptide", "Length"), all=T)

#second add category vector indicating partition peptide belongs to and write out as csv
mergedLists <- makeCategoryVector(mergedLists)
write.csv(mergedLists, file = paste0(c(prefix,"lengths.csv"), collapse = "_"))

print("Preparing score comparison graphs...")
  #score comparison for overlap peptides - using 'filtered' which has peptide reps but only confident ids
#first get max score for each peptide id i.e. get subset dataframes with only max score for each Peptide level
maxScores <- lapply(filtered, getMaxScore)

#second make new dataframe with just peptide and score columns from each set for peptides %IN% overlap
overlapScores <- lapply(maxScores, getPeptidesAndScores, vennPartitions$overlap)

#third rename score column in each dataframe to reflect dataset, then merge dataframes by peptide id
overlapRenamed <- mapply(renameScores, overlapScores, names(overlapScores), SIMPLIFY=FALSE)

#fourth merge dfs by Peptide ie. get max score per overlap peptide in each search
combinedOverlap <- merge(overlapRenamed[[1]],overlapRenamed[[2]], by="Peptide")

#fifth plot it
xname <- names(combinedOverlap)[2]
yname <- names(combinedOverlap)[3]
#get min x/y from thresholds, set as origin!
origin <- min(thresholds)
biggest <- max(c(combinedOverlap[[2]],combinedOverlap[[3]]))

overlapFile <- paste0(c(prefix,"overlapScores.png"),collapse="_")
png(overlapFile)
overlapPlot <- ggplot(combinedOverlap, aes_string(x=xname, y=yname)) + 
  geom_point() + geom_abline(intercept=0,slope=1, linetype = "dotdash", color="red") + 
  coord_fixed(xlim = c(origin, biggest), ylim = c(origin, biggest)) +
  theme_ipsum(base_size = 12, axis_title_size = 12) +
  geom_smooth(method="lm", se=T, fullrange=T, level=0.95,color="purple", fill="grey",show.legend = T) +
  labs(title="Peptides found with both databases") +
  xlab(paste("score with",xname,"database", sep=" ")) +
  ylab(paste("score with",yname,"database", sep=" ")) +
  theme(plot.title = element_text(hjust = 0.5), aspect.ratio=1) +
  annotate(geom="text", label="reference (identical)",
           x=biggest-((biggest - origin)/3)-1, 
           y=(biggest-origin)-5,
           angle=45,hjust=0, 
           family="Arial Narrow", 
           color="red") 
  
print(overlapPlot)
invisible(dev.off())

#score comparison for non-overlap peptides - start with maxScores which has 1 rep of each peptide from both sets
# 1: copy score columns from each set for peptides %IN% appropriate Venn partition to new vector list; 
# name vectors in list by partition;
# convert to dataframe with stack

onlyScores <- mapply(getOnlyScores, maxScores, names(maxScores))
stackedOnlyScores <- stack(onlyScores)
# 2: plot it
onlyScoresFile <- paste0(c(prefix,"onlyScores.png"),collapse="_")
png(onlyScoresFile, width = 760, height = 480)
onlyPlot <- stackedOnlyScores %>% ggplot(aes(x=values, fill=ind)) +
  geom_histogram(color="#e9ecef", alpha=0.75, position = 'identity', binwidth=0.5) +
  scale_fill_manual(values=c("#ff8000", "#4dd2ff"),labels=easyNames) +
  theme_ipsum(base_size = 12, axis_title_size = 12) +
  labs(fill="", title="Peptides found with only one database", x="score", y="number of peptides") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = c(0.8, 0.7))
print(onlyPlot)
invisible(dev.off())

# go back to unfiltered scans for scan comparisons - only if discard option not triggered
if (is.null(discard)) {
  print("Preparing unfiltered scan comparisons...")
  merged <- merge(paredDown[[1]], paredDown[[2]], by=c("Scan", "Peptide"), all=T, suffixes=paste0(".",names(paredDown)))

  scoreNames <- paste0("X.10lgP.", easyNames)
  nNames <- paste(easyNames, "N", sep="_")
  confNames <- paste(easyNames, "conf", sep="_")
  print("Assessing multiples...")
  multiples <- merged %>% 
            group_by(Scan) %>% 
            summarise(number = length(Scan), 
                      conf = any(X.10lgP.cryptic > thresholds[1] || X.10lgP.standard > thresholds[2], 
                      na.rm = T))

  idNumFile <- paste0(c(prefix,"idsperscan.png"), collapse = "_")
  png(idNumFile)
  idNumPlot <- ggplot(multiples) +
    geom_bar(aes(x=as.factor(number), fill=conf), position = position_dodge2(preserve = "single")) +
    theme_ipsum(base_size = 12, axis_title_size = 12) +
    scale_y_continuous(sec.axis = sec_axis(~ . / nrow(multiples), labels=scales::percent, name = "% total scans")) +
    labs(title="Total identifications per scan") +
    xlab("ids per scan") +
    ylab("scan count")
  print(idNumPlot)
  invisible(dev.off)

  print("Doing scan comparison...")
  scanComparison <- merged %>%
                group_by(Scan) %>%
                group_modify(assessScanSet) %>%
                droplevels()
  names(scanComparison)[2:3] <- scoreNames
  names(scanComparison)[4:5] <- nNames
  allScoreFile <- paste0(c(prefix,"all.png"), collapse = "_")
  png(allScoreFile)
  allScorePlot <- ggplot(scanComparison, aes_string(x=scoreNames[1], y=scoreNames[2])) +
    coord_fixed() +
    geom_vline(xintercept = thresholds[1], linetype = "dotdash") +
    geom_hline(yintercept = thresholds[2], linetype = "dotdash") +
    geom_point(aes(alpha = 0.5, colour = factor(Same,labels = c("id different in different databases", "id same in each database")))) +
    geom_point(data = scanComparison[scanComparison$Multiples,],shape = 1, aes(colour = "multiple ids per scan", fill="white")) +
    scale_colour_manual(name= "", values = c("multiple ids per scan" = "black", "id same in each database" = "#E69F00", "id different in different databases" = "#56B4E9")) +
    guides(alpha = F, fill = F) +
    theme_ipsum(base_size = 12, axis_title_size = 12) +
    labs(title="Score comparisons for each scan") +
    xlab(easyNames[1]) +
    ylab(easyNames[2]) +
    annotate(geom="text", family="Arial Narrow", label=paste("FDR5", thresholds[2]),x=biggest, y=thresholds[2]+4,hjust=1) +
    annotate(geom="text", family="Arial Narrow", label=paste("FDR5", thresholds[1]),x=thresholds[1]+2, y=biggest,hjust=0) +
    theme(plot.title = element_text(hjust = 0.5)) 
 
  print(allScorePlot)
  invisible(dev.off)

#get dataframe with just n counts for each scan
  print("Further assessing multiples...")
  multiplesCnt <- scanComparison %>% 
              group_by(Scan) %>% 
              filter(row_number()==1) %>% 
              droplevels()

#rename columns to be sure graph legend will be right
  names(multiplesCnt)[grep(nNames[1], names(multiplesCnt))] <- easyNames[1]
  names(multiplesCnt)[grep(nNames[2], names(multiplesCnt))] <- easyNames[2]

#replace scores with booleans for confidence (over thresholds)
  multiplesCnt$conf1 <- multiplesCnt[[scoreNames[1]]] >= thresholds[1]
  multiplesCnt$conf2 <- multiplesCnt[[scoreNames[2]]] >= thresholds[2]
  names(multiplesCnt)[names(multiplesCnt) == "conf1"] <- confNames[1]
  names(multiplesCnt)[names(multiplesCnt) == "conf2"] <- confNames[2]

  totalScans <- nrow(multiplesCnt)

  multiplesCnt <- melt(multiplesCnt[ ,c("Scan",easyNames,confNames)], id.vars = c("Scan",confNames))

  var <- multiplesCnt$variable == easyNames[1]
  conf <- logical(nrow(multiplesCnt))
  for (i in 1:nrow(multiplesCnt)) {
    conf[i] <- ifelse(var[i]==T, multiplesCnt[[confNames[1]]][i], multiplesCnt[[confNames[2]]][i])
  }
  multiplesCnt$conf <- conf

  labeldf <- as.data.frame(table(as.factor(multiplesCnt$value)))

  multiplesBarchartFile <- paste0(c(prefix,"multiples.png"), collapse = "_")
  png(multiplesBarchartFile, width = 760)
  multiplesBarchart <- ggplot(multiplesCnt) +
    geom_bar(aes(x=as.factor(value), fill=conf), position = position_dodge2(preserve = "single")) +
    theme_ipsum(base_size = 12, axis_title_size = 12) +
    scale_y_continuous(sec.axis = sec_axis(~ . / totalScans, labels=scales::percent, name = "% total scans")) +
    labs(title="Peptides suggested per identified spectrum", fill="confident id") +
    facet_wrap(~variable) +
    xlab("Peptides per scan") +
    ylab("Scan count")
  print(multiplesBarchart) 
  invisible(dev.off)

  if (nrow(labeldf) >3) {
    labeldf2 <- as.data.frame(table(as.factor(multiplesCnt$value), multiplesCnt$variable, multiplesCnt$conf))
    multiplesBarchart2File <- paste0(c(prefix,"multiples_withinsets.png"), collapse = "_")
    png(multiplesBarchart2File, width = 760)

    multiplesBarchart2 <- ggplot(labeldf2) +
      geom_col(aes(x=as.factor(Var1), y=Freq,fill=Var3), position = position_dodge2(preserve = "single")) +
      theme_ipsum(base_size = 12, axis_title_size = 12) +
      scale_y_continuous(sec.axis = sec_axis(~ . / totalScans, labels=scales::percent, name = "% total scans")) +
      labs(title="Peptides suggested per identified spectrum", fill="confident id") +
      facet_wrap(~Var2) +
      xlab("peptides per scan") +
      ylab("scan count")
      
    ymin <- round(max(labeldf2$Freq)/2)
    ymax <- max(labeldf2$Freq)
    xmin <- as.numeric(labeldf[4,1]) -1
    xmax <- as.numeric(labeldf[nrow(labeldf),1]) +0.5
    
    keep <- levels(labeldf2$Var1)[-c(1,2,3)]
    labeldf_zoom <- labeldf2[labeldf2$Var1 %in% keep, ]
    insets2 <- labeldf_zoom %>% 
      split(f = .$Var2) %>%
      purrr::map(~annotation_custom2(
        grob = ggplotGrob(get_inset2(.)), 
        data = data.frame(Var2=unique(.$Var2)),
        ymin = ymin, ymax=ymax, xmin=xmin, xmax=xmax)
      )
    multiplesBarchart2 <- multiplesBarchart2 + insets2
    print(multiplesBarchart2) 
    invisible(dev.off)
  }
}
print("Finished")