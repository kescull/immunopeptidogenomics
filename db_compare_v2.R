library(VennDiagram)
library(hrbrthemes)
library(optparse)
library(tidyverse)
library(UpSetR)
sessionInfo()
options(warn=1,hrbrthemes.loadfonts=TRUE)

pareDown <- function(dat) {
  dat <- dat %>%
    select(Peptide,X.10lgP,Length,Scan,Accession,PTM)
  dat$Scan <- as.factor(dat$Scan)
  dat
}

find.contam.scans <- function(dat) {
  scans <- dat %>%
    filter(grepl("*#CONTAM#*", dat$Accession) | 
             grepl("*iRT-pep*", dat$Accession)) %>%
    pull(Scan) %>%
    droplevels()
}
get.contam.rows <- function(dat) {
  dat <- subset(dat, (dat$Scan %in% all.contam.scans)) %>%
    droplevels()
}
get.noncontam.rows <- function(dat) {
  dat <- subset(dat, !(dat$Scan %in% all.contam.scans)) %>%
    droplevels()
}
discard.artificial <- function(dat) {
  dat %>% 
    filter(!(Peptide %in% discard.list)) %>% 
    droplevels() %>%
    as.data.frame()
}
getMaxChimera <- function(dat,y) {
  dat[which.max(dat$X.10lgP), ]
}
remove.chimeras <- function(dat) {
  dat %>% 
    group_by(Scan,Peptide) %>%
    group_modify(getMaxChimera) %>%
    ungroup() %>%
    droplevels() %>%
    as.data.frame()
}
filterThreshold <- function(csv, threshold) {
  newCsv <- csv %>% filter(X.10lgP >= threshold) %>% droplevels() 
}
filterLength <- function(csv) {
  newCsv <- csv %>% filter(Length < 16, Length > 6) %>% droplevels()
}
makePartitionFn <- function(x) {
  textFile <- paste0(paste0(c(prefix,x), collapse="_"),".txt")
}
rename.columns <- function(dat,prefix) {
  dat <- dat %>%
    rename_at(vars(-Scan), ~ paste(prefix, ., sep="."))
}

mark.unco <- function(dat) {
  dat %>% mutate(unco = (Peptide %in% Nlist))
}
count.ids.per.scan <- function(dat) {
  ids.per.scan <- dat %>% 
    group_by(Scan) %>%
    summarise(ids = n())
}
rename.id.summary <- function(dat,colname) {
  dat <- dat %>% rename(!!colname := n)
}
# # # start 'main' script here
#get data

option_list <- list(
  make_option(c("-c", "--cryptic"), type="character", default=NULL, 
              help="cryptic db PEAKS results filename", metavar="character"),
  make_option(c("-n", "--normal"), type="character", default=NULL, 
              help="normal (e.g. uniprot) PEAKS results filename", metavar="character"),
  make_option(c("-p", "--prefix"), type="character", default="out",
              help="prefix for all output files (e.g. myCellsRep1)", metavar="character"),
  make_option(c("-d", "--cryptic_threshold"), type="double", default=0, 
              help="threshold score for cryptic search results", metavar="number"),
  make_option(c("-m", "--norm_threshold"), type="double", default=0,
              help="threshold score for normal (uniprot) search results", metavar="number"),
  make_option(c("-j", "--junction"), type="character", default=NULL, 
              help="Optional: artificial junction peptide list for discard (txt)", metavar="character"),
  make_option(c("-u", "--unconventional"), type="character", default=NULL, 
              help="Optional: list of unconventional peptides (txt)", metavar="character")
); 

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

if (is.null(opt$cryptic) || 
    is.null(opt$normal) || 
    opt$cryptic_threshold == 0 || 
    opt$norm_threshold == 0) {
  print_help(opt_parser)
  stop("Options -c, -d, -n, -m are required (input files and thresholds)", call.=FALSE)
}

filenames <- c(opt$cryptic,opt$normal)
thresholds <- c(opt$cryptic_threshold,opt$norm_threshold)
message("Files:")
message(filenames)
message("Thresholds:")
thresholds
prefix <- opt$prefix
message(paste("prefix:", prefix))

if ((!is.null(opt$junction) & is.null(opt$unconventional)) |
    (is.null(opt$junction) & !is.null(opt$unconventional))) {
      print_help(opt_parser)
      stop("Options -j and -u must be used in conjunction (lists of artificial and unconventional peptides, respectively)", 
         call.=FALSE)
}
discard.list <- NULL
Nlist <- NULL
if (!is.null(opt$junction)) {
  message(paste("discarding peptides in file:", opt$junction))
  prefix <- paste0(c(prefix,"phase2"), collapse = "_")
  message(paste("prefix now:", prefix))
  discard.list <- scan(opt$junction, what = character())
  str(discard.list)
  message(paste("list of unconventional peptides provided:", opt$unconventional))
  Nlist <- scan(opt$unconventional,what = character())
  str(Nlist)
}
easyNames <- (c("cryptic","standard"))
allDat <- lapply(filenames, read.csv)
names(allDat) <- easyNames

#Next, get rid of unwanted columns -
#i.e. keep Peptide, length, X.10lgP, Scan, Accession, PTM
paredDown <- lapply(allDat, pareDown)
                 
message("Dealing with any iRT peptides and contaminants...")
#get scan numbers for contam peptides and write to file, then delete rows from all sets

contam.scans <- lapply(paredDown, find.contam.scans)
all.contam.scans <- as.factor(unlist(lapply(contam.scans, levels)))
contam.rows <- lapply(paredDown, get.contam.rows)

renamed.contam <- mapply(rename.columns,contam.rows, easyNames, SIMPLIFY = F)
merged.contam <- reduce(renamed.contam, full_join, by="Scan") %>%
  relocate(Scan)

contam.fn <- paste0(c(prefix,"contam.csv"), collapse = "_")
write.csv(merged.contam, file = contam.fn, na="")
message(paste("Output: ", contam.fn))

#continue with only non-contamination scans
no.contam <- lapply(paredDown, get.noncontam.rows)

#remove chimeras (when they cause multiple rows the same Peptide for same scan. 
#By this strategy, chimeric spectra that show 2 different sequences count as 'multiple' and get put in ambiguous list)
no.chimeras <- lapply(no.contam, remove.chimeras)

#remove discard lists if present
if (!is.null(discard.list)) {
  message("Removing peptides on 'discard' list...")
  no.chimeras <- lapply(no.chimeras, discard.artificial)
}

#filter for those over threshold
message("Filtering for confident identifications (over score threshold)...")
filteredThresh <- mapply(filterThreshold, no.chimeras, thresholds, SIMPLIFY = FALSE)

#filter for length:
filtered <- lapply(filteredThresh, filterLength)

message("Preparing Venn diagram and partition lists (for those that pass length filter)...")
#venn diagram
peptideLists <- lapply(filtered, function(dat) dat[,"Peptide"])
vennFile <- paste0(c(prefix,"Venn.png"),collapse="_")
invisible(futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger"))

invisible(venn.diagram(x = peptideLists, 
                       category.names = names(peptideLists),
                       filename = vennFile,
                       imagetype="png",
                       cat.prompts=TRUE,
                       margin=0.1,
                       cat.dist=rep(0.05,length(peptideLists)),
                       output = TRUE))
partitions <- get.venn.partitions(peptideLists)

#for each list of values in 'partitions', write a txt file named appropriately and containing the list of peptides found in that row
partitionNames <- c("overlap",
                    "standard_only",
                    "cryptic_only")

partitionFn <- unlist(lapply(partitionNames, makePartitionFn))

invisible(mapply(writeLines,partitions$..values..,partitionFn,sep = "\n"))
message("Output files:")
print(partitionFn)

###Phase 2 analysis:
if (!is.null(Nlist)) {
  message("\nPhase 2: Analysing levels of ambiguity and getting unambiguous and unconventional ids:")
  #add column to dats filtered by threshold but not by length - is it unconventional?
  unconv.marked <-lapply(filteredThresh, mark.unco)
  
  #plot how many ids per scan
  message("Assessing number of ids for each spectrum...")
  ids.per.scan <- lapply(unconv.marked, count.ids.per.scan)
  id.summary <- lapply(ids.per.scan, function(dat) dat %>% count(ids))
  renamed.id.summary <- mapply(rename.id.summary, id.summary, easyNames, SIMPLIFY = F) 
  merged.id.summary <- reduce(renamed.id.summary, full_join, by="ids")
  id.summary.long <- merged.id.summary %>% 
    pivot_longer(all_of(easyNames)) %>%
    replace(is.na(.), 0) %>%
    group_by(name) %>%
    mutate(total = sum(value),
           pc = value/total*100) %>%
    ungroup()
  
  idsperscan.fn <- paste0(c(prefix,"idsPerScan.png"), collapse = "_")
  
  idsPerScanPlot <- ggplot(id.summary.long,aes(x=ids, fill=name, y=pc)) +
    geom_bar(position = "dodge", stat="identity") +
    geom_text(aes(label=value), position=position_dodge(width=0.9), vjust=-0.25, size=3) +
    theme_ipsum(base_size = 12, axis_title_size = 12) +
    labs(title="Confident ids per spectrum") +
    theme(legend.title = element_blank()) +
    xlab("Number of ids") +
    ylab("Percentage of spectra")
  ggsave(idsperscan.fn, width = 20, height = 15, units = "cm")
  message(paste("Output: ", idsperscan.fn))
  
  message("Annotating and merging data from different searches...")
  #mark ids/scan in each row - adds column "ids" to unconv.marked dataframes
  with.multiples <- mapply(merge,unconv.marked,ids.per.scan, SIMPLIFY=F)
  
  #rename columns in preparation for merging
  renamed <- mapply(rename.columns, with.multiples, easyNames, SIMPLIFY=F)
  
  #merge results from different searches into one table by scan
  merged <- reduce(renamed, full_join, by="Scan")
  
  #3) make new column to indicate if searches reached different (confident) conclusions
  merged <- merged %>%
    mutate(different.per.search = if_else((!is.na(standard.Peptide) & !is.na(cryptic.Peptide) & standard.Peptide != cryptic.Peptide),
                                          T,F)) %>%
    mutate(multiple.per.search = if_else((!is.na(standard.ids) & standard.ids > 1) | 
                                           (!is.na(cryptic.ids) & cryptic.ids > 1)
                                         ,T,F)) %>%
    mutate(has.unconventional = if_else((!is.na(standard.unco) & standard.unco == T) | 
                                          (!is.na(cryptic.unco) & cryptic.unco == T),
                                        T,F)) %>%
    relocate(Scan)
  fn <- paste0(c(prefix,"all_merged.csv"), collapse = "_")
  write.csv(merged, file = fn, na="")
  message(paste("Output: ", fn))
  
  ##########################################################
  #Trying to simplify - use UpSet to classify each Scan by Different per search, Multiple per search, Unconventional.
  #For Different By search - want F if same results in all searches, 
  #   so F if scan copies diff= F; or all scan copies have same multiple id count & count(diff) = (copy number - multiple id count).
  message("Assessing ambiguity among identifications for each spectrum...")
  one.row.per.scan <- merged %>%
    replace(is.na(.), 0) %>%
    group_by(Scan) %>%
    summarise(different = if_else(sum(different.per.search) == 0 | 
                                    (sum(multiple.per.search > 0) & 
                                       (sum(standard.ids)/n() == sum(cryptic.ids)/n()) & 
                                       (sum(different.per.search) == (n() - sum(standard.ids)/n()))),
                                  F,T),
              multiple = if_else(sum(multiple.per.search) == 0,F,T),
              unambiguous = if_else(multiple + different == 0, T,F),
              unconventional = if_else(sum(has.unconventional) == 0,F,T)) %>%
    ungroup()
  #write.csv(one.row.per.scan, file = paste0(c(prefix,"check_for_upset.csv"), collapse = "_"), na="")
  
  for.Upset2 <- one.row.per.scan %>% select(-Scan) %>%
    mutate(across(c(different,multiple,unambiguous),as.integer)) %>%
    as.data.frame()
  fn2 <- paste(prefix,"spectra_upset.png",sep="_")
  png(file = fn2, width = 760, height = 480)
  upsetPlot <- upset(for.Upset2, 
                     order.by = "freq",
                     sets = c("different","multiple","unambiguous"),
                     mainbar.y.label = "Spectra", sets.x.label = "Category", 
                     text.scale=c(2,2,2,1,2,2),
                     query.legend = "top",
                     queries = list(
                       list(
                         query = elements,
                         params = list("unconventional", T),
                         color = "#Df5286", 
                         active = T,
                         query.name = "Unconventional peptides"
                       )
                     ))
  upsetPlot$legend$children$GRID.cellGrob.67$children$GRID.text.66$gp$fontsize <- 12
  print(upsetPlot)
  invisible(dev.off())
  message(paste("Output: ", fn2))
  
  
  ##########################################################
  
  # Need to get lists, overlap, 
  # Then apply length filter (need to apply after ambiguity test but before overlap)
  
  #add one.row.per.scan categories to merged
  with.categories <- merge(merged, one.row.per.scan, by="Scan") 
  #print ambiguous scans 
  ambiguous <- with.categories %>%
    filter(unambiguous == F) %>%
    select(-c(different.per.search,multiple.per.search,has.unconventional))
  fn3 <- paste0(c(prefix,"ambiguous.csv"), collapse = "_")
  write.csv(ambiguous, file = fn3, na="")
  message(paste("Output: ", fn3))
  
  message("Proceeding with spectra which received an unambiguous identification, filtered by length...")
  #make unambiguous list of scans - by definition one row per scan
  simple <- with.categories %>% 
    filter(unambiguous == T) %>%
    select(-c(different.per.search,
	      multiple.per.search,
	      has.unconventional,
	      different,
	      multiple)) %>%
    select(-ends_with("ids")) %>%
    select(-ends_with("unco")) %>%
    select(-unambiguous) %>%
    mutate(Peptide = if_else(!is.na(standard.Peptide), standard.Peptide, cryptic.Peptide),
           Length = if_else(!is.na(standard.Length), standard.Length, cryptic.Length),
           PTM = if_else(!is.na(standard.PTM), standard.PTM, cryptic.PTM)) %>%
    rename(standard.score = standard.X.10lgP,
           cryptic.score = cryptic.X.10lgP) %>%
    mutate(standard.score = ifelse(is.na(standard.score), 0, standard.score),
           cryptic.score = ifelse(is.na(cryptic.score), 0, cryptic.score)) %>%
    select(Scan,Peptide,unconventional,Length,PTM,ends_with("score"),
	   ends_with("Accession")) %>%
    filter(Length < 16, Length > 6) %>%
    droplevels()
  #str(simple)
  #write.csv(simple, file = paste0(c(prefix,"simple.csv"), collapse = "_"), na="")
  
  #make table showing which peptides found with which db and whether unconventional, 
  #showing max score for that peptide with that database (not necessarily for the same scan)
  peptides.by.db <- simple %>%
    select(Peptide,PTM,ends_with("score"),unconventional,standard.Accession) %>%
    group_by(Peptide) %>%
    summarise(standard = max(standard.score),
              cryptic = max(cryptic.score),
	      Accession = paste0(standard.Accession,collapse=";"),
              unconventional = ifelse(sum(unconventional) > 0,T,F)) %>%
    as.data.frame()  
  fn4 <- paste0(c(prefix,"unambiguous_peptides_by_db.csv"), collapse = "_")
  write.csv(peptides.by.db, file = fn4, na="")
  message(paste("Output: ", fn4))
  
  #make equivalent df but for upset visualisation
  overlap <- simple %>%
    select(Peptide,ends_with("score"),unconventional) %>%
    group_by(Peptide) %>%
    summarise(standard = ifelse(sum(standard.score) > 0,1,0),
              cryptic = ifelse(sum(cryptic.score) > 0,1,0),
              unconventional = ifelse(sum(unconventional) > 0,T,F)) %>%
    select(-ends_with("score"),Peptide) %>%
    as.data.frame()  
  
  fn5 <- paste(prefix,"db_upset.png",sep="_")
  png(file = fn5, width = 760, height = 480)
  query.title <- paste0("Unconventional peptides (n = ", sum(overlap$unconventional), ")")
  upsetPlotDb <- upset(overlap,
                       sets = easyNames,
                       mainbar.y.label = "Peptides", 
		       sets.x.label = "Number found with each database", 
                       order.by = "freq",
                       query.legend = "top",
                       text.scale=c(2,2,1,2,2,2),
                       queries = list(
                         list(
                           query = elements,
                           params = list("unconventional", T),
                           color = "#Df5286", 
                           active = T,
                           query.name = query.title
                         )
                       )) 
  upsetPlotDb$legend$children$GRID.cellGrob.164$children$GRID.text.163$gp$fontsize <- 12
  print(upsetPlotDb)
  invisible(dev.off())
  message(paste("Output: ", fn5))
  
  #get list of unambiguous unconventional (ie to help filter origins for unambiguous unconventional peptides)
  unambiguous.unconventional <- overlap %>%
    filter(unconventional == T) %>%
    pull(Peptide)
  fn6 <- paste(prefix,"unambiguous_unconventional.txt",sep="_")
  writeLines(unambiguous.unconventional,fn6,sep="\n")
  message(paste("Output: ", fn6))

  #get txt lists of unambiguous peptides
  standard.unambiguous <- overlap %>%
    filter(standard == 1, cryptic == 0) %>%
    pull(Peptide)
  cryptic.unambiguous <- overlap %>%
    filter(cryptic == 1, standard == 0) %>%
    pull(Peptide)
  overlap.unambiguous <- overlap %>%
    filter(cryptic == 1, standard == 1) %>%
    pull(Peptide)

  fn7 <- paste(prefix,"standard_unambiguous.txt",sep="_")
  writeLines(standard.unambiguous,
             fn7,
             sep = "\n")
  message(paste("Output: ", fn7))
  fn8 <- paste(prefix,"cryptic_unambiguous.txt",sep="_")
  writeLines(cryptic.unambiguous,
             fn8,
             sep = "\n")
  message(paste("Output: ", fn8))
  fn9 <- paste(prefix,"overlap_unambiguous.txt",sep="_")
  writeLines(overlap.unambiguous,
             fn9,
             sep = "\n")
  message(paste("Output: ", fn9))
}

message("\nFinished")
