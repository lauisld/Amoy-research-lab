## This file will find tprK donor sites around the tprD locus
## based on the tprK variable region nucleotide sequences.
## It creates a GFF file for annotation and creates figures of the lengths, locations, and steps of
## the variable regions.

## The following step must be done before using!
## Set path to folder with final_data_dna.csv files. This can be found in the Illumina_frequencies 
## after running the main pipeline tprk_pipeline.py. Copy these csv files into a new folder that also contains the file contents
## of the "donorsite_figures" folder in tprk_paper2-master (eg. tprDlocus.fasta).
path <- "/Users/uwvirongs/Documents/Michelle/tprk_pipeline/2.01_donorsite_redo/"

list.of.packages <- c("ggplot2","stringr","Biostrings","dplyr","data.table","gggenes","seqinr","phylotools","cowplot")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
suppressMessages(invisible(lapply(list.of.packages,library,character.only=T)))

setwd(path)

colorscheme <- c("V1" = "#e41a1c", "V2" = "#377eb8", "V3" = "#4daf4a", "V4" = "gold3",
                 "V5" = "#ff7f00","V6" = "#984ea3","V7" = "#a65628")

##Read in all of the over5*.csv files to get the nucleotide sequences of the variable regions
filenames <- Sys.glob("over5*")
for (i in c(1:length(filenames))){
  seq <- NULL
  sample<-as.character(strsplit(strsplit(filenames[i],"_")[[1]][2],"\\.")[[1]][1])
  mycsv <- read.table(filenames[i],sep=',',header = TRUE,check.names=FALSE) # Modified by Elaine
  newnames <- paste(sample,mycsv$Region,rownames(mycsv),sep='_')
  class(sample)
  mycsv <- mutate(mycsv,newnames=paste(sample,mycsv$Region,as.character(Count),rownames(mycsv),sep='_'))
  mycsv <- mycsv[mycsv$Count > 10,]
  dna <- DNAStringSet(mycsv$Read,use.names=TRUE)
  names(dna) <- mycsv$newnames
  filename <- paste0(sample,"_variable.fasta")
  writeXStringSet(dna, filename)
}
syscommand <- "cat *_variable.fasta > all-variable.fasta"
system(syscommand)
## Add a step to trim off 14bp from end of V2 (AGTATGGATTGGGG); from end of V5 (TGCTGCCTATATT); 
## and from beginning of V3 (TCATACTCACCTTAGCCCCGACA, edit distance of one)

## Filter variable regions based on prior knowledge of conserved sequences
allvariable <- read.fasta("all-variable.fasta")
allvariable$region <- str_split_fixed(allvariable$seq.name,"_",4)[,2]
allvariable[allvariable$region=="V3",]$seq.text <- str_sub(allvariable[allvariable$region=="V3",]$seq.text, end=-24)
allvariable[allvariable$region=="V2",]$seq.text <- str_sub(allvariable[allvariable$region=="V2",]$seq.text, end=-15)
allvariable[allvariable$region=="V5",]$seq.text <- str_sub(allvariable[allvariable$region=="V5",]$seq.text, end=-14)
allvariable[grep("N",allvariable$seq.text),]$seq.text <- NA
allvariable <- allvariable[complete.cases(allvariable), ]
allvariable <- mutate(allvariable,seq.length = nchar(allvariable$seq.text))
dnaString <- BStringSet(allvariable$seq.text)
names(dnaString) = paste0(allvariable$seq.name)
writeXStringSet(dnaString, "all-variable_edited.fasta",append=FALSE,format="fasta")

### Blast the edited allvariable fasta against our 17.2kb locus
syscommand_amin <- "blastn -db tprDlocus.fasta -dust no -soft_masking false -query all-variable_edited.fasta -evalue 1e1 -word_size 10 -ungapped -perc_identity 100 -penalty -10 -outfmt 6 -max_hsps 3 -subject_besthit -out Step3_search.blast"
system(syscommand_amin)

### Set up the blasthits file
blasthits <- read.table("Step3_search.blast",
col.names = c("qseqid","sseqid","pident","length","mismatch",
              "gapopen","qstart","qend","sstart","send","evalue","bitscore"))

count <- as.numeric(str_split_fixed(blasthits$qseqid,"_",4)[,3])
sample <- str_split_fixed(blasthits$qseqid,"_",4)[,1]
region <- str_split_fixed(blasthits$qseqid,"_",4)[,2]
interval<-paste(blasthits$sstart,blasthits$send,sep="-")
blasthits <- dplyr::mutate(blasthits,count=count,sample=sample,interval=interval,region=region)
blasthits <- blasthits %>% dplyr::group_by(sample, region)  %>% dplyr::mutate(regionsamplecounts = sum(count))
blasthits <- ungroup(blasthits)
blasthits <- dplyr::mutate(blasthits, percentage = round(count / regionsamplecounts * 100,3))

##Filter blasthits above an arbitrary count, in this case 50 reads and within-sample percentage of 0.2%
blasthits_100 <- blasthits[blasthits$count > 50 & blasthits$percentage > 0.2,]

##Put the longest blasthits at the top to ease the sort below.
blasthits_100 <- blasthits_100[order(-blasthits_100$length),]

## Function to find the largest locus for each donor site
UniqueDonorSites <- function(df) {
  variableReg <- data.frame(df$sstart[1],df$send[1],df$region[1],df$count[1],df$sample[1],df$percentage[1])
  variableReg$morethanone <- 1
  colnames(variableReg)<- c("start","end","region","count","sample","percentage","morethanone")
  counter <- 1
  for (i in c(2:length(df$sstart)))
  {
    indf <- FALSE
    for (j in c(1:length(variableReg$start)))
    {
      if (df$region[i] != variableReg$region[j]) next()
      if (inrange(df$sstart[i], variableReg$start[j], variableReg$end[j], incbounds=TRUE) & 
          inrange(df$send[i], variableReg$start[j], variableReg$end[j], incbounds=TRUE)){
        indf <- TRUE
        variableReg$count[j] <- variableReg$count[j] + df$count[i]
        variableReg$percentage[j] <- variableReg$percentage[j] + df$percentage[i]
        if (df$sample[i]!=variableReg$sample[j]) variableReg$morethanone[j] <- 2
      } else if (inrange(df$sstart[i], variableReg$end[j],variableReg$start[j], incbounds=TRUE) & 
                 inrange(df$send[i], variableReg$end[j],variableReg$start[j], incbounds=TRUE)){
        indf <- TRUE
        variableReg$count[j] <- variableReg$count[j] + df$count[i]
        variableReg$percentage[j] <- variableReg$percentage[j] + df$percentage[i]
        if (df$sample[i]!=variableReg$sample[j]) variableReg$morethanone[j] <- 2
      } 
    }
    if (!indf) {
      counter <- counter+1
      print(c(counter,i))
      addition <- data.frame(df$sstart[i],df$send[i],df$region[i],df$count[i],df$sample[1],df$percentage[i])
      addition <- mutate(addition,morethanone=1)
      if ((inrange(df$sstart[i], variableReg$end[j],variableReg$start[j], incbounds=TRUE) | 
          inrange(df$send[i], variableReg$end[j],variableReg$start[j], incbounds=TRUE)) | 
        (inrange(df$sstart[i], variableReg$start[j], variableReg$end[j], incbounds=TRUE) | 
           inrange(df$send[i], variableReg$start[j], variableReg$end[j], incbounds=TRUE))) addition$morethanone <- 2
      variableReg <- rbind(variableReg,setNames(addition,names(variableReg)))
      indf <- FALSE
    }
  }
  return(variableReg)
}

##Get unique donor sites (can take a minute)
variableRegions <- UniqueDonorSites(blasthits_100)

### Read in the fasta locus sequence to store in variableRegions dataframe
sequence <- phylotools::read.fasta("tprDlocus.fasta")
tprDlocusSeq <- sequence$seq.text
seqname <- sequence$seq.name
variableRegions$region <- as.character(variableRegions$region)
variableRegions$region <- factor(variableRegions$region, levels = c("V1","V2","V3","V4","V5","V6","V7"))
variableRegions <- mutate(variableRegions,length=abs(start-end)+1,
                          seqname=seqname,source=".",feature="variable",score=".",frame=0,
                          strand=ifelse(start-end>0,"-","+"))

##Function to get a list of sequences of the variable regions by parsing the original fasta file
getVariableSequences <- function(variableRegions){
loci <- NULL
for (i in c(1:length(variableRegions$start))){
  if (variableRegions$start[i] > variableRegions$end[i]) { 
    start <- (variableRegions$end[i])
    end <- (variableRegions$start[i])
    loci <- getSequence(rev(comp(getSequence(substr(tprDlocusSeq,start,end)
              ,split=""),forceToLower = FALSE)),as.string = TRUE)[[1]] } 
  else {
    start <- variableRegions$start[i]
    end <- variableRegions$end[i]
    loci <-substring(tprDlocusSeq,start,end) }
  variableRegions$sequence[i] <- loci }
return(variableRegions)
}

variableRegions <- getVariableSequences(variableRegions)


## Get a count of number of times a putative donor site sequence is represented in the variable region dataset
for (i in c(1:length(variableRegions$start))){
  variableRegions$sequence_count[i] <- sum(str_count(variableRegions[variableRegions$region==variableRegions$region[i],]$sequence,
                                                     variableRegions$sequence[i])) }

## Keep just the variable regions that are unique
variableRegions <- variableRegions[variableRegions$sequence_count==1,]
variableRegions$sequence_count <- NULL

### Expand donor sites that are overlapping to get full-length of potential donor site
### And count how many samples they are in
variableReg <- variableRegions
variableReg <- variableReg[order(-variableReg$length),]
for (i in c(1:(length(variableReg$feature)-1))){
    for (j in c(i:length(variableReg$feature))){
      if (variableReg$region[i] != variableReg$region[j]) next()
      if (i==j) next()
      if ((inrange(variableReg$start[i], variableReg$start[j], variableReg$end[j], incbounds=TRUE)) | 
        (inrange(variableReg$start[i], variableReg$end[j], variableReg$start[j], incbounds=TRUE))){
        variableReg$start[i] <- variableReg$start[j]
        variableReg$count[i] <- variableReg$count[j] + variableReg$count[i]
        variableReg$morethanone[i] <- 2
        variableReg$count[j] <- 0
        variableReg$start[j] <- 0
        variableReg$end[j] <- 0
      } else if ((inrange(variableReg$end[i], variableReg$start[j], variableReg$end[j], incbounds=TRUE)) | 
                 (inrange(variableReg$end[i], variableReg$end[j], variableReg$start[j], incbounds=TRUE))){
        variableReg$end[i] <- variableReg$end[j]
        variableReg$count[i] <- variableReg$count[j] + variableReg$count[i]
        variableReg$morethanone[i] <- 2
        variableReg$start[j] <- 0
        variableReg$end[j] <- 0
        variableReg$count[j] <- 0
      }
    }
}

### Get rid of donor site sequences that are now redundant or within another site
variableRegions <- variableReg[!variableReg$start==0,]

## For specificity sake, we are only calling variable regions if they are present in more than one sample.
## Cull all variable regions not present in more than one sample
variableRegions <- variableRegions[variableRegions$morethanone==2,]
variableRegions <- mutate(variableRegions,length=abs(start-end)+1)


## Get a count of number of times a sequence is represented in the variable region dataset
## This in case two variable regions became one but had not yet been culled.
variableRegions <- getVariableSequences(variableRegions)
for (i in c(1:length(variableRegions$start))){
  variableRegions$sequence_count[i] <- sum(str_count(variableRegions[variableRegions$region==variableRegions$region[i],]$sequence,
                                                     variableRegions$sequence[i])) }

variableRegions <- variableRegions[variableRegions$sequence_count==1,]
variableRegions$sequence_count <- NULL

## Name the variable regions
variableRegions_ogg <- variableRegions
variableRegions <- variableRegions_ogg 
oldannotations <- read.table(paste0(path,"tprDlocus_oldAnnotations2.csv"),sep=",",header=TRUE,check.names=FALSE) # Modified by Elaine
oldannotations <- oldannotations[order(oldannotations$Minimum),]
oldannotations$region <- str_split_fixed(oldannotations$Name,"-",2)[,1]
variableRegions$name <- NA
for (i in c(1:length(oldannotations$Minimum))){
  for (j in  c(1:length(variableRegions$start))){
    if (oldannotations$region[i] != variableRegions$region[j]) next()
    if ((inrange(variableRegions$start[j], oldannotations$Minimum[i], oldannotations$Maximum[i], incbounds=TRUE)) | 
        (inrange(variableRegions$end[j], oldannotations$Minimum[i], oldannotations$Maximum[i], incbounds=TRUE)) |
        (inrange(oldannotations$Minimum[i], variableRegions$start[j], variableRegions$end[j], incbounds=TRUE)) | 
        (inrange(oldannotations$Minimum[i], variableRegions$end[j], variableRegions$start[j], incbounds=TRUE))) {
        variableRegions$name[j] <- oldannotations$Name[i]
    }
  }
}
newregions <- c("V5-DS56","V6-DS55","V6-DS54","V6-DS53","V3-DS52","V6-DS51","V3-DS50","V3-DS49","V4-DS48")
variableRegions[which(is.na(variableRegions$name)),]$name <- newregions


## Export variable region dataframe as a GFF3 format file
col_order <- c("seqname", "source", "region","start","end","score","strand","frame","feature","name")
variableRegionsGFF<- variableRegions[, col_order]
write.table(variableRegionsGFF,"variableRegions_merge3_50_0.2_edited_oldnames.gff",sep='\t',row.names=FALSE,quote=FALSE,col.names=FALSE)

## Average lengths of the donor sites by variable region
aggregate(variableRegions$length, list(variableRegions$region), mean)

## Number of donor sites by variable region
aggregate(variableRegions$length, list(variableRegions$region), length)

##Plot of the donor site lengths by usage (percentage summed across samples)
ggplot(variableRegions,aes(x=length,y=percentage,fill=region)) + geom_bar(stat="identity",position = "stack") + 
  scale_y_continuous() + theme_classic() + scale_fill_manual(values=colorscheme) +
  ylab("Usage (sum percentage)") + xlab("Length (bp)") + theme(legend.position = "none")
ggsave(filename="Figure3C_DonorSites_Length-Percentage_classic.pdf",height = 4, width=6)

## Function to count the number of occurrences of qseqid in a blast file. 
## This is a proxy for the number of donor sites used to make a tprK variable region
qseqid_count <- function(blasthit_df) {
  newsummary <- as.data.frame(table(blasthit_df$qseqid))
  numsteps <- newsummary[newsummary$Freq != 0, ]
  colnames(numsteps)<- c("qseqid","Freq")
  blasthit_df <- left_join(blasthit_df,numsteps, by = "qseqid")
  return(blasthit_df)
}

myhits <- qseqid_count(blasthits_100)

## Count and Calculate Frequency for How Steps are Put Together Based on Percentage
sumtable2 <- myhits %>% dplyr::group_by(region,Freq) %>% dplyr::mutate(percentsum = sum(percentage))
sumtable3 <- unique(data.frame(sumtable2$region,sumtable2$Freq,sumtable2$percentsum))
names(sumtable3) <- c("region","steps","percentsum")
fillin <- data.frame(region=c("V3","V7"),steps=c(1,1),percentsum=c(0,0))
sumtable3 <- rbind(sumtable3,fillin)
sumtable3 <- as.data.frame(setDT(sumtable3)[, percentage := round((percentsum/sum(percentsum))*100,3), by = region])

sumtable3$steps <- factor(sumtable3$steps,levels=c("3","2","1"))
Fig4A <- ggplot(sumtable3[order(sumtable3$steps, decreasing = T),],aes(x=region,y=percentage,fill=steps)) +
  geom_bar(stat="identity",position = "stack") + theme_classic() +  xlab("Region") + ylab("Percentage") +
  ylim(0, 101) + scale_fill_brewer(palette="Set1",name= "Events",breaks=c("1","2","3")) + theme(text=element_text(size=8))
ggsave(filename="Figure4A_DonorSites_Events_Region-Percentage_Stacked.pdf",height = 4, width=6)

sumtable3$steps <- factor(sumtable3$steps,levels=c("1","2","3"))
ggplot(sumtable3,aes(x=region,y=percentage,fill=as.factor(steps))) +
  geom_bar(stat="identity",position = "dodge",mapping=aes(group=as.factor(steps))) + theme_classic() + 
  ylim(0, 100) + scale_fill_brewer(palette="Set1",name= "Steps") + xlab("Region") + ylab("Percentage")
ggsave(filename="DonorSites_Events_Region-Percentage_Dodge.pdf",height = 4, width=6)

## Count and Calculate Frequency for How Steps are Put Together

sumtable <- as.data.frame(table(myhits$region,myhits$Freq))
names(sumtable) <- c("region","steps","count")
sumtable <- as.data.frame(setDT(sumtable)[, percentage := round((count/sum(count))*100,3), by = region])
#gather(sumtable,key="region",value="count")

region_sites <- aggregate(variableRegions$length, list(variableRegions$region), length)
colnames(region_sites) <- c("region","sites")
region_sites <- mutate(region_sites,total_w_replacement=sites+sites^2+sites^3)
region_sites <- mutate(region_sites,total_wo_replacement=sites+sites*(sites-1)+sites*(sites-1)*(sites-2))

## No single-step V3/V7 events, remove them from the possibility
region_sites$total_wo_replacement[region_sites$region=="V3"] <- region_sites$total_wo_replacement[region_sites$region=="V3"]-region_sites$sites[region_sites$region=="V3"]
region_sites$total_wo_replacement[region_sites$region=="V7"] <- region_sites$total_wo_replacement[region_sites$region=="V7"]-region_sites$sites[region_sites$region=="V7"]

#Simple model - one step/event for all
prod(region_sites$sites)
#Simple model - three steps/events for all - all independence, with replacement
prod(region_sites$total_w_replacement)
#Simple model - three steps all - all independence, without replacement
prod(region_sites$total_wo_replacement)



sumtable <- as.data.frame(table(myhits$region,myhits$Freq,myhits$sample))
names(sumtable) <- c("region","steps","sample","count")
sumtable <- as.data.frame(setDT(sumtable)[, percentage := round((count/sum(count))*100,3), by = c("sample","region")])

##Format the GFF for the gggenes plot
tprDlocusGFF <- read.table("tprDlocus_annot.gff",sep='\t',stringsAsFactors = FALSE)
tprDlocusGFF$V2 <- NULL
tprDlocusGFF$V8 <- NULL
tprDlocusGFF$V6 <- NULL
names(tprDlocusGFF) <- c("molecule","type","start","end","direction","gene")
tprDlocusGFF <- tprDlocusGFF[-c(1,2),]
tprDlocusGFF <- tprDlocusGFF[order(tprDlocusGFF$start),]
str_sub(tprDlocusGFF$gene,-4,-1) <- ""
str_sub(tprDlocusGFF$gene,1,5) <- ""
tprDlocusGFF[tprDlocusGFF$gene=="hypothetical protein",]$gene <- "hp"
#tprDlocusGFF$gene <- c("","","","","","","","tprD","","","","","","","","","ntpJ","dat","")
#tprDlocusGFF$gene <- c("hp","hp","hp","hp","hp","tprD","hp","hp","hp","hp","hp","hp","hp","hp","ntpJ","dat","hp")
tprDlocusGFF[tprDlocusGFF$direction=="-",]$direction <- FALSE
tprDlocusGFF[tprDlocusGFF$direction=="+",]$direction <- TRUE

## Make a plot of the number of counts by location.
variableRegions2 <- variableRegions
variableRegions2$sequence <- NULL
variableRegions2 <- mutate(variableRegions2,type="variable",direction=TRUE,gene=NA,molecule=NA)
tprDlocus2 <- tprDlocusGFF
tprDlocus2$type <- NULL
tprDlocus2 <- mutate(tprDlocus2,region=NA,count=NA,length=abs(end-start),type="gene",percentage=NA)
col_order <- colnames(tprDlocus2)
#col_order <- c("seqname", "source", "region","start","end","score","strand","frame","feature")
variableRegions2 <- variableRegions2[, col_order]
total <- rbind(variableRegions2,tprDlocus2)
total$molecule <- -100
variableRegions2 %>% dplyr::group_by(region) %>% summarise(sum(percentage))
## Make percentage/count and location plot of the variable regions by region with gggenes tag.

ggplot(total,aes(xmin = start, xmax = end,y=molecule,label=gene)) + 
  geom_bar(data=subset(total,type=="variable"),aes(x=start,y=percentage,fill=region),
           stat="identity",position = "stack",width = 25) + scale_y_continuous() + xlim(0,17200) +
  theme(legend.position = "none") + theme_classic() + xlab ("Nucleotide") + ylab("Percentage") + 
  theme(axis.line.x = element_line(size=0)) + 
  scale_fill_manual(values = c("V1" = "#e41a1c", "V2" = "#377eb8", "V3" = "#4daf4a", "V4" = "gold3",
                               "V5" = "#ff7f00","V6" = "#984ea3","V7" = "#a65628", 
                               "hp" = "#fbf6bd", "tprD" = "#e41a1c", "ktrA" ="#999999","ktrB" = "#999999",
                               "potA" ="#999999",  "ogt" ="#999999","phnU" ="#999999"),
                    breaks = region, labels=region, name = "Region") + theme(legend.position = c(0.9,0.6)) + 
  geom_gene_arrow(data=subset(total,type=="gene"),
                  aes(xmin = start, xmax = end,y=molecule, fill=gene, alpha=0.7, forward = direction)) + 
  geom_gene_label(align = "centre") + scale_alpha(guide = 'none') + ylab("Usage (sum percentage)")
ggsave(filename="Figure3A_tprDLocus_DonorSites_Percentage_GeneLabeled2.pdf",height = 3, width=8)

total$molecule <- -20000


####
groupblast <- blasthits_100 %>% group_by(qseqid) %>% summarise(Max_End=max(qend),Min_Start=min(qstart))
groupblast <- mutate(groupblast,blastlength=Max_End-Min_Start+1)

allvariable <- read.fasta("all-variable_edited.fasta")
allvariable$seqlength <- nchar(allvariable$seq.text)
allvariable <- allvariable %>% rename(qseqid = seq.name)
allvariable <- left_join(groupblast,allvariable)
allvariable$diff <- as.integer(allvariable$blastlength-allvariable$seqlength)
allvariable$region <- str_split_fixed(allvariable$qseqid,"_",4)[,2]
allvariable$count <- as.integer(str_split_fixed(allvariable$qseqid,"_",4)[,3])

allvariable_og <- read.fasta("all-variable.fasta")
allvariable_og$seqlength <- nchar(allvariable_og$seq.text)
allvariable_og <- allvariable_og %>% rename(qseqid = seq.name)
allvariable_og <- left_join(groupblast,allvariable_og)
allvariable_og$diff <- as.integer(allvariable_og$blastlength-allvariable_og$seqlength)
allvariable_og$region <- str_split_fixed(allvariable_og$qseqid,"_",4)[,2]
allvariable_og$count <- as.integer(str_split_fixed(allvariable_og$qseqid,"_",4)[,3])

colorscheme <- c("V1" = "#e41a1c", "V2" = "#377eb8", "V3" = "#4daf4a", "V4" = "gold3",
                 "V5" = "#ff7f00","V6" = "#984ea3","V7" = "#a65628")

blastseqlength_scatter <- ggplot(allvariable, aes(x=seqlength,y=blastlength,color=region)) + 
  scale_color_manual(values=colorscheme) + geom_point(size=1) + theme_minimal() + 
  geom_abline(slope=1,intercept=0) + xlim(15,125) + ylim(15,125)   + 
  xlab("Sequence Length (bp)") + ylab("Donor Site Alignment Length (bp)")

ggsave("Figure_blastlength_seqlength_scatter_122319.pdf",width=4,height=3)

blastseqlength_scatter_og <- ggplot(allvariable_og, aes(x=seqlength,y=blastlength,color=region)) + 
  scale_color_manual(values=colorscheme) + geom_point(size=1) + theme_minimal() + 
  geom_abline(slope=1,intercept=0) + xlim(15,125) + ylim(15,125)   + 
  xlab("Sequence Length (bp)") + ylab("Donor Site Alignment Length (bp)")

#blastseqlength_histo <- ggplot(allvariable,aes(x=abs(diff),fill=region)) + geom_histogram() + 
#  scale_fill_manual(values=colorscheme) + theme_minimal() + labs(fill='Region') + 
#  xlab("Sequence-Alignment Length Difference") + ylab("Count")

blastseqlength_histo <- ggplot(allvariable,aes(x=abs(diff),y=count,fill=region)) + geom_bar(stat="identity") + 
  scale_fill_manual(values=colorscheme) + theme_minimal() + labs(fill='Region') + 
  xlab("Sequence-Alignment Length Difference") + ylab("Count")

blastseqlength_histo_og <- ggplot(allvariable_og,aes(x=abs(diff),y=count,fill=region)) + geom_bar(stat="identity") + 
  scale_fill_manual(values=colorscheme) + theme_minimal() + labs(fill='Region') + 
  xlab("Sequence-Alignment Length Difference") + ylab("Count")

plot_grid(blastseqlength_scatter + theme(legend.position="none"), blastseqlength_histo, labels = c('A', 'B'), label_size = 10)

ggsave("Figure_blastlength_seqlength_123119.pdf",width=7,height=3)

plot_grid(blastseqlength_scatter_og + theme(legend.position="none"), blastseqlength_histo_og,
          blastseqlength_scatter + theme(legend.position="none"), blastseqlength_histo,
          labels = c('A','B','C', 'D'), label_size = 10, nrow=2, ncol=2)

ggsave("FigureS3_blastlength_seqlength_A-D.pdf",width=8,height=6)


### Intermediate code that helped identify conserved sequences in V2, V3, and V5.
### V3 check
V3_seq12 <- str_sub(allvariable[allvariable$region == "V3",]$seq.text,-24,-1)
V3_seq12counts <- as.integer(allvariable[allvariable$region == "V3",]$count)
mylist <- as.data.frame(V3_seq12,V3_seq12counts)
xtabs(V3_seq12counts ~ V3_seq12, mylist)
sum(V3_seq12counts)
mylist %>% group_by(V3_seq12 ) %>% summarise(V3_seq12counts=sum(V3_seq12counts))
V3list <- xtabs(V3_seq12counts ~ V3_seq12, mylist)
V3list <- as.data.frame(V3list)
V3list <- mutate(V3list,error=round(Freq/773239*100,3))
adist(V3list$V3_seq12,"TGTCGGGGCTAAGGTGAGTATGA")

ggplot(subset(allvariable,region == "V3"),aes(x=diff)) + geom_histogram()

##V5 check 
V5_seq13 <- str_sub(allvariable[allvariable$region == "V5",]$seq.text,-14,-1)
V5_seq13counts <- as.integer(allvariable[allvariable$region == "V5",]$count)
mylist <- as.data.frame(V5_seq13,V5_seq13counts)
xtabs(V5_seq13counts ~ V5_seq13, mylist)
sum(V5_seq13counts)
mylist %>% group_by(V5_seq13) %>% summarise(V5_seq13counts=sum(V5_seq13counts))

##V2 check 
V2_seq7 <- str_sub(allvariable[allvariable$region == "V2",]$seq.text,-19,-1)
V2_seq7counts <- as.integer(allvariable[allvariable$region == "V2",]$count)
mylist <- as.data.frame(V2_seq7,V2_seq7counts)
xtabs(V2_seq7counts ~ V2_seq7, mylist)
sum(V2_seq7counts)
mylist %>% group_by(V2_seq7) %>% summarise(V2_seq7counts=sum(V2_seq7counts))
allvariable[allvariable$region == "V2",]$diff


allvariable <- allvariable[allvariable$region == "V5",]$seqlength - 13

#######

### Are all donor sites found in each position across the variable region?
### Determine if donor sites show a bias for different positions in the variable region

blasthits2 <- blasthits_100 %>% dplyr::group_by(qseqid) %>% dplyr::mutate(rank=dense_rank(qstart))
variableRegions_order <- variableRegions

variableRegions_order$firstcount <- 0
variableRegions_order$firstpercentage <- 0
variableRegions_order$secondcount <- 0
variableRegions_order$secondpercentage <- 0
variableRegions_order$thirdcount <- 0
variableRegions_order$thirdpercentage <- 0

### Nested for loops to account for where each donor site ends up in variable region event list 
### Not tidy and couldn't figure out a vectorized one-liner alternative.
for (i in c(1:length(blasthits2$qend))){
  for (j in c(1:length(variableRegions_order$start))){
    if (blasthits2$region[i] != as.character(variableRegions_order$region[j])) next()
    if ((inrange(blasthits2$sstart[i], variableRegions_order$start[j], variableRegions_order$end[j], incbounds=TRUE)) |
        (inrange(blasthits2$sstart[i], variableRegions_order$end[j],variableRegions_order$start[j], incbounds=TRUE))){
      variable_match <- variableRegions_order[j,]
      myindex <- j
      print(i)
      print(paste0("start ",blasthits2$sstart[i]))
      print(paste0("var start ",variableRegions_order$start[j]))
      print(paste0("var end ",variableRegions_order$end[j]))
      blasthits2$name[i] <- variableRegions_order$name[j]
      }
  }
  if (blasthits2$rank[i]==1){
    variableRegions_order$firstcount[myindex] <- variableRegions_order$firstcount[myindex] + variable_match$count
    variableRegions_order$firstpercentage[myindex] <- variableRegions_order$firstpercentage[myindex] + variable_match$percentage
  }
  if (blasthits2$rank[i]==2){
    variableRegions_order$secondcount[myindex] <- variableRegions_order$secondcount[myindex] + variable_match$count
    variableRegions_order$secondpercentage[myindex] <- variableRegions_order$secondpercentage[myindex] + variable_match$percentage
  }
  if (blasthits2$rank[i]==3){
    variableRegions_order$thirdcount[myindex] <- variableRegions_order$thirdcount[myindex] + variable_match$count
    variableRegions_order$thirdpercentage[myindex] <- variableRegions_order$thirdpercentage[myindex] + variable_match$percentage
  }
}

variableRegions_order2 <- tidyr::gather(variableRegions_order,key="rankset",value="rankpercentage","firstpercentage","secondpercentage","thirdpercentage")

variableRegions_order3 <- variableRegions_order2 %>% dplyr::group_by(sequence) %>% mutate(rankpercentpercent = rankpercentage/sum(rankpercentage)*100)
variableRegions_order3 <- ungroup(variableRegions_order3)
variableRegions_order3$rankset <- factor(variableRegions_order3$rankset,levels=c("thirdpercentage","secondpercentage","firstpercentage"))
variableRegions_order3$name <- as.character(variableRegions_order3$name)
variableRegions_order3$name <- factor(variableRegions_order3$name,levels=str_sort(unique(variableRegions_order3$name),numeric = TRUE))

#### Figure usage by site.
Fig4B <- ggplot(variableRegions_order3,aes(x=name,y=rankpercentpercent,fill=rankset)) + 
  geom_bar(stat="identity",position="stack") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90,size = 4,hjust=1)) +
  xlab("Variable Region Donor Site") + 
  ylab("Percentage") + 
  theme(text=element_text(size=8)) +
  facet_grid(~region, scales="free",space="free") +   
  scale_fill_brewer(palette="Set1",name="Position",
                    breaks=c("firstpercentage","secondpercentage","thirdpercentage"),
                    labels=c("First","Second","Third"))

plot_grid(Fig4A+ theme(legend.key.size = unit(0.5, "cm")), 
          Fig4B+ theme(legend.key.size = unit(0.5, "cm")),
          labels = c('A','B'),rel_widths = c(1, 2.2))
ggsave(filename="Figure4AB_DonorSites_Events-Positions-Percentage_Stacked.pdf",height = 3, width=8)


ggplot(variableRegions_order2,aes(x=start,y=rankpercentage,fill=rankset)) + geom_bar(stat="identity",width=25) + 
  theme_minimal() + labs(fill='Region') + xlim(0,5000) +
  xlab("Sequence-Alignment Length Difference") + ylab("percentage")

ggplot(variableRegions_order3,aes(x=as.factor(start),y=rankpercentpercent,fill=rankset)) + 
  geom_bar(stat="identity",position="stack") + 
  theme_minimal() + labs(fill='Region') + theme(axis.text.x = element_text(angle = 90)) +
  xlab("Sequence-Alignment Length Difference") + ylab("percentage") + facet_wrap(~region, scales="free", nrow=1)

variableRegions_order4 <- variableRegions_order2[variableRegions_order2$rankpercentage > 0,]
variable4 <- data.frame(variableRegions_order2$start, variableRegions_order2$end,
                        variableRegions_order2$region,variableRegions_order2$rankset,variableRegions_order2$rankpercentage)
allowed <- aggregate(variableRegions_order4$region, list(variableRegions_order4$region,variableRegions_order4$rankset), length)
######

blasthits4 <- blasthits2 %>% group_by(qseqid) %>% mutate(varregions=n(),distinctregions=n_distinct(name))
blasthits_overlap <- ungroup(blasthits4[which(blasthits4$varregions != blasthits4$distinctregions),])
length(unique(blasthits_overlap$qseqid))
blasthits_overlap %>% group_by(qseqid) %>% group_by(region) %>% summarise(n())
blasthits_overlap2 <- data.frame(blasthits_overlap$qseqid,blasthits_overlap$interval,blasthits_overlap$name)

###### Let's figure out if there's replacement or not.  There is not really much replacement.

blasthits3 <- blasthits2
blasthits2 <- blasthits3
blasthits2$replacement <- FALSE
j<-1
blasthits2 <- blasthits2[order(blasthits2$qseqid),]
for (i in c(1:length(blasthits2$qend))){
  j<-i+1
  while (blasthits2$qseqid[i]==blasthits2$qseqid[j]) {
    if ((inrange(blasthits2$sstart[i], blasthits2$sstart[j], blasthits2$send[j], incbounds=TRUE)) |
        (inrange(blasthits2$sstart[i], blasthits2$send[j], blasthits2$sstart[j], incbounds=TRUE))){
        blasthits2$replacement[i] <- TRUE
        blasthits2$replacement[j] <- TRUE
  }
  j<-j+1
  }
}
length(unique(blasthits2$qseqid))
replacementset <- blasthits2[blasthits2$replacement==TRUE,]

######

blasthits2
blasthits2$replacement <- FALSE
j<-1
blasthits2 <- blasthits2[order(blasthits2$qseqid),]
for (i in c(1:length(blasthits2$qend))){
  j<-i+1
  while (blasthits2$qseqid[i]==blasthits2$qseqid[j]) {
    if ((inrange(blasthits2$sstart[i], blasthits2$sstart[j], blasthits2$send[j], incbounds=TRUE)) |
        (inrange(blasthits2$sstart[i], blasthits2$send[j], blasthits2$sstart[j], incbounds=TRUE))){
      blasthits2$replacement[i] <- TRUE
      blasthits2$replacement[j] <- TRUE
    }
    j<-j+1
  }
}

#########

blasthits2 <- blasthits_100 %>% dplyr::group_by(qseqid) %>% dplyr::mutate(rank=dense_rank(qstart))
variableRegions_order1 <- variableRegions
variableRegions_order2 <- variableRegions
variableRegions_order3 <- variableRegions
variableRegions_order1$rank <- 1
variableRegions_order2$rank <- 2
variableRegions_order3$rank <- 3
variableRegions_order <- rbind(variableRegions_order1,variableRegions_order2,variableRegions_order3)
variableRegions_order$rankcount <- 0
variableRegions_order$rankpercentage <- 0
for (i in c(1:length(blasthits2$qend))){
  for (j in c(1:length(variableRegions$start))){
    if (blasthits2$region[i] != as.character(variableRegions$region[j])) next()
    if ((inrange(blasthits2$sstart[i], variableRegions$start[j], variableRegions$end[j], incbounds=TRUE)) |
        (inrange(blasthits2$sstart[i], variableRegions$end[j],variableRegions$start[j], incbounds=TRUE))){
      variable_match <- variableRegions[j,]
      myindex <- j
      print(i)
      print(paste0("start ",blasthits2$sstart[i]))
      print(paste0("var start ",variableRegions$start[j]))
      print(paste0("var end ",variableRegions$end[j]))
    }
  }
  if (blasthits2$rank[i]==1){
    variableRegions_order$firstcount[myindex] <- variableRegions_order$firstcount[myindex] + variable_match$count
    variableRegions_order$firstpercentage[myindex] <- variableRegions_order$firstpercentage[myindex] + variable_match$percentage
  }
  if (blasthits2$rank[i]==2){
    variableRegions_order$secondcount[myindex] <- variableRegions_order$secondcount[myindex] + variable_match$count
    variableRegions_order$secondpercentage[myindex] <- variableRegions_order$secondpercentage[myindex] + variable_match$percentage
  }
  if (blasthits2$rank[i]==3){
    variableRegions_order$thirdcount[myindex] <- variableRegions_order$thirdcount[myindex] + variable_match$count
    variableRegions_order$thirdpercentage[myindex] <- variableRegions_order$thirdpercentage[myindex] + variable_match$percentage
  }
}

ggplot(variableRegions_order,aes(x=start,y=firstpercentage,fill=region)) + geom_bar(stat="identity",position = "dodge") + 
  scale_fill_manual(values=colorscheme) + theme_minimal() + labs(fill='Region') + 
  xlab("Sequence-Alignment Length Difference") + ylab("Count")
