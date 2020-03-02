#!/usr/bin/R

#Initialize with libraries
library(rBLAST)
library(ShortRead)
library(ggplot2)
library(taxonomizr)
library(dplyr)
library(forcats)
require(scales)


#set paths
Sys.setenv(PATH=paste(Sys.getenv("PATH"), "/share/pkg.7/blast+/2.7.1/install/bin", sep=":"))
Sys.setenv(PATH=paste(Sys.getenv("PATH"), "/share/pkg.7/sratoolkit/2.9.2/install/bin/", sep=":"))

#set up fastq-dump and accession numbers
#srr=c('SRS6120502') #largest dataset
srr=c('SRS6120496') #largest dataset
system(paste('fastq-dump', srr, sep=' ')) #load fastq

dna = readFastq('.', pattern=srr) #read DNA sequences
reads = sread(dna, id=id(dna)) # parse DNA sequences
qscores = quality(dna) # parse quality scores

widths = as.data.frame(reads@ranges@width)
numqscores = as(qscores, "matrix") # converts to numeric scores automatically
avgscores = apply(numqscores, 1, mean, na.rm=TRUE) #apply the function mean() across all rows (argument '1') of the matrix and ignore "NA" values
avgscores = as.data.frame(avgscores)
mdata = cbind(widths, avgscores)

##BLAST
bl <- blast(db="/projectnb/ct-shbioinf/blast/nt.fa") #access local BLAST database
cl <- predict(bl, reads, BLAST_args = '-num_threads 12 -evalue 1e-100')

## post-BLAST analysis:
accid = as.character(cl$SubjectID) # accession IDs of BLAST hits

#load the taxonomy database files 'nodes' and 'names'
taxaNodes<-read.nodes.sql("/projectnb/ct-shbioinf/taxonomy/data/nodes.dmp")
taxaNames<-read.names.sql("/projectnb/ct-shbioinf/taxonomy/data/names.dmp")
# Search the taxonomy by accession ID #

# takes accession number and gets the taxonomic ID
ids<-accessionToTaxa(accid,'/projectnb/ct-shbioinf/taxonomy/data/accessionTaxa.sql')

# displays the taxonomic names from each ID #
taxlist=getTaxonomy(ids, taxaNodes, taxaNames)

cltax=cbind(cl,taxlist) # merge taxonomy and blast

cltop = cltax %>% 
  group_by(QueryID) %>% 
  top_n(n=1, wt=Bits)

# plots and figures
(ggplot(data=cltop) +
    geom_bar(aes(x=fct_infreq(genus))) +
    theme_minimal() +
    theme(    
      axis.text.x  = element_text(angle = 45, hjust=1)
    ) +
    xlab('')
  
)


my_breaks = c(1,2,4,8,16,32,64,128,256,512,1024)

ggplot(mdata,aes(x=reads@ranges@width,y=avgscores)) +
  geom_bin2d(bins = 48) +
  scale_fill_gradient(low="lightblue2",high="darkblue", name = "count", trans = "log",
                      breaks = my_breaks, labels = my_breaks) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
  theme_linedraw() +
  ggtitle('2D Density') +
  labs( y='Average Scores', x=expression('Read Length (log'[10]* ' bp)'))

###Fix the following:
sele.hits = select(blasthits, QueryID, SubjectID, Perc.Ident, Alignment.Length, E, Bits, family, genus)
#colnames(sele.hits)
top.hits = filter(blasthits, Bits>350)
#nrow(top.hits)
#nrow(blasthits)

groupgen = group_by(blasthits, genus) #ignore error
topgen = summarize(groupgen, count = n())
head(arrange(topgen, desc(count))) # PRINT JUST THE HIGHEST COUNTS

groupbyread = group_by(blasthits, QueryID)
topbygroup = filter(groupbyread, any(Bits>350))
head(arrange(topbygroup, desc(QueryID))) # PRINT JUST THE HIGHEST COUNTS

topgenera = blasthits %>% 
  group_by(QueryID) %>% 
  top_n(n=1, wt=(-E)) %>% #top by E value this time
  group_by(genus) %>% 
  summarize(count=n()) %>% 
  arrange(desc(count))
head(arrange(topgenera, desc(QueryID))) # PRINT JUST THE HIGHEST COUNTS

topgenera = blasthits %>% 
  group_by(QueryID) %>% 
  filter(any(Bits>500)) %>%
  group_by(genus) %>% 
  summarize(count=n()) %>% 
  arrange(desc(count))


