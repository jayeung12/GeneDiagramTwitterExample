#This script visualizes nextclade .tsv data with detailed gene annotations for nucleotide and AA changes.

library(tidyverse)
library(gggenes)
library(ggpubr)
library(ggrepel)
library(Biostrings)
library(stringi)

refName="TVP23352_1_BA_2"
queryname = "ba2_last2weeks"
#Change file path to the dictory with file
genes = read.delim(paste0("C:/Users/jasonyeung/Documents/Projects/Shi Lab/BA1vsBA2/nextcladedrop/",refName,"genemap.gff"), header=F)
nsp = read.delim(paste0("C:/Users/jasonyeung/Documents/Projects/Shi Lab/BA1vsBA2/nextcladedrop/",refName,"ORF1abgenemap.gff"), header=F)

#Order by start
genes = genes[order(genes$V4),]
#Clean up gene/domain labels
genes$V9 = str_remove(genes$V9 %>% str_trim(),"gene_name=")

#Specifies what row to plot genes on
#Need to write code to make sure no overlaps if automating
geneOrder = c("row3","row3","row1","row3","row1","row1","row3","row3","row3","row3","row1","row3")

#Colors for the different genes
colours = c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9","#BC80BD","#CCEBC5","#FFED6F")
ggDF=cbind(molecule = geneOrder, gene = genes$V9, start = genes$V4, end =genes$V5, strand = "forward", orientation = 1, colour = colours) %>% as.data.frame()

#Make a separate df for nsps and put them on the same row
ggDF1=cbind(molecule = "row2", gene = nsp$V9, start = nsp$V4, end = nsp$V5, strand = "foward", orientation = 1, colour = "#ADD8E6")%>% as.data.frame()
ggDF = rbind(ggDF, ggDF1)
ggDF$start = ggDF$start %>% as.numeric()
ggDF$end = ggDF$end %>% as.numeric()

breaks = c(ggDF$colour[which(ggDF$molecule=="row3")],"#ADD8E6",ggDF$colour[which(ggDF$molecule=="row1")])
glabels = c(ggDF$gene[which(ggDF$molecule=="row3")],"nsp",ggDF$gene[which(ggDF$molecule=="row1")])

#Make the gene diagram
p2=ggplot(ggDF, aes(xmin = start %>% as.numeric(), xmax = end %>% as.numeric(), y = molecule, fill = colour, label=gene)) +
  geom_gene_arrow() + xlim(0,30000)+
  geom_gene_label(align = "left") +scale_fill_identity(guide = 'legend', breaks =breaks, labels = glabels)+
  theme_genes() + ylab("")+ labs(fill="") +
  theme(plot.margin = unit(c(-2,0,0,0), "lines"),axis.text.y =element_blank(),axis.text.x =element_blank(), axis.ticks.x  = element_blank())

#Nucleotide
tsv=data.table::fread(file = paste0("C:/Users/jasonyeung/Documents/Projects/Shi Lab/BA1vsBA2/nextcladedrop/TVP23352vs23328/","nextclade.tsv"))
subs = tsv$substitutions %>% str_split(pattern = ",") %>% unlist() %>% gsub(pattern="[A-Z]", replacement = "") %>% as.numeric %>% table()
freq=subs/nrow(tsv)
names(subs)
sublabs=tsv$substitutions %>% str_split(pattern = ",") %>% unlist()
labpoints = freq[which(freq*100>5)]
labelsNuc=sapply(1:length(labpoints), function(i){  grep(sublabs,pattern =paste0("[A-Z]",names(labpoints)[i],"[A-Z]"), value = T) %>% unique() %>% paste(collapse = ", ")})


#Make lollipop plot 
p1=ggplot() + 
  geom_segment(aes(x = names(subs) %>% as.numeric(), xend = names(subs) %>% as.numeric(), y = 0, yend = freq*100)) +
  geom_point(aes( x = names(subs) %>% as.numeric(), y = freq*100), colour = "black", size = 1) +xlim(0,30000)+
  xlab("")+ theme_classic()+ylab("Frequency (%)")+
  geom_text_repel(size=2.5, point.padding = 0.2,
                  nudge_x = .15,
                  nudge_y = .5,
                  aes(x=labpoints %>% names() %>% as.numeric(),y=labpoints*100,label= labelsNuc))


if(nrow(tsv)==1)
{
  
  freq = runif(subs %>% length(), min=.10, max=.90)
  
  
  p1=ggplot() + 
    geom_segment(aes(x = names(subs) %>% as.numeric(), xend = names(subs) %>% as.numeric(), y = 0, yend = freq*100)) +
    geom_point(aes( x = names(subs) %>% as.numeric(), y = freq*100), colour = "black", size = 1) +xlim(0,30000)+
    xlab("")+ theme_classic()+ylab("")+
    theme(axis.line.y = element_blank(),axis.text.y =element_blank(), axis.ticks.y  = element_blank())+
    geom_text_repel(size=2, max.overlaps = Inf,
                    nudge_x = .15,
                    nudge_y = .5,
                    aes(x=labpoints %>% names() %>% as.numeric(),y=freq*100,label= labelsNuc))
  
}

#Plotting the lollipop plot and gene diagram together
ggarrange(p1,p2,ncol=1, nrow = 2,heights=c(9,3),align = "v", labels = c("",""),common.legend = TRUE, legend="bottom")



#AA only
aasubs = tsv$aaSubstitutions %>% str_split(pattern = ",") %>% unlist()
subgene=aasubs%>% str_extract(pattern = "[^:]+")
subloc=aasubs%>% str_extract(pattern = "(?<=:).*")%>% gsub(pattern="[A-Z]", replacement = "")%>% gsub(pattern="\\*", replacement = "") %>% as.numeric
aasubsDF = cbind(aasubs,subgene,subloc) %>%  as.data.frame()
aasubsDF$subloc=aasubsDF$subloc %>% as.numeric()
aasubsDF$nucloc = NA

for(indexGene in unique(subgene))
{
  start = genes$V4[which(genes$V9==indexGene)]
  
  for(indexLoc in which(aasubsDF$subgene == indexGene))
  {
    aasubsDF$nucloc[indexLoc] = start + aasubsDF$subloc[indexLoc]*3
  }
}

freq=aasubsDF$nucloc %>% table()/nrow(tsv)
labpoints = freq[which(freq*100>5)]
labelsAA=sapply(labpoints %>% names, function(i){  aasubsDF$aasubs[which(aasubsDF$nucloc == i)] %>% unique() %>% paste(collapse = ", ")})



p1.1=ggplot() + 
  geom_segment(aes(x = aasubsDF$nucloc %>% table() %>% names() %>% as.numeric(), xend = aasubsDF$nucloc %>% table() %>% names() %>% as.numeric(), y = 0, yend = freq*100)) +
  geom_point(aes( x = aasubsDF$nucloc %>% table() %>% names() %>% as.numeric(), y = freq*100), colour = "black", size = 1) +xlim(0,30000)+
  xlab("")+ theme_classic()+ylab("Frequency (%)")+
  geom_text_repel(size=2.5, point.padding = 0.2,
                  nudge_x = .15,
                  nudge_y = .5,
                  aes(x=labpoints %>% names() %>% as.numeric(),y=labpoints*100,label= labelsAA))


if(nrow(tsv)==1)
{
  
  
  freq = runif(aasubsDF$nucloc %>% length(), min=.20, max=.50)
  
  
  p1.1=ggplot() + 
    geom_segment(aes(x = aasubsDF$nucloc %>% table() %>% names() %>% as.numeric(), xend = aasubsDF$nucloc %>% table() %>% names() %>% as.numeric(), y = 0, yend = freq*100)) +
    geom_point(aes( x = aasubsDF$nucloc %>% table() %>% names() %>% as.numeric(), y = freq*100), colour = "black", size = 1) +xlim(0,30000)+
    xlab("")+ theme_classic()+ylab("")+
    theme(axis.line.y = element_blank(),axis.text.y =element_blank(), axis.ticks.y  = element_blank())+
    geom_text_repel(size=2, max.overlaps = Inf,
                    nudge_x = .15,
                    nudge_y = .5,
                    aes(x=labpoints %>% names() %>% as.numeric(),y=freq*100,label= labelsAA))
  
}

#Plotting the lollipop plot and gene diagram together
ggarrange(p1.1,p2,ncol=1, nrow = 2,heights=c(9,3),align = "v", labels = c("A",""),common.legend = TRUE, legend="bottom")

