library(taxize)
library(rentrez)
library(Biostrings)
library(UpSetR)
library(drawProteins)
library(dplyr)
library(ggplot2)
library(biomaRt)
library(tidyverse)

df <- read.csv("2_Species_class.csv")

species <- df$Species
df1 <- tax_name(sci = species, get = c("family", "genus"), db = "ncbi")
family <- df1$family
genus <- df1$genus
df_new <- cbind(df, genus, family)
write.csv(df_new, "3_species_class_tax.csv", row.names = FALSE)
df <- read.csv("1_R_genes_all_classes.csv")
class <- c("RLK", "LECRK", "LYK", "LYP", "RLP", "CNL", "TNL", "Other", "N", "CN", "NL", "T")

for (i in 1:length(class)) {
  tmp <- read.csv("1_R_genes_all_classes.csv")
  tmp <- subset(tmp, tmp$R_class==class[i])
  out_file <- paste0(class[i], ".csv")
  write.csv(tmp, out_file, row.names = FALSE)
} 
df <- read.csv("R_class_Other.csv")
nrow(df)  ## 13 rows
df1 <- subset(df,df$Uniprot != "-" | df$Genbank_accession != "-")
nrow(df1)
write.csv(df1, "R_class_Other.csv", row.names = FALSE)
df <- read.csv("R_class_LYK.csv")
accession <- df$Genbank_accession
seq <- entrez_fetch(db="protein", id=accession, rettype="fasta")
write(seq, "LYK.fasta")
fastafiles <- list.files(pattern="*.fasta")
fastafiles

for (i in 1:length(fastafiles)) {
  file.copy(from = paste0("C:/Users/vn81649/Documents/SNAP_Bean/R_genes/1_R_Classes/uniprot/",fastafiles[i]),
            to = "C:/Users/vn81649/Documents/SNAP_Bean/R_genes/1_R_Classes/fasta")
  file.remove(from=paste0("C:/Users/vn81649/Documents/SNAP_Bean/R_genes/1_R_Classes/uniprot/",fastafiles[i]))
} 

file.copy(from = "C:/Users/vn81649/Documents/SNAP_Bean/R_genes/1_R_Classes/ncbi/LYK.fasta",
          to = "C:/Users/vn81649/Documents/SNAP_Bean/R_genes/1_R_Classes/fasta")
file.remove(from="C:/Users/vn81649/Documents/SNAP_Bean/R_genes/1_R_Classes/ncbi/LYK.fasta")
for (i in 1:length(fastafiles)) {
  fastaFile <- readDNAStringSet(fastafiles[i], format="fasta")
  seq_name = names(fastaFile)
  sequence = paste(fastaFile)
  df <- data.frame(seq_name, sequence)
  out_file <- gsub(".fasta", "_fasta.csv", fastafiles[i])
  write.csv(df, out_file, row.names = FALSE)
}

setwd("C:/Users/vn81649/Documents/SNAP_Bean/R_genes/1_R_Classes/fasta")
fastafiles <- list.files(pattern="*.fasta")

df <- read.csv("Class_IPRDomain.csv")
## Select class names
class <- colnames(df)[2:13]
upset(df, sets = class, sets.bar.color = "#56B4E9",
      order.by = "freq", text.scale = 1.5)

upset(df, sets = class, sets.bar.color = "#56B4E9",
      order.by = "freq", empty.intersections = "on")

files <- list.files(pattern="*.csv")
files

classnames <- gsub(".csv", "", files) 
classnames

for (i in 1:length(files)) {
  df <- read.csv(files[i])
  df1 <- subset(df,df$Uniprot != "-")
  outfile <- paste0("UID_", classnames[i], ".csv")
  write.csv(df1, outfile, row.names = FALSE)
}

files <- list.files(pattern="UID_.*")
files

for (i in 1:length(files)) {
  df <- read.csv(files[i])
  ID <- df$Uniprot
  ID <- paste(ID,collapse = " ")
  drawProteins::get_features(ID) -> prot_json
  drawProteins::feature_to_dataframe(prot_json) -> prot_data
  outfile <- paste0("Prot_Sch_", classnames[i], ".csv")
  write.csv(prot_data, outfile, row.names = FALSE)
}

files <- list.files(pattern="Prot_Sch_.*")
files

for (i in 1:length(files)) {
  prot_data <- read.csv(files[i])
  
  tiff(file=paste0("PS_", classnames[i], ".tiff", sep=""))
  p <- draw_canvas(prot_data)
  p <- draw_chains(p, prot_data)
  p <- draw_domains(p, prot_data)
  p <- draw_regions(p, prot_data)
  p <- draw_repeat(p, prot_data)
  p <- draw_motif(p, prot_data)
  
  p <- p + theme_bw(base_size = 12) + # white background & change text size
    theme(panel.grid.minor=element_blank())+
    theme(axis.ticks = element_blank(),
          axis.text.y = element_blank()) +
    theme(panel.border = element_blank())
  print(p)
  dev.off()
}

ensemblp <- useMart(biomart="plants_mart",host="https://plants.ensembl.org")
Pl_Datasets <- listDatasets(ensemblp)
write.csv(Pl_Datasets, "Pl_Datasets.csv", row.names=FALSE) 

IDs <- c("IPR000157", "IPR000858", "IPR001220", 
         "IPR001480", "IPR002182", "IPR008808", 
         "IPR009743", "IPR019825", "IPR024788", 
         "IPR018392", "IPR038005", "IPR041118", 
         "IPR001611", "IPR011713", "IPR013210", 
         "IPR025875")

pv = useDataset("pvulgaris_eg_gene", mart = ensemblp)
attributes <- listAttributes(pv)
write.csv(attributes, "Attributes.csv", row.names = FALSE)


df <- getBM(attributes =c("ensembl_peptide_id", "uniprotswissprot", "interpro", "peptide", "refseq_peptide", "ensembl_gene_id"),
            filters = ("interpro"), 
            values = IDs, 
            mart = pv)
write.csv(df, "R_genes_pv.csv", row.names = FALSE)
df <- read.csv("R_genes_pv.csv")
df1 <- df[!duplicated(df[2]),]  
df1$peptide <- gsub("\\*", "", as.character(df1$peptide)) 
write.csv(df1, "R_genes_U_pv.csv", row.names = FALSE)
fasta <- read.csv("R_genes_U_pv.csv")
fasta.fas<-paste0(">",fasta$ensembl_peptide_id, 
                  "\n",fasta$peptide,"\n")
writeLines(fasta.fas,"R_genes_U_pv.fasta")

df <- read.csv("33_domain_set.csv")
domain <- df$Domain.Sets

for (i in 1:length(domain)){
  df <- read.csv("R_genes_U_pv.csv")
  df1 <- df[!duplicated(df[2]),]  ## Unique sequences only based on column 2
  tmp<-subset(df1,df1$interpro==domain[i])
  out_file <- paste0("Pv_", domain[i], ".csv")
  write.csv(tmp, out_file, row.names=FALSE)
}  

filenames <- list.files(pattern="Pv*") 
filenames
domain_set <- gsub(".csv", "", filenames) 
domain_set
domain_set_1 <- gsub("Pv_", "", domain_set) 
domain_set_1
domains <- gsub(";", "_", domain_set_1) 
domains

for (i in 1:length(filenames)) {
  fasta <- read.csv(filenames[i])
  fasta.fas<-paste0(">",fasta$ensembl_peptide_id, 
                    "\n",fasta$peptide,"\n")
  out_file <- paste0(domains[i], ".fasta")
  writeLines(fasta.fas, out_file)
}
for (i in 1:ncol(df)){
  ID <- df[[i]][df[[i]]!=""]
  pv = useDataset("pvulgaris_eg_gene", mart = ensemblp)
  tmp <- getBM(attributes =c("ensembl_peptide_id", "interpro", "peptide", "refseq_peptide"),
               filters = ("interpro"), 
               values = ID, 
               mart = pv)   
  out_file <- paste0(R_class[i], ".csv")
  write.csv(tmp, out_file, row.names = FALSE)
}

species <- c("G. max", "L. angustifolius", "M. truncatula", "P. sativum", "P. vulgaris", "T. pratense", "V. angularis", "V. radiata", "V. unguiculata")
for (i in 1:length(species)) {
  nlr_data <- read_csv(paste0(species[i], "_nlr.csv"))
  chr_length <- read_csv(paste0(species[i], "_chr.csv"))
  nlr_density <- nlr_data %>%
    distinct(ensembl_peptide_id, .keep_all = TRUE) %>%  # remove duplicate genes
    group_by(chromosome_name, bin = floor((start_position - 1) / 1000000) + 1) %>%
    summarize(nlr_density = n() / 1000, .groups = "keep") %>%
    ungroup() %>%
    left_join(chr_length, by = "chromosome_name") %>%
    mutate(chromosome_length_mb = chromosome_length / 1000000)
  plot_title <- bquote(italic(.(species[i])))  
  plot_file <- paste("NLR_density_", gsub(" ", "_", species[i]), ".tiff", sep="")
  
  nlr_density %>%
    ggplot(aes(x = bin, y = as.numeric(chromosome_name), fill = nlr_density)) +
    geom_tile() +
    geom_rect(aes(xmin = 0, xmax = chromosome_length_mb + 1, 
                  ymin = as.numeric(chromosome_name) - 0.5, ymax = as.numeric(chromosome_name) + 0.5),
              fill = "transparent", color = "black", size = 0.1) +
    scale_fill_gradient(low = "white", high = "red") +
    theme(plot.title = element_text(hjust = 0.5, size = 8),
          legend.position = "right",
          legend.direction = "vertical",
          legend.text = element_text(size = 6), # set legend text size
          legend.title = element_text(size = 6), # set legend title size
          axis.line = element_line(color = "black", size = 0.2),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          axis.title.x = element_text(size = 6),
          axis.title.y = element_text(size = 6),
          panel.background = element_rect(fill = "white"))+
    xlab("1 MB bins along the chromosome") +
    ylab("Chromosomes") +
    ggtitle(plot_title) +
    scale_y_continuous(breaks = seq(1, max(nlr_density$chromosome_name), 1),
                       labels = seq(1, max(nlr_density$chromosome_name), 1))
  ggsave(plot_file, width = 2.5, height = 2.5, dpi = 300)
}



