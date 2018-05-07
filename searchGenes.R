source("http://bioconductor.org/biocLite.R")
biocLite("biomaRt")

# Dados de HipoMetilação.
### Reuni as informações direto da base de dados de genes.
ensembl=useMart("ensembl")

### Realiza o download dos dados referentes ao genoma solicitado.
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

### Cria um vetor com os grupos em que os cromossomos estão divididos.
chrom <- c(1:22, "M", "X","Y")

### Organiza as informações em um data frame.
gene.location <- getBM(attributes = c("chromosome_name",
                                      "start_position",
                                      "end_position", "strand",
                                      "external_gene_name",
                                      "entrezgene"),
                       filters = c("chromosome_name"),
                       values = list(chrom), mart = ensembl)