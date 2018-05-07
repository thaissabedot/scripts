library(GenomicRanges)
library(org.Mm.eg.db)
k <- keys(org.Mm.eg.db, keytype="SYMBOL")
gene <- select(org.Mm.eg.db, keys = k, columns = c("CHRLOC", "CHRLOCEND"), keytype="SYMBOL") #seleciona as colunas cromossomo, inicio, final e symbol para cada gene conhecido de rato
gene <- na.omit(gene) #retirar os genes que nao temos informação das coordenadas genomicas

#Identifica em qual fita do DNA os genes estao localizados
strand <- regmatches(gene$CHRLOC, regexpr("^[^[:digit:]]*", gene$CHRLOC))
for (i in 1:length(strand)){
    if(strand[i] == "")
        strand[i] <- "+"
}
gene$strand <- as.factor(strand)
gene$CHRLOC <- abs(gene$CHRLOC) #retira o 'sinal de menos' da frente das coordenadas genomicas
gene$CHRLOCEND <- abs(gene$CHRLOCEND)
gene <- subset(gene, CHRLOCCHR %in% c(1:19,"X","Y") #Tirar os cromossomos tipo "random"
Mm.genes <- GRanges(seqnames = paste("chr", as.character(gene$CHRLOCCHR), sep=""), ranges = IRanges(start = gene$CHRLOC, end = gene$CHRLOCEND), strand = gene$strand, symbol = gene$SYMBOL) #Transformando o objeto 'gene' em um objeto do tipo GRanges


Act_AD2 <- read.delim("CsiBones_Act_AD2.rmdup.q20.noinput_peaks.txt")
Act_AD3 <- read.delim("CsiBones_Act_AD3.rmdup.q20.noinput_peaks.txt")
Act_D24.2 <- read.delim("CsiBones_Act_D24_2.rmdup.q20.noinput_peaks.txt")
Act_D24.3 <- read.delim("CsiBones_Act_D24_3.rmdup.q20.noinput_peaks.txt")
Dnm_AD1 <- read.delim("CsiBones_Dnm_AD1.rmdup.q20.noinput_peaks.txt")
Dnm_AD2 <- read.delim("CsiBones_Dnm_AD2.rmdup.q20.noinput_peaks.txt")
Dnm_D24.1 <- read.delim("CsiBones_Dnm_D24_1.rmdup.q20.noinput_peaks.txt")
Dnm_D24.2 <- read.delim("CsiBones_Dnm_D24_2.rmdup.q20.noinput_peaks.txt")

### Exemplo para o primeiro (Act_AD2) .. para fazer os outros, só copiar essas 4 próximas linhas e ir mudando o nome dos objetos
Act_AD2.gr <- GRanges(seqnames = as.character(Act_AD2$chr), ranges = IRanges(start = Act_AD2$start, end = Act_AD2$end), peak = as.character(Act_AD2$name)) 
vizinho.act_ad2 <- nearest(Act_AD2.gr, Mm.genes, ignore.strand=TRUE) #identifica o gene mais proximo de cada pico 
Act_AD2$gene_mais_prox <- as.character(gene[vizinho.act_ad2,"SYMBOL"]) #adiciona o nome do gene mais do proximo ao pico correspondente
head(Act_AD2) #Visualizar o resultado

#save(Act_AD2,Act_AD3, Act_D24.2, Act_D24.3, Dnm_AD1, Dnm_AD2, Dnm_D24.1, Dnm_D24.2, file="Genes_Picos.rda") #Para salvar o resultado


