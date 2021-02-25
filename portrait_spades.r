### to do - выбирать на бласт не края, но серединку контига!


library("seqinr")
library("dbscan")
library("data.table")

args <- commandArgs(trailingOnly = TRUE)

assembly_dir = args[1]

blast_length = 10000
blast_cpus=40
blast_path="/srv/common/bin/"
blastdb_path="/data2/bio/blast_shared_storage/"
each_cluster_blast_count=10
dbscan_eps = 0.05


asm =read.fasta(paste0(assembly_dir,"/contigs.fasta"))
names = names(asm)
con <- textConnection(as.character(names))
names.tab <- read.table(con, sep="_")
close(con)

contig.name = paste0(names.tab$V1,"_",names.tab$V2) 
contig.len  = names.tab$V4
contig.cov  = names.tab$V6
contig.gc   = sapply(asm, GC)

T = data.frame(gc = contig.gc, cov = (contig.cov))
pdf(paste0(assembly_dir,"_portrait.pdf"))
plot(T$gc, T$cov,pch=16, col=rgb(0.2, 0.2, 0.2,0.4), cex.lab=1.5, ylab="depth", xlab="GC", log="y")

# with contig names
plot(T$gc, T$cov,pch=16, type='n', cex.lab=1.5,  ylab="depth", xlab="GC", log="y")
text(T$gc, T$cov, contig.name, cex = 0.4)

### cluster
C =dbscan(T,dbscan_eps)
plot(T$gc, T$cov,pch=16,col=C$cluster)
clusters.count = max(C$cluster)
i.to.blast=c()
for(i in 0:clusters.count){
  indexes=which(C$cluster==i)
  to.blast=sample(indexes, min(each_cluster_blast_count, length(indexes)))
  i.to.blast=c(i.to.blast, to.blast)
  }

### blast them

file.for.blast="portrait_tmp_for_blast.fasta"
write.fasta(sequences=lapply(asm[i.to.blast], function(x) x[seq(from = 1, to = min(length(x),blast_length))]), names=names(asm[i.to.blast]), file.out=file.for.blast)

file.blast.out="portrait_tmp_blast.out"
blast.cmd = paste0("export BLASTDB=", blastdb_path,"; ",blast_path,"/blastn -db ", blastdb_path,  "/nt -query ", file.for.blast, " -outfmt \"6 qseqid sseqid slen evalue bitscore qcovs sscinames scomnames stitle\" -num_threads ", blast_cpus," -evalue 1e-10 -out ", file.blast.out)
system(blast.cmd)

blast_result <- fread(file.blast.out, sep="\t",col.names=c("contig", "subject", "alignment_length", "e-value", "bit_score", "qcov","scinames", "scomnames", "stitle"))
blast_result <- blast_result[, .SD[which(qcov==max(qcov))], by = contig] 
blast_result <- blast_result[, total_score:=sum(bit_score), by = c('subject', 'contig')]
blast_result <- blast_result[, .SD[which.max(total_score)], by = contig]
write.table(blast_result, paste0(assembly_dir,"_blast.txt"), sep='\t',quote=FALSE)


### plot blast results
T$name=rownames(T)
B = merge(blast_result,T,by.x="contig",by.y="name")

plot(T$gc, T$cov,pch=16, col=rgb(0, 0.4, 0.4,0.1), cex.lab=1.5, ylab="depth", xlab="GC", log="y", main = assembly_dir)
text(B$gc, B$cov, B$scomnames, cex = 0.4)

dev.off()

png(paste0(assembly_dir,"_portrait.png"), width=1200, height=800)
plot(T$gc, T$cov,pch=16, col=rgb(0, 0.4, 0.4,0.2), cex.lab=1.5,  ylab="depth", xlab="GC", log="y", main = assembly_dir)
text(B$gc, B$cov, B$scomnames, cex = 0.9)
dev.off()
