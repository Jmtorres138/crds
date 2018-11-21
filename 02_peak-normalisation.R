"%&%" <- function(a,b) paste0(a,b)
library("dplyr")
library("data.table")

count.dir <- "/well/mccarthy/users/jason/projects/atac_analyses/evaluate_peaks/eLife2018/"
print(count.dir)

out.file1 <- count.dir %&% "combined-counts.txt"
out.file2 <- count.dir %&% "combined-cpm.txt"
out.file3 <- count.dir %&% "combined-counts_filtered.txt"
out.file4 <- count.dir %&% "combined-cpm_filtered.txt"

count.files <- list.files(count.dir)
count.files <- count.files[grepl("merged_peaks.counts-",count.files)]
count.files <- count.files[!(grepl("summary",count.files))]

f <- count.files[1]
samp <- (strsplit(f,"merged_peaks.counts-")[[1]][2] %>% strsplit(.,".txt"))[[1]][1]
fpath <- count.dir %&% f
out.df <- fread(fpath)
names(out.df)[7] <- samp

for (f in count.files[2:length(count.files)]){
  samp <- (strsplit(f,"merged_peaks.counts-")[[1]][2] %>% strsplit(.,".txt"))[[1]][1]
  fpath <- count.dir %&% f
  df <- fread(fpath)
  names(df)[7] <- samp
  build.df <- df[,c(1,7)]
  out.df <- inner_join(out.df,build.df,by="Geneid")
}

write.table(x=out.df,file=out.file1,quote=F,sep="\t",row.names=F,col.names=T)


# CPM normalization
cpm.df <- out.df
for (i in 7:dim(out.df)[2]){
  vec <- (out.df[,i] / sum(out.df[,i])) * 1e6
  cpm.df[,i] <- vec
}

write.table(x=cpm.df,file=out.file2,quote=F,sep="\t",row.names=F,col.names=T)


# CPM thresholding;
# https://github.com/broadinstitute/gtex-pipeline/blob/master/qtl/README.md
# Using CPM rather than RPKM
#Genes were selected based on the following exression thresholds:
#>0.1 RPKM in ≥10 samples AND
#≥6 reads (unnormalized) in ≥10 samples

id.vec <- out.df$Geneid
keep.vec <- c()

pb <- txtProgressBar(min=0,max=length(id.vec),style=3)
for (i in 1:length(id.vec)){
  setTxtProgressBar(pb,i)
  id <- id.vec[i]
  eval1 <- ((filter(out.df,Geneid==id)[7:dim(out.df)[2]] >= 6) %>% sum(.)) >= 10
  eval2 <- ((filter(cpm.df,Geneid==id)[7:dim(out.df)[2]] >= 0.1) %>% sum(.)) >= 10
  if (eval1 == TRUE & eval2 == TRUE){
    keep.vec <- append(keep.vec,id)
  }
}

write.table(x=filter(out.df,Geneid%in%keep.vec),file=out.file3,quote=F,sep="\t",row.names=F,col.names=T)
write.table(x=filter(cpm.df,Geneid%in%keep.vec),file=out.file4,quote=F,sep="\t",row.names=F,col.names=T)



#keep.index1 <- apply(out.df,1,function(vec){((vec[7:length(vec)] %>% as.integer(.) >= 6) %>% sum(.))>=10})
#keep.vec1 <- id.vec[keep.index1]

#keep.index2 <- apply(cpm.df,1,function(vec){((vec[7:length(vec)] %>% as.integer(.) >= 0.1) %>% sum(.))>=10})
#keep.vec2 <- id.vec[keep.index2]

#keep.vec <- set.intersect(keep.vec1,keep.vec2) %>% as.character(.)
