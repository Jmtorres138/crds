"%&%" <- function(a,b) paste0(a,b)
library("dplyr")
library("data.table")

count.dir <- "/well/mccarthy/users/jason/projects/crds/peak_counts/"
print(count.dir)

out.file1 <- count.dir %&% "combined-counts.txt"
out.file2 <- count.dir %&% "combined-cpm.txt"

count.files <- list.files(count.dir)
count.files <- count.files[grepl("atac_peaks.counts-",count.files)]
count.files <- count.files[!(grepl("summary",count.files))]

f <- count.files[1]
samp <- (strsplit(f,"atac_peaks.counts-")[[1]][2] %>% strsplit(.,".txt"))[[1]][1]
fpath <- count.dir %&% f
out.df <- fread(fpath)
names(out.df)[7] <- samp

for (f in count.files[2:length(count.files)]){
  samp <- (strsplit(f,"atac_peaks.counts-")[[1]][2] %>% strsplit(.,".txt"))[[1]][1]
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

