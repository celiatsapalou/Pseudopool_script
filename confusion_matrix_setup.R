library(data.table)
library(dplyr)
library(ggplot2)

#############################
inversions<-read.table("/home/celia/Downloads/Mosaicatcher/variants_callset_inversions.tsv", header = TRUE, sep = '\t')
inversions<- subset(inversions, select = -c(misorient, CALLERSET_LIST, lowConfCount, noreadsCount, categ) )

inv_mosaicatcher<-read.table("/home/celia/Downloads/Mosaicatcher/lenient_filterFALSE.tsv", header = TRUE, sep = '\t')
inv_mosaicatcher$cell=sub("HG02059102.*", "HG02059x102", inv_mosaicatcher$cell)
inv_mosaicatcher$cell=sub("_.*", "", inv_mosaicatcher$cell)
inv_mosaicatcher$cell=sub("x.*", "", inv_mosaicatcher$cell)
inv_mosaic <-inv_mosaicatcher[ !( inv_mosaicatcher$sv_call_name%in% c("dup_h1", "complex", 
                                                                      "dup_h2", "idup_h1", "idup_h2", "del_h1", 
                                                                      "del_h2", "dup_hom", "del_hom")), ]


##filter 200kb
inversions <- filter(inversions, inversions$len>400)


#g <- ggplot(table_inv_new, aes(stack))
#g + geom_bar(aes(fill=sv_call_name), width = 0.5) +
 # theme(axis.text.x = element_text(angle=90, vjust=0.6)) +
 # labs(title="Samples with Contradiction on Inversion calls")+ scale_fill_brewer(palette = "Set2")
###make the one for bedtools
#mosaic_new<-read.table("/home/celia/Downloads/Mosaicatcher/output_101022.csv", header = TRUE, sep = '\t')

#renaming rows
##create list of names
no_inv <- c("0|0", "0|0_lowconf")
inv_hom<- c("1|1", "1|1_lowconf", "inv_hom", "0101")
inv_het<- c("0|1", "1|0", "0110", "1001", "0|1_lowconf", "1|0_lowconf", "inv_h1", "inv_h2")

#rename row variables
inv_mosaic <- inv_mosaic %>%
  mutate_at(vars(contains("sv_call_name")),
            funs(ifelse(. %in% inv_hom, "inv_hom",
                        ifelse(. %in% inv_het, "inv_het", NA_character_))))


inversions<- inversions %>%
  mutate_at(vars(contains("NA1")),
            funs(ifelse(. %in% no_inv, "no_inv",
                        ifelse(. %in% inv_hom, "inv_hom",
                               ifelse(. %in% inv_het, "inv_het", NA_character_)))))
inversions<- inversions %>%
  mutate_at(vars(contains("NA20")),
            funs(ifelse(. %in% no_inv, "no_inv",
                        ifelse(. %in% inv_hom, "inv_hom",
                               ifelse(. %in% inv_het, "inv_het", NA_character_)))))

inversions<- inversions %>%
  mutate_at(vars(contains("NA24")),
            funs(ifelse(. %in% no_inv, "no_inv",
                        ifelse(. %in% inv_hom, "inv_hom",
                               ifelse(. %in% inv_het, "inv_het", NA_character_)))))


inversions<-inversions%>%
  mutate_at(vars(contains("HG")),
            funs(ifelse(. %in% no_inv, "no_inv",
                        ifelse(. %in% inv_hom, "inv_hom",
                               ifelse(. %in% inv_het, "inv_het", NA_character_)))))

#inversions<-na.omit(inversions)

####BEDTOOLS
write.table(inv_mosaic, "/home/celia/Downloads/Mosaicatcher/mosaic_new.bed", sep = '\t', row.names=FALSE, quote = F, col.names=FALSE)
write.table(inversions, "/home/celia/Downloads/Mosaicatcher/inv_paper_400kb.bed", sep = '\t', row.names=FALSE,quote = F, col.names=FALSE)

###
mosaic_file= "/home/celia/Downloads/Mosaicatcher/mosaic_new.bed"
inversions_all="/home/celia/Downloads/Mosaicatcher/inv_paper_400kb.bed"
output= "/home/celia/Downloads/Mosaicatcher/intersect400.bed"

system(paste0("bedtools intersect -a ", inversions_all,  " -b " , mosaic_file, " -f 0.7 -r -wao > " , output))
output_intersect_4 <-read.table("/home/celia/Downloads/Mosaicatcher/intersect400.bed",  sep = '\t')

output_confusion<- subset(output_intersect_4, select = -c(V4,V5,V7,V8,V9, V10, V11, V12, V13, V14, V62, V65,V66, V64, V68,V69
                                                          ,V70,V71,V72,V73,V74) )

colnames(output_confusion) <- c('chrom','start','end','stack',"NA19434", "HG00096","HG00171","HG00268","HG00512","HG00513","HG00514","HG00731","HG00732","HG00733","HG00864","HG01114","HG01352", "HG01505" ,"HG01573","HG01596","HG02011",
"HG02018","HG02059","HG02106","HG02492","HG02587","HG02818","HG03009","HG03065","HG03125", "HG03371","HG03486" ,"HG03683","HG03732","HG04217","NA12329","NA12878","NA18534","NA18939","NA19036","NA19238","NA19239",
"NA19240", "NA19650","NA19983","NA20509","NA20847","NA24385","chrom_mosaic","start_mosaic", "end_mosaic", "sample","SV")

output_confusion_unique<-output_confusion %>% distinct(stack, sample, SV, .keep_all = TRUE)

write.table(output_confusion_unique, "/home/celia/Downloads/Mosaicatcher/celia-testdir/celia-testdir/output_confusion_400.csv", sep='\t', quote = F, col.names=TRUE)








#rownames(inv_paper)<- case_when(rownames(inv_paper) == "0|0" ~ "no_inv",
                       #   rownames(inv_paper) == "1|1" ~ "inv_hom",
                       #   TRUE ~ rownames(inv_paper))

#rownames(inv_paper)[which(rownames(inv_paper) %in% c("0|0", "1|1"))] <- c("no_inv","inv_hom")

#renaming rows
##create list of names
no_inv <- c("0|0", "0|0_lowconf")
inv_hom<- c("1|1", "1|1_lowconf", "inv_hom")
inv_het<- c("0|1", "1|0", "0|1_lowconf", "1|0_lowconf", "inv_h1", "inv_h2")

#rename row variables
inv_mosaicatcher <- inv_mosaicatcher %>%
  mutate_at(vars(contains("sv_call_name")),
            funs(ifelse(. %in% inv_hom, "inv_hom",
              ifelse(. %in% inv_het, "inv_het", NA_character_))))


confusion<- confusion %>%
  mutate_at(vars(contains("NA1")),
            funs(ifelse(. %in% no_inv, "no_inv",
                        ifelse(. %in% inv_hom, "inv_hom",
                               ifelse(. %in% inv_het, "inv_het", NA_character_)))))

confusion<- confusion %>%
  mutate_at(vars(contains("HG")),
            funs(ifelse(. %in% no_inv, "no_inv",
                        ifelse(. %in% inv_hom, "inv_hom",
                               ifelse(. %in% inv_het, "inv_het", NA_character_)))))


#prepare bed files for bedtools intersect
#names(inv_paper)[names(inv_paper) == 'seqnames'] <- 'chrom'
#inv_intersect<-inv_paper%>%select(chrom, start, end)
#inv_mosaic_short<- inv_mosaic %>%select(chrom, start, end)
#inv_mosaic_all<- inv_mosaic %>%select(chrom, start, end, cell)

#output="/home/celia/Downloads/Mosaicatcher/intersect_all.bed"
#output1="/home/celia/Downloads/Mosaicatcher/testing.bed"


write.table(inv_paper, "/home/celia/Downloads/Mosaicatcher/inv_paper.bed", sep = '\t', row.names=FALSE, quote = F, col.names=FALSE)
write.table(inv_mosaic_short, "/home/celia/Downloads/Mosaicatcher/inv_mosaic.bed", sep = '\t', row.names=FALSE,quote = F, col.names=FALSE)
write.table(inv_intersect, "/home/celia/Downloads/Mosaicatcher/inv_intersect.bed", sep = '\t', row.names=FALSE)
write.table(inv_mosaic_all, "/home/celia/Downloads/Mosaicatcher/inv_4.bed", sep = '\t', row.names=FALSE)


mosaic_file= "/home/celia/Downloads/Mosaicatcher/mosaic_int.bed"
inversions_paper="/home/celia/Downloads/Mosaicatcher/test_int.bed"

system(paste0("bedtools intersect -a ", mosaic_file,  " -b " , inversions_paper, " -f 0.75 -r > " , output))



inv_intersect_all <-read.table("/home/celia/Downloads/Mosaicatcher/intersect_all.bed", header = TRUE, sep = '\t')
inversions_testing<-read.table("/home/celia/Downloads/Mosaicatcher/testing.bed", header = TRUE, sep = '\t')
colnames(inversions_testing) <- c("chrom","start","end","cell")

as.data.table(inv_intersect_all)

merged2<-merge(inversions_testing, inv_paper, by=c('chrom', 'start','end'), all.x=TRUE)



melt.inversions <-melt(inv_paper,id=c("seqnames","start","end"),
                  measure.vars= c("NA19434","HG00096", "HG00171","HG00268","HG00512","HG00513", "HG00514","HG00731","HG00732","HG00733",
                 "HG00864","HG01114","HG01352" ,"HG01505","HG01573","HG01596","HG02011","HG02018" ,"HG02059","HG02106","HG02492","HG02587",
                 "HG02818","HG03009","HG03065","HG03125","HG03371","HG03486" ,"HG03683","HG03732","HG04217","NA12329","NA12878","NA18534",
                 "NA18939","NA19036","NA19238","NA19239","NA19240","NA19650","NA19983", "NA20509","NA20847","NA24385"),
                 variable.name="sv_call_name",
                 value.name= "cell")

