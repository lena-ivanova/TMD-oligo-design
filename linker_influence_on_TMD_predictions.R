# code written by : Elena Ivanova
# Lab : Khemlinskii, Institute of Molecular Biology
# Project: Testing the influence of different linkers on the TMD predictions

setwd("S:/linker fusions universal/GAvsGS")

mydir = "S:/linker fusions universal/GAvsGS"

#load packages
library("dplyr")
library("stringi")
library("htmltab")
library("XML")
library("ggplot2")

yeast_TA <- read.csv("yeast_TA.csv", header=TRUE) #read in starting table with the TA protein AA sequences and Phobious start and end predictions for the TMDs
mNG_linker <- read.csv("mNG_linker.csv", header=TRUE) #read in starting table with the fluorescent reporter sequence and linker sequence

#####################################################################################################################################################################################################
# "make_fusions" function creates the fusion of the fluorescent reporter sequence, linker and TMD sequence (up_TMD_extension - flanking residues from the N-terminal side of the TMD, 
#                                                                                              down_TMD_extension - flanking residues from the C-terminal side of the TMD,
#                                                                                              stopseq - sequence flanking the TMD at the C-terminus
#                                                                                              linker - sequence between the fluorescent reporter and N-terminal part of the TMD)
#####################################################################################################################################################################################################

make_fusions <- function(up_TMD_extension, down_TMD_extension, stopseq, linker) {
  
  yeast_TA$TMD_seq <- substr(yeast_TA$Sequence, yeast_TA$Phobius_start, yeast_TA$Phobius_end) #get TMD sequence based on start and end
  yeast_TA$TMD_length <-stri_length(yeast_TA$TMD_seq)
  
  yeast_TA$C_term_seq <- substr(yeast_TA$Sequence, yeast_TA$Phobius_end+1, yeast_TA$Length) #get the entire C-terminal seq after TMD
  
  yeast_TA$up_TMD <- substr(yeast_TA$Sequence, yeast_TA$Phobius_start-up_TMD_extension, yeast_TA$Phobius_start-1) #get sequence in front of TMD (define the amount of flanking residues a "up_TMD_extension")
  yeast_TA$up_TMD_length <-stri_length(yeast_TA$up_TMD)
  
  yeast_TA$down_TMD <- substr(yeast_TA$Sequence, yeast_TA$Phobius_end+1, yeast_TA$Phobius_end+down_TMD_extension) #get sequence after the TMD (define the amount of flanking residues a "down_TMD_extension")
  yeast_TA$down_TMD_length <-stri_length(yeast_TA$down_TMD)
  
  yeast_TA$stopseq <- stopseq # to exchange natural C-term to different defined sequence, e.g. "GPGG"
  
  
  mNG_linker$linker_seq <- linker # e.g., "GAGAGAGAGA" (GA)x5 linker, "GGPG" stopseq, "EQKLISEEDLGSGSGSGSGS" Flag-(GS)x3 linker, "EQKLISEEDLGSGSGSGSGS" Myc-(GS)x3 linker
  yeast_TA$mNG_seq <- mNG_linker$mNG_seq
  yeast_TA$mNG_length <-stri_length(yeast_TA$mNG_seq)
  yeast_TA$linker_seq <- mNG_linker$linker_seq
  yeast_TA$linker_length <-stri_length(yeast_TA$linker_seq)
  
  yeast_TA$mNG_linker_start <- yeast_TA$mNG_length + yeast_TA$linker_length+up_TMD_extension+1 #get expected TMD start
  yeast_TA$mNG_linker_end <- yeast_TA$mNG_length + yeast_TA$linker_length + yeast_TA$up_TMD_length + yeast_TA$TMD_length #get expected TMD end
  
  yeast_TA$mNG_linker_TMD <- paste(yeast_TA$mNG_seq, yeast_TA$linker_seq, yeast_TA$up_TMD, yeast_TA$TMD_seq, yeast_TA$down_TMD, yeast_TA$stopseq, sep = "", collapse = NULL) # create combinations of mNG_linker_TMD_C_term_seq
  #yeast_TA$mNG_linker_TMD_length <- stri_length(yeast_TA$mNG_linker_TMD)
  #yeast_TA$mNG_linker_TMD_TMDseq <- substr(yeast_TA$mNG_linker_TMD, yeast_TA$mNG_linker_start, yeast_TA$mNG_linker_end)
  #yeast_TA$mNG_linker_TMD_TMDseq_length <- stri_length(yeast_TA$mNG_linker_TMD_TMDseq)
  
  yeast_TA$mNG_ORF <- paste("mNG_", yeast_TA$ORF, sep = "", collapse = NULL)
  yeast_TA$mNG_ORF_FASTA <- paste(">", yeast_TA$mNG_ORF, sep = "", collapse = NULL)
  yeast_TA$FASTA <- paste(yeast_TA$mNG_ORF_FASTA, yeast_TA$mNG_linker_TMD, sep = "\n", collapse = NULL) #final sequences in FASTA format
  
  yeast_TA$up_TMD_extension <- up_TMD_extension
  yeast_TA$down_TMD_extension <- down_TMD_extension
  yeast_TA$stopseq <- as.character(stopseq)
  yeast_TA$linker <- linker
  yeast_TA$ID <- paste(up_TMD_extension, down_TMD_extension, stopseq, linker, sep="_", collapse = NULL)
  
  
  fusionsname1 <- paste0("fusions_FASTA", "_", up_TMD_extension, "_", down_TMD_extension, "_", stopseq, "_", linker, sep="", collapse = NULL)
  fusionsname2 <- paste0(fusionsname1, ".txt", sep="", collapse = NULL)
  fusionsnametxt <- paste0("S:/linker fusions universal/GAvsGS/", fusionsname2, sep="", collapse = NULL)
  write.table(select(yeast_TA, FASTA), fusionsnametxt, quote = FALSE, row.names = FALSE, col.names = FALSE) #export FASTA file
  
  fusionsname3 <- paste0(fusionsname1, ".csv", sep="", collapse = NULL)
  fusionsnamecsv <- paste0("S:/linker fusions universal/GAvsGS/", fusionsname3, sep="", collapse = NULL)
  write.csv(yeast_TA, fusionsnamecsv) #export csv file
  
  
  return(yeast_TA)
}



########
### feed table top PHOBIUS server
### output in "short" format, saved as html
### open output with notepad, remove useless header and footnote
########

#read in PHOBIUS output from notepad

parsing_prediction <- function(filename) {
  
  mycol = c("mNG_ORF", "TM", "SP", "prediction")
  pred = read.csv (file = filename, header = FALSE, sep = "", col.names = mycol)
  
  #parsing prediction to find TMD start and end
  k=nrow(pred)
  
  for (i in 1:k) {
    
    if(pred$TM[i] == 1) {
      
      pred$new_start[i] = as.numeric(sub("^[i,o]([[:digit:]]+)[-][[:digit:]]+[i,o]", "\\1", 
                                         pred$prediction[i], perl = TRUE))
      pred$new_end[i] = as.numeric(sub("^[i,o][[:digit:]]+[-]([[:digit:]]+)[i,o]", "\\1", pred$prediction[i],
                                       perl = TRUE))
      
    } else {
      
      if(pred$TM[i] == 2) {
        
        pred$new_start[i] = as.numeric(sub("^[i,o]([[:digit:]]+)[-][[:digit:]]+[i,o][[:digit:]]+[-][[:digit:]]+[i,o]", "\\1", 
                                           pred$prediction[i], perl = TRUE))
        
        pred$new_end[i] = as.numeric(sub("^[i,o][[:digit:]]+[-]([[:digit:]]+)[i,o][[:digit:]]+[-][[:digit:]]+[i,o]", "\\1", pred$prediction[i], 
                                         perl = TRUE))
        
      } else {
        
        pred$new_start[i] = NA
        pred$new_end[i] = NA
        
      }
    }
  }
  return(pred)
}

#wt GA linker

fusions_0_1_0_GAGAGAGAGA <- make_fusions(0, 1, "", "GAGAGAGAGA")
pred_0_1_0_GAGAGAGAGA <- parsing_prediction("pred_0_1_0_GAGAGAGAGA.txt")
df_combined_0_1_0_GAGAGAGAGA <- inner_join (fusions_0_1_0_GAGAGAGAGA, pred_0_1_0_GAGAGAGAGA, by="mNG_ORF")
#df_combined_0_1_0_GAGAGAGAGA$lost_tmds <-count_(NA, df_combined_0_1_0_GAGAGAGAGA$new_start)

#5_1_0_GA_linker

fusions_5_1_0_GAGAGAGAGA <- make_fusions(5, 1, "", "GAGAGAGAGA")
pred_5_1_0_GAGAGAGAGA <- parsing_prediction("pred_5_1_0_GAGAGAGAGA.txt")
df_combined_5_1_0_GAGAGAGAGA <- inner_join (fusions_5_1_0_GAGAGAGAGA, pred_5_1_0_GAGAGAGAGA, by="mNG_ORF")

# 0_1_0_EQKLISEEDLGSGSGSGSGS  Myc-(GS)x3 linker
fusions_0_1_0_EQKLISEEDLGSGSGS <- make_fusions(0, 1, "", "EQKLISEEDLGSGSGS")
pred_0_1_0_EQKLISEEDLGSGSGS <- parsing_prediction("pred_0_1_0_EQKLISEEDLGSGSGS.txt")
df_combined_0_1_0_EQKLISEEDLGSGSGS <- inner_join (fusions_0_1_0_EQKLISEEDLGSGSGS, pred_0_1_0_EQKLISEEDLGSGSGS, by="mNG_ORF")

# 5_1_0_EQKLISEEDLGSGSGSGSGS  Myc-(GS)x3 linker
fusions_5_1_0_EQKLISEEDLGSGSGS <- make_fusions(5, 1, "", "EQKLISEEDLGSGSGS")
pred_5_1_0_EQKLISEEDLGSGSGS <- parsing_prediction("pred_5_1_0_EQKLISEEDLGSGSGS.txt")
df_combined_5_1_0_EQKLISEEDLGSGSGS <- inner_join (fusions_5_1_0_EQKLISEEDLGSGSGS, pred_5_1_0_EQKLISEEDLGSGSGS, by="mNG_ORF")

# 0_1_0_EQKLISEEDLGSGSGSGSGS  Myc-(GS)x3 linker with GGPG
fusions_0_1_0_EQKLISEEDLGSGSGSGGPG <- make_fusions(0, 1, "", "EQKLISEEDLGSGSGSGGPG")
pred_0_1_0_EQKLISEEDLGSGSGSGGPG <- parsing_prediction("pred_0_1_0_EQKLISEEDLGSGSGSGGPG.txt")
df_combined_0_1_0_EQKLISEEDLGSGSGSGGPG <- inner_join (fusions_0_1_0_EQKLISEEDLGSGSGSGGPG, pred_0_1_0_EQKLISEEDLGSGSGSGGPG, by="mNG_ORF")


df_combined_for_boxplot <- bind_rows(df_combined_0_1_0_GAGAGAGAGA, df_combined_0_1_0_EQKLISEEDLGSGSGS, df_combined_0_1_0_EQKLISEEDLGSGSGSGGPG)

df_combined_for_boxplot$start_change <- df_combined_for_boxplot$mNG_linker_start-df_combined_for_boxplot$new_start

df_combined_for_boxplot$end_change <- df_combined_for_boxplot$mNG_linker_end-df_combined_for_boxplot$new_end

df_combined_for_boxplot$length_change <- (df_combined_for_boxplot$new_end-df_combined_for_boxplot$new_start+1)-df_combined_for_boxplot$TMD_length

df_final_for_boxplot <- select(df_combined_for_boxplot, ID, mNG_ORF, start_change, end_change, length_change)

write.csv(df_final_for_boxplot, "S:/linker fusions universal/GAvsGS/df_final_for_boxplot_0_1_0_.csv")

### Plotting the TMD start, end and length difference

#boxplot_start

boxplot_start <- #ggplot(data=df_final_for_boxplot)+
  #geom_boxplot(aes(ID, start_change))+
  #labs(x="condition", y="start difference")+
  ggplot(data=df_final_for_boxplot, aes(ID, start_change, fill=ID))+
  stat_boxplot( aes(ID, start_change), 
                geom='errorbar', linetype=1, width=0.5)+  #whiskers
  geom_boxplot( aes(ID, start_change)) +    
  #stat_summary(fun.y=mean, geom="point", size=2) + 
  #stat_summary(fun.data = mean_se, geom = "errorbar")+
  scale_fill_manual(values=c("hotpink1", "turquoise4","darkgoldenrod2"))+
  labs(x="condition", y="TMD start difference, AA")+
  
  scale_y_continuous(limits=c(-10,10), breaks=seq(-10,10, by = 1))

boxplot_start

boxplot_start + coord_flip()

pdf("start_change.pdf", width=12, height=10)

print(boxplot_start + coord_flip())

dev.off()

boxplot_end <- #ggplot(data=df_final_for_boxplot)+
  #geom_boxplot(aes(ID, end_change))+
  #labs(x="condition", y="end difference")+
  ggplot(data=df_final_for_boxplot, aes(ID, end_change, fill=ID))+
  stat_boxplot( aes(ID, start_change), 
                geom='errorbar', linetype=1, width=0.5)+  #whiskers
  geom_boxplot( aes(ID, start_change)) +    
  #stat_summary(fun.y=mean, geom="point", size=2) + 
  #stat_summary(fun.data = mean_se, geom = "errorbar")+
  scale_fill_manual(values=c("hotpink1","turquoise4","darkgoldenrod2"))+
  labs(x="condition", y="TMD end difference, AA")+
  
  scale_y_continuous(limits=c(-10,10), breaks=seq(-10,10, by = 1))

boxplot_end

boxplot_end + coord_flip()

pdf("end_change.pdf", width=12, height=10)

print(boxplot_end + coord_flip())

dev.off()

boxplot_length <- #ggplot(data=df_final_for_boxplot)+
  #geom_boxplot(aes(ID, end_change))+
  #labs(x="condition", y="length difference")+
  ggplot(data=df_final_for_boxplot, aes(ID, length_change, fill=ID))+
  stat_boxplot( aes(ID, length_change), 
                geom='errorbar', linetype=1, width=0.5)+  #whiskers
  geom_boxplot( aes(ID, length_change)) +    
  #stat_summary(fun.y=mean, geom="point", size=2) + 
  #stat_summary(fun.data = mean_se, geom = "errorbar")+
  scale_fill_manual(values=c("hotpink1","turquoise4","darkgoldenrod2"))+
  labs(x="condition", y="TMD length difference, AA")+
  
  scale_y_continuous(limits=c(-10,10), breaks=seq(-10,10, by = 1))

boxplot_length

boxplot_length + coord_flip()

pdf("length_change.pdf", width=12, height=10)

print(boxplot_length + coord_flip())

dev.off()
