################################################################
#
#   Jacques Corbeil <jacques.corbeil@crchul.ulaval.ca>
#    Frédéric Raymond <frederic.raymond@crchul.ulaval.ca>
#   Sébastien Boisvert <Sébastien.Boisvert@UShebrooke.ca> (T1)
#
# author : Sébastien Boisvert
# creation date 2006-%m-%d
#
# this file is part of the leishmania project
#
# all the source code is in a subversion repository
#
# ls3.genome.ulaval.ca:/svn/leishmania
#
###############################################################

#
#
# source code by Frédéric Raymond
#
# modified by Sébastien Boisvert
#
#
# red -green data structure
#

################################################################
  ## Function that gets results from a RG variable for a probe list

            getgeneRG <- local( function(RG = NULL, x1 = "") {
                    listedegenescherchable <- cbind(RG$genes$Name, RG$R, RG$G)
                    nbgenesinlist <- length(x1)
                    nbgenes <- length(listedegenescherchable[,1])
                    tempgene <- NULL;
                    for (j in (1:nbgenesinlist)) {
                            for (i in (1:nbgenes)) {
                                    if (listedegenescherchable[i,1] == x1[j]) {
                                            tempgene <- rbind(tempgene,
                                             listedegenescherchable[i,])
                                    }
                            }
                    }
                    tempgene
           })

################################################################
            ### Function that gets results from a MA variable for a probe list

            getgeneMA <-local( function(MA = NULL, x1 = "") {
                    listedegenescherchable <- cbind(MA$genes$Name, MA$M, MA$A)
                    nbgenesinlist <- length(x1)
                    nbgenes <- length(listedegenescherchable[,1])
                    tempgene <- NULL
                    for (j in (1:nbgenesinlist)) {
                            for (i in (1:nbgenes)) {
                                    if (listedegenescherchable[i,1] == x1[j]) {
                                            tempgene <- rbind(tempgene,
                                             listedegenescherchable[i,])
                                    }
                            }
                    }
                    tempgene
            })

################################################################
leishmania_generate_files_for_ma <- local(function(RG, MA.between, 
                              where_to_put_the_files, controllist){

        #### MA Section

        ### Generates a variable that
        #contains the names to use for each graph or table in the MA section

        namestogive <- NULL
        resulttype <- c("M","A")
        yaxis <- NULL
        for (i in (1:2)) {
                for (j in (1:length(rownames(RG$targets)))) {
                        namestogive <- rbind(namestogive,
                         paste(colnames(RG$weights)[j], resulttype[i]))
                        yaxis <- c(yaxis, resulttype[i])
                }
        }

        ### Generate graphs for quality control of
        #mismatch probes for MA results (normalized M and A results)
        ## Writes graphs in a PDF

        pdf(file = paste( where_to_put_the_files,"PlotMismatchMA.pdf",sep=""),
        width = 6, height = 6, onefile = TRUE, family = "Helvetica",
        title = "Quality Control Spots Analysis and Graphs", version = "1.1",
        paper = "letter")

        suffixlistmismatch <- c("", "-1m", "-2m", "-3m", "-5m", "-7m", "-10m")

        for (nbcontrol in (1:length(controllist))){

                genelistz <- c(	paste(controllist[nbcontrol],
                 suffixlistmismatch[1], sep =""),
                                paste(controllist[nbcontrol], 
                                      suffixlistmismatch[2], sep =""),
                                paste(controllist[nbcontrol], 
                                      suffixlistmismatch[3], sep =""),
                                paste(controllist[nbcontrol], 
                                      suffixlistmismatch[4], sep =""),
                                paste(controllist[nbcontrol], 
                                      suffixlistmismatch[5], sep =""),
                                paste(controllist[nbcontrol], 
                                      suffixlistmismatch[6], sep =""),
                                paste(controllist[nbcontrol], 
                                      suffixlistmismatch[7], sep =""))

                MAcontrol <- getgeneMA(MA.between, genelistz)

                nbcolumnforgene <- length(MAcontrol[1,])

                if(nbcolumnforgene>0) {
                        par(mfrow=c(2,2), oma = c(0,0,2,0))
                        for (i in (2:nbcolumnforgene)) {

                                if (yaxis[i-1] == c("M")) {
                                        limit <- c(-2,2)
                                }
                                if (yaxis[i-1] == c("A")) {
                                        limit <- c(0,14)
                                }
                                plot(c(0,0,0,0,1,1,2,2,3,3,5,5,7,7,10,10), 
                                     MAcontrol[,i],
                                 type="p", main = namestogive[i-1],
                                 xlab = "Number of mismatch", 
                                     ylab = yaxis[i-1],
                                 ylim = limit)
                                abline(0,0, lty=3)
                                title(paste(c("Mismatch plot for"),
                                 controllist[nbcontrol],c("- page"),
                                  ceiling((i-1)/4),c("of"),
                                            ceiling((nbcolumnforgene-1)/4),
                                   sep = " "), outer = TRUE)
                        }
                }

                if (nbcolumnforgene <1) {
                        cat(controllist[nbcontrol])
                }
        }

        dev.off()
        ## End of writing to PDF


        #### Calculates the statistics for probes
        #resynthesis for MA results (normalized M and A results)
        ### t test for each pair of resynthesis results assuming equal variance
        # for both probes
        ## Writes the results in a file using a serie of write.table
        #commands that append new lines to the file

        write.table(NULL, file= paste ( where_to_put_the_files,
         "resynthesisMA.txt",sep=""), append = FALSE, row.names=FALSE,
          col.names=FALSE, quote=FALSE)

        for (nbcontrol in (1:length(controllist))){

                totalvariable <- NULL
                genelistz <- c(controllist[nbcontrol])
                MAcontrol <- getgeneMA(MA.between, genelistz)
                nbcolumnforgene <- length(MAcontrol[1,])-1

                if(nbcolumnforgene>0) {
                        for (k in (1:length(namestogive))) {
                                x <- as.numeric(c(MAcontrol[1, k+1],
                                                  MAcontrol[2, k+1]))
                                y <- as.numeric(c(MAcontrol[3, k+1], 
                                                  MAcontrol[4, k+1]))

                                the_t_test_result <- t.test(x, y,
                                  var.equal = TRUE)
                                variable <- c(mean(x), sd(x), mean(y), sd(y), 
                                              the_t_test_result$statistic, 
                                              the_t_test_result$p.value)
                                totalvariable <- rbind(totalvariable, variable)
                        }

                        colnames(totalvariable) <-
                         c("Synt1.Mean","Synt1.SD","Synt2.Mean","Synt2.SD",
                           "t",                         "p.value")
                        rownames(totalvariable) <- namestogive

                        write.table(controllist[nbcontrol], file=paste(where_to_put_the_files, "resynthesisMA.txt", sep=""),
                         append = TRUE, row.names=FALSE, col.names=FALSE,
                          quote=FALSE)
                        write.table(cbind("Slide","Synt1.Mean","Synt1.SD",
                        "Synt2.Mean","Synt2.SD", "t", "p.value"),
                         file=paste(where_to_put_the_files,"resynthesisMA.txt",sep=""), append = TRUE, row.names=FALSE,
                          col.names=FALSE, quote=FALSE,sep=";")
                        write.table(totalvariable, file=paste (where_to_put_the_files,"resynthesisMA.txt",sep=""),
                         append = TRUE, col.names = FALSE, sep=";")
                }

                if (nbcolumnforgene <1) {
                        cat(controllist[nbcontrol])
                }
        }
        # End of writing to TXT


        ### Generate graphs for quality control
        # of positions of CDS probes for MA results (normalized M and A results)
        ## Writes graphs in a PDF

        pdf(file = paste(where_to_put_the_files,"PlotPositionMA.pdf",sep=""),
        width = 6, height = 6, onefile = TRUE, family = "Helvetica",
         title = "Quality Control Spots Analysis and Graphs",
          version = "1.1", paper = "letter")

        suffixposition <- c("","-p1", "-p2", "-p3")

        for (nbcontrol in (1:length(controllist))){

                genelistz <- c(	paste(controllist[nbcontrol], suffixposition[1],
                sep =""),
                                paste(controllist[nbcontrol], suffixposition[2], sep =""),
                                paste(controllist[nbcontrol], suffixposition[3], sep =""),
                                paste(controllist[nbcontrol], suffixposition[4], sep =""))

                MAcontrol <- getgeneMA(MA.between, genelistz)

                nbcolumnforgene <- length(MAcontrol[1,])

                if(nbcolumnforgene>0) {
                        par(mfrow=c(2,2), oma = c(0,0,2,0))
                        for (i in (2:nbcolumnforgene)) {
                                if (yaxis[i-1] == c("M")) {
                                        limit <- c(-2,2)
                                }
                                if (yaxis[i-1] == c("A")) {
                                        limit <- c(0,14)
                                }
                                plot(c(0,0,0,0,1,1,2,2,3,3), MAcontrol[,i], type="p",
                                 main = namestogive[i-1], xlab = "Position on CDS",
                                 ylab = yaxis[i-1], ylim = limit)
                                abline(0,0, lty=3)
                                title(paste(c("Position on CDS plot for"),
                                 controllist[nbcontrol],c("- page"),
                                 ceiling((i-1)/4),c("of"),ceiling((nbcolumnforgene-1)/4),
                                 sep = " "), outer = TRUE)
                        }
                }

                if (nbcolumnforgene <1) {
                        cat(controllist[nbcontrol])
                }
        }

        dev.off()
        ## End of writing to PDF

})


################################################################
#
# this function will call the other which generate MA stuff
#
#

leishmania_generate_files_for_rg <- local(function(RG, RGb, 
                              where_to_put_the_files, controllist ){
        #### RG Section

        ### Generates a variable that contains
        #the names to use for each graph or table in the RG section

        namestogive <- NULL
        for (i in (3:4)) {
                for (j in (1:length(rownames(RG$targets)))) {
                        namestogive <- rbind(namestogive,
                         paste(colnames(RG$weights)[j], colnames(RG$targets)[i]))
                }
        }

        ### Generate graphs for quality control of
        # mismatch probes for RG results (non normalized intensities)
        ## Writes graphs in a PDF

        pdf(file = paste(where_to_put_the_files,"PlotMismatchRG.pdf",sep=""),
        width = 6,
        height = 6, onefile = TRUE, family = "Helvetica",
         title = "Quality Control Spots Analysis and Graphs",
         version = "1.1", paper = "letter")

        suffixlistmismatch <- c("", "-1m", "-2m", "-3m", "-5m", "-7m", "-10m")

          for (nbcontrol in (1:length(controllist))){

                genelistz <- c(	paste(controllist[nbcontrol],
                suffixlistmismatch[1], sep =""),
                                paste(controllist[nbcontrol], suffixlistmismatch[2], sep =""),
                                paste(controllist[nbcontrol], suffixlistmismatch[3], sep =""),
                                paste(controllist[nbcontrol], suffixlistmismatch[4], sep =""),
                                paste(controllist[nbcontrol], suffixlistmismatch[5], sep =""),
                                paste(controllist[nbcontrol], suffixlistmismatch[6], sep =""),
                                paste(controllist[nbcontrol], suffixlistmismatch[7], sep =""))

                RGcontrol <- getgeneRG(RGb, genelistz)

                nbcolumnforgene <- length(RGcontrol[1,])
                maxyaxis <-
                 max(as.numeric(as.vector(RGcontrol)[17:length(RGcontrol)]))+
                 max(as.numeric(as.vector(RGcontrol)[17:length(RGcontrol)]))/10

                if(nbcolumnforgene>0) {
                        par(mfrow=c(2,2), oma = c(0,0,2,0))
                        for (i in (2:nbcolumnforgene)) {
                                plot(c(0,0,0,0,1,1,2,2,3,3,5,5,7,7,10,10), RGcontrol[,i],
                                 type="p", main = namestogive[i-1],
                                 xlab = "Number of mismatch",
                                 ylab = "Fluorescence intensity", ylim = c(0, maxyaxis))
                                title(paste(c("Mismatch plot for"),
                                 controllist[nbcontrol],c("- page"),
                                  ceiling((i-1)/4),c("of"),ceiling((nbcolumnforgene-1)/4),
                                   sep = " "), outer = TRUE)
                        }
                }

                if (nbcolumnforgene <1) {
                        cat(controllist[nbcontrol])
                }
        }

        dev.off()
        ## End of writting to PDF

        #### Calculates the statistics for probes
        #resynthesis  for RG results (non normalized intensities)
        ### t test for each pair of resynthesis results
        # assuming equal variance for both probes
        ## Writes the results in a file using a serie of
        #write.table commands that append new lines to the file

        write.table(NULL, file= paste ( where_to_put_the_files,
         "resynthesisRG.txt", sep=""), append = FALSE, row.names=FALSE,
          col.names=FALSE, quote=FALSE)

        for (nbcontrol in (1:length(controllist))){

                totalvariable <- NULL
                genelistz <- c(controllist[nbcontrol])
                RGcontrol <- getgeneRG(RGb, genelistz)
                nbcolumnforgene <- length(RGcontrol[1,])-1

                if(nbcolumnforgene>0) {
                        for (k in (1:length(namestogive))) {
                                x <- as.numeric(c(RGcontrol[1, k+1],RGcontrol[2, k+1]))
                                y <- as.numeric(c(RGcontrol[3, k+1], RGcontrol[4, k+1]))
                                variable <- c(mean(x), sd(x), mean(y), sd(y), t.test(x, y,
                                 var.equal = TRUE)$statistic, t.test(x, y,
                                 var.equal = TRUE)$p.value)
                                totalvariable <- rbind(totalvariable, variable)
                        }

                        colnames(totalvariable) <-
                         c("Synt1.Mean","Synt1.SD","Synt2.Mean","Synt2.SD", "t",
                          "p.value")
                        rownames(totalvariable) <- namestogive

                        write.table(controllist[nbcontrol], file= paste  (
                         where_to_put_the_files,"resynthesisRG.txt", sep=""), append = TRUE,
                          row.names=FALSE, col.names=FALSE, quote=FALSE)
                        write.table(cbind("Slide","Synt1.Mean","Synt1.SD",
                        "Synt2.Mean","Synt2.SD", "t", "p.value"), file= paste (
                         where_to_put_the_files,"resynthesisRG.txt",sep=""), append = TRUE,
                          row.names=FALSE, col.names=FALSE, quote=FALSE,sep=";")
                        write.table(totalvariable, file= paste ( where_to_put_the_files,
                         "resynthesisRG.txt",sep=""), append = TRUE, col.names = FALSE,sep=";")
                }

                if (nbcolumnforgene <1) {
                        cat(controllist[nbcontrol])
                }
        }
        ## End of writing to TXT

        ### Generate graphs for quality control
        #of positions of CDS probes for RG results (non normalized intensities)
        ## Writes graphs in a PDF

        pdf(file = paste ( where_to_put_the_files, "PlotPositionRG.pdf",sep=""),
        width = 6, height = 6, onefile = TRUE, family = "Helvetica",
        title = "Quality Control Spots Analysis and Graphs", version = "1.1",
        paper = "letter")

        suffixposition <- c("","-p1", "-p2", "-p3")

        for (nbcontrol in (1:length(controllist))){

                genelistz <- c(	paste(controllist[nbcontrol], suffixposition[1],
                 sep =""),
                                paste(controllist[nbcontrol], suffixposition[2], sep =""),
                                paste(controllist[nbcontrol], suffixposition[3], sep =""),
                                paste(controllist[nbcontrol], suffixposition[4], sep =""))

                RGcontrol <- getgeneRG(RGb, genelistz)
                nbcolumnforgene <- length(RGcontrol[1,])
                maxyaxis <-
                 max(as.numeric(as.vector(RGcontrol)[11:length(RGcontrol)]))+
                max(as.numeric(as.vector(RGcontrol)[11:length(RGcontrol)]))/10

                if(nbcolumnforgene>0) {
                        par(mfrow=c(2,2), oma = c(0,0,2,0))
                        for (i in (2:nbcolumnforgene)) {
                                plot(c(0,0,0,0,1,1,2,2,3,3), RGcontrol[,i], type="p",
                                 main = namestogive[i-1], xlab = "Location ID",
                                 ylab = "Fluorescence intensity", ylim = c(0, maxyaxis))
                                title(paste(c("Position on CDS plot for"),
                                 controllist[nbcontrol],c("- page"),
                                  ceiling((i-1)/4),c("of"),ceiling((nbcolumnforgene-1)/4),
                                  sep = " "), outer = TRUE)
                        }
                }

                if (nbcolumnforgene <1) {
                        cat(controllist[nbcontrol])
                }
        }

        dev.off()
        # End of writing to PDF
})

################################################################
control_analysis <- local(function( RG, RGb, MA.between){

        #
        #like  a sym link
        #

        where_to_put_the_files <- "control_analysis/"

        ### Start Analysis of control probes

        ## Variable containing list of control probes

        controllist <- c(
        "LinJ03.0010:f:601712-603886",
        "LinJ04.1250:r:1388553-1389683",
        "LinJ06.0080:r:1896416-1898389",
        "LinJ10.0780:r:4380886-4382679",
        "LinJ11.0040:f:4797802-4798707",
        "LinJ11.1250:f:5307734-5313088",
        "LinJ15.0790:r:7500544-7504266",
        "LinJ19.0630:f:10294671-10296509",
        "LinJ21.0610:r:12126638-12128617",
        "LinJ23.0240:f:13443645-13448333",
        "LinJ23.0290:f:13460442-13465151",
        "LinJ23.0420:f:13509239-13512979",
        "LinJ25.0540:r:15166757-15168727",
        "LinJ26.2650:f:17100851-17104654",
        "LinJ27.0380:f:17240623-17242689",
        "LinJ27.0550:r:17414952-17420483",
        "LinJ29.1960:r:20338849-20340240",
        "LinJ30.1700:f:21536374-21537699",
        "LinJ30.2950:r:22198913-22199998",
        "LinJ31.1210:r:23221742-23224267",
        "LinJ31.1650:r:23414010-23418623",
        "LinJ31.1660:r:23420411-23425072",
        "LinJ32.3640:f:26171279-26173393",
        "LinJ33.0340:r:26648879-26650885",
        "LinJ33.1860:f:27256747-27258669",
        "LinJ33.2840:f:27936614-27940627",
        "LinJ33.3050:f:28015462-28018311",
        "LinJ34.0690:f:28310216-28316485",
        "LinJ34.0860:f:28383964-28384473",
        "LinJ34.0950:r:28424389-28428414",
        "LinJ36.3580:r:33378067-33380118"
        )

        leishmania_generate_files_for_rg(RG, RGb, where_to_put_the_files, 
                                         controllist)

        leishmania_generate_files_for_ma(RG, MA.between, where_to_put_the_files, 
                                         controllist)
        # The End

}
)
