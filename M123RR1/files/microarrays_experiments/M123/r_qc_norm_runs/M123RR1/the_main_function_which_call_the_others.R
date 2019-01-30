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

draw_weights_image <- local(function(file_path,
                                       width,
                                       height,
                                     data_structure) {
    # data_Structure
    
    already_found <- c()
    z_values <- data_structure@maW
    
    for(the_index in 1:length(z_values)) {
        the_current_value <- z_values [the_index]
          if(!Array.includes(already_found, the_current_value))   {
              already_found <- 
                cbind(already_found, c(the_current_value))
            }
    }

    #
    # we actually need at least 2 colors (usually 1 and 0)
    #
    
    length_of_already_found <- length(already_found)

    print (length_of_already_found)

    if(length_of_already_found <2) {
        print (file_path)
        print ("CANNOT DRAW IMAGE : x:ok, y:ok, z: not ok")
        return (FALSE)
    } else {

    }

    print ("AFTER")

    png(filename =file_path, width = width, height =height)

    # Picture of the distribution of the flags in
    # the third hybridization
                    
    image(data_structure, xvar = "maW", col = rainbow(10))
    dev.off()
})

the_main_function_which_call_the_others <- local(function() {

        targets=readTargets(parameters.readTargets.targets_file)
        RG<-NULL

        gpr_columns = leish_get_gpr_columns (targets)

        #
        # Read the grp files
        #
        if (is.null(parameters.read.maimages.wt.fun)) {

        RG <- read.maimages(targets,
        columns=gpr_columns,
                annotation = c("Block","Row","Column","ID","Name"),
        source=parameters.read.maimages.source)

        }    else {

        RG <- read.maimages(targets,
        columns=gpr_columns,
                annotation = c("Block","Row","Column","ID","Name"),
        source=parameters.read.maimages.source,
        wt.fun=parameters.read.maimages.wt.fun)
        }

        # READ THE GAL FILE
        RG$genes<-readGAL(parameters.readGAL.gal_file)
        RG$printer<-getLayout(RG$genes)

        #
        # Load the spottypes
        #

        spottypes<-readSpotTypes(parameters.readSpotTypes.spottypes_file)
        RG$genes$Status<-controlStatus(spottypes,RG)

        # create the weights
        w<-modifyWeights(RG$weights, RG$genes$Status,
        parameters.modifyWeights.w_keys, parameters.modifyWeights.w_values)

        #
        # The signals with a background correction
        #
        RGb<-backgroundCorrect(RG,method=parameters.backgroundCorrect.method)

        p_weights=NULL
        if (parameters.normalizeWithinArrays.weights == "rg"){
            p_weights=RG$weights
        }else if( parameters.normalizeWithinArrays.weights =="w"){
            p_weights=w
        }

        # normalize within array (printtiploess works well)
        #

        me1=parameters.normalizeWithinArrays.method

        if(me1=="loeess" || me1== "printtiploess" || me1=="composite"){

                if(me1=="composite") {
                                MA<-normalizeWithinArrays(RGb,method=
        parameters.normalizeWithinArrays.method,weights=p_weights,
                                span= parameters.normalizeWithinArrays.span,
                                iterations=
        parameters.normalizeWithinArrays.iterations,
                            controlspots=
        parameters.normalizeWithinArrays.controlspots)
                }else{

                        # "loess"
                            # "printtiploess"
                                    MA<-normalizeWithinArrays(RGb,
        method=parameters.normalizeWithinArrays.method,weights=p_weights,
                                span= parameters.normalizeWithinArrays.span,
                                iterations=
                parameters.normalizeWithinArrays.iterations)
                }
        }else if (me1=="robustspline") {# else
                    MA<-normalizeWithinArrays(RGb,
        method=parameters.normalizeWithinArrays.method,weights=p_weights,
                    df=parameters.normalizeWithinArrays.df)
        }else{
        # "none"
        # "median"
        # "control"

                MA<-normalizeWithinArrays(RGb,
        method=parameters.normalizeWithinArrays.method,weights=p_weights)
        }

        #
        # Now we have to normalize the data between arrays
        #
        #
        MA.between<-normalizeBetweenArrays(MA,
        method=parameters.normalizeBetweenArrays.method)

        #
        # Write the intensities files
        #
        #intensities
        RG.normalized <- RG.MA(MA.between)
        nbcol<-length(RG.normalized$R[1,])
        RG.normalized.matrix<-NULL
        namecolumn<-c("Genes")

        for (i in 1:nbcol) {
                RG.normalized.matrix<-cbind(RG.normalized.matrix,
                        RG.normalized$R[,i],RG.normalized$G[,i])
                tempnameR<-paste(RG.normalized$targets$SlideID[i],"-R", sep="")
                tempnameG<-paste(RG.normalized$targets$SlideID[i],"-G", sep="")
                namecolumn<-cbind(namecolumn,tempnameR,tempnameG)
        }

 #       quote=parameters.quote,sep="\t",eol="\n",na="NA",
 #               dec=parameters.dec,row.names=FALSE,col.names=FALSE,
 #               qmethod=c("double"))
        #
        # #write content
        #


        nb_slides = length(RG.normalized$targets$SlideID)

        nbcol<-length(MA.between$M[1,])
        namecolumn<-c("Genes")
        for(i in 1:nbcol) {
        namecolumn<-cbind(namecolumn,MA.between$targets$SlideID[i])

        }

        t_file = paste ("normalize/",parameters.normalization.unique_name,
        "-M.txt",sep="")

        # Define colors for plots
        Gcol=maPalette(low="white",high="green",k=50)
        Rcol=maPalette(low="white",high="red",k=50)
        RGcol=maPalette(low="green",high="red",mid="white",k=50)

        # Loess normalization (within-array, red vs green)

        # Normalize all arrays at the same time
        # Algorithm = Printtip Loess3

        marray_data <-as(RG, "marrayRaw")
        leishmania.norm<-as(MA.between,"marrayNorm")

        #
        #
        # PNG pictures OUTPUT
        #
        #####################
        png(filename = "qc/maBoxplot.png", width = parameters.image_width,
                height = parameters.image_height)

        # Boxplot of the logratios of
        #all the hybridizations after normalization:

        maBoxplot(leishmania.norm,y="maM",main="arrays after normalization")

        dev.off()

        #
        # densities

        png(filename ="qc/densities.png", width = parameters.image_width,
                height = parameters.image_height)
        par(mfrow=c(2,2))
        # Picture of the Red background of the first hybridization
        plotDensities.with.main(RG,main="RAW DATA")
        plotDensities.with.main(RGb,main="RAW DATA (background corrected)")
        plotDensities.with.main(MA,main="Normalized data within arrays")

        # Picture of the Red background of the first hybridization
        plotDensities.with.main(MA.between,main="Normalized data between arrays")
        dev.off()

        png(filename =paste("qc/","maBoxplot-maM.png",sep=""),
width = parameters.image_width, height = parameters.image_height)

      # Boxplot of the logratios of all the hybridizations (save picture)
        maBoxplot(marray_data,y="maM",main="arrays pre-normalization")
        dev.off()

        for (i in 1:nb_slides) {

              slide_use_code = targets["SlideID"][,1][i]

                    png(filename =paste("qc/",slide_use_code,"image-maRb.png",
        sep=""), width = parameters.image_width,
        height = parameters.image_height)

                    # Picture of the Red background of the first hybridization
                    image(marray_data[, i], xvar = "maRb", subset = TRUE, col = Rcol, contours = FALSE, bar = TRUE)
                    dev.off()

                    png(filename =paste("qc/",slide_use_code,"image-maRb.png",sep=""), width = parameters.image_width, height = parameters.image_height)

                    # Picture of the Red background of the first hybridization
                    image(marray_data[, i], xvar = "maRb", subset = TRUE, col = Rcol, contours = FALSE, bar = TRUE)
                    dev.off()

                    png(filename =paste("qc/",slide_use_code,"image-maGb.png",sep=""), width = parameters.image_width, height = parameters.image_height)

                    # Picture of the Green background of the second hybridization
                    image(marray_data[, i], xvar = "maGb", subset = TRUE, col = Gcol, contours = FALSE, bar = TRUE)
                    dev.off()
                    
                    draw_weights_image(
                      paste("qc/",slide_use_code,"image-maW.png",sep=""),
                      parameters.image_width,
                      parameters.image_height,
                      marray_data[, i]
                                       )

                    print("After the drawing")

                    # Log ratio of the two colors in one hybridization, nr 5 here
                    # main =  titel of the picture
                    # bar = TRUE gives a bar on the side with the scale
                    png(filename = paste("qc/",slide_use_code,"maImage-maM.png",sep=""), width = parameters.image_width, height = parameters.image_height)

                    tmp=maImage(marray_data[,i],x="maM",bar=TRUE,main="Image of prenormalization M")
                    dev.off()

                    # Boxplot of the logratios per block (= printtip) of one hybridization
                    # nr 4 here (copy paste this picture to Powerpoint)
                    png(filename =paste("qc/",slide_use_code,"maBoxplot-maPrintTip.png",sep=""), width = parameters.image_width, height = parameters.image_height)

                    maBoxplot(marray_data[,i],x="maPrintTip",y="maM",main="pre-normalization")
                    dev.off()

                    png(filename = paste("qc/",slide_use_code,"maPlot.png",sep=""), width = parameters.image_width*1, height = (parameters.image_height*1))

                    # MA plot with means per printtip, here for hybridization nr. 5
                    maPlot(marray_data[,i], legend.func=NULL)
                    dev.off()

                    # Boxplot of logratios per block (= printtip) of one hybridization
                    # after normalization, nr. 4 (save picture and compare to raw data)
                    png(filename = paste("qc/",slide_use_code,"maBoxplot-norm.png",sep=""), width = parameters.image_width, height = (parameters.image_height*1))

                    maBoxplot(leishmania.norm[,i],x="maPrintTip",y="maM",main=" arrays after normalization")
                    dev.off()

                    # MA plot with means per printtip, after normalization nr. 5
                    png(filename =paste("qc/",slide_use_code,"maPlot-norm.png",sep=""), width = parameters.image_width*1, height = parameters.image_height)

                    maPlot(leishmania.norm[,i], legend.func=NULL)
                    dev.off()
        }
        #

        #
        #
        #  code snippet from Frédéric Raymond
        # this code calculate the M ratio for each condition. It can later be used to check the results
        #

        reference = parameters.reference #condition de r��ence. Va pr�enter rations Pro/Ama

        design <- modelMatrix(targets, ref = reference)
        corfit <- duplicateCorrelation(MA.between, design, ndups =2, spacing =200)
        fit <- lmFit(MA.between, design, correlation=corfit$consensus, ndups =2, spacing =200, weights = w)
        fitbayes <- eBayes(fit)

        # Liste de g�es diff�entiellement exprim�

        #restable <- topTable(fitbayes, genelist = fitbayes$genes, adjust = "BH", number = 100, sort.by="p", resort.by="M")

        #Diagramme de venne des genes significativement exprim�

        # method = "separate", "global", "heirarchical", "nestedF"
        # adjust.method = "none", "BH", "fdr" (equivalent to "BH"), "BY" and "holm"
        # p = valeur entre 0 et 1

        results <- decideTests(fitbayes,method=parameters.decideTests.method,
        adjust.method =parameters.decideTests.ajust.method, p=parameters.decideTests.p)
        a <- vennCounts(results)
        #print(a)

        png(filename =paste("results/","venn.png",sep=""), width = parameters.image_width,
        height = parameters.image_height)
        vennDiagram(a)
                    dev.off()

        #
        #
        # Finally, what we were looking for is created, !!!
        #
        write.fit(fitbayes, results=results, file="results/results.txt", digits=2, adjust="none", sep="\t")

        control_analysis ( RG, RGb, MA.between )

}
)

