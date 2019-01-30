################################################################
#
#   Jacques Corbeil <jacques.corbeil@crchul.ulaval.ca>
#    Frédéric Raymond <frederic.raymond@crchul.ulaval.ca>
#   Sébastien Boisvert <Sébastien.Boisvert@USherbrooke.ca> (T1)
#
# author : Sébastien Boisvert                  ..
# creation date 2006-%m-%d                       ....
#                                                   ... .
# this file is part of the leishmania project           .....
#                                                           ..
# all the source code is in a subversion repository          .
#                                                            .
# ls3.genome.ulaval.ca:/svn/leishmania
#
###############################################################

############################################################

leish_debugger <-local( function(message, verbosity = FALSE) {

        if ( verbosity == TRUE) {
                print (message)
        }
})

############################################################
leish_get_columns_name <- local(function (first_file) {
        leish_debugger("----> leish_get_columns_name")
        conn=file(first_file,"r")
        content <- readLines ( conn ,n=2)
        line_with_amount_of_header_lines = content[2]
        columns=strsplit(line_with_amount_of_header_lines,split="\t",
                fixed=TRUE)

        nb_lines_in_header =as.numeric(columns[[1]][1])

#        /*skip the header*/
        readLines(conn, n= nb_lines_in_header)

# can I do do-while in R ???????
#        do {
#                line_with_columns = readLines (conn,n=1)
#        } while ( length( grep ("Block", line_with_columns)) == 0)

        line_with_columns = readLines (conn,n=1)
        while (length( grep ("Block", line_with_columns)) == 0) {
                line_with_columns = readLines (conn,n=1)
        }

        close(conn)

        columns <- strsplit(line_with_columns, split="\t")
        columns_to_process <- columns [[1]]
        columns_to_return <- c()

        for (a_column in columns_to_process)    {
            a_clean_column = gsub("\"", "", a_column)
            columns_to_return <- cbind(c(a_clean_column), 
                                       columns_to_return)
        }

        columns_to_return
})

############################################################
#
# Description: Extract columns which match a pattern
#
#       Arguments :
#
#            columns_name : columns
#               pattern : regexp
#
leish_get_columns_using_a_pattern <- local(function ( columns_name, pattern ) {
        leish_debugger ("------> leish_get_columns_using_a_pattern")

        columns <- c()

        leish_debugger("BEGIN;")
        leish_debugger(columns_name)
        leish_debugger("COMMIT;")

        for ( column_name in columns_name ) {
                leish_debugger(paste("COLUMN : ", column_name))
                leish_debugger (paste("Pattern : ", pattern))
                if (length( grep ( pattern, column_name)) > 0) {
                        columns <- cbind (columns, c(column_name))
                }
        }

        columns
})

############################################################
leish_get_mean_columns <- local( function ( columns_name ) {
        leish_debugger ("--------> leish_get_mean_columns")

        leish_get_columns_using_a_pattern (columns_name, "^F.+Mean$")
})

############################################################
leish_get_median_columns <-local( function ( columns_name ) {
        leish_debugger ("-----------> leish_get_median_columns")
        leish_get_columns_using_a_pattern (columns_name, "^B.+Median$")
})

############################################################
leish_get_wave_lengths <- local(function (mean_columns) {

        leish_debugger ("-----------> leish_get_wave_lengths")
        wave_lengths <- c()

        leish_debugger (paste("nb_means : ", length(mean_columns)))

        for (mean_column in mean_columns ){
                wave_length <- gsub ( "[^0-9]","",mean_column)
                leish_debugger ("wave length  : ")
                leish_debugger (wave_length)
                wave_lengths <- cbind (wave_lengths, c(wave_length))
        }

        leish_debugger ("wave_lengths : ")
        leish_debugger (wave_lengths)

        wave_lengths
})

############################################################
leish_get_columns_list <- local(function (mean_columns, median_columns,
                red_wave_length, green_wave_length){

        leish_debugger("----------------> leish_get_columns_list")
        red_mean <- NULL # "F647 Mean"
        green_mean <- NULL # "F555 Mean"
        red_median <- NULL # "B647 Median"
        green_median <- NULL # "B555 Median"

        leish_debugger (paste("red pattern : ", red_wave_length))
        leish_debugger (paste("green pattern : " , green_wave_length))

        for( mean_column in mean_columns){
                if (length(grep (red_wave_length, mean_column))>0){
                        red_mean = mean_column
                } else if (length(grep (green_wave_length, mean_column))>0) {
                        green_mean = mean_column
                }
        }

        for (median_column in median_columns) {
                if (length(grep (red_wave_length, median_column))>0) {
                        red_median = median_column
                } else if (length(grep ( green_wave_length, median_column))>0) {
                        green_median = median_column
                }
        }

        list(R= red_mean, G= green_mean,
                Rb=red_median, Gb= green_median)
})

############################################################
leish_get_gpr_columns<- local( function(targets){
        leish_debugger("---------> leish_get_gpr_columns")

        first_file <- targets$FileName[1]

        columns_name <- leish_get_columns_name(first_file)

        leish_debugger ("column names ")
        leish_debugger (columns_name)

        mean_columns <- leish_get_mean_columns (columns_name)
        median_columns <- leish_get_median_columns (columns_name)

        leish_debugger("mean columns : ")
        leish_debugger (mean_columns)



        wave_lengths <- leish_get_wave_lengths (mean_columns)

        leish_debugger("wave_lengths : " )
        leish_debugger ( wave_lengths)

        red_wave_length <- NULL
        green_wave_length <- NULL

        if (wave_lengths[1] < wave_lengths[2]){
                red_wave_length <- wave_lengths[2]
                green_wave_length <- wave_lengths[1]
        }else if (wave_lengths[2] < wave_lengths[1]){
                red_wave_length <- wave_lengths[1]
                green_wave_length <- wave_lengths[2]
        }

        leish_get_columns_list (mean_columns, median_columns, red_wave_length, green_wave_length)

})

# EOF
