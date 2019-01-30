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
#  researcher  : Jacques Corbeil <jacques.corbeil@crchul.ulaval.ca>
#  Trainee Supervisor :  Frédéric Raymond <frederic.raymond.1@ulaval.ca>
#  Trainee: Sébastien Boisvert <Sébastien.Boisvert@UShebrooke.ca> (T1)
#
#
# Centre de recherche
# Pavillon CHUL
# Université Laval
#
#
# Centre de recherche en endocrinologie moléculaire et oncologique (CREMO) -- Projet Leishmania
#
#
# DATE 2006-07-11
#
#
# HEADER REVISION 1
#
#

# projet Leishmania
#
#
######################
#
# This file has been written by Frédéric Raymond and Sébastien Boisvert
#
####################
#

library(limma)
library(marray)
library(convert)

source("plotDensities.with.main.R")
source("l_parameters.R")
source("search.R")
source("leish_get_gpr_columns.R")
source("control_analysis.R")
source("the_main_function_which_call_the_others.R")

the_main_function_which_call_the_others ()
