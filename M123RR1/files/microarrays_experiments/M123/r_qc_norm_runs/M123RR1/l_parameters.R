
parameters.image_width = 800
parameters.image_height=600
parameters.reference="A"
parameters.read.maimages.source="genepix"

parameters.modifyWeights.w_keys = c("Spike1", "Oligo", "BLANK")
parameters.modifyWeights.w_values = c(0, 1, 0)
    parameters.read.maimages.wt.fun=wtflags(0.0)
parameters.readTargets.targets_file="targets.txt"
parameters.readGAL.gal_file="control-corrected_060210_corbeil_gal.gal"
parameters.readSpotTypes.spottypes_file="spottypes.txt"
parameters.backgroundCorrect.method="minimum"
parameters.normalizeWithinArrays.method="printtiploess"
parameters.normalizeWithinArrays.weights="w"
parameters.normalizeWithinArrays.span=0.3
    parameters.normalizeWithinArrays.iterations=4
    parameters.normalizeBetweenArrays.method="Aquantile"
parameters.normalization.unique_name="M123RR1"
parameters.quote=FALSE
parameters.dec="."

parameters.decideTests.method="global"
parameters.decideTests.ajust.method="holm"
parameters.decideTests.p=0.01
