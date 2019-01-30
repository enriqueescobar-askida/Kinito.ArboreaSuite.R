#Bloc1: Calculs

#load required libraries
library(marray)									;
library(limma)									;
library(convert)								;
library(nnNorm)									;


myTarget		<- readTargets("in_Files/targets.txt")	;
myDesign		<- modelMatrix(myTarget,ref="C");
#factors<- colnames(myDesign)
print(myDesign)									;


#directories
in_Path			<- "in_Files"					;
galPath			<- "galFiles"					;
rawPath			<- "rawFiles"					;
outPath			<- "outFiles"					;

#Bring data in BioC
datadir			<- paste(getwd(),rawPath,sep="/")	;
files			<- dir(path=	datadir,
						pattern=".txt")			;
slides			<- length(files)				;
s				<- slides						;
cat(paste("There are ",slides, "files to analyse!"))	;
RG				<- read.maimages(files,
								path=	datadir,
								"quantarray")	;
aGalFile		<- paste(galPath,"/","SpruceDir8_MNC-11k.gal",sep="")
RG$genes		<- readGAL(aGalFile)			;
RG$printer		<- getLayout(RG$genes)			;

#produce boxplots of log2-ratio for unnormalized slides
MA				<-normalizeWithinArrays(RG,
										method=	"none")	;
png(	file=	paste(outPath,"Boxplot_notnorm_raw.png",sep="/"),
		width=	800,
		height=	600,
		res=	600)							;

	boxplot(MA$M~col(MA$M),
			names=	colnames(MA$M))				;
			
dev.off()										;

#convert data to marray format
#raw data
dt				<- as(RG,"marrayRaw")			;

#obtain data with zeroed backgrounds
dtz				<- dt							;
dtz@maRb		<- dtz@maGb	<- matrix(0,
										dim(dtz)[1],
										dim(dtz)[2])	;

#perform normalizations with background subtraction
dtL				<- maNorm(dt)					;
dtC				<- maNormMain(dt,
							f.loc=	list(maNorm2D(),
										maNormLoess(x="maA",
													y="maM",
													z="maPrintTip")),
							a.loc=	maCompNormEq()
							)					;

#perform normalizations without background subtraction
dtzL			<- maNorm(dtz)					;
dtzC			<- maNormMain(dtz,
							f.loc=	list(maNorm2D(),
										maNormLoess(x="maA",
													y="maM",
													z="maPrintTip")),
							a.loc=	maCompNormEq()
							)					;

#convert to Limma format
dtL				<- as(dtL,	"MAList")			;
dtC				<- as(dtC,	"MAList")			;
dtzL			<- as(dtzL,	"MAList")			;
dtzC			<- as(dtzC,	"MAList")			;

dtL_b			<- dtL							;
dtC_b			<- dtC							;
dtzL_b			<- dtzL							;
dtzC_b			<- dtzC							;

#perform between slide normalization
dtL				<- normalizeBetweenArrays(dtL,	method="scale")	;
dtC				<- normalizeBetweenArrays(dtC,	method="scale")	;
dtzL			<- normalizeBetweenArrays(dtzL,	method="scale")	;
dtzC			<- normalizeBetweenArrays(dtzC,	method="scale")	;

#produce boxplots of log2-ratio of normalized slides
png(	file=	paste(outPath,"Boxplot_intraslides.png",sep="/"),
		width=	800,
		height=	600,
		res=	600)							;

	boxplot(dtC_b$M~col(dtC_b$M),
			names=	colnames(dtC_b$M))			;

dev.off()										;

png(	file=	paste(outPath,"Boxplot_interslides.png",sep="/"),
		width=	800,
		height=	600,
		res=	600)							;

	boxplot(dtC$M~col(dtC$M),
			names=	colnames(dtC$M))			;

dev.off()										;

dtL				<- as(dtL,	"marrayNorm") #raw data
dtC				<- as(dtC,	"marrayNorm") #raw data
dtzL			<- as(dtzL,	"marrayNorm") #raw data
dtzC			<- as(dtzC,	"marrayNorm") #raw data

#get the positions of duplicate spots
b				<- maGnames(dt)@maInfo$Name		;
nspt			<- dim(dt)[1]					;
dim(dt)[1]		->	nspt						;

almn			<- grep("MNC",b)				;
h1				<- almn[seq(1,
							length(almn),
							by=2)]				;
h2				<- almn[seq(2,
							length(almn),
							by=2)]				;

#define the function to compute the inside corelation
inslidecor		<- function(batch){
	incor		<- function(batch,i){
		a		<- na.omit(cbind(maM(batch)[h1,i],
							maM(batch)[h2,i]))	;
		cor(a[,1],
			a[,2])								;
	}
	incorval	<- NULL							;
	for( i in 1:s){
		incorval[i]	<- incor(batch,i)			;
	}
	incorval									;
}

#compute withinslide correlations 
inc				<- cbind(inslidecor(dt),
						inslidecor(dtL),
						inslidecor(dtC),
						inslidecor(dtz),
						inslidecor(dtzL),
						inslidecor(dtzC))		;
rownames(inc)	<- files						;
colnames(inc)	<- c("None_B",
					"PTLOess_B",
					"Composite_B",
					"None_NB",
					"PTLOess_NB",
					"Composite_NB")				;
write.table(c("Within slides correlation using all spots"),
 			file=	paste(outPath,"correlations.csv",sep="/"),
			sep=	"\t",
			row.names=FALSE,
			col.names=FALSE)					;
			
write.table(round(inc,3),
			file=	paste(outPath,"correlations.csv",sep="/"),
			sep=	"\t",
			append=TRUE)						;

#define a function to compute the between slides correlation
corm			<- function(batch){
	if (slides>1){
		mat		<- (maM(batch)[h1,]+maM(batch)[h2,])/2	;
		mat		<- na.omit(mat)					;
		cor(mat)								;
	}else{
		1										;
	}
}

#define a function to compute the average absolute between slides correlation
av				<- function(crs){
	acrs		<- abs(crs)						;
	(apply(acrs,2,sum)-1)/(slides-1)			;
}


#compute bewteen slides correlation
cdtz			<- corm(dtz)					;
cdtzL			<- corm(dtzL)					;
cdtzC			<- corm(dtzC)					;
cdt				<- corm(dt)						;
cdtL			<- corm(dtL)					;
cdtC			<- corm(dtC)					;


#compute bewteen slides correlation
acdtz			<- av(cdtz)						;
acdtzL			<- av(cdtzL)					;
acdtzC			<- av(cdtzC)					;

acdt			<- av(cdt)						;
acdtL			<- av(cdtL)						;
acdtC			<- av(cdtC)						;



write.table(c(" "),
			file=		paste(outPath,"correlations.csv",sep="/"),
			sep=		"\t",
			append=		TRUE,
			row.names=	FALSE,
			col.names=	FALSE)					;
			
write.table(c("Mean Absolute between slides correlation for every slide."),
			file=		paste(outPath,"correlations.csv",sep="/"),
			sep=		"\t",
			append=		TRUE,
			row.names=	FALSE,
			col.names=	FALSE)					;
			
inc				<- cbind(acdt,
						acdtL,
						acdtC,
						acdtz,
						acdtzL,
						acdtzC)					;
						
colnames(inc)	<- c("None_B",
					"PTLOess_B",
					"Composite_B",
					"None_NB",
					"PTLOess_NB",
					"Composite_NB")				;
					
rownames(inc)	<- files						;

write.table(round(inc,3),
			file=		paste(outPath,"correlations.csv",sep="/"),
			sep=		"\t",
			append=		TRUE)					;


write.table(c(" "),
			file=		paste(outPath,"correlations.csv",sep="/"),
			sep=		"\t",
			append=		TRUE,
			row.names=	FALSE,
			col.names=	FALSE)					;
			
write.table(c("Correlation matrix for Backg. corr. composite norm with back. corr."),
			file=		paste(outPath,"correlations.csv",sep="/"),
			sep=		"\t",
			append=		TRUE,
			row.names=	FALSE,
			col.names=	FALSE)					;
			
write.table(round(cdtC,3),
			file=		paste(outPath,"correlations.csv",sep="/"),
			sep=		"\t",
			append=		TRUE)					;

write.table(c(" "),
			file=		paste(outPath,"correlations.csv",sep="/"),
			sep=		"\t",
			append=		TRUE,
			row.names=	FALSE,
			col.names=	FALSE)					;
			
write.table(c("Control Values comp norm: GUS."),
			file=		paste(outPath,"correlations.csv",sep="/"),
			sep=		"\t",
			append=		TRUE,
			row.names=	FALSE,
			col.names=	FALSE)					;
			
write.table(round(maM(dtC)[grep("GUS",b),],3),
			file=		paste(outPath,"correlations.csv",sep="/"),
			sep=		"\t",
			append=		TRUE)					;

write.table(c(" "),
			file=		paste(outPath,"correlations.csv",sep="/"),
			sep=		"\t",
			append=		TRUE,
			row.names=	FALSE,
			col.names=	FALSE)					;
			
write.table(c("Control Values comp norm: ptMYB5."),
			file=		paste(outPath,"correlations.csv",sep="/"),,
			sep=		"\t",
			append=		TRUE,
			row.names=	FALSE,
			col.names=	FALSE)					;
			
write.table(round(maM(dtC)[grep("ptMYB5",b),],3),
			file=		paste(outPath,"correlations.csv",sep="/"),
			sep=		"\t",
			append=		TRUE)					;

write.table(c(" "),
			file=		paste(outPath,"correlations.csv",sep="/"),
			sep=		"\t",
			append=		TRUE,
			row.names=	FALSE,
			col.names=	FALSE)					;
			
write.table(c("Control Values comp norm: ptHB4."),
			file=		paste(outPath,"correlations.csv",sep="/"),
			sep=		"\t",
			append=		TRUE,
			row.names=	FALSE,
			col.names=	FALSE)					;
			
write.table(round(maM(dtC)[grep("ptHB4",b),],3),
			file=		paste(outPath,"correlations.csv",sep="/"),
			sep=		"\t",
			append=		TRUE)					;


ms				<- maM(dtC)[grep("BAR|BT|Globin|pECE-1",b),]	;
ci				<- matrix(0,
						s,
						2)						;
						
for (i in 1:s){
	ci[i,]		<- c((t.test(ms[,i]))$conf.int[1],
				(t.test(ms[,i]))$conf.int[2])	;
}

write.table(c(" "),
			file=		paste(outPath,"correlations.csv",sep="/"),
			sep=		"\t",
			append=		TRUE,
			row.names=	FALSE,
			col.names=	FALSE)					;
			
write.table(c("Confidence interval on controls BAR,BT,Globin,pECE-1. Comp norm with back. corr"),
			file=		paste(outPath,"correlations.csv",sep="/"),
			sep=		"\t",
			append=		TRUE,
			row.names=	FALSE,
			col.names=	FALSE)					;
			
colnames(ci)	<- c("Lower","Upper")			;
rownames(ci)	<- files						;
write.table(round(ci,4),
			file=		paste(outPath,"correlations.csv",sep="/"),
			sep=		"\t",
			append=		TRUE)					;


isok			<- function(x,l,h){
	(x>=h)|(x<=l)								;
}

mat				<- (maM(dtC)[h1,]+maM(dtC)[h2,])/2	;
mat_WithNA		<- (maM(dtC)[h1,]+maM(dtC)[h2,])/2	;
mat				<- na.omit(mat)					;
cs				<- NULL							;

for (j in 1:dim(mat)[1]){
	sums		<- 0							;
	for (i in 1:s){
		sums	<- sums+ isok(mat[j,i],
							ci[i,1],
							ci[i,2])			;
	}
	if(sums==slides){
		cs		<- c(cs,TRUE)					;
	}else{
		cs		<- c(cs,FALSE)					;
	}
}

cormT			<- function(batch){
	if (slides>1){
		mat		<- (maM(batch)[h1,]+maM(batch)[h2,])/2	;
		mat		<- na.omit(mat)					;
		cor(mat[cs,])							;
	}else{
		1										;
	}
}

write.table(c(" "),
			file=		paste(outPath,"correlations.csv",sep="/"),
			sep=		"\t",
			append=		TRUE,
			row.names=	FALSE,
			col.names=	FALSE)					;
			
write.table(c("Correlation matrix for Backg. corr. comp. norm. using only spot with M values outside confid. int."),
			file=		paste(outPath,"correlations.csv",sep="/"),
			sep=		"\t",
			append=		TRUE,
			row.names=	FALSE,
			col.names=	FALSE)					;
			
write.table(round(cormT(dtC),3),
			file=		paste(outPath,"correlations.csv",sep="/"),
			sep=		"\t",
			append=		TRUE)					;


#define the function to compute the inside corelation
inslidecorT		<- function(batch){
	incor		<- function(batch,i){
		a		<- na.omit(cbind(maM(batch)[h1,i],
							maM(batch)[h2,i]))	;
		am		<- apply(a,
						1,
						mean)					;
		okp		<- NULL							;
		for (k in 1:length(am)){
			okp[k]	<- isok(am[k],
							ci[i,1],
							ci[i,2])			;
		}
		cor(a[okp,1],
			a[okp,2])							;
	}
	incorval	<- NULL							;
	for( i in 1:s){
		incorval[i]	<- incor(batch,
							i)					;
	}
	incorval									;
}


#compute withinslide correlations 
inc				<- cbind(inslidecorT(dtC))		;
rownames(inc)	<- files						;
colnames(inc)	<- c("Composite_B_TRunc")		;
write.table(c(" "),
			file=		paste(outPath,"correlations.csv",sep="/"),
			sep=		"\t",
			append=		TRUE,
			row.names=	FALSE,
			col.names=	FALSE)					;
			
write.table(c("Within slides correlation using spots outside conf. int."),
			file=		paste(outPath,"correlations.csv",sep="/"),
			sep=		"\t",
			row.names=	FALSE,
			col.names=	FALSE,
			append=		TRUE)					;
			
write.table(round(inc,3),
			file=		paste(outPath,"correlations.csv",sep="/"),
			sep=		"\t",
			append=		TRUE)					;


ResNomComp		<- as(dtC,
					"MAList")					;

#Write MAplots
for (s in 1:slides){
# 	print(paste(outPath,"/","MAplots_",s,"_",colnames(RG)[s],".png",sep=""))	;
	
	png(file=	paste(outPath,"/","MAplots_",s,"_",colnames(RG)[s],"__raw",".png",sep=""),
		width=	800,
		height=	800,
		res=	600	)							;
		
		plotMA(RG,
				array	= s,
				xlab    = "A:Intensity(log)",
				ylab    = "M:Fold(log)",
				main	= paste(colnames(RG)[s],"__raw",sep=""),
				xlim    = c(6,16),
				ylim    = c(-4,4)
				)								;

		abline(	h=	ci[s,1],
				col="green")					;

		abline(	h=	ci[s,2],
				col="green")					;

		abline(	h=	0,
				col="blue")						;

	dev.off()									;
	
	png(file=	paste(outPath,"/","MAplots_",s,"_",colnames(RG)[s],"_norm",".png",sep=""),
		width=	800,
		height=	800,
		res=	600	)							;
		
		plotMA(ResNomComp,
				array	= s,
				xlab    = "A:Fold(log)",
				ylab    = "M:Intensity(log)",
				main	= paste(colnames(RG)[s],"_norm", sep=""),
				xlim    = c(6,16),
				ylim    = c(-4,4)
				)								;
				
		abline(	h=	ci[s,1],
				col="red")						;

		abline(	h=	ci[s,2],
				col="red")						;

		abline(	h=	0,
				col="blue")						;

	dev.off()									;
}


#Write_Limma_ANN est une fonction qui annote directement les r�ultats limma
# 
#   Format: Write_Limma_ANN<- function(batch, annName, AnnFile)  
#
#   Ex: Write_Limma_ANN(tt1, "Analyse_Coeff1.csv", "R-SpruceDir7ann.txt")
#
#   batch: objet r�ultat de limma. le nom peut varier si vous faites sortir plusieurs analyses limma en m�e temps
#   annName : repr�ente le fichier excel dans lequel les r�ultats seront �ris
#   AnnFile: le fichier EASE qui servira de r��ence pour l'annotation


Write_Limma_ANN	<- function(batch, annName, AnnFile){
	x			<- read.csv(AnnFile,
							header=		TRUE,
							sep=		"\t",
							fill=		TRUE)	;
	x			<- read.table(AnnFile,
							header=		FALSE,
							col.names=	c("MNCID", "ANNOTATION"),
							sep=		"\t",
							fill=		TRUE)	;
	
	MNS			<- batch$ID						;
	nbMNS		<- length(MNS)					;
	
	MNSANN		<- matrix(0,
						nbMNS,
						2,
						dimnames=	list(c(1:nbMNS),c("MNCID", "ANNOTATION")))	;
	
	ci			<- matrix(0,
							length(batch$M),
							7)					;
	for (i in 1:length(batch$M)){
		ci[i,]	<- c(batch$ID[i],
					as.character(x[x$MNCID==MNS[i],]$ANNOTATION[1]),
					batch$M[i],
					batch$A[i],
					batch$t[i],
					batch$P.Value[i],
					batch$B[i])					;
	}
	entete		<- c("ID",
					"ANNOTATION",
					"M",
					"A",
					"t",
					"P.Value",
					"B")						;
	
	annName		<- paste(annName,".csv",sep="")	;
	cat(paste(annName, " created !\n"))			;
	write.table(t(entete),
				file=		annName,
				sep=		"\t",
				row.names=	FALSE,
				col.names=	FALSE)				;
				
	write.table(ci,
				file=		annName,
				sep=		"\t",
				append=		TRUE,
				row.names=	FALSE,
				col.names=	FALSE)				;
}

#take input from user
# 1 quand c'est du bon c�� -1 quand c'est un dye swap
#designS		<- c(1,-1,1,1,-1,1,-1,1,-1,1,-1,-1)	;
#designS		<- c(-1,-1,1,-1,1,1,1,-1)		;
#designS		<- c(1,1,1,-1,-1,-1,1,1,1,-1,-1,-1,1,1,1,-1,-1,-1)	;
#designS		<- c(1,1,1,-1,-1,-1,1,1,1,-1,-1,-1)	;
designS			<- c(1,1,1,-1,-1,-1,-1,1,-1,1)	;
#designS		<- c(1,-1-1,1,1,-1,-1,-1,-1,1,-1,1,1,-1)
#designS		<- c(-1,-1,1,-1,1,1,1,1,1,-1,1,-1,-1,-1,1,-1)
print(designS)									;

#give the name and adjust for signs and write standford
rownames(mat_WithNA)	<- b[h1]				;
mm				<- mat_WithNA*
					t(matrix(rep(designS,dim(mat_WithNA)[1]),slides,dim(mat_WithNA)[1]))	;
write.table(mm,
			file=	paste(outPath,"standford.txt",sep="/"),
			sep=	"\t")						;


#Bloc2: Analyse LIMMA
LigneeRegardee	<- 1							;
targets			<- readTargets(paste(in_Path,"/","targets.txt",sep=""))	;
design			<- modelMatrix(targets,
								ref="C")		;
factors			<- colnames(design)				;

##########Analyse des g�es diff�entiellement exprim�
#analyse par lign�s
fit				<- lmFit(mat_WithNA,
						design,
						correlation=	mean(inslidecor(dtC)))	;
fit2			<- eBayes(fit,
						proportion=	0.01)		;

for (m in 1:length(factors)){
	tt			<- topTable(fit2,
							coef=			m,
							n=				11036,
							adjust.method=	"fdr",
							sort.by=		"P")	;
	write.table(tt,
				file=	paste(outPath,"/",factors[m],".csv",sep=""),
				sep=	"\t")					;
	liste		<- tt[tt$P.Value<0.01,]			;
		if (nrow(liste)>0){
			Write_Limma_ANN(liste,
							paste(outPath,"/",factors[m],"_Pfam",			sep=""),
							paste(in_Path,"/","SpruceDir8_Pfam.dat",		sep=""))	;
			Write_Limma_ANN(liste,
							paste(outPath,"/",factors[m],"_nr",				sep=""),
							paste(in_Path,"/","SpruceDir8_nr.dat",			sep=""))	;
			Write_Limma_ANN(liste,
							paste(outPath,"/",factors[m],"_PGI5",			sep=""),
							paste(in_Path,"/","SpruceDir8_PGI5.dat",		sep=""))	;
			Write_Limma_ANN(liste,
							paste(outPath,"/",factors[m],"_SGI1",			sep=""),
							paste(in_Path,"/","SpruceDir8_SGI1.dat",		sep=""))	;
			Write_Limma_ANN(liste,
							paste(outPath,"/",factors[m],"_UnigeneID",		sep=""),
							paste(in_Path,"/","SpruceDir8_UnigeneID.tab",	sep=""))	;
			Write_Limma_ANN(liste,
							paste(outPath,"/",factors[m],"_ContigDir8",		sep=""),
							paste(in_Path,"/","SpruceDir8_ContigDir8.dat",	sep=""))	;
		}
}

#analyse de toutes les lign�s ensemble
#EL005targets <-  readTargets("targets3.txt")
#design <- modelMatrix(EL005targets,ref="C")
#factors<- colnames(design)

#fit3 <- lmFit(mat_WithNA,design,correlation=mean(inslidecor(dtC)))
#fit4<- eBayes(fit3,proportion=0.01)

#for (m in 1:length(factors)){
#tt<- topTable(fit4,coef=m ,n=9053,adjust.method="fdr",sort.by="P")
#write.table(tt,file=paste(factors[m],".csv",sep=""),sep="\t")
#liste<- tt[tt$P.Value<0.01,]
#    if (nrow(liste)>0){
#    Write_Limma_ANN(liste, paste(outPath,"/",factors[m],"_pfam"), "SpruceDir7_Pfam_EASE.dat")
#    Write_Limma_ANN(liste, paste(outPath,"/",factors[m],"_Dir7"), "SpruceDir7_ContigDir7_EASE.dat")
#    Write_Limma_ANN(liste, paste(outPath,"/",factors[m],"_nr"), "SpruceDir7_nr_EASE.dat")
#    Write_Limma_ANN(liste, paste(outPath,"/",factors[m],"_PGI5"), "SpruceDir7_PGI5_EASE.dat")
#    Write_Limma_ANN(liste, paste(outPath,"/",factors[m],"_SGI1"), "SpruceDir7_SGI1_EASE.dat")
#    }
#}

#Bloc3: Annotation

#Cette fonction annotera vos r�ultats limma avec un fichier d'annotation EASE
# 
#   Ex: Write_Limma_ANN(tt1, "Analyse_Coeff1.csv", "R-SpruceDir7ann.txt")
#
#   tt1: objet r�ultat de limma. le nom peut varier si vous faites sortir plusieurs analyses limma en m�e temps
#   "Analyse_Coeff1.csv" : le fichier excel dans lequel les r�ultats seront �ris
#   "R-SpruceDir7ann.txt": le fichier EASE qui servira de r��ence pour l'annotation
#
#


paste("Analysis Done!")												;
