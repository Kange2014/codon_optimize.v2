##############################################################
# This script implemets the approach described in Nature paper
# "Codon influence on protein expression in E. coli correlates 
# with mRNA levels" by Grégory Boël, et al. in 2016.
# 
# It could be used to calculate the probability of obtaining
# the high expression level vs. no expression, just based on
# a given DNA sequence. Therefore, we could use it to help 
# redesign heterogeous genes for their expression in E.Coli.
#
# 2016/09/07, LUZE@novonordisk.com

library(seqinr)

codon_influence_cal <- function (dnaseq="",utr="TTAATACGACTCACTATAGGGGAATTGTGAGCGGATAACAATTCCCCTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACAT"){

	dnaseq_obj = s2c(dnaseq) ;
	
	# codon slopes from linear logistic-regression models
	# they provide a useful metric to quantify the influence of individual codons on protein expression
	codon_slopes <- list(
		AAA = 4.3048, AAC = 1.53155, AAG	= 6.70028, AAT	= -5.2547, ACA	= 3.6303,
		ACC	= 9.69558, ACG	= 3.1174, ACT	= -0.62752, AGA	= 4.7637, AGC	= 6.76485,
		AGG	= 11.59071, AGT	= 2.11559, ATA	= -31.62539, ATC	= -6.70654, ATT	= -0.19064,
		CAA	= 8.85073, CAC	= 2.31177, CAG	= 1.7476, CAT	= 10.49089, CCA	= -4.98274,
		CCC	= 12.35937, CCG	= 6.54786, CCT	= -6.86108, CGA	= -8.45732,CGC	= 2.54922,
		CGG	= -13.00096, CGT	= -4.11513, CTA	= -10.82582, CTC	= 3.53066, CTG	= 0.99496,
		CTT	= -5.45645, GAA	= 22.51082, GAC	= 16.18596, GAG	= 15.51138, GAT	= 23.84976,
		GCA	= 1.45082, GCC	= 0.97019, GCG	= 9.61392, GCT	= 2.187, GGA	= 10.6617,
		GGC	= 5.05263, GGG	= 7.83171, GGT	= 10.83417, GTA	= 11.32142, GTC	= 0.30884,
		GTG	= 3.79534, GTT	= 4.91921, TAC	= -5.45385, TAT	= 3.0145, TCA	= 8.17619,
		TCC	= -3.25351, TCG	= -9.66795, TCT	= -6.23922, TGC	= -10.7029, TGG	= -7.48736,
		TGT	= -12.47396, TTA	= -5.24472, TTC	= -3.94608, TTG	= -5.3648, TTT	= 1.97578
	);
	
	# RNAstructure predicted free energy of folding of 5'UTR + first 48bp (16 codons)  
	
	utr_head <- toupper(paste(utr,substr(dnaseq,1,48),sep=""));
	utr_head <- gsub("T","U",utr_head); ## convert to rna for RNAstructure
	utr_head_obj <- s2c(utr_head);
	fasta <- list(utr_head_obj);
	write.fasta(fasta,"tmp01",file.out="~/tmpseq.fasta");
	#shell(paste("RNAfold.exe", "-p <tmpseq.fasta>tmpseq.rnafold.txt"));
	system("partition-smp ~/tmpseq.fasta ~/tmpseq.pfs")
	system("EnsembleEnergy ~/tmpseq.pfs > ~/tmpseq.txt")
	
	#raw_fold_data <- readLines("tmpseq.rnafold.txt");
	raw_fold_data <- readLines("~/tmpseq.txt");
	efes <- strsplit(raw_fold_data[4],": "); ## ensemble free energy
	efes <- sub(" kcal/mol","",efes[[1]][2]);
	deltaG <- as.numeric(efes);

	system("rm ~/tmpseq*");
	
	# frequency of A and G in codons 2-6
	freq <- seqinr::count(tolower(dnaseq_obj[4:18]),wordsize=1,freq=TRUE);
	A_freq <- freq[1][[1]];
	G_freq <- freq[3][[1]];
	
	# GC content of codons 2-6
	GC = freq[2][[1]] + freq[3][[1]];
	
	# binary indicator I that is 1 only if deltaG < -39 kcal/mol and GC content of codons 2-6 is greater than 62%
	if(deltaG < -39 && GC > 0.62 ){ I = 1; }
	else { I = 0; }
	
	# T(U) frequency at the third position in codons 2-6
	freq <- seqinr::count(tolower(dnaseq_obj[4:18]),wordsize=1,start=2,by=3,freq=TRUE);
	T3_freq <- freq[4][[1]];
	
	# frequency of each non-termination codon (TAG, TAA, TGA), excluding "ATG"
	exc_codons <- c("atg","tag","taa","tga");
	last3_chars <- tolower(substr(dnaseq,length(dnaseq_obj) - 2, length(dnaseq_obj)));
	if(any(exc_codons[2:4] == last3_chars)){freq <- seqinr::count(tolower(dnaseq_obj)[1:(length(dnaseq_obj) - 3)],wordsize=3,by=3,freq=TRUE);}
	else{ freq <- seqinr::count(tolower(dnaseq_obj),wordsize=3,by=3,freq=TRUE); }
	codon_freq <- freq[!(names(freq) %in% exc_codons)];
	
	# mean slopes for codons 7-16 and 17-32
	codons1 <- splitseq(tolower(dnaseq_obj[19:48]));
	codons2 <- splitseq(tolower(dnaseq_obj[49:96]));
	
	#try(if codons1 %in% exc_codons) stop("There is/are initial or stop codon/s in codons 7-16."))
	#try(if(codons2 %in% exc_codons) stop("There is/are initial or stop codon/s in codons 17-32."))
	codons1 <- codons1[!(codons1 %in% exc_codons)]
	codons2 <- codons2[!(codons2 %in% exc_codons)]
	
	count1 <- table(codons1);
	count2 <- table(codons2);
	
	av_slope1 <- sum(as.numeric(codon_slopes[tolower(names(codon_slopes)) %in% names(count1)])*as.vector(count1))/sum(as.vector(count1));
	av_slope2 <- sum(as.numeric(codon_slopes[tolower(names(codon_slopes)) %in% names(count2)])*as.vector(count2))/sum(as.vector(count2));
	
	# AUA-AUA di-codon match
	codons = splitseq(tolower(dnaseq_obj));
	di_codons = paste(codons[1:(length(codons)-1)],codons[2:length(codons)],sep="")
	if(any(di_codons == "ataata")) { d_AUA = 1; }
	else{ d_AUA = 0; }
	
	# amino acid repetition rate r
	aa <- translate(dnaseq_obj);
	d <- c();
	for(i in 1:(length(aa)-1)){
		if( any(aa[(i+1):length(aa)] == aa[i]) ){ locs <- which(aa[(i+1):length(aa)] == aa[i]); d <-c(d, 1/locs[1]); }
		else{ d <- c(d,0); }
	}
	d <- c(d,0);
	r <- mean(d);
	
	# multi-parameter logistic-regression model of codon influence
	sAll = as.numeric(codon_slopes)*as.numeric(codon_freq);
	y = 2.0 + 0.054*deltaG - 1.5*I + 6.6*A_freq - 6.3*A_freq^2 - 1.8*G_freq^2 + 0.80*T3_freq + 
		sum(0.86*sAll) + 0.078*av_slope1 + 0.063*av_slope2 -
		1.8*d_AUA - 16*r - 0.0012*length(dnaseq_obj) - 520/length(dnaseq_obj);

	# calculate the probability of obtaining high expression vs no expression
	p <- exp(y)/(1+exp(y));
	
	parameters <- list(deltaG=deltaG,GC=GC,A_freq=A_freq,A_freq2=A_freq^2,G_freq2=G_freq^2,T3_freq=T3_freq,ATAATA=d_AUA,sAll=sum(sAll));
	codons <- as.list(sAll);
	names(codons) <- names(codon_slopes);
	parameters <- c(parameters,codons);
	parameters <- c(parameters,list(Av_slope7_16=av_slope1,Av_slope17_32=av_slope2,r=r,theta=y,Probability=p))
	
	return(parameters);
}
