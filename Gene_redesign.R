##############################################################
# This script implemets the approach described in Nature paper
# "Codon influence on protein expression in E. coli correlates 
# with mRNA levels" by Grégory Boël, et al. in 2016.
# 
# It could be used to re-design genes by substuiting codons 
# with synonymous ones to improve the recombinant expression 
# levels in E.Coli. Generally, it tries to reduce the native
# degeneracy of the genetic code to eliminate codons that correlates
# with reduced protein expression in single-variable logistic-
# regression analysis of the large-scale data set.
#
# 2016/09/07, LUZE@novonordisk.com

library(seqinr)

# This funtion works as 6AA method in the paper
# In this method, a single codon was used for 6 amino aids, while codons for other 14 codons were not changed. 

sixAA <- function (dnaseq = ""){
	dnaseq_obj = s2c(dnaseq) ;
	aa <- aaa(translate(dnaseq_obj));
	codons <- splitseq(tolower(dnaseq_obj));
	six_aa <- c("Arg","Asp","Gln","Glu","His","Ile");
	sixCodons <- c("cgt","gat","caa","gaa","cat","att");
	
	for(i in 1:length(six_aa)){
		if(any(aa == six_aa[i])) { loc <- which(aa == six_aa[i]); codons[loc] = sixCodons[i]; }
	}
	return(c2s(codons));
}

# These two functions work as 31C method in the paper
# In the method, the free energy was optimizied using only the indicated subset of codons for each amino acid

AA_codons <- list(
					Ala=c("GCT","GCA"),Arg=c("CGT","CGA"),Asn="AAT",Asp="GAT",
					Cys="TGT",Gln=c("CAA","CAG"),Glu="GAA",Gly="GGT",His=c("CAT","CAC"),
					Ile=c("ATT","ATC"),Leu=c("TTA","TTG","CTA"),Lys="AAA",Met="ATG",
					Phe="TTT",Pro=c("CCT","CCA"),Ser=c("AGT","TCA"),Thr=c("ACA","ACT"),
					Trp="TGG",Tyr="TAT",Val=c("GTT","GTA")
					);

AA_codons_names <- names(AA_codons);

f <- function(x,y) paste(x,y,sep="");

# optimize the head (16 codons) + 5'UTR with free energy of folding: maximizing
head_folding_optimize <- function(dnaseq = "",utr="TTAATACGACTCACTATAGGGGAATTGTGAGCGGATAACAATTCCCCTCTAGAAATAATTTTGTTTAACTTTAAGAAGGAGATATACAT"){
	dnaseq_obj = s2c(dnaseq) ;
	aa <- aaa(translate(dnaseq_obj));
	codons <- splitseq(tolower(dnaseq_obj));
	seq_len <- length(aa);
	
	syn_codons <- list();
	for(i in 1:seq_len){
		if(any(AA_codons_names == aa[i])){ 
			loc <- which(AA_codons_names == aa[i]);
			tmp <- list(AA_codons[[loc]])
			syn_codons <- c(syn_codons,tmp); 
		}
	}
	
	combn <- "";
	for(i in 1:seq_len){
			combn <- outer(combn,syn_codons[[i]],f)
			combn <- as.vector(combn)
	}
	combn <- paste(utr,toupper(combn),sep="");
	deltaG <- efe_folding_cal(combn);
	loc = which(deltaG[,2] == max(deltaG[,2]));
	subseqs <- deltaG[loc,];
	subseqs[,1] <- substr(subseqs[,1],nchar(utr)+1,nchar(utr)+nchar(dnaseq));
	return(subseqs);
}

# optimize the tail sequence with the free energy of folding near to -10 kcal/mol in windows of 48bp
tail_folding_optimize <- function (dnaseq = ""){
	dnaseq_obj = s2c(dnaseq) ;
	aa <- aaa(translate(dnaseq_obj));
	codons <- splitseq(tolower(dnaseq_obj));
	seq_len <- length(aa);

	syn_codons <- list();
	for(i in 1:seq_len){
		if(any(AA_codons_names == aa[i])){ 
			loc <- which(AA_codons_names == aa[i]);
			tmp <- list(AA_codons[[loc]]);
			syn_codons <- c(syn_codons,tmp); 
		}
	}
	
	if(seq_len >= 16){
		window_num <- seq_len%/%16;
		remainder <- seq_len%%16;
		if(remainder != 0) { window_num = window_num + 1; }
		
		for(i in 1:window_num){
			combn <- "";
			# for the last window whose size is less than 48bp
			if(i == window_num){
				for(j in ((i-1)*16+1):seq_len){
					combn <- outer(combn,syn_codons[[j]],f);
					combn <- as.vector(combn);
				}
				## extract specific length of nucleotides to make a window size of 48bp with the last window
				print(combn[1])
				preseqs <- substr(subseqs[,1],((i-2)*16+remainder)*3 + 1,(i-1)*16*3)
				print(preseqs[1])
				combn <- outer(preseqs,combn,f);
				combn <- as.vector(combn);
				print(combn[1])
				deltaG <- efe_folding_cal(combn);
					
				max_value <- max(deltaG[,2])
				min_value <- min(deltaG[,2])
				if(max_value < -10){ loc = which(deltaG[,2] == max_value); }
				else if(min_value > -10) { loc = which(deltaG[,2] == min_value); }
				else{ loc = which(deltaG[,2] <= -9.5 && deltaG[,2] >= -10.5); }
				
				lastseqs <- substr(deltaG[loc,1],48- remainder*3 + 1,48);
				index <- duplicated(lastseqs);
				lastseqs <- unique(lastseqs);
				values <- deltaG[loc,2];
				values <- values[!index];
					
				tmp_seqs <- outer(subseqs[,1],lastseqs,f);
				tmp_seqs <- as.vector(tmp_seqs);
				tmp_values <- outer(subseqs[,2],values,paste);
				tmp_values <-as.vector(tmp_values);
				subseqs <- cbind(tmp_seqs,tmp_values);
			}
			else{
				for(j in ((i-1)*16+1):(i*16)){
					combn <- outer(combn,syn_codons[[j]],f);
					combn <- as.vector(combn);
				}
				deltaG <- efe_folding_cal(combn);
				
				max_value <- max(deltaG[,2])
				min_value <- min(deltaG[,2])
				if(max_value < -10){ loc = which(deltaG[,2] == max_value); }
				else if(min_value > -10) { loc = which(deltaG[,2] == min_value); }
				else{ loc = which(deltaG[,2] <= -9.5 && deltaG[,2] >= -10.5); }
				
				if(i == 1){ subseqs <- deltaG[loc,]; }
				else{
					tmp_seqs <- outer(subseqs[,1],deltaG[loc,1],f);
					tmp_seqs <- as.vector(tmp_seqs);
					tmp_values <- outer(subseqs[,2],deltaG[loc,2],paste);
					tmp_values <-as.vector(tmp_values);
					subseqs <- cbind(tmp_seqs,tmp_values);
				}
			}	
		}
	}
	else{
		combn <- "";
		for(i in 1:seq_len){
			combn <- outer(combn,syn_codons[[i]],f)
			combn <- as.vector(combn)
		}
		deltaG <- efe_folding_cal(combn);
		max_value <- max(deltaG[,2])
		min_value <- min(deltaG[,2])
		if(max_value < -10){ loc = which(deltaG[,2] == max_value); }
		else if(min_value > -10) { loc = which(deltaG[,2] == min_value); }
		else{ loc = which(deltaG[,2] <= -9.5 && deltaG[,2] >= -10.5); }
		subseqs <- deltaG[loc,];	
	}
	
	return(subseqs);
}

remove_motif <- function(dnaseq = ""){
	# remove Shine-Dalgarno sequence robosome binds to initiate transcription: AGGAGG
	index1 <- grepl("AGGAGG",toupper(dnaseq));
	dnaseq <- dnaseq[!index1];
	
	# remove RNA destabilizing motif: ATTTA
	index2 <- grepl("ATTTA",toupper(dnaseq));
	dnaseq2 <- dnaseq[!index2];
	
	# remove undesired restriction enzyme sites
	sites <- c("ACCGGT","AGATCT","GGATCC","GATATC","GAATTC","AAGCTT","GGTACC","ACGCGT","TTAATTAA","CCATGG","CATATG","GCTAGC",
				"GAGCTC","ACTAGT","TCTAGA","CTCGAG","CTGCAG","GTCGAC","AGCGCT","CAATTG","GCGGCCGC"
				);
	
	for(i in sites){
		index <- grepl(i, dnaseq);
		dnaseq <- dnaseq[!index];
	}
	
	return(dnaseq);
}

## ViennaRNA RNAfold
efe_folding_cal <- function(all_seqs = c()){
	num = length(all_seqs)
	ids <- paste(">tmp",1:num,sep="")
	new_seqs <- c()
	new_seqs[seq(1,2*num,2)] = ids
	new_seqs[seq(2,2*num,2)] = all_seqs
	write(new_seqs,"~/tmpseq.fasta")
	
	system("RNAfold -p < ~/tmpseq.fasta > ~/tmpseq.rnafold.txt" );
	
	raw_fold_data <- readLines("~/tmpseq.rnafold.txt");
	
	efe_vec <- c()
	for(i in seq(1,length(raw_fold_data),6)){
		efes <- strsplit(raw_fold_data[i+3]," ")
		## ensemble free energy
		len <- length(efes[[1]]);
		efes <- sub("[[]","",efes[[1]][len]);
		efes <- sub("[]]","",efes);
		efe_vec <- c(efe_vec,efes)
	}
	
	deltaG <- data.frame(Seq=all_seqs,EFE=as.numeric(efe_vec));

	system("rm ~/tmpseq.rnafold.txt");
	system("rm ~/tmpseq.fasta");
	system("rm ~/*.ps");


	return(deltaG);
}
