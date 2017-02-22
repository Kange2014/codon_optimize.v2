library(shiny)
library(shinydashboard)
library(seqinr)
library(ggplot2)
library(dplyr)
library(polycor)
library(AlgDesign)
source("Codon_influence_calculate.R")
source("Gene_redesign.R")

options(shiny.maxRequestSize=30*1024^2)

readInput <- function(inText,inFile){
		if ((inText == "" | is.null(inText)) & is.null(inFile)) return(NULL)
		if(inText != "" & (!is.null(inText)) ){
			if(grepl("^>",inText)){
				inText <- unlist(strsplit(inText,">"))
				SeqId <- c()
				inSeq <- c()
				for(i in 1:length(inText)){
					if(nchar(inText[i]) == 0 || inText[i] == "" || inText[i] == " " || inText[i] == "	"){ next }
					tmp_seq <- unlist(strsplit(inText[i],"\r?\n"))
					
					if(grepl("\\s+",tmp_seq[1])) { ids <- unlist(strsplit(tmp_seq[1],"\\s+")) }
					else {ids <- tmp_seq }
					
					SeqId <- c(SeqId,ids[1])
					inSeq <- c(inSeq,paste(tmp_seq[2:length(tmp_seq)],collapse=""))
				}
			}
			else{ 
				inText <- unlist(strsplit(inText,"\r?\n"))
				inSeq = inText
				SeqId = paste("tmp",1:length(inText),sep="")
			}
		}
		
		seq_strings <- c()
		if (! is.null(inFile)){
			seqs <- read.fasta(file = inFile$datapath, as.string = FALSE, seqtype = "DNA");
			for(i in 1:length(seqs)){
				seq_strings <- c(seq_strings,c2s(seqs[[i]]))
			}
			inSeq <- seq_strings
			SeqId = names(seqs)
		}
		
		return(list(Seq=inSeq,Id=SeqId))
}
	
# combinations
f <- function(x,y) paste(x,y,sep="")

shinyServer(function(input, output, session) {
	
	##################################################################################################

	source("copt_funs.R",local = TRUE)
	source("doe.R",local = TRUE)
	source("doe_ui.R",local = TRUE)
	
	source("regress.R",local = TRUE)
	source("regress_ui.R",local = TRUE)
	
	#################################################################################################
	
	# set up input interface
	output$resetable_input1 <- renderUI({
		times <- input$reset_input1
        div(id=letters[(times %% length(letters)) + 1],
            textAreaInput("text1", label = "Paste your sequence(s) here:", value = "")
			)
    })
	
	output$resetable_input4 <- renderUI({
		times <- input$reset_input4
        div(id=letters[(times %% length(letters)) + 1],
            textAreaInput("text4", label = "Paste your sequence(s) here:", value = "")
			)
    })
	
	# set up input interface
	output$resetable_input2 <- renderUI({
		times <- input$reset_input2
        div(id=letters[(times %% length(letters)) + 1],
            textAreaInput("text2", label = "Paste your sequence(s) here:", value = "")
			)
    })
	
	# set up input interface
	output$resetable_input3 <- renderUI({
		times <- input$reset_input3
        div(id=letters[(times %% length(letters)) + 1],
            textAreaInput("text3", label = "Paste your sequence(s) here:", value = "")
			)
    })
	
	# capture the input value
	v <- reactiveValues(data1 = NULL,data2 = NULL,data3 = NULL, data4 = NULL)
	
	observeEvent(input$calculate1,{
		v$data1 <- input$calculate1
	})
	
	observeEvent(input$reset_input1, {
		v$data1 <- NULL
	})  
	
	observeEvent(input$calculate4,{
		v$data4 <- input$calculate4
	})
	
	observeEvent(input$reset_input4, {
		v$data4 <- NULL
	})  
	
	observeEvent(input$calculate2,{
		v$data2 <- input$calculate2
	})
	
	observeEvent(input$reset_input2, {
		v$data2 <- NULL
	})  
		
	observeEvent(input$calculate3, {
		v$data3 <- input$calculate3
	})
	
	observeEvent(input$reset_input3, {
		v$data3 <- NULL
	})
	
	output$parameters <- renderTable({ 
	
		if (is.null(v$data1)) return()

		isolate(inText <- input$text1)
	
		# input$file1 will be NULL initially. After the user selects and uploads a 
		# file, it will be a data frame with 'name', 'size', 'type', and 'datapath' 
		# columns. The 'datapath' column will contain the local filenames where the 
		# data can be found.

		isolate(inFile <- input$file1)
	
		result <- readInput(inText,inFile)
		validate(
			need(!is.null(result),"Please input sequence(s)")
		)
	
		inSeq <- result$Seq
		SeqId <- result$Id
    
		pVec <- c()
	
		if(length(inSeq) > 50){
			# Create a Progress object
			progress <- shiny::Progress$new()
			# Make sure it closes when we exit this reactive, even if there's an error
			on.exit(progress$close())
		}
		for(i in 1:length(inSeq)){
			#try(if (nchar(inSeq[i]) < 96 ) stop("The length of input sequences must be >= 96 bp."))
			validate(
				need( nchar(inSeq[i]) >= 96, "The length of input sequence(s) must be >= 96 bp.")
			)
		
			if(length(inSeq) > 50){
				progress$set(message = "Calculating", value = 0)
				# Increment the progress bar, and update the detail text.
				progress$inc(1/length(inSeq), detail = paste("Finishing", i))
				# Pause for 0.1 seconds to simulate a long computation.
				Sys.sleep(0.1)
			}
		
			p <- codon_influence_cal(dnaseq = inSeq[i]);
			pVec <- rbind(pVec,unlist(p));
		}
	
		codon_influence <- data.frame(SeqId = SeqId,Seq = inSeq,pVec)

		output$downloadData1 <- downloadHandler(
			filename = function() { paste("codon_influence", '.txt', sep="") },
			content = function(file) { write.table(codon_influence, file,sep="\t",row.names=FALSE) }
		)
	
		num_col = ncol(codon_influence)
		num_row = nrow(codon_influence)
		if(num_row > 21){ num_row = 20; }
	
		codon_influence[1:num_row,c(1,3:10,(num_col-4):num_col)]
	
	})

	#####################################################################################################
	
	output$deltaG <- renderTable({ 
	
		if (is.null(v$data4)) return()
		isolate(inText <- input$text4)
		isolate(inFile <- input$file4)
		
		result <- readInput(inText,inFile)
	
		validate(
			need(!is.null(result),"Please input sequence(s)")
		)
		inSeq <- result$Seq
		SeqId <- result$Id
	
		isolate(len <- 3*input$num2)

		for(i in 1:length(inSeq)){
			validate(
				need( nchar(inSeq[i]) >= len, paste0("The length of input sequence(s) must be >= ", len, " bp." ))
			)
		}
		
		tmp <- folding_cal(substr(inSeq,1,len))
		deltaG <- data.frame(SeqId,inSeq,tmp[,2])
		colnames(deltaG) <- c("SeqID","Seq","Mnimum free energy")
		
		output$downloadData4 <- downloadHandler(
			filename = function() { paste("N_terminal_deltaG", '.txt', sep="") },
			content = function(file) { write.table(deltaG,file,sep="\t",row.names=FALSE) }
		)
	
		tmp <- head(deltaG)
		tmp[,2] <- substr(tmp[,2],1,len)
		tmp
	
	})
	
	#####################################################################################################
	
	output$design <- renderTable({
	
		if (is.null(v$data2)) return()
	
		isolate(inText <- input$text2)
		isolate(inFile <- input$file2)
	
		result <- readInput(inText,inFile)
		validate(
			need(!is.null(result),"Please input sequence(s)")
		)
	
		inSeq <- result$Seq
		SeqId <- result$Id
	
		target_array <- c()
		target_id <- c()
	
		isolate(method1 <- input$method1)
		isolate(method2 <- input$method2)
	
		# Create a Progress object
		progress <- shiny::Progress$new()
		# Make sure it closes when we exit this reactive, even if there's an error
		on.exit(progress$close())
		progress$set(message = "Designing", value = 0)
	
		for(i in 1:length(inSeq)){
			validate(
				need( nchar(inSeq[i]) >= 48, "The length of input sequence(s) must be >= 48 bp.")
			)
			# Increment the progress bar, and update the detail text.
			progress$inc(1/length(inSeq), detail = paste("Finishing", i))

			# Pause for 0.1 seconds to simulate a long computation.
			Sys.sleep(0.1)

			## for head optimization
			if(method1 == 'None'){ target_head <- substr(inSeq[i],1,48) }
			if(method1 == '6AA'){ target_head <- sixAA(substr(inSeq[i],1,48)) }
			if(method1 == '31C'){ target_head <- head_folding_optimize(substr(inSeq[i],1,48)); target_head <- target_head[,1] }
		
			## for tail optimization
			if(method2 == 'None'){ target_tail <- substr(inSeq[i],49,nchar(inSeq[i])) }
			if(method2 == '6AA'){ target_tail <- sixAA(substr(inSeq[i],49,nchar(inSeq[i])))  }
			if(method2 == '31C'){ target_tail <- tail_folding_optimize(substr(inSeq[i],49,nchar(inSeq[i]))) }

			target <- outer(target_head,target_tail,f)
			target <- toupper(as.vector(target))
		
			target_id <- c(target_id,paste(SeqId[i],"_head_",method1,"_tail_",method2,"_",1:length(target),sep=""))
			target_array <- c(target_array,target)
		}
	
		fasta <-list()
		for(j in 1:length(target_array)){
			a <- s2c(target_array[j])
			fasta <- c(fasta,list(a))
		}
	
		output$downloadData2 <- downloadHandler(
			filename = function() { paste("gene_design", '.fasta', sep="") },
			content = function(file) { write.fasta(fasta,target_id, file.out=file) }
		)
	
		num_seq <- length(target_array)
		if(num_seq > 11) {num_seq = 10; }
		tmp <- c()
		for(i in 1:num_seq){
			id <- paste(">",target_id[i],sep="")
			tmp <- c(tmp,id)
			tmp <- c(tmp,target_array[i])
		}
		tmp <- data.frame(tmp)
		colnames(tmp) <- ""
		tmp
	})
	
	###########################################################################################################
	
	output$noptimize <- renderTable({
	
		if (is.null(v$data3)) return()
	
		# Create a Progress object
		progress <- shiny::Progress$new()
		# Make sure it closes when we exit this reactive, even if there's an error
		on.exit(progress$close())
	
		#inputData <- v$data()
		#inText <- inputData$text3
		#inFile <- inputData$file3
	
		isolate(inText <- input$text3)
		isolate(inFile <- input$file3)
		
		result <- readInput(inText,inFile)
	
		validate(
			need(!is.null(result),"Please input sequence(s)")
		)
		inSeq <- result$Seq
		SeqId <- result$Id
	
		isolate(len <- 3*input$num)
	
		deltaG <- data.frame()
		for(i in 1:length(inSeq)){
			validate(
				need( nchar(inSeq[i]) >= 120, "The length of input sequence(s) must be >= 120 bp." )
			)
			dna_string <- substr(inSeq[i],1,len)
			fasta <- s2c(dna_string)
			syn_codons <- syncodons(splitseq(fasta))
			count = 1
			for(j in 1:(len/3)){
				count = count*length(unlist(syn_codons[j]))
			}
			validate(
				need( count <= 300000, "The total combinations seem too large (>300,000). Please input a smaller number")
			)
			combn <-syn_codons[[1]]	
			for(j in 2:(len/3)){
				combn <- outer(combn,syn_codons[[j]],f)
				combn <- as.vector(combn)
			}

			combn <- toupper(paste(combn,substr(inSeq[i],len+1,nchar(inSeq[i])),sep=""))

			progress$set(message = paste("Doing sequence",SeqId[i]), value = 0)
		
			combn1 <- substr(combn,1,18)
			tmp1 <- folding_cal(combn1)
			progress$inc(1/3, detail = paste("Finishing", 1))
			Sys.sleep(0.1)
		
			combn2 <- substr(combn,1,48)
			tmp2 <- folding_cal(combn2)
			progress$inc(1/3, detail = paste("Finishing", 2))
			Sys.sleep(0.1)
		
			combn3 <- substr(combn,1,120)
			tmp3 <- folding_cal(combn3)
			progress$inc(1/3, detail = paste("Finishing", 3))
			Sys.sleep(0.1)
		
			GC1 <- c();
			GC2 <- c();
			GC3 <- c();
			for(j in 1:length(combn)){
				GC1 <- c(GC1,GC(s2c(combn1[j])))
				GC2 <- c(GC2,GC(s2c(combn2[j])))
				GC3 <- c(GC3,GC(s2c(combn3[j])))
			}
		
			ids <- paste(SeqId[i],"_",1:length(combn),sep="")
			deltaG_tmp <- data.frame(ids,combn,tmp1[,2],GC1,tmp2[,2],GC2,tmp3[,2],GC3)
			colnames(deltaG_tmp) <- c("SeqID","Seq","MFE of 6","GC of 6","MFE of 16","GC of 16","MFE of 40","GC of 40")
			deltaG_tmp <- deltaG_tmp[order(-deltaG_tmp[,5],deltaG_tmp[,6]),]
			deltaG <- rbind(deltaG,deltaG_tmp)
		}
		
		output$downloadData3 <- downloadHandler(
			filename = function() { paste("N_terminal_optimize", '.txt', sep="") },
			content = function(file) { write.table(deltaG,file,sep="\t",row.names=FALSE) }
		)
	
		tmp <- head(deltaG)
		tmp[,2] <- substr(tmp[,2],1,len)
		tmp
	})

})
