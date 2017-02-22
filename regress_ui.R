################################################################
# Regression - UI
################################################################

reg_args <- as.list(formals(regress))

## list of function inputs selected by user
reg_inputs <- reactive({
  ## loop needed because reactive values don't allow single bracket indexing
  
  inText = input$string
  tmpseq <- unlist(strsplit(inText,"\r?\n"))
  start_line = 1
  if(grepl("^>",inText)){ start_line = 2 }
  ori_seq <- paste(tmpseq[start_line:length(tmpseq)],collapse="")
  
  inFile <- input$dataset
  sep = ","
  if(input$dataType == "tab") sep = "\t"

  filename <- basename(inFile$name)
  fext <- tools::file_ext(filename) %>% tolower
  
  validate(
			need(!(fext %in% c("xls","xlsx")),"### It does not support excel files. Please save the data as a csv/tab file and try again.")
		)
  if(fext == "csv"){
	validate(need(input$dataType == fext,"### The input file data type is not consistent with the selected data type.\n\n### Please select the same data type file.") )
  }
  
  else{
	validate(need(input$dataType == "tab","### The input file data type is not consistent with the selected data type.\n\n### Please select the same data type file."))
  }
  
  rawdata <- read.table(inFile$datapath, header = input$man_header, sep = sep, comment.char = "", quote = "\"", fill = TRUE, stringsAsFactors = FALSE)			
  dat <- getData(rawdata,ori_seq, input$man_met)
  
  for (i in r_drop(names(reg_args)))
    reg_args[[i]] <- input[[paste0("reg_",i)]]

  reg_args$trainingdata <- dat$trainingdata
  reg_args$testingdata <- dat$testingdata
  
  reg_args
})

reg_sum_args <- as.list(formals(summary.regress))

## list of function inputs selected by user
reg_sum_inputs <- reactive({
  ## loop needed because reactive values don't allow single bracket indexing
  for (i in names(reg_sum_args))
    reg_sum_args[[i]] <- input[[paste0("reg_",i)]]
  reg_sum_args
})

reg_plot_regress_args <- as.list(formals(plot.regress))

## list of function inputs selected by user
reg_plot_regress_inputs <- reactive({
  ## loop needed because reactive values don't allow single bracket indexing
  for (i in names(reg_plot_regress_args))
    reg_plot_regress_args[[i]] <- input[[paste0("reg_",i)]]
	
  reg_plot_regress_args
})


reg_pred_args <- as.list(formals(predict.regress))

## list of function inputs selected by user
reg_pred_inputs <- reactive({
  ## loop needed because reactive values don't allow single bracket indexing
  for (i in names(reg_pred_args))
    reg_pred_args[[i]] <- input[[paste0("reg_",i)]]
  
  if (is_empty(input$pred_data)) dat <- ""
  else dat <- dataNumeric(input$string,input$pred_data)
  reg_pred_args$pred_data <- dat
  
  reg_pred_args
})

output$ui_fileUpload <- renderUI({
  req(input$dataType)
  
  fileInput('dataset', '', multiple = FALSE,
      accept = c('text/csv','text/comma-separated-values',
                 'text/tab-separated-values', 'text/plain','.csv','.tsv','.tab','.txt'))
 }
)

## data ui and tabs
output$ui_regress <- renderUI({
  #req(input$dataset)
  
  data_types_in <- c("csv" = "csv","tab/txt" = "tab")
  
  tagList(
    wellPanel(
      actionButton("reg_run", "Estimate", width = "100%")
    ),
	
    conditionalPanel(condition = "input.tabs_regress == 'Predict'",
      wellPanel(
          checkboxInput("reg_pred_new", "Predict new codons", state_init("reg_pred_new", FALSE)),
          conditionalPanel("input.reg_pred_new == true",
            textInput("pred_data", label = "Paste your new codons here including initial ATG", value = "")
          )
      )
    ),
    
    wellPanel(
	  textAreaInput("string", label = "Paste your original coding sequence here:", value = "")
    ),
	
	wellPanel(
      selectInput("dataType", label = "Load data of type:", data_types_in, selected = "tab"),
      with(tags, table(td(checkboxInput('man_header', 'Header', TRUE)),
		td(HTML("&nbsp;&nbsp;")),
		td(checkboxInput("man_met", "Including initial Met",FALSE)))
		),
      uiOutput("ui_fileUpload")
    )
  )
})

# output is called from the main radiant ui.R
output$regress <- renderUI({

    register_print_output("summary_regress", ".summary_regress")
	register_plot_output("plot_regress", ".plot_regress")
    register_print_output("predict_regress", ".predict_print_regress")

    # two separate tabs
    reg_output_panels <- tabsetPanel(
      id = "tabs_regress",
      tabPanel("Summary",
		downloadLink("dl_reg_training", "", class = "fa fa-download alignright"), br(),
		verbatimTextOutput("summary_regress"),
		plotOutput("plot_regress", width = "100%", height = "100%")
		),
      tabPanel("Predict",
        downloadLink("dl_reg_pred", "", class = "fa fa-download alignright"), br(),
        verbatimTextOutput("predict_regress")
      )
    )
    stat_tab_panel(
                  tool = "Linear regression (LM)",
				  data = NULL,
                  tool_ui = "ui_regress",
                  output_panels = reg_output_panels)
})

.regress <- eventReactive(input$reg_run, {
  do.call(regress, reg_inputs())
})

.summary_regress <- reactive({
  if (not_pressed(input$reg_run)) return("** Press the Estimate button to estimate the model **\n\n** Upload file should with similar format to the following example: **\n
No.	A	S	C	P	Expression
1	GCG	TCG	TGT	CCT	0.284097793
2	GCG	TCT	TGT	CCC	0.598907223
3	GCG	AGT	TGC	CCA	0.115092922
4	GCA	TCG	TGC	CCA	0.269014562
5	GCA	TCA	TGT	CCC	2.486855177
\n")
  do.call(summary.regress, c(list(object = .regress()), reg_sum_inputs()))
})

.predict_regress <- reactive({
  if (not_pressed(input$reg_run)) return("** Press the Estimate button to estimate the model **")
  
  withProgress(message = "Generating predictions", value = 1, {
    do.call(predict, c(list(object = .regress()), reg_pred_inputs()))
  })
})

.predict_print_regress <- reactive({
  .predict_regress() %>% {if (is.character(.)) cat(.,"\n") else print(.)}
})

.plot_regress <- reactive({
  if (not_pressed(input$reg_run)) return(invisible())

  do.call(plot.regress, c(list(object = .regress()), reg_plot_regress_inputs()))
})


## reset prediction settings when the dataset changes
##observeEvent(input$dataset, {
 
## updateSelectInput(session = session, inputId = "reg_predict", selected = "none")
## })


output$dl_reg_training <- downloadHandler(
  filename = function() { "reg_training.csv" },
  content = function(file) {
    if (pressed(input$reg_run)) {
      .regress()$model$training_pred %>%
        write.csv(file = file, row.names = FALSE)
    } else {
      cat("No output available. Press the Estimate button to generate results", file = file)
    }
  }
)

output$dl_reg_pred <- downloadHandler(
  filename = function() { "reg_predictions.csv" },
  content = function(file) {
    if (pressed(input$reg_run)) {
      .predict_regress() %>%
        write.csv(file = file, row.names = FALSE)
    } else {
      cat("No output available. Press the Estimate button to generate results", file = file)
    }
  }
)