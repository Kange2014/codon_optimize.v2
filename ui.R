#-------------------------------------------------------------------------
#  This application is free to anyone. 
#  You can  use, modify and/ or redistribute these codes
#
#  LUZE@novonordisk
#  Feb 20, 2017
#-------------------------------------------------------------------------

library(shiny)
library(shinydashboard)

header <- dashboardHeader(title = "Codon optimization")

sidebar <- dashboardSidebar(
    hr(),
	sidebarMenu(id = 'tabs',
	  menuItem("DoE + Machine learning",  icon = icon("line-chart"), 
                       menuSubItem("Design of Experiments", tabName = "doe", icon = icon("angle-right")), selected=TRUE,
                       menuSubItem("Linear regression model", tabName = "lrm", icon = icon("angle-right"))
              ),
      menuItem("N-terminal optimize", tabName = "nterminal",icon = icon("table")),
	  menuItem("Preferred codon replace", tabName = "genedesign",icon = icon("file-text-o")),
	  menuItem("Codon influence calculate", tabName = "codoninflu",icon = icon("list")),
	  menuItem("Minimum free energy calculate", tabName = "mfe", icon = icon("calculator")),
	  menuItem("Back translate",icon = icon("external-link"),href = "Taihe.v2/rev_trans.html"), 
	  menuItem("ReadMe", tabName = "readme", icon = icon("mortar-board")),
      menuItem("About", tabName = "about", icon = icon("question"))
    )
)

body <- dashboardBody(
	tabItems(
	  tabItem("doe",
	   uiOutput('doe')
	  ),
	  
	  tabItem("lrm",
	   uiOutput('regress')
	  ),
	  
	  tabItem("codoninflu",
		# Copy the line below to make a text input box
		#textAreaInput("text1", label = "Paste your sequence(s) here:", value = ""),

		fluidRow(style = "padding-bottom: 20px;",
			column(12, uiOutput('resetable_input1')),
			column(12, actionButton("reset_input1", "Clear"))
		),
		
		#HTML("<button id='reset_input1' class='action-button clearButton'>Clear</button>"),
		
		fileInput('file1', 'Or upload a file in FASTA format'),
		tags$script('$( "#file1" ).on( "click", function() { this.value = null; });'),
		hr(),
		#submitButton("Calculate codon influence"),
		fluidRow(
			column(4,actionButton("calculate1", "Codon influence calculate",icon("paper-plane"), style="color: #fff; background-color: #337ab7; border-color: #2e6da4")),
			column(4,downloadButton('downloadData1', 'Download'))
		),
		hr(),
		tableOutput('parameters')
	  ),
	  
	  tabItem("mfe",
		fluidRow(style = "padding-bottom: 20px;",
			column(12, uiOutput('resetable_input4')),
			column(12, actionButton("reset_input4", "Clear"))
		),
		
		#HTML("<button id='reset_input1' class='action-button clearButton'>Clear</button>"),
		
		fileInput('file4', 'Or upload a file in FASTA format'),
		tags$script('$( "#file4" ).on( "click", function() { this.value = null; });'),
		
		numericInput("num2", label = ("N-terminal amino acid length for delta G calculating:"), value = 16),
		
		hr(),
		fluidRow(
			column(4,actionButton("calculate4", "Delta G calculate",icon("paper-plane"), style="color: #fff; background-color: #337ab7; border-color: #2e6da4")),
			column(4,downloadButton('downloadData4', 'Download'))
		),
		hr(),
		tableOutput('deltaG')
	  ),
	  
	  tabItem("genedesign",
	 		fluidRow(style = "padding-bottom: 20px;",
			column(12, uiOutput('resetable_input2')),
			column(12, actionButton("reset_input2", "Clear"))
		),
			
		fileInput('file2', 'Or upload a file in FASTA format'),
		tags$script('$( "#file2" ).on( "click", function() { this.value = null; });'),
		
		fluidRow(
			box(
				width = 4, status = "primary",
				radioButtons("method1", "Head optimize method:",
							list("None" = "None",
								"6AA Method" = "6AA",
								"31C Method" = "31C"
							), selected = "31C")
			),
			box(
				width = 4, status = "success",
				radioButtons("method2", "Tail optimize method:",
							list("None" = "None",
								"6AA Method" = "6AA",
								"31C Method" = "31C"
								), selected = "6AA")
			)
		),
		#submitButton("Gene design"),
		fluidRow(
			column(4,actionButton("calculate2", "Gene design",icon("paper-plane"), style="color: #fff; background-color: #337ab7; border-color: #2e6da4")),
			column(4,downloadButton('downloadData2', 'Download'))
		),
		
		htmlOutput("design")
		#verbatimTextOutput('design')
		#textOutput('design')
	  ),
	  
	  tabItem(tabName = "nterminal",
		fluidRow(style = "padding-bottom: 20px;",
			column(12, uiOutput('resetable_input3')),
			column(12, actionButton("reset_input3", "Clear"))
		),

		fileInput('file3', 'Or upload a file in FASTA format'),
		tags$script('$( "#file3" ).on( "click", function() { this.value = null; });'),
		
		numericInput("num", label = ("N-terminal amino acid length:"), value = 5),
		
		#submitButton("N-terminal optimize"),
		fluidRow(
			column(4,actionButton("calculate3", "N-terminal optimize",icon("paper-plane"), style="color: #fff; background-color: #337ab7; border-color: #2e6da4")),
			column(4,downloadButton('downloadData3', 'Download'))
		),
		tags$hr(),
		tableOutput('noptimize')
	  ),
	  tabItem(tabName = "readme",
            includeMarkdown("readMe.Rmd")
		),
	  tabItem(tabName = "about",
            includeMarkdown("about.Rmd")
		)
    )
)

dashboardPage(
  skin = "purple",
  header,
  sidebar,
  body
)