#
# UI for lambdaPrimeR Shiny app
#

library(shiny)

# Define UI:
shinyUI(navbarPage("lambdaPrimeR", fluid = T,
  

  # Primer design tab -------------------------------------------------------
  tabPanel("Design primers",
    
    # Sidebar with input and output options ----
    sidebarLayout(
      
      # Panel in sidebar with sequence inputs ----
      sidebarPanel(
        
        # Sidebar panel title ----
        titlePanel("Sequence input"),
        
        # Input: input sequence type selection ----
        selectInput('sequence_type', label = 'Sequence type', 
                    choices = list('Vector' = 'vector', 'Target' = 'target'),
                    selected = 'vector'),
        
        # Input: Paste sequence ----
        textAreaInput("input_seq", label = "Paste input sequence below:", 
                      width = '100%'),
        
        # Input: Upload file ----
        fileInput("input_file", "Or upload FASTA file with sequence:", width = '100%'), 
        
        # Input: load sequences button ---- 
        actionButton('load_seq_button', 'Load sequence')
      ), 
    
      # Main panel displaying outputs ----
      mainPanel(
        fluidPage(
        
          # Ouput: Input sequence print ----
          verbatimTextOutput("pasted_seq") 
        )
      )
    )
  )
))
