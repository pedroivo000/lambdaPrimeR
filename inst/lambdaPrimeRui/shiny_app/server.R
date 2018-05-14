#
# Server logic of the lambdaPrimeR Shiny app
#
# 
# library(shiny)
# library(lambdaPrimeR)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  # Global variables ----
  globals <- reactiveValues(
    # Dataframe with sequence inputs
    sequence_inputs = tibble::tibble(
      id  = integer(),
      path = character(), 
      header = character(),
      name = character(), 
      seq = character(), 
      type = character(), 
      length = integer()
    )
  )
  
  # Reactive object: input sequence ---- 
  input_sequence <- eventReactive(input$load_seq_button, {
    read_sequences_fasta(input$input_seq, input_type = input$sequence_type)
  })
  
  # Reactive event: "Load sequence" button click ----
  observeEvent(input$load_seq_button, {
    # Update globals$sequence_inputs:
    globals$sequence_inputs <- rbind(input_sequence(), globals$sequence_inputs)
    output$sequence_table <- DT::renderDataTable({
      globals$sequence_inputs %>%
        select(-path, -seq)
    })
  })
  
  # Reactive event: "design primers" button click:
  observeEvent(input$design_primers_button, {
    #Create run object from sequence inputs:
    run <- read_sequences_csv(globals$sequence_inputs)
  })
  

})
