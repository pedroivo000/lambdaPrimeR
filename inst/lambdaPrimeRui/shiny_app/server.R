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
    sequence_inputs = tibble::tibble()
  )
  
  # Reactive object: input sequence ---- 
  input_sequence <- eventReactive(input$load_seq_button, {
    read_sequences(input$input_seq, input_type = input$sequence_type)
  })
  
  # Reactive event: "Load sequence" button click ----
  observeEvent(input$load_seq_button, {
    # Update globals$sequence_inputs:
    globals$sequence_inputs <- collect_input_shiny(input_sequence(), globals$sequence_inputs)
    output$sequence_table <- DT::renderDataTable({
      globals$sequence_inputs
    })
  })

})
