#
# Server logic of the lambdaPrimeR Shiny app
#

library(shiny)
library(lambdaPrimeR)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  # # Reactive object: input sequence type ----
  # sequence_type <- reactive({
  #   
  # })
  # Reactive object: input sequence ---- 
  input_sequence <- eventReactive(input$load_seq_button, {
    read_sequences(input$input_seq, input_type = input$sequence_type)
  })
  
  observeEvent(input$input_seq, {
    output$sequence_table <- DT::renderDataTable({
      input_sequence()
    })
  })
  # output$pasted_seq <- renderText({input$input_seq})

})
