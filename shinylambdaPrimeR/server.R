#
# Server logic of the lambdaPrimeR Shiny app
#

library(shiny)
library(lambdaPrimeR)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  # Reactive object: input sequence ---- 
  input_sequence <- eventReactive(input$load_seq_button, {
    read_sequences(input$input_seq)
  })
  
  output$pasted_seq <- renderText({input$input_seq})

})
