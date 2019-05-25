library(shiny)

# Add elements to app as arguments to fluidPage
ui <- fluidPage("Hello, World!"
    # *Input() functions,
    # *Output() functions
    # numeric input for taxid to jump without interactive parsing
    # todo get root taxid
    numericInput(inputID = "taxid", label = "Choose a taxid", value = 1),
    # slider for coverage
    sliderInput(inputID = "coverage", label = "Choose a coverage rate", value = 0.1, min = 0.0, max = 1.0),
    # Select box for forward primer sequence, todo: retrieve list and set first to default
    selectInput(inputID = "primer_fwd", label = "Choose a forward primer"),
    # Select box for reverse primer sequence
    selectInput(inputID = "primer_rev", label = "Choose a reverse primer")

    # output Info box with primerID, DNA sequence, Tm (Celsius)
)

server <- function(input, output) {}

shinyApp(ui = ui, server = server)
