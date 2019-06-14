library(shiny)

# load taxonomy from /library/*.tax io_cfg.tax_file()
#taxonomy = read.csv("<tax_file>", header = TRUE)
taxonomy = read.csv("/Users/troja/priset/library/root_335928.tax", header = TRUE)

results = read.csv("/Users/troja/priset/")

primer_fwd <- c("<None>", "primer 1", "primer 2", "primer 3")
primer_rev <- c("<None>", "primer 2", "primer 4", "primer 5")

# variable name will be column header
primer = c("primer 1", "primer 2", "primer 3", "primer 4", "primer 5")
sequence = c("AAAAAAA", "CCCCCCC", "TTTTTTTT", "ATATATATAT", "GCATATATAA")
Tm = c(43.4, 39.0, 40.1, 38.4, 42.7)

df = data.frame(primer, sequence, Tm)

root_taxid = 1
#head(df)




# Add elements to app as arguments to fluidPage
ui <- fluidPage(titlePanel("PriSeT"),
        fluidRow(
            column(6,
                # numeric input for taxid to jump without interactive parsing
                # todo get root taxid
                numericInput(inputId = "taxid", label = "Choose a taxid", value = root_taxid),
                # slider for coverage
                sliderInput(inputId = "coverage", label = "Choose a coverage rate", value = 0.1, min = 0.0, max = 1.0),
                # Select box for forward primer sequence, todo: retrieve list and set first to default
                selectInput(inputId = "primer_fwd", label = "Choose a forward primer", c("Choose one" = "", primer_fwd) ),
                # Select box for reverse primer sequence
                selectInput(inputId = "primer_rev", label = "Choose a forward primer", c("Choose one" = "", primer_rev) ),

                # output Info box with primerId, DNA sequence, Tm (Celsius)
                tableOutput(outputId = "primer_info")
            ),
            column(6,
                plotOutput("taxonomy")
            )
        )
)

server <- function(input, output)
{

    output$primer_info <- renderTable(return(df), bordered = TRUE, spacing = c("m"), width = "auto")
}

shinyApp(ui = ui, server = server)
