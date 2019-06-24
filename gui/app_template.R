library(DT)
library(shiny)
library(treemap)
library(d3treeR)

##############################  Tree built ####################################
# load taxonomy from /library/*.tax io_cfg.tax_file()
#tax = read.csv("<tax_file>")
tax = read.csv("../src/tests/library/root_1234.tax")
# initial root node
root = tax$p_taxid[1]
# direct descendants of root
desc = tax$taxid[tax$p_taxid == root]
# augment data frame by subtree size (number of taxids in clade including root = tax$taxid)
ctrs <- table(tax$p_taxid)
group_ctrs <- rep(1, length(taxid))
for (i in 1:length(taxid)){
    ctr = ctrs[names(ctrs) == taxid[i]]
    if (length(ctr) > 0){
            group_ctrs[i] <- ctr
    }
}

tax$clade_size = group_ctrs
tax

clade <- subset(tax, tax$p_taxid %in% desc, select=c(p_taxid, taxid, clade_size))

# build treemap
tree = treemap(clade, index = c("p_taxid", "taxid"), vSize = "clade_size", type = "index")

# make it interactive ("rootname" becomes the title of the plot):
#tree_inter = d3tree2(tree,  rootname = paste("root", root))

##############################  Color map ####################################
#taxid,fwd,rev,matches,coverage,ID_list
#data = read.csv("<result_file>", row.names=NULL)
data = read.csv("../src/tests/work/table/results.csv", row.names=NULL)
# assign manually headers
names(data) <- c("taxid", "fwd", "rev", "matches", "coverage", "ID_list")
head(data)

#primer_info = read.csv("<primer_info>")
primer_info = read.csv("../src/tests/work/table/primer_info.csv")
names(primer_info) <- c("kmer_ID", "kmer_sequence", "Tm")
head(primer_info)

# set of unique IDs 2nd column are fwdIDs
primer_fwd <- sort(unique(data$fwd))
names(primer_fwd) <- sprintf("primer %d", primer_fwd)

primer_rev <- sort(unique(data$rev))
stopifnot(primer_rev[0] != 0)
length(primer_rev[-1])
names(primer_rev) <- c("NA", sprintf("primer %d", primer_rev[-1]))

heat <- rep(0, length(clade))
for (i in 1:length(clade)){
    row <- data[data$taxid == clade$taxid[i] && data$fwd == primer_fwd[1] && data$rev == 0]
    heat[i] <- row$matches / row$coverage
}

clade$heat = heat

############################## Build UI ####################################
ui <- fluidPage(titlePanel("PriSeT"),
        fluidRow(
            column(6, # see https://shiny.rstudio.com/tutorial/written-tutorial/lesson3/
                # numeric input for taxid to jump without interactive parsing
                # todo get root taxid
                textInput(inputId = "root_taxid", label = "Choose a taxid", value = init_taxid),
                # slider for coverage
                sliderInput(inputId = "coverage", label = "Choose a coverage rate", value = 0.1, min = 0.0, max = 1.0),
                # Select box for forward primer sequence, todo: retrieve list and set first to default

                selectInput(inputId = "primer_fwd", label = "Choose a forward primer", choices = primer_fwd, selected = primer_fwd[1]),
                # Select box for reverse primer sequence
                selectInput(inputId = "primer_rev", label = "Choose a forward primer", choices = primer_rev, selected = primer_rev[1]),
                # output Info box with primerId, DNA sequence, Tm (Celsius)
                #renderTable(expr, striped = FALSE, hover = FALSE, bordered = FALSE, rownames = FALSE, colnames = TRUE, digits = #NULL, outputArgs = list())
                #tableOutput(outputId = "primer_info", primer_info$kmer_ID, primer_info$kmer_sequence, primer_info$Tm
                DT::datatable(primer_info, options = list(
                    order = list(list(1, 'asc')),
                    autoWidth = TRUE,
                    columnDefs = list(list(width = '25px', targets = c(0, 1, 2, 3)), list(className = 'dt-center', targets = c(0,1,3))),
                    pageLength = 10,
                    lengthMenu = c(10,25,50)))
            ),
            column(6,
                plotOutput(outputId = "treecircle")
            )
        )
)

server <- function(input, output)
{
    output$primer_fwd <- renderUI({
        primer_fwd <- input$primer_fwd
    })
    output$primer_rev <- renderUI({
        primer_rev <- input$primer_rev
    })
    output$primer_info <- DT::renderDataTable(return(primer_info), bordered = FALSE, spacing = c("s"), rownames = TRUE)
    output$treecircle <- renderUI({ # or renderPlot
        root_taxid <- input$root_taxid
    })
}

shinyApp(ui = ui, server = server)
