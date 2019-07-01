library(colorRamps)
library(DT)
library(shiny)
library(treemap)
library(d3treeR)


# File names to be replaced by gui helper --------------------------------------

primer_file = "<primer_info>"
tax_file = "<tax_file>"
result_file = "<result_file>"
## for debugging:
#primer_file = "../PriSeT/src/tests/work/table/primer_info.csv"
#tax_file = "../PriSet/src/tests/library/root_123.tax"
#result_file = "../PriSeT/src/tests/work/table/results.csv"

# Load taxonomy file ------------------------------------
# load taxonomy from /library/*.tax io_cfg.tax_file()
tax = read.csv(tax_file)

# initial root node
root <- tax$p_taxid[1]
# augment data frame by subtree size (number of taxids in clade including root = tax$taxid)
p_taxids <- c(tax$p_taxid)
group_ctrs <- rep(1, nrow(tax))
group_ctrs

p_taxids
for (i in 1:length(tax$taxid)){
    ctr = length(p_taxids[p_taxids == tax$taxid[i]])
    group_ctrs[i] = group_ctrs[i] + ctr
}
group_ctrs
tax$clade_size <- group_ctrs
tax

# Load result columns taxid,fwd,rev,matches,coverage,ID_list -------------------
results <- read.csv(result_file, header = TRUE, sep = ",", row.names = NULL, col.names = c("taxid","fwd","rev","matches","coverage","ID_list"))

# Load and build (in results listed) primers legend ----------------------------
primer_info <- read.csv(primer_file)

# set of unique IDs 2nd column are fwdIDs
primer_fwd <- sort(unique(results$fwd))
names(primer_fwd) <- sprintf("primer %d", primer_fwd)

primer_rev <- sort(unique(results$rev))
stopifnot(primer_rev[0] != 0)
length(primer_rev[-1])
names(primer_rev) <- c("NA", sprintf("primer %d", primer_rev[-1]))

# Update data subset selection on GUI change -----------------------------------
update <- function(root_, coverage_, fwd_, rev_)
{
    print("Enter update ...")
    print("tax: ")
    print(tax)
    # TODO: make (or not coverage reactive)
    # direct descendants of root
    desc <- tax$taxid[tax$p_taxid == root_]
    print("input root: ")
    print(root_)
    print(desc)
    clade <- subset(tax, tax$p_taxid %in% desc, select=c(p_taxid, taxid))
    # add clade sizes of p_taxid
    clade_size_p <- tax$clade_size[tax$taxid == clade$p_taxid]

    print(clade_size_p)

    clade$clade_size <- clade_size_p
    names(clade) <- c("taxid_L1", "taxid_L2", "clade_size_L1")

    print("clade:")
    print(clade)

    # filter result table by clade resident taxids and selected primer combination
    print("result_sub")
    data_sub = results[(results$taxid %in% clade$taxid_L1) & (results$fwd == fwd_) & (results$rev == rev_), c("taxid", "matches", "coverage")]
    print(data_sub)

    # join clade and selected matching stats based on clade$taxid_L1 == data_sub$taxid
    data_sub_c <- merge(clade, data_sub, by.clade = c("taxid_L1"), by.data_sub = c("taxid"))[,c("taxid_L1", "taxid_L2", "clade_size_L1", "matches", "coverage")]
    print(data_sub_c)

    data_sub_c$heat <- data_sub_c$matches / data_sub_c$coverage

    data_sub_c <- data_sub_c[,!(names(data_sub_c) %in% c("matches", "coverage"))]

    data_sub_c$heat
    names(data_sub_c) <- c("taxid_L1", "taxid_L2", "clade_size_L1", "heat")
    return(data_sub_c)
}

# UI Layout --------------------------------------------------------------------
ui <- fluidPage(titlePanel("PriSeT"),
    sidebarLayout(
        sidebarPanel(
            # numeric input for taxid to jump without interactive parsing
            textInput(inputId = "root_taxid", label = "Choose a taxid", value = root),
            # slider for coverage
            #  sliderInput("range", "Range:", min = 1, max = 1000, value = c(200,500)),
            sliderInput(inputId = "coverage", label = "Lower bound of coverage rate (%)", value = c(20, 100), min = 0, max = 100),
            # Select box for forward primer sequence, todo: retrieve list and set first to default
            selectInput(inputId = "primer_fwd", label = "Choose a forward primer", choices = primer_fwd, selected = primer_fwd[1]),
            # Select box for reverse primer sequence
            selectInput(inputId = "primer_rev", label = "Choose a forward primer", choices = primer_rev, selected = primer_rev[1]),
            # output Info box with primerId, DNA sequence, Tm (Celsius)
            DT::datatable(primer_info,  width = '400px', options =
                list(   order = list(list(1, 'asc')),
                        autoWidth = TRUE,
                        columnDefs = #list(  list( width = '10%', targets = "_all"),
                                     #list(list(className = 'dt-center', targets = c(1,3))),
                                     list(list(targets=c(0), visible=TRUE, width='90'),
                                             list(targets=c(1), visible=TRUE, width='145')
                                         ),
                        pageLength = 5, lengthMenu = c(5,10,25,50),
                        deferRender = TRUE,
                        scrollX = TRUE,
                        scrollY = 400,
                        scrollCollapse = TRUE
                ),
                rownames = FALSE # suppress display of row IDs
            )
        ),
        mainPanel(
            plotOutput(outputId = "tree_plot")
        )
    )
)

# Server function --------------------------------------------------------------
server <- function(input, output)
{
    # Reactive elements --------------------------------------------------------
    sliderValues <- reactive({
        data.frame(
            Name = c("Root Taxid", "Coverage", "Forward Primer", "Reverse Primer"),
            Value = as.character(c( input$root_taxid,
                                    input$coverage,
                                    input$primer_fwd,
                                    input$primer_rev)
                                ),
            stringsAsFactors = FALSE)
    })

    # Output elements ----------------------------------------------------------
    output$primer_info <- DT::renderDataTable({
                                datatable(primer_info) %>% formatStyle(columns = c(1,2,3), width='10%')
                            })

    output$tree_plot <- renderPlot({
        data_show <- update(input$root_taxid, input$coverage, input$primer_fwd, input$primer_rev)
        tree <- treemap(
            data_show, index = c("taxid_L1", "taxid_L2"),
            vSize = "clade_size_L1",
            type = "dens",
            vColor = "heat",
            title = "Clade",
            palette = matlab.like2(20),
            range = c(0,1),
            title.legend = "primer coverage (%)"
        )
        d3tree2(tree,  rootname = paste("root", root))
    })
}

shinyApp(ui = ui, server = server)
