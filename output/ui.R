library(shiny)

#library(BiocManager)
#options(repos = BiocManager::repositories())
#library(rsconnect)
#deployApp(appName = 'MS_broad',account = 'malhotralab')

# Define UI for application
ui <- fluidPage(
    titlePanel("Brain cell type cis-eQTLs"),
    p("This shinyApp displays cis-eQTLs, cell type specific cis-eQTLs, and colocalization results in 8 cell types from the central nervous system"),
    tabsetPanel(
      tabPanel("cis-eQTLs", fluid = TRUE,

               sidebarLayout(
                 sidebarPanel(
                     selectInput("select_cell_type", label = h3("Cell type"),
                                                    choices = list("All" = 'All',
                                                                    "Astrocytes" = 'Astrocytes',
                                                                   "Endothelial cells" = 'Endothelial cells',
                                                                   "Excitatory neurons" = 'Excitatory neurons',
                                                                   "Inhibitory neurons" = 'Inhibitory neurons',
                                                                   "Microglia" = 'Microglia',
                                                                   "Oligodendrocytes" = 'Oligodendrocytes',
                                                                   "OPCs / COPs" = 'OPCs / COPs',
                                                                   "Pericytes" = 'Pericytes'
                                                    ),
                                                    selected = 'All'),
                     numericInput("num",h3("FDR threshold"),value = 0.05),
                     DT::dataTableOutput("mytable")),
                 mainPanel(
                     plotOutput('plot_eqtl',width="70%")
                 )
               )
      ),
      tabPanel("Cell type specific cis-eQTLs", fluid = TRUE,

               sidebarLayout(
                   sidebarPanel(
                       selectInput("select_model", label = h3("Statistical model"),
                                   choices = list("Aggregate" = 'Aggregate',
                                                  "At least one" = 'At least one',
                                                  "All" = 'All'
                                   ),
                                   selected = 'Aggregate'),
                       numericInput("num2",h3("FDR threshold"),value = 0.05),
                       DT::dataTableOutput("mytable2")
                   ),
                   mainPanel(
                       plotOutput('plot_eqtl_spe',height=700)
                   )
               )
      ),
      tabPanel("Colocalization", fluid = TRUE,

               sidebarLayout(
                   sidebarPanel(
                       selectInput("select_trait", label = h3("Trait"),
                                   choices = list("Alzheimer's disease" = 'ad',
                                                  "Parkinson's disease" = 'pd',
                                                  "Schizophrenia" = 'scz',
                                                  "Multiple sclerosis" = 'ms',
                                                  "Multiple sclerosis (GTEx & DICE)" = 'ms_gtex_dice'
                                   ),
                                   selected = 'AD'),
                       sliderInput("num3", label = h3("Posterior probability"), min = 0.5, max = 1, value = 0.7),
                   ),
                   mainPanel(
                       plotOutput('plot_coloc',width="70%")
                   )
               )
      )
    )
)
