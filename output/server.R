# You can run the application by clicking 'Run App' above (RStudio)
# Author: Julien Bryois
# Date: 4.08.2021

#Upload to Shinyapps.io
#library(BiocManager)
#options(repos = BiocManager::repositories())
#library(rsconnect)
#setwd("/Volumes/GoogleDrive/My Drive/Projects/Brain_cell_type_eQTL")
#deployApp(account='malhotralab')

library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(rhdf5)
library(DT)
library(forcats)
library(ComplexHeatmap)
library(viridis)

# Option to increase max size of the file to be loaded
options(shiny.maxRequestSize=90*1024^2)

# Define server logic
shinyServer(function(input, output) {

    #Set path to hdf5 data
    hdF5_file_path <- 'data.h5'

    #Read eQTL list
    eqtl_df <- h5read(hdF5_file_path,'/eqtl_results/eqtl_results_all')

    #Read eQTL list
    eqtl_spe_df <- h5read(hdF5_file_path,'/eqtl_results/eqtl_results_specific')

    #Get eQTL df
    eqtl_list <- reactive({
        req(input$select_cell_type)
        req(input$num)
        if(input$select_cell_type!='All'){
          d <- eqtl_df %>% filter(cell_type%in%input$select_cell_type)
        } else{
          d <- eqtl_df
        }
        d <- filter(d,adj_p<=input$num)
        return(d)
    })

    #Get eQTL spe df
    eqtl_spe_list <- reactive({
        req(input$select_model)
        req(input$num2)
        d <- eqtl_spe_df
        if (input$select_model=='Aggregate'){
            d <- d[c(1,2,3,4)] %>%
                setNames(c('cell_type','gene','SNP','p_adj')) %>%
                arrange(p_adj)
        }
        if (input$select_model=='At least one'){
            d <- d[c(1,2,3,5)] %>%
                setNames(c('cell_type','gene','SNP','p_adj')) %>%
                arrange(p_adj)
        }
        if (input$select_model=='All'){
            d <- d[c(1,2,3,6)] %>%
                setNames(c('cell_type','gene','SNP','p_adj')) %>%
                arrange(p_adj)
        }
        d <- filter(d,p_adj<=input$num2)
        return(d)
    })

    #Display list of DE genes
    output$mytable = renderDT({
        eqtl_list() %>%
            separate(gene,into=c('symbol','ensembl'),sep='_') %>%
            dplyr::select(-ensembl)
    }, selection = list(mode='single',selected=1L))

    #Display list of eqtl spe genes
    output$mytable2 = renderDT({
        eqtl_spe_list() %>%
            separate(gene,into=c('symbol','ensembl'),sep='_') %>%
            dplyr::select(-ensembl)
    }, selection = list(mode='single',selected=1L))

   #Ploting function cis-eQTLs
   plot_eqtl <- function(cell_type_name,gene_name,snp_name,p_adj){

        gene_short <- gsub('_.+','',gene_name)

        genotype <- h5read(hdF5_file_path, paste0("genotype/",snp_name)) %>% filter(!is.na(genotype))
        expression <- h5read(hdF5_file_path, paste0("expression/",gene_name)) %>%
          dplyr::rename(individual=individual_id)

        d <- inner_join(genotype,expression,by='individual') %>%
            mutate(genotype_label = factor(case_when(
                genotype == 0 ~ paste0(REF,REF),
                genotype == 1 ~ paste0(REF,ALT),
                genotype == 2 ~ paste0(ALT,ALT)
            ),levels=c(unique(paste0(REF,REF)),unique(paste0(REF,ALT)),unique(paste0(ALT,ALT))))) %>%
            filter(cell_type%in%cell_type_name)

        p <- ggplot(d, aes(genotype_label,log2_cpm,col=genotype_label)) +
            ggbeeswarm::geom_quasirandom(size=0.5) +
            geom_boxplot(alpha=0.05,aes(group=genotype_label),outlier.shape = NA) +
            theme_classic() +
            theme(legend.position='none',text=element_text(size=16), plot.title = element_text(size=16)) +
            xlab('Genotype') +
            ylab('Expression (cpm) (log2+1)') +
            ggtitle(paste0(gene_short,' - ',snp_name,' - ',cell_type_name,' - ',p_adj))
        return(p)
   }

   #Ploting functions specific cis-eQTLs
   plot_eqtl_specific <- function(cell_type_name,gene_name,snp_name,p_adj){

       gene_short <- gsub('_.+','',gene_name)

       genotype <- h5read(hdF5_file_path, paste0("genotype/",snp_name)) %>% filter(!is.na(genotype))
       expression <- h5read(hdF5_file_path, paste0("expression/",gene_name)) %>%
         dplyr::rename(individual=individual_id)

       d <- inner_join(genotype,expression,by='individual') %>%
           mutate(genotype_label = factor(case_when(
               genotype == 0 ~ paste0(REF,REF),
               genotype == 1 ~ paste0(REF,ALT),
               genotype == 2 ~ paste0(ALT,ALT)
           ),levels=c(unique(paste0(REF,REF)),unique(paste0(REF,ALT)),unique(paste0(ALT,ALT))))) %>%
           mutate(cell_type=fct_relevel(cell_type,cell_type_name))


       p <- ggplot(d, aes(genotype_label,log2_cpm,col=genotype_label)) +
           ggbeeswarm::geom_quasirandom(size=0.5) +
           geom_boxplot(alpha=0.05,aes(group=genotype_label),outlier.shape = NA) +
           facet_wrap(~cell_type)+
           theme_classic() +
           theme(legend.position='none',text=element_text(size=16), plot.title = element_text(size=16)) +
           xlab('Genotype') +
           ylab('Expression (cpm) (log2+1)') +
           ggtitle(paste0(gene_short,' - ',snp_name,' - ',cell_type_name,' - ',p_adj))
       return(p)
   }

   #Plot function for coloc

   coloc_heatmap <- function(d,threshold){

       coloc_all <- d

       selected_trait = ifelse(input$select_trait=='ms_gtex_dice','ms',input$select_trait)

       metabrain_s1 <- readxl::read_xlsx('media-1.xlsx',skip=1,sheet = 2) %>%
         mutate(disease=case_when(
           outcome=="Alzheimer’s disease" ~ 'ad',
           outcome=="Multiple sclerosis" ~ 'ms',
           outcome=="Parkinson’s disease" ~ 'pd',
           outcome=="Schizophrenia" ~ 'scz',
           TRUE ~ outcome
         )) %>%
         dplyr::select(gene,disease) %>%
         filter(disease==selected_trait)

       coloc_all <- coloc_all %>% mutate(symbol=ifelse(symbol%in%metabrain_s1$gene,paste0(symbol,'*'),symbol))


       #Get absolute value of effect size of the GWAS
       coloc_all$beta_top_GWAS <- abs(coloc_all$beta_top_GWAS)

       coloc_all_matrix <- dplyr::select(coloc_all,tissue,symbol,PP.H4.abf) %>%
           unique() %>%
           spread(tissue,PP.H4.abf,fill=0) %>%
           column_to_rownames('symbol')

       if(threshold>max(coloc_all$PP.H4.abf)){
         threshold=max(coloc_all$PP.H4.abf)
       }

       coloc_all_format_per_gene <- coloc_all %>%
           group_by(ensembl) %>%
           filter(PP.H4.abf==max(PP.H4.abf),PP.H4.abf>=threshold) %>%
           ungroup() %>%
           dplyr::rename(symbol_coloc=symbol) %>%
           filter(!is.na(symbol_coloc)) %>%
           arrange(-PP.H4.abf) %>%
           mutate(risk=case_when(
               direction ==1 ~ 'Up',
               direction ==-1 ~ 'Down',
               direction ==0 ~ 'Isoform'
           )) %>%
           column_to_rownames('symbol_coloc') %>%
           dplyr::rename(LOEUF=oe_lof_upper_bin)

       coloc_all_matrix_subset <- coloc_all_matrix[apply(coloc_all_matrix,1,function(x) any(x>=threshold)),] %>% as.matrix()

       coloc_all_format_per_gene_direction <- coloc_all_format_per_gene[rownames(coloc_all_matrix_subset),] %>%
           dplyr::select(closest_gene,risk,LOEUF,`Beta GWAS`=beta_top_GWAS)

       ha = HeatmapAnnotation(df = coloc_all_format_per_gene_direction[,2:3,drop=FALSE],
                              which ='row',
                              col = list('risk' = c("Up"="#E31A1C","Down"="#1F78B4","Isoform"="darkorange"),
                                         'LOEUF' = circlize::colorRamp2(c(0,10),c('white','springgreen4'))),
                              show_annotation_name = c('risk' = TRUE,'LOEUF'=TRUE),
                              annotation_name_rot = 90
       )

       hb = HeatmapAnnotation(df = coloc_all_format_per_gene_direction[,4,drop=FALSE],
                              which ='row',
                              col = list(circlize::colorRamp2(c(ifelse(nrow(coloc_all_format_per_gene_direction)==1,0,min(coloc_all_format_per_gene_direction[[4]],na.rm=T)),
                                                                max(coloc_all_format_per_gene_direction[[4]],na.rm=T)),c('white','palevioletred3'))) %>%
                                  setNames(colnames(coloc_all_format_per_gene_direction)[4]),
                              annotation_name_rot = 45
       )

       ht <- Heatmap(coloc_all_matrix_subset,col=viridis(100),
                     right_annotation=ha,column_names_rot = 45,
                     left_annotation = hb,
                     row_split = coloc_all_format_per_gene_direction[, 1],
                     row_title_rot = 0,
                     cluster_rows = FALSE,
                     layer_fun = function(j, i, x, y, width, height, fill) {
                         grid.text(sprintf("%.2f", pindex(coloc_all_matrix_subset, i, j)),
                                   x, y, gp = gpar(fontsize = 10,col="black"))
                     },
                     heatmap_legend_param = list(title="PP",at = c(0,0.5,1),
                                                 lables = c(0,0.5,1)))
       draw(ht,heatmap_legend_side = "right")
   }

   #Do the plotting for cis-eQTLs
   output$plot_eqtl <- renderPlot({
       req(input$mytable_rows_selected)
       gene_name <- eqtl_list()$gene[input$mytable_rows_selected]
       cell_type_name <- eqtl_list()$cell_type[input$mytable_rows_selected]
       snp_name <- eqtl_list()$SNP[input$mytable_rows_selected]
       p_adj <- eqtl_list()$adj_p[input$mytable_rows_selected]
       if(!is.na(snp_name)){
        plot_eqtl(cell_type_name,gene_name,snp_name,p_adj)
       }
   })

   #Do the plotting for specific cis-eQTLs
   output$plot_eqtl_spe <- renderPlot({
       req(input$mytable2_rows_selected)
       gene_name <- eqtl_spe_list()$gene[input$mytable2_rows_selected]
       cell_type_name <- eqtl_spe_list()$cell_type[input$mytable2_rows_selected]
       snp_name <- eqtl_spe_list()$SNP[input$mytable2_rows_selected]
       p_adj <- eqtl_spe_list()$p_adj[input$mytable2_rows_selected]
       if(!is.na(snp_name)){
           plot_eqtl_specific(cell_type_name,gene_name,snp_name,p_adj)
       }
  })

   #Get coloc_df
   coloc_df <- reactive({
     req(input$num3)
     req(input$select_trait)
     threshold <- input$num3
     trait <- input$select_trait

     coloc_all <- h5read(hdF5_file_path, paste0("coloc/",trait))

     #Get best coloc if same gene in same cell type was tested in two different loci
     coloc_all <- coloc_all %>% group_by(ensembl,tissue) %>%
       slice_max(n=1,order_by=PP.H4.abf,with_ties=FALSE) %>%
       ungroup() %>%
       filter(!is.na(symbol),symbol!='NA')

   })

   #Do the plotting for coloc
   output$plot_coloc <- renderPlot(height = function(){
     n_unique_genes <- filter(coloc_df(),PP.H4.abf>=input$num3) %>% pull(symbol) %>% unique() %>% length()
     return(200+(n_unique_genes*15))
   },{
       req(input$num3)
       req(input$select_trait)
       threshold <- input$num3
       trait <- input$select_trait
       coloc_heatmap(coloc_df(),threshold)
   })

})
