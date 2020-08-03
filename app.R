library(shiny)
library(tidyverse)
library(cowplot)

# things to figure out
# 1. add SNP cline plts > added, currently show subset of clines with high (>0.8 or 0.9 Pearson correlation), hard to show all the correlations without dividing into subsets 
# 2. add gene function plts > changed to tables for now, but would like to link back to what genes/SNPs go with the pathways
# 3. add admixture K plts?? >> maybe not as useful, compared to static image? OR, just add the static image(s) with the PCA plots? but then a lot going on in one tab...

# ---------------------------------------
# UI function
# ---------------------------------------
ui <- fluidPage(
  title='Hbird data playground',
  titlePanel('Hummingbird pop gen data playground'),
  
  tabsetPanel(
    # ---------------------------------------
    # home page - info about this shiny app
    # ---------------------------------------
    tabPanel(
      title='About',
      h2('Further exploration of data from', em('in revision'), 'manuscript:'),
      h3(em('Pervasive genomic signatures of local adaptation to altitude across highland specialist Andean hummingbird populations')),
      h4('Marisa C.W. Lim, Ke Bi, Christopher C. Witt, Catherine H. Graham, Liliana M. Davalos')
    ),
    
    # ---------------------------------------
    # PCA plot tab
    # ---------------------------------------
    tabPanel(
      title='PCA plots',
      sidebarLayout(
        position='right',
        sidebarPanel('Select data to plot',
                     selectInput(inputId = 'xvarcv',
                                 label=h4('x-axis variable'),
                                 choices = c('PC1', 'PC2', 'PC3', 'PC4', 'PC5')),
                     selectInput(inputId ='yvarcv',
                                 label=h4('y-axis variable'),
                                 choices = c('PC2', 'PC3', 'PC4', 'PC5'))),
        
        mainPanel(
          h3('PCA plot for', em('Coeligena violifer'), ':'),
          plotOutput('PCAcviol')
        )
    ),
    sidebarLayout(
      position='right',
      sidebarPanel('Select data to plot',
                   selectInput(inputId = 'xvar',
                               label=h4('x-axis variable'),
                               choices = c('PC1', 'PC2', 'PC3', 'PC4', 'PC5')),
                   selectInput(inputId ='yvar',
                               label=h4('y-axis variable'),
                               choices = c('PC2', 'PC3', 'PC4', 'PC5'))),
      
      mainPanel(
        h3('PCA plot for', em('Colibri coruscans'), ':'),
        plotOutput('PCAccoru')
      )
    )),
    
    # ---------------------------------------
    # SNP cline plot tab
    # ---------------------------------------
    tabPanel(
      title='SNP clines',
      p('... maybe the boxplots if i can figure out getting them on same tab'),
      h2(em('Coeligena violifer'), 'plots:'),
      sidebarLayout(
        position='right',
        sidebarPanel('Select Pearson correlation range to plot',
                     sliderInput('cvslider1',
                                 label=h3('Select range'),
                                 min=0.9, max=1, value=0.9)),
        mainPanel(
          h3('Upward trending SNP clines'),
          plotOutput('SNPcviol_up')
        )
      ),
      sidebarLayout(
        position='right',
        sidebarPanel('Select Pearson correlation range to plot',
                     sliderInput('cvslider2',
                                 label=h3('Select range'),
                                 min=-1, max=-0.9, value=-0.9)),
        mainPanel(
          h3('Downward trending SNP clines'),
          plotOutput('SNPcviol_down')
        )
      ),
      h2(em('Colibri coruscans'), 'plots:'),
      sidebarLayout(
        position='right',
        sidebarPanel('Select Pearson correlation range to plot',
                     sliderInput('ccslider1',
                                 label=h3('Select range'),
                                 min=0.9, max=1, value=0.9)),
        mainPanel(
          h3('Upward trending SNP clines'),
          plotOutput('SNPccoru_up')
        )
      ),
      sidebarLayout(
        position='right',
        sidebarPanel('Select Pearson correlation range to plot',
                     sliderInput('ccslider2',
                                 label=h3('Select range'),
                                 min=-1, max=-0.8, value=-0.9)),
        mainPanel(
          h3('Downward trending SNP clines'),
          plotOutput('SNPccoru_down')
        )
      )
    ),
    
    # ---------------------------------------
    # candidate gene pathways tab
    # ---------------------------------------
    tabPanel(
      title='gene function',
      sidebarLayout(
        position='right',
        sidebarPanel('Select number of genes with x function',
                     selectInput('genecounter1',
                                 label=h3('Select number of genes'),
                                 choices=c(seq(1,8)))),
                     
        mainPanel(
          h3('Pathway functions of candidate genes with clinal SNPs for', em('Coeligena violifer'), ':'),
          tableOutput('candgenepaths1')
        )
      ),
      sidebarLayout(
        position='right',
        sidebarPanel('Select number of genes with x function',
                     selectInput('genecounter2',
                                 label=h3('Select number of genes'),
                                 choices=c(seq(1,10)))),
        
        mainPanel(
          h3('Pathway functions of candidate genes with clinal SNPs for', em('Colibri coruscans'), ':'),
          tableOutput('candgenepaths2')
        )
      )
      
    )
  )
)

# ---------------------------------------
# server function
# ---------------------------------------
server <- function(input, output){
  output$PCAcviol <- renderPlot({
    # read in data
    cviol_pca_dat <- read.csv('./PCA_plotdata_cviol.csv', header=T)
    cviol_eigen <- read.csv('./PCA_eigenvals_cviol.csv', header=T)
    cviol_eigen_t <- as.data.frame(t(cviol_eigen))
    names(cviol_eigen_t) <- names(cviol_pca_dat[1:59])
    
    #set up sample colors by geo subpop
    cviolcols <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#c51b7d")
    
    # plot
    cviol_pca_dat %>%
      ggplot() +
      geom_point(aes_string(x = input$xvarcv, y = input$yvarcv, fill = 'GeoPop'), size=6, alpha=0.7, pch=21) +
      scale_fill_manual(values=cviolcols) +
      xlab(paste(input$xvarcv, ' (', signif(cviol_eigen_t[[input$xvarcv]], digits=3)*100, '%)', sep='')) +
      ylab(paste(input$yvarcv, ' (', signif(cviol_eigen_t[[input$yvarcv]], digits=3)*100, '%)', sep='')) +
      theme_set(theme_cowplot())
  })
  
  output$PCAccoru <- renderPlot({
    # read in data
    ccoru_pca_dat <- read.csv('./PCA_plotdata_ccoru.csv', header=T)
    ccoru_eigen <- read.csv('./PCA_eigenvals_ccoru.csv', header=T)
    ccoru_eigen_t <- as.data.frame(t(ccoru_eigen))
    names(ccoru_eigen_t) <- names(ccoru_pca_dat[1:97])
    
    #set up sample colors by geo subpop
    ccorucols <- c("#a6cee3", "#fdbf6f", "#ff7f00", "#cab2d6", "#33a02c", "#c51b7d", "#6a3d9a")
    
    # plot
    ccoru_pca_dat %>%
      ggplot(aes_string(x = input$xvar, y = input$yvar, fill = 'GeoPop')) + 
      geom_point(size=6, alpha=0.7, pch=21) +
      scale_fill_manual(values=ccorucols) + 
      xlab(paste(input$xvar, ' (', signif(ccoru_eigen_t[[input$xvar]], digits=3)*100, '%)', sep='')) +
      ylab(paste(input$yvar, ' (', signif(ccoru_eigen_t[[input$yvar]], digits=3)*100, '%)', sep='')) +
      theme_set(theme_cowplot())
  })
  
  # input data with Pearson corrs and allele frequencies per SNP by Gene
  snp_df_cviol <- read.csv('./SNPclinedata_cviol.csv', header=T)
  snp_df_ccoru <- read.csv('./SNPclinedata_ccoru.csv', header=T)
  output$SNPcviol_up <- renderPlot({
    # positive correlations of minor allele frequency and elevation bin
    dattoplot_sub_trendup <- snp_df_cviol[snp_df_cviol$Pearson_coef_littlea >= input$cvslider1, ]
    ggplot(data=dattoplot_sub_trendup, aes(x=thebin, y=littlea_freq, group=SNP, col=Gene)) +
      # facet_wrap(~Gene) +
      geom_line(alpha=0.7) +
      geom_point() +
      theme(axis.text.x=element_text(angle=45, hjust=1)) +
      ylab('Minor allele frequency') + xlab('Elevation bin (meters)')
    })
  output$SNPcviol_down <- renderPlot({
    # negative correlations of minor allele frequency and elevation bin
    dattoplot_sub_trenddown <- snp_df_cviol[snp_df_cviol$Pearson_coef_littlea <= input$cvslider2, ]
    ggplot(data=dattoplot_sub_trenddown, aes(x=thebin, y=littlea_freq, group=SNP, col=Gene)) +
      # facet_wrap(~Gene) +
      geom_line(alpha=0.7) +
      geom_point() +
      theme(axis.text.x=element_text(angle=45, hjust=1)) +
      ylab('Minor allele frequency') + xlab('Elevation bin (meters)')
  })
  output$SNPccoru_up <- renderPlot({
    # positive correlations of minor allele frequency and elevation bin
    dattoplot_sub_trendup <- snp_df_ccoru[snp_df_ccoru$Pearson_coef_littlea >= input$ccslider1, ]
    ggplot(data=dattoplot_sub_trendup, aes(x=thebin, y=littlea_freq, group=SNP, col=Gene)) +
      # facet_wrap(~Gene) +
      geom_line(alpha=0.7) +
      geom_point() +
      theme(axis.text.x=element_text(angle=45, hjust=1)) +
      ylab('Minor allele frequency') + xlab('Elevation bin (meters)')
  })
  output$SNPccoru_down <- renderPlot({
    # negative correlations of minor allele frequency and elevation bin
    dattoplot_sub_trenddown <- snp_df_ccoru[snp_df_ccoru$Pearson_coef_littlea <= input$ccslider2, ]
    ggplot(data=dattoplot_sub_trenddown, aes(x=thebin, y=littlea_freq, group=SNP, col=Gene)) +
      # facet_wrap(~Gene) +
      geom_line(alpha=0.7) +
      geom_point() +
      theme(axis.text.x=element_text(angle=45, hjust=1)) +
      ylab('Minor allele frequency') + xlab('Elevation bin (meters)')
  })
  
  output$candgenepaths1 <- renderTable({
    # Previously ID'd candidate genes
    cviol_Wilcoxcand_pathway <- read.table('cviol_wilcox_pathway.txt', sep='\t', header=F)
    cvioltab <- cviol_Wilcoxcand_pathway[cviol_Wilcoxcand_pathway$V3 == input$genecounter1, 2]
  })
  
  output$candgenepaths2 <- renderTable({
    # Previously ID'd candidate genes
    ccoru_Wilcoxcand_pathway <- read.table('ccoru_wilcox_pathway.txt', sep='\t', header=F)
    ccorutab <- ccoru_Wilcoxcand_pathway[ccoru_Wilcoxcand_pathway$V3 == input$genecounter2, 2]
  })
  

  
}

shinyApp(ui = ui, server = server)
