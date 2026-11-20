#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(ggplot2)
library(umap)
library(NMF)
library(mclust)

gene_expr <- as.matrix(TBnanostring[, -1])
gene_scaled <- scale(gene_expr)
col_mins <- apply(gene_expr, 2, min)
expr_shift <- sweep(gene_expr, 2, col_mins, FUN = "-")
TBStatus <- factor(TBnanostring$TB_Status)

# Define UI
ui <- fluidPage(
  
  # Put a titlePanel here
  titlePanel("k-means clustering of TBNanostring data"),
  
  sidebarLayout(
    # Sidebar. Put your inputs inside the sidebarPanel
    sidebarPanel(
      selectInput("method", "Dimension reduction method", 
                  choices = c("PCA", "UMAP", "NMF")),
      numericInput("clusters", "Cluster count (k)", 2,
                   min = 1, max = 9)
    ),
    
    # Main panel. put your output plot here
    mainPanel(
      plotOutput("plot1"),
      h4("Clusters vs. TB Status"),
      tableOutput("crosstab"),
      h4("Adjusted Rand index"),
      textOutput("ARI")
    )
  )
)

server <- function(input, output, session) {
  
  set.seed(123)
  embedding <- reactive({
    emb <- switch(input$method,
                  PCA = prcomp(gene_scaled)$x[, 1:2],
                  UMAP = umap(gene_scaled)$layout,
                  NMF = basis(nmf(expr_shift, rank = 2, nrun = 10))
    )
    emb <- as.matrix(emb)
    colnames(emb) <- c("Dim1", "Dim2")
    emb
  }) %>% bindCache(input$method)
  
  axis_labels <- reactive({
    switch(input$method,
           PCA = c("PC1", "PC2"),
           UMAP = c("UMAP1", "UMAP2"),
           NMF = c("NMF1", "NMF2"))
  })
  
  clusters <- reactive({
    set.seed(123)
    kmeans(embedding(), input$clusters, nstart = 25)
  })
  
  output$plot1 <- renderPlot({
    emb <- embedding()
    km <- clusters()
    
    plot_df <- data.frame(
      Dim1 = emb[, 1],
      Dim2 = emb[, 2],
      Cluster = factor(km$cluster),
      TB_Status = TBStatus
    )
    centers_df <- as.data.frame(km$centers)
    colnames(centers_df) <- c("Dim1", "Dim2")
    
    ggplot() +
      geom_point(data = plot_df,
                 aes(x = Dim1, y = Dim2, color = Cluster, shape = TB_Status),
                 size = 3, alpha = 0.8) +
      geom_point(data = centers_df, aes(x = Dim1, y = Dim2),
                 shape = 4, size = 6, stroke = 2, color = "black") +
      scale_color_brewer(palette = "Set1") +
      labs(x = axis_labels()[1], y = axis_labels()[2],
           color = "k-means cluster", shape = "TB Status",
           title = paste(input$method, "with", input$clusters, "clusters")) +
      theme_minimal(base_size = 14)
  })
  
  output$crosstab <- renderTable({
    tab <- table(Cluster = clusters()$cluster, TB_Status = TBStatus)
    as.data.frame.matrix(tab)
  }, rownames = TRUE)
  
  output$ari <- renderText({
    ari <- adjustedRandIndex(clusters()$cluster, TBStatus)
    paste0("ARI = ", round(ari, 3))
  })
}

shinyApp(ui = ui, server = server)