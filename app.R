suppressWarnings(suppressPackageStartupMessages(library(shiny)))
suppressWarnings(suppressPackageStartupMessages(library(readxl)))
suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))
suppressWarnings(suppressPackageStartupMessages(library(scales)))
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(DT)))
suppressWarnings(suppressPackageStartupMessages(library(patchwork)))
suppressWarnings(suppressPackageStartupMessages(library(ggpmisc)))
suppressWarnings(suppressPackageStartupMessages(library(lemon)))


options(shiny.maxRequestSize = 50 * 1024^2)
# Path to data
get_plot_dat <- function(gene, MANE, BIC, parameters, lasso.coefs, ccle.tpms, boot_data) {
  
  if (str_detect(gene, "ENSG")) {
    genesym <- MANE$GeneSymbol[MANE$Ensembl_Gene == gene]
  } else {
    genesym <- gene
  }
  
  best_model <- BIC[genesym, "BestModel"]
  
  file <- list.files(path = "csvfiles/") %>% as.data.frame() %>%
    filter(. == paste0(genesym, ".csv")) %>%
    pull()
  
  ccle_dat <- fread(paste0("csvfiles/", file)) %>%
    select(one_of("CellLine", "TPM", "Protein"))
  
  params <- parameters %>%
    filter(GeneSymbol == genesym)
  
  if (best_model == "Median") {
    
    ybest <- rep(boot_data$Median[boot_data$GeneSymbol == genesym], length.out = length(ccle_dat$TPM))
    yupper <- rep(boot_data$Upper_boot[boot_data$GeneSymbol == genesym], length.out = length(ccle_dat$TPM))
    ylower <- rep(boot_data$Lower_boot[boot_data$GeneSymbol == genesym], length.out = length(ccle_dat$TPM))
    
  } else {
    
    ybest <- params$A*(params$B*ccle_dat$TPM/(params$C+ccle_dat$TPM)+ccle_dat$TPM)
    ybest[ybest < 0] <- 0
    
    yupper <- data.frame(
      A = params$A_upper*(params$B*ccle_dat$TPM / (params$C + ccle_dat$TPM) + ccle_dat$TPM),
      B = params$A*(params$B_upper*ccle_dat$TPM / (params$C + ccle_dat$TPM) + ccle_dat$TPM),
      C = params$A*(params$B*ccle_dat$TPM / (params$C_upper + ccle_dat$TPM) + ccle_dat$TPM)
    ) %>% 
      mutate(yupper = {apply(., 1, max)}) %>%
      mutate(yupper = replace(yupper, which(yupper < 0), 0)) %>%
      pull()
    
    ylower <- data.frame(
      A = params$A_lower*(params$B*ccle_dat$TPM / (params$C + ccle_dat$TPM) + ccle_dat$TPM),
      B = params$A*(params$B_lower*ccle_dat$TPM / (params$C + ccle_dat$TPM) + ccle_dat$TPM),
      C = params$A*(params$B*ccle_dat$TPM / (params$C_lower + ccle_dat$TPM) + ccle_dat$TPM)
    ) %>%
      mutate(ylower = {apply(., 1, min)}) %>%
      mutate(ylower = replace(ylower, which(ylower < 0), 0)) %>%
      pull()
    
    if (best_model == "LASSO") {
      
      coefs <- lasso.coefs[[genesym]] %>% as.data.frame()
      lasso.tpm <- ccle.tpms %>%
        column_to_rownames(var = "GeneSymbol") %>%
        select(all_of(ccle_dat$CellLine)) %>%
        t() %>% as.data.frame() %>%
        select(one_of(colnames(coefs)))
      
      coefs <- coefs %>%
        select(colnames(lasso.tpm))
      
      mod <- data.matrix(lasso.tpm) %*% t(coefs) %>%
        as.data.frame() %>%
        rename(LASSO_mod = "V1") %>%
        arrange({factor(rownames(.), levels = ccle_dat$CellLine)})
      
      ybest <- ybest - mod$LASSO_mod
      yupper <- yupper - mod$LASSO_mod
      ylower <- ylower - mod$LASSO_mod
      
    }
    
  }
  
  Regression <- tibble(
    TPM = ccle_dat$TPM,
    yfit = ybest,
    yupper = yupper,
    ylower = ylower
  ) %>%
    mutate(ymax = pmax(yupper, ylower),
           ymin = pmin(yupper, ylower),
           fill = yupper >= ylower) %>%
    select(one_of("TPM", "yfit", "ymax", "ymin", "fill"))
  Regression[Regression < 0] <- 0
  tbl <- params %>% 
    column_to_rownames(var = "GeneSymbol") %>% 
    select("A", "B", "C") %>% 
    t() %>% as.data.frame() %>%
    `colnames<-`(value = "Best fit") %>%
    rownames_to_column(var = "Parameter")
  tbl$`Best fit` <- round(tbl$`Best fit`, digits = 2)
  plot_dat <- list(ccle_dat, Regression, best_model, tbl, genesym)
  names(plot_dat) <- c("ccle_dat", "Regression", "BestModel", "param_table", "gene")
  return(plot_dat)
}

make_predictions <- function(dat, BIC, parameters, lasso.coefs, meds, lasso_tpms) {
  # Select the genes for which predictions are able to be made.
  dat <- dat %>%
    filter(GeneSymbol %in% parameters$GeneSymbol)
  
  preds <- apply(dat, 1, function(x) {
    best_model <- BIC[x["GeneSymbol"], "BestModel"]
    df <- as.data.frame(x) %>%
      filter({rownames(.)} != "GeneSymbol") %>%
      `colnames<-`(value = "TPM") %>%
      mutate(Pred = NA,
             Upper = NA,
             Lower = NA)
    df$TPM <- as.numeric(df$TPM)
    
    if (best_model == "Median") {
      df$Pred <- rep(meds[meds$GeneSymbol == x["GeneSymbol"], "Median"], length.out = length(df$Pred))
      df$Upper <- rep(meds[meds$GeneSymbol ==x["GeneSymbol"], "Upper_boot"], length.out = length(df$Upper))
      df$Lower <- rep(meds[meds$GeneSymbol ==x["GeneSymbol"], "Lower_boot"], length.out = length(df$Lower))
      df <- df %>% rownames_to_column(var = "Sample")
    } else {
      A <- parameters$A[parameters$GeneSymbol == x["GeneSymbol"]]
      B <- parameters$B[parameters$GeneSymbol == x["GeneSymbol"]]
      C <- parameters$C[parameters$GeneSymbol == x["GeneSymbol"]]
      df$Pred <- A*(B*df$TPM/(C+df$TPM)+df$TPM)
      
      A_upper <- parameters$A_upper[parameters$GeneSymbol == x["GeneSymbol"]]
      A_lower <- parameters$A_lower[parameters$GeneSymbol == x["GeneSymbol"]]
      B_upper <- parameters$B_upper[parameters$GeneSymbol == x["GeneSymbol"]]
      B_lower <- parameters$B_lower[parameters$GeneSymbol == x["GeneSymbol"]]
      C_upper <- parameters$C_upper[parameters$GeneSymbol == x["GeneSymbol"]]
      C_lower <- parameters$C_lower[parameters$GeneSymbol == x["GeneSymbol"]]
      
      yupper <- data.frame(
        A = A_upper*(B*df$TPM / (C + df$TPM) + df$TPM),
        B = A*(B_upper*df$TPM / (C + df$TPM) + df$TPM),
        C = A*(B*df$TPM / (C_upper + df$TPM) + df$TPM)
      ) %>%
        mutate(yupper = {apply(., 1, max)}) %>%
        mutate(yupper = replace(yupper, which(yupper < 0), 0))
      df$Upper <- yupper$yupper
      
      ylower <- data.frame(
        A = A_lower*(B*df$TPM / (C + df$TPM) + df$TPM),
        B = A*(B_lower*df$TPM / (C + df$TPM) + df$TPM),
        C = A*(B*df$TPM / (C_lower + df$TPM) + df$TPM)
      ) %>% 
        mutate(ylower = {apply(., 1, min)}) %>%
        mutate(ylower = replace(ylower, which(ylower < 0), 0))
      df$Lower <- ylower$ylower
      
      df <- df %>% rownames_to_column(var = "Sample")
      
      if (best_model == "LASSO" & !is.null(lasso_tpms)) {
        sig.feats <- lasso.coefs[[x["GeneSymbol"]]] %>% as.data.frame()
        sig.tpm <- lasso_tpms %>%
          filter(GeneSymbol %in% colnames(sig.feats)) %>%
          column_to_rownames(var = "GeneSymbol")
        sig.feats <- sig.feats %>%
          select(rownames(sig.tpm))
        
        mod <- t(as.matrix(sig.tpm)) %*% t(as.matrix(sig.feats)) %>% as.data.frame() %>%
          `rownames<-`(value = df$Sample) %>%
          `colnames<-`(value = "mod") %>%
          rownames_to_column(var = "Sample")
        df <- df %>%
          left_join(mod, by = "Sample") %>%
          mutate(Pred = Pred - mod,
                 Upper = Upper - mod,
                 Lower = Lower - mod)
      }
    }
    df[,c("Pred", "Upper", "Lower")][df[,c("Pred", "Upper", "Lower")] < 0] <- 0
    return(df)
    
  }) 
  names(preds) <- dat$GeneSymbol
  
  return(preds)
  
}

table_data <- function(dat, index) {
  sample_names <- sapply(dat, "[[", "Sample")
  sample_names <- sample_names[,1] %>% as.character()
  dat_tbl <- round(sapply(dat, "[[", index)) %>%
    `rownames<-`(value = sample_names) %>%
    t() %>% as.data.frame() %>%
    rownames_to_column(var = "Gene Symbol")
  return(dat_tbl)
}

# Non-user input data needed for calculations.
MANE <- fread("IDs.csv", data.table = F) %>%
  select("GeneSymbol", "Ensembl_Gene") %>%
  mutate(Ensembl_Gene = sub("\\..*", "", Ensembl_Gene))
BIC <- fread("BIC.csv", data.table = F) %>%
  column_to_rownames(var = "GeneSymbol") %>%
  select("BestModel")
parameters <- fread("RegressionParameters.csv", data.table = F)
lasso.coefs <- readRDS("LASSO_coefficients.RData")
ccle.tpms <- fread("CCLE_expression.csv", data.table = F)
medians <- fread("StDevs_Meds.csv", data.table = F)

# Define UI ----
# Define UI ----
ui <- fluidPage(
  titlePanel("Pinferna"),
  sidebarLayout(
    sidebarPanel(
      helpText("To use Pinferna: 1) Upload a table of genes with mRNA normalized as transcripts per million (TPM). 2) Optional: Upload a table of the full transcriptomics to allow for HL+LASSO predictions. You may also select the checkbox to use the same data as step 1. 3) Click Go!"),
      helpText("Upload genes of interest normalized as TPM in either .csv, .tsv, or .xls(x) file format. Genes should be in rows; samples should be in columns. Gene identifiers should be in a column containing gene symbols matching those in MANE Select (1). The column name should be “GeneSymbol.” A conversion table and data entry template are provided with the publication associated with Pinferna (2). The maximum file size is 50 MB."),
      fileInput("upload", "Upload a file", buttonLabel = "Upload...", multiple = F,
                accept = c(".csv", ".tsv", ".xls", ".xlsx")),
      helpText("To incorporate LASSO, upload a second file of the full transcriptomic output from the RNA sequencing experiment normalized as TPM. The file format should be as above, and all (gene and sample) identifiers should match between the two tables. Or, check the box to use the data already uploaded above."),
      checkboxInput("lasso_data", "Use the table I uploaded above.", value = F),
      fileInput("upload_lasso", "Upload a file", buttonLabel = "Upload for LASSO...", multiple = F,
                accept = c(".csv", ".tsv", ".xls", ".xlsx")),
      actionButton("go", "Go!"),
      selectizeInput("Gene", "Select gene:", choices = NULL, multiple = F),
      checkboxGroupInput("plotOptions", "Plot Adjustments:",
                         c("Show CCLE training data" = "showCCLE",
                           "Show input data" = "showInput",
                           "Show regression" = "showRegression",
                           "Show 95% confidence interval" = "showCI"),
                         selected = c("showCCLE",
                                      "showInput",
                                      "showRegression",
                                      "showCI")),
      downloadButton("currentPlot", "Export Current Plot"),
      # downloadButton("allPlots", "Export All Plots"),
      helpText("References: (1) Morales et al., Nature 604:310-315 (2022) \n (2) Sweatt et al., citation pending.")
    ),
    
    mainPanel(
      plotOutput("regression"),
      tabsetPanel(id = "Predictions",
                  tabPanel("Predicted protein abundance", DTOutput("best_table")),
                  tabPanel("Upper 95% CI", DTOutput("upper_table")),
                  tabPanel("Lower 95% CI", DTOutput("lower_table"))
      )
    )
  )
)


# Define server logic ----
server <- function(input, output, session) {
  
  user_data <- reactive({
    req(input$upload)
    
    ext <- tools::file_ext(input$upload$name)
    switch(ext,
           csv = fread(input$upload$datapath),
           tsv = fread(input$upload$datapath),
           xls = read_excel(input$upload$datapath),
           xlsx = read_excel(input$upload$datapath),
           validate("Error with file extension: Please upload a .csv, .tsv, or .xls(x) file.")
    )
  })
  
  all_tpms <- reactive({
    if (input$lasso_data) {
      return(user_data())
    } else if (is.null(input$upload_lasso)) {
      return(NULL)
    } else {
      ext <- tools::file_ext(input$upload_lasso$name)
      switch(ext,
             csv = fread(input$upload_lasso$datapath),
             tsv = fread(input$upload_lasso$datapath),
             xls = read_excel(input$upload_lasso$datapath),
             xlsx = read_excel(input$upload_lasso$datapath),
             validate("Error with file extension: Please upload a .csv, .tsv, or .xls(x) file.")
      )
    }
  })
  
  
  observeEvent(input$go, {
    choices <- as.character(unique(user_data()$GeneSymbol[user_data()$GeneSymbol %in% rownames(BIC)]))
    updateSelectizeInput(session, inputId = "Gene",
                         choices = choices, server = T)
  })
  
  predict_dat <- eventReactive(input$go, {
    progress <- Progress$new(session, min = 0, max = 1)
    on.exit(progress$close())
    progress$set(message = "Predictions in progress...",
                 detail = "For more than 100 genes, this may take a while.")
    return(make_predictions(dat = user_data(), BIC = BIC, parameters = parameters, lasso.coefs = lasso.coefs, lasso_tpms = all_tpms(), meds = medians))
  })
  
  plot_dat <- reactive({
    req(input$Gene)
    get_plot_dat(gene = input$Gene, MANE = MANE, BIC = BIC, ccle.tpms = ccle.tpms, lasso.coefs = lasso.coefs, parameters = parameters, boot_data = medians)
  })
  
  regression_data <- reactive({
    ccle_dat <- plot_dat()[["ccle_dat"]]
    Regression <- plot_dat()[["Regression"]]
    BestModel <- plot_dat()[["BestModel"]]
    if (BestModel == "LASSO") {BestModel <- "HL+LASSO"}
    param_tbl <- plot_dat()[["param_table"]]
    genesym <- plot_dat()[["gene"]]
    xmin <- 10^(floor(log10(min(ccle_dat$TPM+1))))
    xmax <- 10^(ceiling(log10(max(ccle_dat$TPM+1))))
    ymin <- 10^(floor(log10(min(Regression$ymin+1))))
    ymax <- 10^(ceiling(log10(max(Regression$ymax+1))))
    g <- ggplot() +
      scale_x_continuous(limits = c(xmin, xmax),
                         breaks = c(1 %o% 10^(0:log10(xmax))),
                         trans = "log10",
                         labels = comma) +
      scale_y_continuous(limits = c(ymin, ymax),
                         breaks = c(1 %o% 10^(0:log10(ymax))),
                         trans = "log10",
                         labels = comma) +
      coord_capped_cart(bottom = "right", left = "top", gap = 0) +
      xlab("TPM + 1") + ylab("Protein copies + 1") + ggtitle(genesym) +
      theme(aspect.ratio = 1, plot.title = element_text(h = 0.5),
            panel.background = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(color = "black", linewidth = 0.5),
            axis.ticks = element_line(color = "black", linewidth = 0.5),
            axis.ticks.length = unit(0.1, "inch"),
            axis.text = element_text(family = "sans", color = "black", size = 14),
            axis.title = element_text(family = "sans", color = "black", size = 16),
            title = element_text(family = "sans", color = "black", size = 18))

    if ("showCCLE" %in% input$plotOptions) {
      g <- g +
        geom_point(data = ccle_dat, aes(x = TPM+1, y = Protein+1), color = "black", stroke = NA, alpha = 0.3, size = 2)
    }
    if ("showRegression" %in% input$plotOptions) {
      if (BestModel == "Median") {
        g <- g +
          geom_smooth(data = Regression, aes(x = TPM+1, y = yfit+1), color = "#28c742", linewidth = 0.75, se =F) +
          geom_text(aes(x = xmax, y = ymin, label = BestModel), color = "#28c742", hjust = "inward", vjust = "inward", size = 6)
      } else if (BestModel == "HL") {
        g <- g +
          geom_smooth(data = Regression, aes(x = TPM+1, y = yfit+1), color = "#aa3481", linewidth = 0.75, se =F) +
          geom_text(aes(x = xmax, y = ymin, label = BestModel), color = "#aa3481", hjust = "inward", vjust = "inward", size = 6)
      } else {
        g <- g +
          geom_smooth(data = Regression, aes(x = TPM+1, y = yfit+1), color = "#e66100", linewidth = 0.75, se =F) +
          geom_text(aes(x = xmax, y = ymin, label = BestModel), color = "#e66100", hjust = "inward", vjust = "inward", size = 6)
      }

    }
    if ("showCI" %in% input$plotOptions) {
      Regression$ymin_rib <- 10^predict(loess(log10(Regression$ymin+1)~log10(Regression$TPM+1)))-1
      Regression$ymax_rib <- 10^predict(loess(log10(Regression$ymax+1)~log10(Regression$TPM+1)))-1
      if (BestModel == "Median") {
        g <- g +
          geom_ribbon(data = Regression,
                      aes(x = TPM+1, ymin = ymin_rib+1, ymax = ymax_rib+1),
                      fill = "#28c742", alpha = 0.2, inherit.aes = F)
      } else if (BestModel == "HL") {
        g <- g +
          geom_ribbon(data = Regression,
                      aes(x = TPM+1, ymin = ymin_rib+1, ymax = ymax_rib+1),
                      fill = "#aa3481", alpha = 0.2, inherit.aes = F)
      } else {
        g <- g +
          geom_ribbon(data = Regression,
                      aes(x = TPM+1, ymin = ymin_rib+1, ymax = ymax_rib+1),
                      fill = "#e66100", alpha = 0.2, inherit.aes = F)
      }
    }
    if ("showInput" %in% input$plotOptions) {
      point_dat <- predict_dat()[[genesym]] %>% as.data.frame()
      g <- g +
        geom_point(data = point_dat, aes(x = TPM + 1, y = Pred + 1),
                   shape = 23, fill = "#e7298a", color = "black", size = 4)

      if (min(point_dat$TPM) < min(ccle_dat$TPM)) {
        g <- g +
          geom_text(aes(x = xmin, y = ymax, label = "Warning: Some input data extrapolated."),
                    color = "red", hjust = "inward", vjust = "inward", size = 4)
      }
    }
    if (BestModel %in% c("HL", "HL+LASSO")) {
      ggtbl <- ggplot() +
        theme_void() +
        annotate(geom = "table", x = 1, y = 1, label = list(param_tbl), fontface = "sans", size = 6)
    } else {
      ggtbl <- ggplot() + theme_void()
    }
    regression_data <- list(g, ggtbl)
    names(regression_data) <- c("g", "ggtbl")
    return(regression_data)
  })
  
  output$regression <- renderPlot({
    g <- regression_data()[["g"]]
    ggtbl <- regression_data()[["ggtbl"]]
    g + ggtbl
  })
  
  output$best_table <- renderDataTable(server = F,
    {table_data(predict_dat(), "Pred")},
    extensions = "Buttons",
    options = list(
      "dom" = 'T<"clear">lBfrtip',
      buttons = list('copy', list(
        extend = 'collection',
        buttons = list(
          list(extend = 'csv', filename = "prediction"),
          list(extend = 'excel', filename = "prediction")),
        text = "Download"
      )),
      scrollY = TRUE
    )
  )
  
  output$upper_table <- renderDataTable(server = F,
    {table_data(predict_dat(), "Upper")},
    extensions = "Buttons",
    options = list(
      "dom" = 'T<"clear">lBfrtip',
      buttons = list('copy', list(
        extend = 'collection',
        buttons = list(
          list(extend = 'csv', filename = "Upper95CI"),
          list(extend = 'excel', filename = "Upper95CI")),
        text = "Download"
      )),
      scrollY = TRUE
    )
  )
  output$lower_table <- renderDataTable(server = F,
    {table_data(predict_dat(), "Lower")},
    extensions = "Buttons",
    options = list(
      "dom" = 'T<"clear">lBrftip',
      buttons = list('copy', list(
        extend = 'collection',
        buttons = list(
          list(extend = 'csv', filename = "Lower95CI"),
          list(extend = 'excel', filename = "Lower95CI")),
        text = "Download"
      )),
      scrollY = TRUE
    )
  )
  
  output$currentPlot <- downloadHandler(
    filename = function() { paste0(input$Gene, ".pdf") },
    content = function(file) {
      ggsave(file, plot = regression_data()[["g"]]+regression_data()[["ggtbl"]], device = "pdf", width = 8, height = 4, units = "in", dpi = 300)
    }
  )
  
}

# Run the app ----
shinyApp(ui = ui, server = server)
