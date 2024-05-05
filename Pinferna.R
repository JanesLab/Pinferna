## Code for Pinferna Shinyapp.
## Andrew Sweatt
## April 29, 2024

suppressWarnings(suppressPackageStartupMessages(library(shiny)))
suppressWarnings(suppressPackageStartupMessages(library(shinyBS)))
suppressWarnings(suppressPackageStartupMessages(library(readxl)))
suppressWarnings(suppressPackageStartupMessages(library(tidyverse)))
suppressWarnings(suppressPackageStartupMessages(library(scales)))
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
suppressWarnings(suppressPackageStartupMessages(library(DT)))
suppressWarnings(suppressPackageStartupMessages(library(patchwork)))
suppressWarnings(suppressPackageStartupMessages(library(ggpmisc)))
suppressWarnings(suppressPackageStartupMessages(library(lemon)))

options(shiny.maxRequestSize = 50 * 1024^2)

set_identifiers <- function(data) {
  id_cols <- colnames(data)[colnames(data) %in% c('GeneSymbol','Ensembl_Gene')]
  if (length(id_cols) == 2) {
    data <- data %>%
      select(setdiff({colnames(.)}, 'Ensembl_Gene')) %>%
      relocate(GeneSymbol)
  } else if (id_cols == 'Ensembl_Gene') {
    if (any(str_detect(data$Ensembl_Gene, '.'))) {
      data$Ensembl_Gene <- sub('\\..*', '', data$Ensembl_Gene)
    }
    data <- data %>%
      left_join(MANE, by = 'Ensembl_Gene') %>%
      select(setdiff({colnames(.)}, 'Ensembl_Gene')) %>%
      relocate(GeneSymbol)
  } else {
    data <- data %>% relocate(GeneSymbol)
  }
  return(data)
}

get_calibration<-function(gene, calibration){
  
  if (data_type$Absolute[data_type$Gene == gene] == 'yes') {
    ccle_dat <- swath.prot %>%
      filter(GeneSymbol == gene) %>%
      column_to_rownames('GeneSymbol') %>%
      t() %>% as.data.frame() %>%
      rename(Protein = as.character(gene)) %>%
      rownames_to_column('CellLine') %>%
      left_join(
        ccle.tpms %>%
          filter(GeneSymbol == gene) %>%
          column_to_rownames('GeneSymbol') %>%
          t() %>% as.data.frame() %>%
          rename(TPM = as.character(gene)) %>%
          rownames_to_column('CellLine'),
        by = 'CellLine'
      )
  } else {
    ccle_dat <- tmt.prot %>%
      filter(GeneSymbol == gene) %>%
      column_to_rownames('GeneSymbol') %>%
      t() %>% as.data.frame() %>%
      rename(Protein = as.character(gene)) %>%
      {2^.} %>%
      rownames_to_column('CellLine') %>%
      left_join(
        ccle.tpms %>%
          filter(GeneSymbol == gene) %>%
          column_to_rownames('GeneSymbol') %>%
          t() %>% as.data.frame() %>%
          rename(TPM = as.character(gene)) %>%
          rownames_to_column('CellLine'),
        by = 'CellLine'
      )
  }
  
  calibration <- calibration %>%filter(GeneSymbol == gene)
  
  if (length(calibration$GeneSymbol) == 0){
    scaling_factor = 1
    return(scaling_factor)
  } else {
    calibration_celllines <- ccle_dat %>% filter(CellLine %in% calibration$CellLine) %>% pull(CellLine)
    scaling_factor <-ccle_dat %>% 
      filter(CellLine %in% calibration$CellLine) %>% 
      left_join(calibration, by = 'CellLine', suffix = c(".rel", ".abs")) %>% 
      mutate(scaling = Protein.abs/Protein.rel)
    
    scaling_factor<-mean(scaling_factor$scaling)
    calibration_data <- list(scaling_factor, calibration_celllines)
    names(calibration_data) <- c("scaling_factor", "calibration_celllines")
    return(scaling_factor)
  }
}


get_plot_dat <- function(gene, MANE, BIC, parameters, lasso.coefs, ccle.tpms, boot_data) {
  
  best_model <- BIC[gene, "BestModel"]
  
  if (data_type$Absolute[data_type$Gene == gene] == 'yes') {
    ccle_dat <- swath.prot %>%
      filter(GeneSymbol == gene) %>%
      column_to_rownames('GeneSymbol') %>%
      t() %>% as.data.frame() %>%
      rename(Protein = as.character(gene)) %>%
      rownames_to_column('CellLine') %>%
      left_join(
        ccle.tpms %>%
          filter(GeneSymbol == gene) %>%
          column_to_rownames('GeneSymbol') %>%
          t() %>% as.data.frame() %>%
          rename(TPM = as.character(gene)) %>%
          rownames_to_column('CellLine'),
        by = 'CellLine'
      )
  } else {
    ccle_dat <- tmt.prot %>%
      filter(GeneSymbol == gene) %>%
      column_to_rownames('GeneSymbol') %>%
      t() %>% as.data.frame() %>%
      rename(Protein = as.character(gene)) %>%
      {2^.} %>%
      rownames_to_column('CellLine') %>%
      left_join(
        ccle.tpms %>%
          filter(GeneSymbol == gene) %>%
          column_to_rownames('GeneSymbol') %>%
          t() %>% as.data.frame() %>%
          rename(TPM = as.character(gene)) %>%
          rownames_to_column('CellLine'),
        by = 'CellLine'
      )
  }
  
  na_dat <- ccle_dat %>%
    filter(is.na(TPM) | is.na(Protein))
  
  ccle_dat <- ccle_dat %>%
    filter(!(is.na(TPM) | is.na(Protein)))
  
  params <- parameters %>%
    filter(GeneSymbol == gene)
  
  if (best_model == "Median") {
    
    ybest <- rep(boot_data$Median[boot_data$GeneSymbol == gene], length.out = length(ccle_dat$TPM))
    yupper <- rep(boot_data$Upper_boot[boot_data$GeneSymbol == gene], length.out = length(ccle_dat$TPM))
    ylower <- rep(boot_data$Lower_boot[boot_data$GeneSymbol == gene], length.out = length(ccle_dat$TPM))
    
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
      
      coefs <- lasso.coefs[[gene]] %>% as.data.frame()
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
      ybest[ybest < 0] <- 0
      yupper <- yupper - mod$LASSO_mod
      yupper[yupper < 0] <- 0
      ylower <- ylower - mod$LASSO_mod
      ylower[ylower < 0] <- 0
      
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
    select(one_of("TPM", "yfit", "ymax", "ymin", "fill")) %>%
    arrange(TPM)
  
  Regression[Regression < 0] <- 0
  tbl <- params %>% 
    column_to_rownames(var = "GeneSymbol") %>% 
    select("A", "B", "C") %>% 
    t() %>% as.data.frame() %>%
    `colnames<-`(value = "Best fit") %>%
    rownames_to_column(var = "Parameter")
  tbl$`Best fit` <- round(tbl$`Best fit`, digits = 2)
  if (data_type$Absolute[data_type$Gene == gene] == 'yes') {type = 'Absolute'} else { type = 'Relative'}
  plot_dat <- list(ccle_dat, Regression, best_model, tbl, gene, na_dat, type)
  names(plot_dat) <- c("ccle_dat", "Regression", "BestModel", "param_table", "gene", "na_dat", "type")
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
          filter(GeneSymbol %in% colnames(sig.feats)) 
        
        if (length(which(duplicated(sig.tpm$GeneSymbol))) > 0) {
          dups <- sig.tpm$GeneSymbol[which(duplicated(sig.tpm$GeneSymbol))]
          for (dup in dups) {
            tmp <- sig.tpm[sig.tpm$GeneSymbol == dup,] %>%
              select(-one_of("GeneSymbol"))
            max_dup <- apply(tmp, 2, function(x) return(max(x, na.rm = T))) %>% t() %>% as.data.frame()
            df_tmp <- data.frame(GeneSymbol = dup) %>% bind_cols(max_dup)
            sig.tpm <- sig.tpm %>%
              filter(!(GeneSymbol %in% dup)) %>%
              bind_rows(df_tmp)
          }
        }
        
        sig.tpm <- sig.tpm %>%
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
  if (is.null(dat)){
    dat_tbl <-NULL
    return(dat_tbl)
  } else{
    sample_names <- sapply(dat, "[[", "Sample")
    if (length(unique(sample_names)) > 1) {
      sample_names <- sample_names[,1] %>% as.character()
      dat_tbl <- round(sapply(dat, "[[", index)) %>%
        `rownames<-`(value = sample_names) %>%
        t() %>% as.data.frame() %>%
        rownames_to_column(var = "Gene Symbol")
    } else {
      sample_names <- sample_names[1] %>% as.character()
      dat_tbl <- round(sapply(dat, "[[", index)) %>%
        t() %>% as.data.frame() %>%
        `rownames<-`(value = sample_names) %>%
        t() %>% as.data.frame() %>%
        rownames_to_column(var = "Gene Symbol")
    }
    return(dat_tbl)
  }
}

# Non-user input data needed for calculations.
MANE <- fread("IDs.csv", data.table = F) %>%
  select("GeneSymbol", "Ensembl_Gene") %>%
  mutate(Ensembl_Gene = sub("\\..*", "", Ensembl_Gene))
BIC <- fread("BIC_ALL.csv", data.table = F) %>%
  column_to_rownames(var = "GeneSymbol") %>%
  select("BestModel")
parameters <- fread("RegressionParameters_ALL.csv", data.table = F)
lasso.coefs <- readRDS("LASSO_ALL.RData")
ccle.tpms <- fread("CCLE_expression.csv", data.table = F)
swath.prot <- fread("Proteomics.csv", data.table = F)
tmt.prot <- fread("Uncalibrated_Proteomics.csv", data.table = F) %>%
  filter(!(GeneSymbol %in% swath.prot$GeneSymbol)) %>%
  select(colnames(swath.prot))
medians <- fread("StDevs_Meds_ALL.csv", data.table = F)
data_type <-fread("swath_vs_tmt_list.csv", header = T, drop = 'V1')

# Define UI ----
# Define UI ----
ui <- fluidPage(
  titlePanel("Pinferna"),
  sidebarLayout(
    sidebarPanel(
      helpText("Welcome to Pinferna! You may use this application to peruse relationships between a gene and its protein and upload transcriptomic data to receive protein abundance predictions."),
      selectizeInput("Gene", label = list(
        "Select gene:", 
        bsButton('select-info', label = '', icon = icon('info', lib = 'font-awesome'), style = 'default', size = 'extra-small')
      ), choices = rownames(BIC), multiple = F, selected = "A2M"),
      bsPopover(
        id = 'select-info',
        title = 'More information',
        content = 'The scroll bar shows a limited number, but you may type a gene of interest.',
        placement = 'right',
        trigger = 'hover',
        options = list(container = 'body')
      ),
      checkboxGroupInput("plotOptions", label = list(
                            "Plot Adjustments:",
                            bsButton('plot-info', label = '', icon = icon('info', lib = 'font-awesome'), style = 'default', size = 'extra-small')),
                         c("Show CCLE training data" = "showCCLE",
                           "Show regression" = "showRegression",
                           "Show 95% confidence interval" = "showCI",
                           "Show input data" = "showInput"),
                         selected = c("showCCLE",
                                      "showRegression",
                                      "showCI")),
      bsPopover(
        id = 'plot-info',
        title = 'More information',
        content = 'Show input data will only work after you have made predictions in a new dataset. Otherwise, selecting will display an error.',
        placement = 'right',
        trigger = 'hover',
        options = list(container = 'body')
      ),
      downloadButton("currentPlot", "Export Current Plot"),
      helpText("To make predictions, upload a table of RNA-seq data below (required). To incorporate LASSO, upload a table of the full transcriptomics (optional). Otherwise, only HL predictions are made. Click Go!"),
      fileInput("upload", label = list(
                  "Upload a file (required)", 
                  bsButton('file-info', label = '', icon = icon('info', lib = 'font-awesome'), style = 'default', size = 'extra-small')
                  ), buttonLabel = "Upload...", multiple = F,
                accept = c(".csv", ".tsv", ".xls", ".xlsx")),
      bsPopover(
        id = 'file-info',
        title = 'Input requirements',
        content = 'Upload RNA-seq data normalized as transcripts per million in either .csv, .tsv, or .xls(x) file format. Genes should be in rows; samples should be in columns. Gene identifiers should be in a column containing gene symbols and/or ensembl gene identifiers matching those in MANE Select (1). The column name should be “GeneSymbol” or “Ensembl_Gene”. An identifer conversion table and data entry template are provided with the publication associated with Pinferna (2). The maximum file size is 50 MB.',
        placement = 'right',
        trigger = 'hover',
        options = list(container = 'body')
      ),
      fileInput("upload_lasso", label = list(
                  "Upload a file (optional)", 
                  bsButton('lasso-info', label = '', icon = icon('info', lib = 'font-awesome'), style = 'default', size = 'extra-small')
                ), buttonLabel = "Upload for LASSO...", multiple = F,
                accept = c(".csv", ".tsv", ".xls", ".xlsx")),
      bsPopover(
        id = 'lasso-info',
        title = 'Input requirements',
        content = 'To incorporate LASSO, either 1) select the checkbox below to use the data already uploaded or 2) upload a second file of the full transcriptome from the RNA sequencing experiment normalized as TPM. The file format should be as above, and all (gene and sample) identifiers should match between the two tables.',
        placement = 'right',
        trigger = 'hover',
        options = list(container = 'body')
      ),
      checkboxInput("lasso_data", "Use the table I uploaded above.", value = F),
      actionButton("go", "Go!"),
      fileInput("upload_calibration", label = list(
                  "Upload a calibration file (optional)", 
                  bsButton('calibration-info', label = '', icon = icon('info', lib = 'font-awesome'), style = 'default', size = 'extra-small')
                ), buttonLabel = "Upload a file...", multiple = F,
                accept = c(".csv", ".tsv", ".xls", ".xlsx")),
      bsPopover(
        id = 'calibration-info',
        title = 'Input requirements',
        content = 'This is a table with absolute measurements for a specific protein (gene selected above). The format should have three columns: CellLine, GeneSymbol (or Ensembl_Gene), and Protein (as the absolute abundance in the cell line). The file may contain one or more samples and one or more genes.',
        placement = 'right',
        trigger = 'hover',
        options = list(container = 'body')
      ),
      actionButton("calibrate", "Calibrate"),
      
      # downloadButton("allPlots", "Export All Plots"),
      helpText("References: (1) Morales et al., Nature 604:310-315 (2022) \n (2) Sweatt et al., bioRxiv preprint doi: 10.1101/2023.07.10.548432.")
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
           csv  = set_identifiers(fread(input$upload$datapath, data.table = F)),
           tsv  = set_identifiers(fread(input$upload$datapath, data.table = F)),
           xls  = set_identifiers(read_excel(input$upload$datapath)),
           xlsx = set_identifiers(read_excel(input$upload$datapath)),
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
             csv  = set_identifiers(fread(input$upload_lasso$datapath, data.table = F)),
             tsv  = set_identifiers(fread(input$upload_lasso$datapath,data.table = F)),
             xls  = set_identifiers(read_excel(input$upload_lasso$datapath)),
             xlsx = set_identifiers(read_excel(input$upload_lasso$datapath)),
             validate("Error with file extension: Please upload a .csv, .tsv, or .xls(x) file.")
      )
    }
  })
  
  observeEvent(input$go, {
    choices <- as.character(unique(user_data()$GeneSymbol[user_data()$GeneSymbol %in% rownames(BIC)]))
    updateSelectizeInput(session, inputId = "Gene",
                         choices = choices, server = T)
  })
  
  calibration <- reactive({
    # req(input$upload_calibration)
    
    ext <- tools::file_ext(input$upload_calibration$name)
    switch(ext,
           csv  = set_identifiers(fread(input$upload_calibration$datapath, data.table = F)),
           tsv  = set_identifiers(fread(input$upload_calibration$datapath, data.table = F)),
           xls  = set_identifiers(read_excel(input$upload_calibration$datapath)),
           xlsx = set_identifiers(read_excel(input$upload_calibration$datapath)),
           validate("Error with file extension: Please upload a .csv, .tsv, or .xls(x) file.")
    )
  })
  
  predict_dat <- reactiveVal(NULL)
  
  observeEvent(input$go, {
    choices <- as.character(unique(user_data()$GeneSymbol[user_data()$GeneSymbol %in% rownames(BIC)]))
    updateSelectizeInput(session, inputId = "Gene",
                         choices = choices, server = T)
    progress <- Progress$new(session, min = 0, max = 1)
    on.exit(progress$close())
    progress$set(message = "Predictions in progress...",
                 detail = "For more than 100 genes, this may take a while.")
    p <- make_predictions(dat = user_data(), BIC = BIC, parameters = parameters, lasso.coefs = lasso.coefs, lasso_tpms = all_tpms(), meds = medians)
    predict_dat(p)
  })
  
  
  observeEvent(input$calibrate, {
    dat<-predict_dat()
    for (g in  names(dat)){
      scaling_factor<- get_calibration(g, calibration())
      dat[[g]] <-
        dat[[g]]  %>%
        mutate(Pred = Pred*scaling_factor)
    }
    predict_dat(dat)
    
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
    na_dat <- plot_dat()[["na_dat"]]
    type <- paste(plot_dat()[["type"]], 'copy numbers')
    
    if (input$calibrate){
      scaling_factor<- get_calibration(genesym, calibration())
      if (scaling_factor != 1){param_tbl<- rbind(param_tbl,c("Scaling_Factor",round(scaling_factor,2)))}
    } else {
      scaling_factor<- 1
    }
    xmin <- 10^(floor(log10(min(ccle_dat$TPM+1))))
    xmax <- 10^(ceiling(log10(max(ccle_dat$TPM+1))))
    ymin <- 10^(floor(log10(min(scaling_factor*(Regression$ymin+1)))))
    ymax <- 10^(ceiling(log10(max(scaling_factor*(Regression$ymax+1)))))
    
    if (nrow(na_dat > 0)) {
      na_dat$Protein <- ymin
      na_lab <- 'No protein data'
    } else {
      na_lab <- ''
    }
    
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
      xlab("TPM + 1")  + ggtitle(genesym, subtitle = type) +
      theme(aspect.ratio = 1, plot.title = element_text(h = 0.5),
            panel.background = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(color = "black", linewidth = 0.5),
            axis.ticks = element_line(color = "black", linewidth = 0.5),
            axis.ticks.length = unit(0.1, "inch"),
            axis.text = element_text(family = "sans", color = "black", size = 14),
            axis.title = element_text(family = "sans", color = "black", size = 16),
            title = element_text(family = "sans", color = "black", size = 18),
            plot.subtitle = element_text(family = "sans", color = "black", size = 16, hjust = 0.5))
    
    if(data_type[data_type$Gene == input$Gene,]$Absolute == 'yes'){
      g <- g +
        ylab("Protein copies + 1")
    } else if(input$calibrate){
      g <- g +
        ylab("Calibrated Protein copies + 1")
    }else{
      g <- g +
        ylab("Relative protein copies + 1")
    }
    
    
    if ("showCCLE" %in% input$plotOptions) {
      g <- g +
        geom_point(data = ccle_dat, aes(x = TPM+1, y =scaling_factor*(Protein+1)), color = "black", stroke = NA, alpha = 0.6, size = 2) +
        geom_point(data = na_dat, aes(x = TPM+1, y = Protein), color = "gray", shape = 4, size = 2) + 
        geom_text(aes(x = xmax, y = ymin, label = "CCLE data"), color = "black", hjust = "inward", vjust = "center", size = 6, check_overlap = T) +
        geom_text(aes(x = xmax, y = ymin, label = na_lab), color = 'gray', hjust = "inward", vjust = -1, size = 6, check_overlap = T)
    }
    if ("showRegression" %in% input$plotOptions) {
      if (BestModel == "Median") {
        g <- g +
          geom_smooth(data = Regression, aes(x = TPM+1, y = scaling_factor*(yfit+1)), color = "#28c742", linewidth = 0.75, se =F) +
          geom_text(aes(x = xmax, y = ymin, label = BestModel), color = "#28c742", hjust = "inward", vjust = -2.5, size = 6, check_overlap = T)
      } else if (BestModel == "HL") {
        g <- g +
          geom_smooth(data = Regression, aes(x = TPM+1, y = scaling_factor*(yfit+1)), color = "#aa3481", linewidth = 0.75, se =F) +
          geom_text(aes(x = xmax, y = ymin, label = BestModel), color = "#aa3481", hjust = "inward", vjust = -2.5, size = 6, check_overlap = T)
      } else {
        g <- g +
          geom_smooth(data = Regression, aes(x = TPM+1, y = scaling_factor*(yfit+1)), color = "#e66100", linewidth = 0.75, se =F) +
          geom_text(aes(x = xmax, y = ymin, label = BestModel), color = "#e66100", hjust = "inward", vjust = -2.5, size = 6, check_overlap = T)
      }
    }
    if ("showCI" %in% input$plotOptions) {
      Regression$ymin_rib <- 10^predict(loess(log10(scaling_factor*(Regression$ymin+1))~log10(Regression$TPM+1)))-1
      Regression$ymax_rib <- 10^predict(loess(log10(scaling_factor*(Regression$ymax+1))~log10(Regression$TPM+1)))-1
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
        geom_point(data = point_dat, aes(x = TPM + 1, y = (Pred + 1)), #####scaling_factor*
                   shape = 23, fill = "#e7298a", color = "black", size = 4) +
        geom_text(aes(x = xmax, y = ymin, label = 'User data'), color = "#e7298a", hjust = 'inward', vjust = -4, size = 6, check_overlap = T)
      
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
          list(extend = 'excel', filename = "prediction")
        ),
        text = "Download")
      ), scrollY = TRUE
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
          list(extend = 'excel', filename = "Upper95CI")
        ),
        text = "Download")
      ), scrollY = TRUE
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
          list(extend = 'excel', filename = "Lower95CI")
        ),
        text = "Download")
      ), scrollY = TRUE
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
