#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

require(shiny)
require(ggplot2)
require(shinyWidgets)
require(rpart)
require(e1071)
require(plotly)
require(plyr)
require(jsonlite)
require(kableExtra)
source('directoryInput.R')

# set the max file upload size to 20Mb
options(shiny.maxRequestSize = 20*1024^2)

# Define UI for application that draws a histogram
ui <- fluidPage(
  # Application title
  titlePanel("Tomato"),
  
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(
        tabPanel("Pre Processing",
                 directoryInput("directory", label = "Choose events director"),
        br(),
        actionButton('process', 'Process')
      ),
      tabPanel("load CSV",
               fileInput('datafile', 'Choose CSV file',
                         accept=c('text/csv', 'text/comma-separated-values,text/plain')),
               uiOutput("Pore_Number"),
               hr(),
               tableOutput("csvdata_copy")
      ),
      tabPanel("Modeling",
               h2("Data sets & Experiment"),
               tableOutput("csvdata"),
               
               hr(),
               h2("Load a pre-build model"),
               selectInput('pre_model', label = "Choose pre-built model", 
                           choice = list.dirs(paste(getwd(),'models',sep = '/'), 
                                              full.names = FALSE,recursive = FALSE)),
               actionButton("Load_Model", "Load Model"),
               hr(),
               h2("Build a model"),
               selectInput('model', label = "Choose classification algorithm", 
                            choice = c('Auto-Q', 'SVM-linear', 'SVM-Gaussian', 'SVM-Poly')),
               hr(),
               prettyCheckboxGroup(inputId = "paras",  label = "Choose parameters for model",
                                   choices = list("dwell" = "dwell", "Median Amplitude" = "medAmp", "Max Amplitude" = "maxAmp", "Area" = "NewArea"),
                                   outline = TRUE,
                                   plain = TRUE, icon = icon("thumbs-up")),
                 
                 hr(),
                 uiOutput("controls"),
                 actionButton('Build', 'Build Model'),
                 downloadButton('Save', 'Save Model')
               ),
      
      tabPanel("Reporter",
               fileInput("report_meta", "Meta Data",
                         accept = c("application/json", ".json")),
               hr(),
               checkboxInput("report_model_information", "Model information", TRUE),
               textAreaInput("report_model_notes", "Notes"),
               hr(),
               checkboxInput("report_scatter_plots", "Scatter Plot(s)", TRUE),
               uiOutput("report_scatter_plots_sets"),
               textAreaInput("report_scatter_plot_notes", "Notes"),
               hr(),
               checkboxInput("report_model_prediction_plots", "Model Prediction Plot(s)", TRUE),
               uiOutput("report_model_prediction_plot_sets"),
               textAreaInput("report_model_prediction_notes", "Notes"),
               hr(),
               checkboxInput("report_q_plots", "Q Plot(s)", TRUE),
               uiOutput("report_q_plot_sets"),
               textAreaInput("report_q_plots_notes", "Notes"),
               hr(),
               actionGroupButtons(inputIds = c("Generate_report","Publish_report"),
                                  labels = list("Generate","Publish"),
                                  status = "primary"),
               downloadButton("Download_report", "Download Report")
      )
      )
    ),

    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(id = "tabs",
        tabPanel("Raw Plot",
                 # 2D plot
                 h2("2D plot"),
                 radioButtons("X", label = h4("X"),
                              choices = list("dwell" = "dwell", "Median Amplitude" = "medAmp", "Max Amplitude" = "maxAmp", "Area" = "NewArea"),
                              selected = "dwell", inline = TRUE),
                 radioButtons("Y", label = h4("Y"),
                              choices = list("dwell" = "dwell", "Median Amplitude" = "medAmp", "Max Amplitude" = "maxAmp", "Area" = "NewArea"), 
                              selected = "maxAmp", inline = TRUE),
                 plotlyOutput("distPlot",height = "600px"),
                 
                 # 3D plot
                 h2("3D plot"),
                 radioButtons("X_3D", label = h4("X"),
                              choices = list("dwell" = "dwell", "Median Amplitude" = "medAmp", "Max Amplitude" = "maxAmp", "Area" = "NewArea"),
                              selected = "dwell", inline = TRUE),
                 radioButtons("Y_3D", label = h4("Y"),
                              choices = list("dwell" = "dwell", "Median Amplitude" = "medAmp", "Max Amplitude" = "maxAmp", "Area" = "NewArea"), 
                              selected = "medAmp", inline = TRUE),
                 radioButtons("Z_3D", label = h4("Z"),
                              choices = list("dwell" = "dwell", "Median Amplitude" = "medAmp", "Max Amplitude" = "maxAmp", "Area" = "NewArea"), 
                              selected = "maxAmp", inline = TRUE),
                 plotlyOutput("distPlot_3D",height = "600px")
                 ),
        
        tabPanel("Class Plot",uiOutput("set_num"), plotlyOutput("classPlot"),plotOutput("barplot")),
        tabPanel("Summary", uiOutput("Summary_Sets") ,tableOutput("table"), downloadButton("downloadData", "Download")),
        tabPanel("Q Plot", selectInput('target_molecule',label = "Target molecule", choices = NA), 
                 selectInput('calibrant_molecule',label = "Calibrant/Negative molecule", choices = NA), uiOutput("Q_Plots_Sets"), plotOutput("Q_Plots"),plotOutput("Q_barplot")),
        tabPanel("FA Calculator",
                 h4('Documentation'),
                 HTML("This is a simple FA Caculator used in single plex molecule test (Such as Monsanto 3-primer). For the inputs: <br/>
                      1. False Positive rate [0-1]: The model prediction of target moleucle's percentage (in decimal) on Negative or Calibrant molecule. Leave the value to 0 if you do not have this run. <br/>
                      2. True Positive rate [0-1]: The model's prediction of target moleucle's percentage (in decimal) on Positive sample. Leave the value to 1 if you do not have this run. <br/>
                      3. 50/50 correction [0-1]: The model's prediction of target moleucle's percentage (in decimal) on the 50/50 sample. Leave the value to 0.5 if you do not have this run. <br/>
                      4. Unknown [0-1]: The model's prediction of target moleucle's percentage (in decimal) on the unknonw sample.<br/>"),
                 hr(),
                 numericInput("FP_Rate", "False Positive rate: ", 0, min = 0, max = 1, step = 0.1),
                 numericInput("TP_Rate", "True Positive rate: ", 1, min = 0, max = 1, step = 0.1),
                 numericInput("correction_50", "50/50 correction: ", 0.5, min = 0, max = 1, step = 0.1),
                 numericInput("Unknown", "Unknown: ", 0.5, min = 0, max = 1, step = 0.1),
                 hr(),
                 h4('Unknown prediction after correction:'),
                 verbatimTextOutput("Unknown_after_correction")
        )
      ))
  )
)

# 
server <- function(input, output, session) {
  #This function is repsonsible for loading in the selected file
  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$directory
    },
    handlerExpr = {
      if (input$directory > 0) {
        # condition prevents handler execution on initial app launch
        
        # launch the directory selection dialog with initial path read from the widget
        path = choose.dir(default = readDirectoryInput(session, 'directory'))
        
        # update the widget value
        updateDirectoryInput(session, 'directory', value = path)
      }
    }
  )
  
  observeEvent(input$process, {
    dir_path <- readDirectoryInput(session, 'directory')
    system(paste('python3 main.py',dir_path, dir_path, sep = " "))
    showNotification("Preprocessing done!")
  })
  
  filedata <- reactive({
    infile <- input$datafile
    if (is.null(infile)) {
      # User has not uploaded a file yet
      return(NULL)
    }
    
    na.omit(read.csv(infile$datapath,row.names = 1))
  })
  
  observeEvent(input$Load_Model, {
    model <<- readRDS(paste(getwd(),'models', input$pre_model, 'model.rds',sep = '/'))
    
    if (grepl("rpart", toString(model$call))) {
      # predict_type <- "class"
      class_choices <- attr(model, "ylevels")
    }
    else {
      # predict_type <- "raw"
      class_choices <- model$levels
    }
    
    updateSelectInput(session, "target_molecule", choices = class_choices)
    updateSelectInput(session, "calibrant_molecule", choices = class_choices)
    
    showNotification("Model loaded!")
  })
  
  observeEvent(input$Build, {
    pore_num <- input$Pore_Number
    model_selected <- input$model
    para_list <- input$paras
    controls <- input$pure_control
    
    df <-filedata()
    # df_training <- df[which(df$PoreNum == as.integer(input$Pore_Number) & df$setNum %in% c(as.integer(pos_controls), as.integer(neg_controls))), c(para_list,'Exp')]
    df_training <- df[which(df$PoreNum == as.integer(input$Pore_Number) & df$setNum %in% as.integer(controls)), c(para_list,'Exp')]
    df_training$Exp <- factor(droplevels(df_training$Exp))
    switch (input$model,
            "Auto-Q" = {model <<- rpart(Exp ~., data = df_training)},
            "SVM-linear" = {model <<- svm(Exp ~., data = df_training, kernel = "linear")},
            "SVM-Gaussian" = {model <<- svm(Exp ~., data = df_training, kernel = "radial")},
            "SVM-Poly" = {model <<- svm(Exp ~., data = df_training, kernel = "polynomial")}
    )
    
    if (grepl("rpart", toString(model$call))) {
      # predict_type <- "class"
      class_choices <- attr(model, "ylevels")
    }
    else {
      # predict_type <- "raw"
      class_choices <- model$levels
    }
    
    updateSelectInput(session, "target_molecule", choices = class_choices)
    updateSelectInput(session, "calibrant_molecule", choices = class_choices)
    
    showNotification("Model building is finished!")
  })
  
  output$Save <- downloadHandler(
    filename = function() {paste("model",".rds",sep = "")},
    content = function(file) {
      saveRDS(model,file)
    }
  )
  
  #The following set of functions populate the column selectors
  output$Pore_Number <- renderUI({
    df <-filedata()
    items=unique(df$PoreNum)
    selectInput("Pore_Number", "Pore Number:",items)
  })
  
  output$set_num <- renderUI({
    df <-filedata()
    items <- unique(df[which(df$PoreNum == as.integer(input$Pore_Number)),]$setNum)
    items <- items[order(strtoi(items))]
    selectInput("Set_Number", "Set Number:",items)
  })
  
  output$distPlot <- renderPlotly({
    # scatter plot
    df <-filedata()
    
    # if(is.null(df)) {
    #   return(NULL)
    # }

    df_selected <- df[which(df$PoreNum == as.integer(input$Pore_Number)),]
    x_para <- input$X
    y_para <- input$Y
    # qplot(df_selected[,x_para], df_selected[,y_para], xlab=x_para,ylab=y_para, geom = 'point')
    # color_dataset <- paste0(df_selected$setNum,'-',df_selected$Exp)
    # print(length(color_dataset))
    p <- plot_ly(x=df_selected[,x_para], y=df_selected[,y_para], mode = "markers", type = "scatter", color = factor(paste(df_selected$setNum,df_selected$Exp,sep = "-")), colors = "Set1") %>%
    layout(xaxis=list(title=x_para), yaxis=list(title=y_para))
    p$elementId <- NULL
    # suppressWarnings(print(p))
    p
  })
  
  # 3D plot
  output$distPlot_3D <- renderPlotly({
    # scatter plot
    df <-filedata()
    
    df_selected <- df[which(df$PoreNum == as.integer(input$Pore_Number)),]
    x_para <- input$X_3D
    y_para <- input$Y_3D
    z_para <- input$Z_3D
    p_3D <- plot_ly(as.data.frame(df_selected), x=df_selected[,x_para], y=df_selected[,y_para], z=df_selected[,z_para], mode = "markers", type = "scatter3d", 
                    color = factor(paste(df_selected$setNum,df_selected$Exp,sep = "-")),colors = "Set1") %>%
      layout(scene = list(xaxis=list(title=x_para), yaxis=list(title=y_para), zaxis=list(title=z_para)))
    p_3D$elementId <- NULL
    # suppressWarnings(print(p_3D))
    p_3D
  })
  
  output$classPlot <- renderPlotly({
    # class plot
    df <-filedata()
    
    df_selected_class <<- df[which(df$PoreNum == as.integer(input$Pore_Number) & df$setNum == input$Set_Number),]

    if(grepl("rpart",toString(model$call))) {
      predict_type <- "class"
    }
    else {
      predict_type <- "raw"
    }
      
    class_type <<- predict(model, df_selected_class, type = predict_type)
    # write.csv(class_type,paste(input$Set_Number,'csv',sep = "."))
    x_para <- input$X
    y_para <- input$Y
    # qplot(df_selected_class[,x_para], df_selected_class[,y_para], color=class_type,xlab=x_para,ylab=y_para,geom = 'point')
    p_class <- plot_ly(as.data.frame(df_selected_class), x=df_selected_class[,x_para], y=df_selected_class[,y_para], mode="markers", type="scatter",color = class_type) %>%
      layout(xaxis=list(title=x_para), yaxis=list(title=y_para))
    p_class$elementId <-NULL
    # suppressWarnings(print(p_class))
    p_class
  })
  
  output$barplot <- renderPlot({
    df <-filedata()
    df_selected_class <- df[which(df$PoreNum == as.integer(input$Pore_Number) & df$setNum == input$Set_Number),]
    
    if(grepl("rpart",toString(model$call))) {
      predict_type <- "class"
    }
    else {
      predict_type <- "raw"
    }
    
    class_type <- predict(model, df_selected_class, type = predict_type)

    b <- barplot(100*table(class_type)/sum(table(class_type)), main="Precentage of each class", xlab = "Molecule",ylab="Percentage (%)",ylim = c(0,110))
    text(x=b, y=100*table(class_type)/sum(table(class_type))+3,
         labels = as.character(round(100*table(class_type)/sum(table(class_type)),digits = 2)))
  })
  
  output$Summary_Sets <- renderUI({
    df <- filedata()
    df_selected_pore <- df[which(df$PoreNum == as.integer(input$Pore_Number)),]
    dset_name <- unique(df_selected_pore$setNum)
    dset_name <- dset_name[order(strtoi(dset_name))]
    prettyCheckboxGroup(inputId = "Summary_sets_group",  label = "Select datasets for Summary",
                        choices = dset_name,
                        outline = TRUE,
                        plain = TRUE, icon = icon("thumbs-up"), inline = TRUE)
  })
  
  output$table <- renderTable({
    df <-filedata()
    
    output_df = data.frame(
      'Dataset' = as.integer(),
      'Molecule' = as.character(),
      'number of events' = as.integer(),
      'percentage(%)' = as.double(),
      'error(95%)' = as.double(),
      check.names = FALSE
    )
    
    for (Set_Number in input$Summary_sets_group) {
      df_selected_class <- df[which(df$PoreNum == as.integer(input$Pore_Number) & df$setNum == Set_Number),]
      
      if(grepl("rpart",toString(model$call))) {
        predict_type <- "class"
      }
      else {
        predict_type <- "raw"
      }
      
      class_type <- predict(model, df_selected_class, type = predict_type)
      class_names = unique(class_type)
      
      temp <- data.frame(
        'Dataset' = as.integer(Set_Number),
        'Molecule' = as.character(''),
        'number of events' = as.integer(''),
        'percentage(%)' = as.double(''),
        'error(95%)' = as.double(''),
        check.names = FALSE
      )
      output_df = rbind(output_df,temp)
      
      for (item in class_names) {
        q_value <- sum(class_type == item)/length(class_type)
        temp <- data.frame('Dataset' =' ',
                           'Molecule' = item,
                           'number of events' = sum(class_type == item),
                           'percentage(%)' = round(100*q_value,digits = 2),
                           'error(95%)' = round(100*1.96*sqrt(q_value*(1-q_value)/length(class_type)),digits = 2),
                           check.names = FALSE)
        output_df = rbind(output_df,temp)
      }
    }
    output_df_buff <<- output_df
    output_df
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("summmary",".csv",sep="")
    },
    content = function(file) {
      write.csv(output_df_buff,file,row.names = FALSE)
    }
  )
  
  output$csvdata <- renderTable({
    df <-filedata()
    df_selected_pore <- df[which(df$PoreNum == as.integer(input$Pore_Number)),]
    dset_name <- unique(df_selected_pore[c("setNum","Exp")])
    dset_name[order(strtoi(dset_name$setNum)),]
  })
  
  output$csvdata_copy <- renderTable({
    df <-filedata()
    df_selected_pore <- df[which(df$PoreNum == as.integer(input$Pore_Number)),]
    dset_name <- unique(df_selected_pore[c("setNum","Exp")])
    dset_name[order(strtoi(dset_name$setNum)),]
  })
  
  output$Q_Plots_Sets <- renderUI({
    df <- filedata()
    df_selected_pore <- df[which(df$PoreNum == as.integer(input$Pore_Number)),]
    dset_name <- unique(df_selected_pore$setNum)
    dset_name <- dset_name[order(strtoi(dset_name))]
    
    prettyCheckboxGroup(inputId = "Q_sets",  label = "Select datasets for Q-Plot",
                            choices = dset_name,
                            outline = TRUE,
                            plain = TRUE, icon = icon("thumbs-up"), inline = TRUE)
  })
  
  Q_df_data <- reactive({
    df <- filedata()
    target_molecule <- input$target_molecule
    calibrant_molecule <- input$calibrant_molecule
    set_list <- input$Q_sets
    
    if(grepl("rpart",toString(model$call))) {
      predict_type <- "class"
    }
    else {
      predict_type <- "raw"
    }
    
    Q_df <- data.frame('index' = as.integer(), 'q_value' = as.double(), 'q_error' = as.double(), 'd_set' = as.integer())
    for (sub_set in set_list) {
      df_selected_subset <- df[which(df$PoreNum == as.integer(input$Pore_Number) & df$setNum == sub_set),]
      class_type_subset <- predict(model, df_selected_subset, type = predict_type)
      q_value <- vector(length = length(class_type_subset), mode = "double")
      q_error <- vector(length = length(class_type_subset), mode = "double")
      
      num_target <- 0
      num_all <- 0
      for (i in 1:length(class_type_subset)) {
        if(i == 1) {
          q_value[1] <- as.integer(class_type_subset[1] == target_molecule)
          q_error[i] <- 0
          num_target <- as.integer(class_type_subset[1] == target_molecule)
          num_all <- as.integer(class_type_subset[1] == target_molecule | class_type_subset[1] == calibrant_molecule)
          if(num_all == 0) {
            q_value[i] <- 0
            q_error[i] <- 0
          }
          else {
            q_value[i] <- num_target/num_all
            q_error[i] <- 1.96*sqrt(q_value[i]*(1-q_value[i])/num_all)
          }
        }
        
        else {
          num_target <- num_target + as.integer(class_type_subset[i] == target_molecule)
          num_all <- num_all + as.integer(class_type_subset[i] == target_molecule | class_type_subset[i] == calibrant_molecule)
          if(num_all == 0) {
            q_value[i] <- 0
            q_error[i] <- 0
          }
          else {
            q_value[i] <- num_target/num_all
            q_error[i] <- 1.96*sqrt(q_value[i]*(1-q_value[i])/num_all)
          }
        }
      }
      Q_df <- rbind(Q_df,data.frame('index' = 1:length(class_type_subset), 
                                    'q_value' = q_value, 
                                    'q_error' = q_error, 
                                    'd_set' = as.vector(rep(sub_set, length(class_type_subset)))))
    }
    Q_df
  })
  
  output$Q_Plots <- renderPlot({
    Q_df <- Q_df_data()
    ggplot(data = Q_df, aes(x=index, y=100*q_value, colour = d_set)) + 
      geom_errorbar(data=Q_df, aes(ymin=100*(q_value-q_error), ymax=100*(q_value+q_error))) + 
      geom_line() +
      labs(x = "Number of Events", y = "Q Value")
  })
  
  output$Q_barplot <- renderPlot({
    Q_df <- Q_df_data()
    if (!empty(Q_df)) {
      Q_df_split <- split(Q_df,f=Q_df$d_set)
      Q_bar <- do.call(rbind.data.frame,lapply(Q_df_split, function(x) tail(x,1)))
      # b_p <- barplot(100*(Q_bar$q_value), main = "Q Values of each data set", xlab = "data set number", ylab = "Precentage(%)", col = Q_bar$d_set, ylim = c(0,110), names.arg = Q_bar$d_set)
      b_p <- ggplot(data = Q_bar, aes(x=d_set,y=100*q_value, fill=d_set)) +
        geom_bar(color=Q_bar$d_set,stat = "identity") +
        labs(x="dataset number",y="Precentage(%)") +
        geom_text(aes(label=paste(round(100*(Q_bar$q_value),digits = 2), round(100*(Q_bar$q_error),digits = 2), sep = " +/- ")), vjust=1.6)
      b_p
    }
  })
  
  output$controls <- renderUI({
    df <-filedata()
    df_selected_pore <- df[which(df$PoreNum == as.integer(input$Pore_Number)),]
    dset_name <- unique(df_selected_pore$setNum)
    dset_name <- dset_name[order(strtoi(dset_name))]
    
    prettyCheckboxGroup(inputId = "pure_control",  label = "Pure Controls",
                        choices = dset_name,
                        outline = TRUE,
                        plain = TRUE, icon = icon("thumbs-up"), inline = TRUE)
  })
  
  output$report_scatter_plots_sets <- renderUI({
    df <- filedata()
    dset_name <- unique(df$setNum)
    dset_name_Exp <- unique(df[,c("setNum","Exp")])
    dset_name_Exp <- paste(dset_name_Exp$setNum, dset_name_Exp$Exp, sep = " - ")
    dset_name_Exp <- dset_name_Exp[order(strtoi(dset_name))]
    prettyCheckboxGroup(inputId = "report_scatter_sets",  label = "Choose plotting dataset(s)",
                        choices = dset_name_Exp,
                        outline = TRUE,
                        plain = TRUE, icon = icon("thumbs-up"), inline = TRUE)
  })
  
  output$report_model_prediction_plot_sets <- renderUI({
    df <- filedata()
    dset_name <- unique(df$setNum)
    dset_name_Exp <- unique(df[,c("setNum","Exp")])
    dset_name_Exp <- paste(dset_name_Exp$setNum, dset_name_Exp$Exp, sep = " - ")
    dset_name_Exp <- dset_name_Exp[order(strtoi(dset_name))]
    prettyCheckboxGroup(inputId = "report_prediction_sets",  label = "Choose plotting dataset(s)",
                        choices = dset_name_Exp,
                        outline = TRUE,
                        plain = TRUE, icon = icon("thumbs-up"), inline = TRUE)
  })
  
  output$report_q_plot_sets <- renderUI({
    df <- filedata()
    dset_name <- unique(df$setNum)
    dset_name_Exp <- unique(df[,c("setNum","Exp")])
    dset_name_Exp <- paste(dset_name_Exp$setNum, dset_name_Exp$Exp, sep = " - ")
    dset_name_Exp <- dset_name_Exp[order(strtoi(dset_name))]
    prettyCheckboxGroup(inputId = "report_Q_sets",  label = "Choose plotting dataset(s)",
                        choices = dset_name_Exp,
                        outline = TRUE,
                        plain = TRUE, icon = icon("thumbs-up"), inline = TRUE)
  })
  
  observeEvent(input$Generate_report, {
    # parameters to be passed to report markdown
    para_config <- vector(mode = "list")
    para_config$Model_Info_Enable = input$report_model_information
    para_config$Model_Info$X = input$X
    para_config$Model_Info$Y = input$Y
    para_config$Model = input$model
    para_config$Model_Para = as.list(input$paras)
    para_config$Model_dsets = as.vector(input$pure_control)
    para_config$Model_Info_Notes = input$report_model_notes
    para_config$Target_Molecule = input$target_molecule
    para_config$Calibrant_Molecule = input$calibrant_molecule
    para_config$Scatter_Plot_Enable = input$report_scatter_plots
    para_config$Scatter_Plot_dsets = as.vector(input$report_scatter_sets)
    para_config$Scatter_Plot_Notes = input$report_scatter_plot_notes
    para_config$Model_Prediction_Plot_Enable = input$report_model_prediction_plots
    para_config$Model_Prediction_Plot_dsets = as.vector(input$report_prediction_sets)
    para_config$Model_Prediction_Plot_Notes = input$report_model_prediction_notes
    para_config$Q_Plot_Enable = input$report_q_plots
    para_config$Q_Plot_dsets = as.vector(input$report_Q_sets)
    para_config$Q_Plot_Notes = input$report_q_plots_notes
    # temp Experiment_Info
    report_meta_data <- fromJSON(input$report_meta$datapath)
    Exp_Info <<- report_meta_data$experiment_metadata
    report.params <- list(Events = filedata(), Model_file = model, Experiment_Info = Exp_Info, config = para_config)
    
    # pass parameters to markdown
    tempReport <- file.path(getwd(), "report.Rmd")
    rmarkdown::render(tempReport, output_file = "report.html", params = report.params)
    
    # insert preview Tab
    insertTab(inputId = "tabs",
      tabPanel("Preview Report", includeHTML("report.html")),
      target = "FA Calculator", select = TRUE)
  })
  
  output$Download_report <- downloadHandler(
    filename = function() {
      return("report.html")
    },
    content = function(con) {
      file.copy("report.html", con)
    }
  )
  
  observeEvent(input$Publish_report, {
    library(RMySQL)
    library(DBI)
    library(rdrop2)
    con <- dbConnect(MySQL(),
                     host = "sql3.freemysqlhosting.net",
                     dbname = "sql3250877",
                     user = "sql3250877",
                     password = "H6MwbuquLr",
                     port = 3306)
    dbWriteTable(con,"2PG_Experiment_DB", data.frame("Project" = Exp_Info$project,
                                                     "Goal" = Exp_Info$goal,
                                                     "Date" = Exp_Info$datetime,
                                                     "Experimenter" = Exp_Info$user,
                                                     "Pore_Number" = Exp_Info$pore_number),
                 append = TRUE,
                 row.names = FALSE
    )
    
    dropbox_file_path <- file.path("2PG_Reports", Exp_Info$project, Exp_Info$goal, Exp_Info$datetime, Exp_Info$user, Exp_Info$pore_number)
    drop_upload("report.html", path = dropbox_file_path)
    showNotification("Report Published!")
  })
  
  output$Unknown_after_correction <- renderText({
    if (input$TP_Rate < input$FP_Rate | input$correction_50 < input$FP_Rate | input$correction_50 > input$TP_Rate) {
      final <- 'ERROR'
    } else if (input$Unknown < input$FP_Rate) {
      final <- 0
    } else if (input$Unknown > input$TP_Rate) {
      final <- 1
    } else {
      Cali_after_TFP <- (input$correction_50 - input$FP_Rate)/(input$TP_Rate - input$FP_Rate)
      unknown_after_TFP = (input$Unknown - input$FP_Rate)/(input$TP_Rate - input$FP_Rate)
      final = 1/(1+(Cali_after_TFP*(1-unknown_after_TFP)/(unknown_after_TFP*(1-Cali_after_TFP))))
    }
    final
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)