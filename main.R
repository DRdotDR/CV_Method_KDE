library(shiny)
library(ggplot2)
library(ggthemes)
library(MASS)
library(future)
library(future.apply)
library(compiler)
library(progressr)
library(Rcpp)

source("helpers.R")
enableJIT(3)

SAVE_CORES <- 0
USE_CORES <- max(1, min(availableCores() - SAVE_CORES, availableCores()))
plan(multisession, workers = USE_CORES)

ui <- fluidPage(
  titlePanel("Оцінка щільності суміші"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      h4("Компонента 1"),
      selectInput("type1", "Тип розподілу:", 
                  choices = c("Нормальний" = "normal", "Нормальний двогорбий" = "multi_normal",  "Рівномірний" = "uniform", "Коші" = "cauchy")),
      uiOutput("params1_ui"),
      
      hr(),
      
      h4("Компонента 2"),
      selectInput("type2", "Тип розподілу:", 
                  choices = c("Нормальний" = "normal", "Нормальний двогорбий" = "multi_normal", "Рівномірний" = "uniform", "Коші" = "cauchy")),
      uiOutput("params2_ui"),
      
      hr(),
      
      h4("Налаштування експерименту"),
      h5("L (Розмір вибірки)"),
      sliderInput("L", "", min = 50, max = 10000, value = 100, step = 50),
      
      h5("T (Кількість ітерацій)"),
      sliderInput("T", "", min = 10, max = 1000, value = 1000, step = 10),

      br(),
      actionButton("runExperiment", "Експеримент", class = "btn-primary", width = "100%"),
      br(),
      downloadButton("downloadResults", "Зберегти результати", width = "100%")
    ),
    
    mainPanel(
      width = 9,
      h4("Компонента 1 — Справжня щільність"),
      plotOutput("plot_density1", height = "250px"),
      
      h4("Компонента 2 — Справжня щільність"),
      plotOutput("plot_density2", height = "250px"),
      
      h4("Щільність суміщі компонентів"),
      plotOutput("plot_mixture", height = "250px"),
      
      hr(),
      h3("Результати експерименту"),
      
      h4("Компонента 1 — Порівняння оцінок та справжньої щільності"),
      plotOutput("plot_result1", height = "300px"),
      
      h4("Компонента 2 — Порівняння оцінок та справжньої щільності"),
      plotOutput("plot_result2", height = "300px"),
      
      h4("Результати моделювання"),
      tableOutput("ise_table")
    )
  )
)

server <- function(input, output, session) {
  experiment_results <- reactiveVal(NULL)
  experiment_running <- reactiveVal(FALSE)
  
  output$params1_ui <- renderUI({
    type <- input$type1
    if (type == "normal") {
      list(
        sliderInput("param1_mean", "Середнє:", min = -5, max = 5, value = 0, step = 0.5),
        sliderInput("param1_sd", "Станд. відх.:", min = 0.1, max = 3, value = 1, step = 0.1)
      )
    } else if (type == "uniform") {
      list(
        sliderInput("param1_min", "Мінімум:", min = -10, max = 5, value = 0, step = 0.5),
        sliderInput("param1_max", "Максимум:", min = -5, max = 10, value = 1, step = 0.5)
      )
    } else if (type == "cauchy") {
      list(
        sliderInput("param1_loc", "Локація:", min = -5, max = 5, value = 0, step = 0.5),
        sliderInput("param1_scale", "Масштаб:", min = 0.1, max = 3, value = 1, step = 0.1)
      )
    } else if (type == "multi_normal") {
      list(
        sliderInput("param1_mean", "Розділення:", min = -5, max = 5, value = 0, step = 0.5),
        sliderInput("param1_sd", "Станд. відх.:", min = 0.1, max = 3, value = 1, step = 0.1)
      )
    }
  })
  
  output$params2_ui <- renderUI({
    type <- input$type2
    if (type == "normal") {
      list(
        sliderInput("param2_mean", "Середнє:", min = -5, max = 5, value = 1, step = 0.5),
        sliderInput("param2_sd", "Станд. відх.:", min = 0.1, max = 3, value = 1, step = 0.1)
      )
    } else if (type == "uniform") {
      list(
        sliderInput("param2_min", "Мінімум:", min = -10, max = 5, value = 1, step = 0.5),
        sliderInput("param2_max", "Максимум:", min = -5, max = 10, value = 2, step = 0.5)
      )
    } else if (type == "cauchy") {
      list(
        sliderInput("param2_loc", "Локація:", min = -5, max = 5, value = 1, step = 0.5),
        sliderInput("param2_scale", "Масштаб:", min = 0.1, max = 3, value = 1, step = 0.1)
      )
    } else if (type == "multi_normal") {
      list(
        sliderInput("param2_mean", "Розділення:", min = -5, max = 5, value = 0, step = 0.5),
        sliderInput("param2_sd", "Станд. відх.:", min = 0.1, max = 3, value = 1, step = 0.1)
      )
    }
  })
  
  get_params1 <- reactive({
    type <- input$type1
    if (type == "normal") {
      list(type = "normal", param1 = input$param1_mean, param2 = input$param1_sd)
    } else if (type == "uniform") {
      list(type = "uniform", param1 = input$param1_min, param2 = input$param1_max)
    } else if (type == "cauchy") {
      list(type = "cauchy", param1 = input$param1_loc, param2 = input$param1_scale)
    } else if (type == "multi_normal") {
      list(type = "multi_normal", param1 = input$param1_mean, param2 = input$param1_sd)
    }
  })
  
  get_params2 <- reactive({
    type <- input$type2
    if (type == "normal") {
      list(type = "normal", param1 = input$param2_mean, param2 = input$param2_sd)
    } else if (type == "uniform") {
      list(type = "uniform", param1 = input$param2_min, param2 = input$param2_max)
    } else if (type == "cauchy") {
      list(type = "cauchy", param1 = input$param2_loc, param2 = input$param2_scale)
    } else if (type == "multi_normal") {
      list(type = "multi_normal", param1 = input$param2_mean, param2 = input$param2_sd)
    }
  })
  
  output$plot_density1 <- renderPlot({
    p1 <- get_params1()
    param1 <- if(p1$type == "uniform") min(p1$param1, p1$param2) else p1$param1
    param2 <- if(p1$type == "uniform") max(p1$param1, p1$param2) else p1$param2
    x_min <- if(p1$type == "uniform") param1 - 1 else param1 - 3 * param2
    x_max <- if(p1$type == "uniform") param2 + 1 else param1 + 3 * param2
    param1_safe <- ifelse(param1 == 0, 1, param1)
    x_min <- if(p1$type == "multi_normal") -3 * param2 * param1_safe else x_min
    x_max <- if(p1$type == "multi_normal") 3 * param2 * param1_safe  else x_max
    
    x <- seq(x_min, x_max, length.out = 500)
    y <- get_density_function(p1$type, param1, param2)(x)
    
    ggplot(data.frame(x=x, y=pmax(y,0)), aes(x=x, y=y)) +
      geom_line(color = "steelblue", size = 1) +
      geom_area(fill = "steelblue", alpha = 0.3) +
      theme_minimal() + labs(x = "x", y = "Щільність")
  })
  
  output$plot_density2 <- renderPlot({
    p2 <- get_params2()
    param1 <- if(p2$type == "uniform") min(p2$param1, p2$param2) else p2$param1
    param2 <- if(p2$type == "uniform") max(p2$param1, p2$param2) else p2$param2
    x_min <- if(p2$type == "uniform") param1 - 1 else param1 - 3 * param2
    x_max <- if(p2$type == "uniform") param2 + 1 else param1 + 3 * param2
    param1_safe <- ifelse(param1 == 0, 1, param1)
    x_min <- if(p2$type == "multi_normal") -3 * param2 * param1_safe else x_min
    x_max <- if(p2$type == "multi_normal") 3 * param2 * param1_safe  else x_max
    
    x <- seq(x_min, x_max, length.out = 500)
    y <- get_density_function(p2$type, param1, param2)(x)
    
    ggplot(data.frame(x=x, y=pmax(y,0)), aes(x=x, y=y)) +
      geom_line(color = "coral", size = 1) +
      geom_area(fill = "coral", alpha = 0.3) +
      theme_minimal() + labs(x = "x", y = "Щільність")
  })
  
  output$plot_mixture <- renderPlot({
    p1 <- get_params1()
    p2 <- get_params2()
    
    p1_1 <- if(p1$type == "uniform") min(p1$param1, p1$param2) else p1$param1
    p1_2 <- if(p1$type == "uniform") max(p1$param1, p1$param2) else p1$param2
    p2_1 <- if(p2$type == "uniform") min(p2$param1, p2$param2) else p2$param1
    p2_2 <- if(p2$type == "uniform") max(p2$param1, p2$param2) else p2$param2
    
    x_min <- -10
    x_max <- 10
    
    x <- seq(x_min, x_max, length.out = 500)
    y <- mixture_density(x, p1$type, p1_1, p1_2, p2$type, p2_1, p2_2, w1 = 0.5)
    
    ggplot(data.frame(x=x, y=pmax(y,0)), aes(x=x, y=y)) +
      geom_line(color = "darkgreen", size = 1) +
      geom_area(fill = "darkgreen", alpha = 0.3) +
      theme_minimal() + labs(x = "x", y = "Щільність")
  })
  
  observeEvent(input$runExperiment, {
    p1 <- get_params1()
    p2 <- get_params2()
    
    type1 <- p1$type
    param1_1 <- p1$param1
    param1_2 <- p1$param2
    type2 <- p2$type
    param2_1 <- p2$param1
    param2_2 <- p2$param2
    
    L <- input$L
    T <- input$T
    sigma <- 1
    
    withProgress(message = "Експеримент триває...", value = NULL, {
    with_progress({
      p <- progressor(along = 1:input$T)
      
      results <- future_lapply(
        1:input$T,
        function(i) {
          res <- run_iteration(
            iter = i,
            type1 = type1,
            param1_1 = param1_1,
            param1_2 = param1_2,
            type2 = type2,
            param2_1 = param2_1,
            param2_2 = param2_2,
            L = L,
            T = T,
            sigma = sigma
          )
          p()
          res
        },
        future.seed = TRUE
      )
    })
    })
    
    experiment_results(list(
      results = results,
      p1 = p1,
      p2 = p2,
      L = L,
      T = T,
      sigma = sigma
    ))
  })
  
  output$plot_result1 <- renderPlot({
    res <- experiment_results()
    if (is.null(res)) return(NULL)
    gen_plot_1(res, second = FALSE)
  })
  
  output$plot_result2 <- renderPlot({
    res <- experiment_results()
    if (is.null(res)) return(NULL)
    gen_plot_1(res, second = TRUE)
  })
  
  output$ise_table <- renderTable({
    res <- experiment_results()
    if (is.null(res)) return(NULL)
    results <- res$results
    
    Min1 <- sapply(results, `[[`, "Min1")
    Min2 <- sapply(results, `[[`, "Min2")
    Sil1 <- sapply(results, `[[`, "Sil1")
    Sil2 <- sapply(results, `[[`, "Sil2")
    ISE1Min <- sapply(results, `[[`, "ISE1Min")
    ISE2Min <- sapply(results, `[[`, "ISE2Min")
    ISE1Sil <- sapply(results, `[[`, "ISE1Sil")
    ISE2Sil <- sapply(results, `[[`, "ISE2Sil")

    df <- data.frame(
      Компонента = c("Комп. 1", "Комп. 2"),
      
      T = c(res$T, res$T),
      L = c(res$L, res$L),
      
      mean_hN_CV = c(mean(Min1), mean(Min2)),
      mean_hN_Silver = c(mean(Sil1), mean(Sil2)),
      
      var_hN_CV = c(var(Min1), var(Min2)),
      var_hN_Sil = c(var(Sil1), var(Sil2)),
      
      mean_ISE_CV = c(mean(ISE1Min), mean(ISE2Min)),
      mean_ISE_Sil = c(mean(ISE1Sil), mean(ISE2Sil)),
      
      var_ISE_CV = c(var(ISE1Min), var(ISE2Min)),
      var_ISE_Sil = c(var(ISE1Sil), var(ISE2Sil))
    )
    
    scientific_cols <- c(
      "var_hN_CV",
      "var_hN_Sil",
      "mean_ISE_CV",
      "mean_ISE_Sil",
      "var_ISE_CV",
      "var_ISE_Sil"
    )
    df[scientific_cols] <- lapply(df[scientific_cols], function(x) {
      sprintf("%.4e", x)
    })
    df$mean_hN_CV <- round(df$mean_hN_CV, 4)
    df$mean_hN_Silver <- round(df$mean_hN_Silver, 4)
    df
  }, striped = TRUE, hover = TRUE, bordered = TRUE)
  
  output$downloadResults <- downloadHandler(
    filename = function() {
      timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
      paste0("experiment_results_", timestamp, ".csv")
    },
    
    content = function(file) {
      res <- experiment_results()
      
      if (is.null(res)) {
        return(NULL)
      }
      
      results <- res$results
      
      Min1 <- sapply(results, `[[`, "Min1")
      Min2 <- sapply(results, `[[`, "Min2")
      Sil1 <- sapply(results, `[[`, "Sil1")
      Sil2 <- sapply(results, `[[`, "Sil2")
      
      ISE1Min <- sapply(results, `[[`, "ISE1Min")
      ISE2Min <- sapply(results, `[[`, "ISE2Min")
      ISE1Sil <- sapply(results, `[[`, "ISE1Sil")
      ISE2Sil <- sapply(results, `[[`, "ISE2Sil")
      
      results_frame <- data.frame(
        Component = c("Component 1", "Component 2"),
        T = c(res$T, res$T),
        L = c(res$L, res$L),
        
        mean_hN_CV = c(mean(Min1), mean(Min2)),
        mean_hN_Silver = c(mean(Sil1), mean(Sil2)),
        
        var_hN_CV = c(var(Min1), var(Min2)),
        var_hN_Silver = c(var(Sil1), var(Sil2)),
        
        mean_ISE_CV = c(mean(ISE1Min), mean(ISE2Min)),
        mean_ISE_Sil = c(mean(ISE1Sil), mean(ISE2Sil)),
        
        var_ISE_CV = c(var(ISE1Min), var(ISE2Min)),
        var_ISE_Sil = c(var(ISE1Sil), var(ISE2Sil))
      )
      
      write.csv(results_frame, file, row.names = FALSE)
    }
  )
}

shinyApp(ui = ui, server = server)