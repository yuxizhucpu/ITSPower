library(shiny)
library(simr)
library(lme4)
library(glmmTMB)
library(ggplot2)
library(foreach)
library(doParallel)
library(DT)
library(reshape2)
library(patchwork)

ui <- fluidPage(
  titlePanel("Power Analysis for Interrupted Time Series Designs with Cluster Effects"),
  tags$head(
    tags$style(HTML(".shiny-notification {position: fixed; top: 50%; left: 50%; transform: translate(-50%, -50%);}"))
  ),
  sidebarLayout(
    sidebarPanel(
      radioButtons("phases", "Select Study Design:", choices = c("Two Phases" = "two", "Three Phases" = "three"), selected = "three"),
      radioButtons("outcome_type", "Outcome Type:", choices = c("Binary" = "binary", "Continuous" = "continuous"), selected = "binary"),
      conditionalPanel(
        condition = "input.phases == 'two'",
        selectInput("twomodel", "Two Phase Model Type:", choices = c("Step Model" = "step", "Slope Model" = "slope"), selected = "step")
      ),
      conditionalPanel(
        condition = "input.phases == 'three'",
        selectInput("threemodel", "Three Phase Model Type:", choices = c("Step Model" = "step", "Slope-Plateau" = "slopeplateau"), selected = "step")
      ),
      numericInput("n_clusters", "Number of Clusters:", 37),
      numericInput("n_time_points", "Number of Time Points:", 42),
      conditionalPanel(
        condition = "input.phases == 'two'",
        textInput("firstphase", "Phase 1 End (comma-separated):", "6")
      ),
      conditionalPanel(
        condition = "input.phases == 'three'",
        numericInput("firstphase", "Phase 1 End:", 6),
        numericInput("secondphase", "Phase 2 End:", 30)
      ),
      textInput("icc_range", "ICC Values (comma-separated):", "0.1,0.3,0.5,0.7"),
      conditionalPanel(condition = "input.outcome_type == 'binary'", 
                       numericInput("p0", "Baseline Proportion (p0):", 0.53), 
                       numericInput("p1", "Target Proportion (p1):", 0.58)),
      conditionalPanel(condition = "input.outcome_type == 'continuous'", 
                       numericInput("y0", "Baseline Mean (y0):", 0.53), 
                       numericInput("y1", "Target Mean (y1):", 0.58), 
                       numericInput("sigma_y", "Intra-cluster variability (σ):", value = 1)),
      numericInput("n_sim", "Number of Simulations:", 500),
      textInput("sample_range", "Sample size per time (comma-separated):", "5,6,7,8,9,10,11,12,13,14,15"),
      actionButton("runSim", "Run Simulation"),
      downloadButton("downloadPower", "Download Power Table")
    ),
    mainPanel(
      verbatimTextOutput("runtime"),
      plotOutput("combinedPlot"),
      h4(textOutput("tableTitle")),
      DTOutput("powerTable")
      
    )
  ),
  div(class = "footer-fixed",
      helpText("⚠️ ICC Guidance:",
               "\nvalues below 0.50 indicate low correlation; 
               values between 0.50 and 0.75 indicate moderate correlation; 
               values between 0.75 and 0.90 suggest strong correlation; 
               and values above 0.90 reflect very strong correlation.")
  )
)

server <- function(input, output) {
  observeEvent(input$runSim, {
    powerf <- function(n_clusters, n_time_points, patients_per_timepoint,
                       firstphase, secondphase = NULL, icc, p0, p1, design = "three",
                       twomodel = "step", threemodel = "step", 
                       outcome_type = "binary", sigma_e = 1) {
      
      n_obs <- n_clusters * n_time_points * patients_per_timepoint
      data_sim <- expand.grid(cluster = 1:n_clusters,
                              time = 1:n_time_points,
                              patient = 1:patients_per_timepoint)
      
      if (design == "three") {
        data_sim$phase <- with(data_sim, cut(time, breaks = c(0, firstphase, secondphase, n_time_points),
                                             labels = c("Phase1", "Phase2", "Phase3"), include.lowest = TRUE))
        data_sim$phase2 <- as.numeric(data_sim$phase == "Phase2")
        data_sim$phase3 <- as.numeric(data_sim$phase == "Phase3")
        if (threemodel == "slopeplateau") {
          data_sim$time_slope <- ifelse(data_sim$phase == "Phase2", data_sim$time - firstphase, 0)
        }
      } else {
        data_sim$phase <- with(data_sim, ifelse(time <= firstphase, "Phase1", "Phase2"))
        if (twomodel == "step") {
          data_sim$phase2 <- as.numeric(data_sim$phase == "Phase2")
        } else if (twomodel == "slope") {
          data_sim$time_slope <- pmax(0, data_sim$time - firstphase)
        }
      }
      
      if (outcome_type== "binary"){
        varsigma <- sqrt(pi^2 / (3 * ((1 / icc) - 1)))
        random_intercepts <- rnorm(n_clusters, mean = 0, sd = varsigma)
        data_sim$random_effect <- random_intercepts[data_sim$cluster]
        
        base_log_odds <- -log((1 / p0) - 1)
        if (design == "three") {
          if (threemodel == "step") {
            phase2_effect <- -log((2 / (p0 + p1)) - 1) - base_log_odds
            phase3_effect <- -log((1 / p1) - 1) - (-log((2 / (p0 + p1)) - 1))
            log_odds <- with(data_sim, base_log_odds + phase2_effect * phase2 + phase3_effect * phase3 + random_effect)
          } else if (threemodel == "slopeplateau") {
            slope_effect <- (-log((1 / p1) - 1) - base_log_odds) / (secondphase - firstphase)
            log_odds <- with(data_sim, base_log_odds + slope_effect * time_slope + random_effect)
          }
        } else if (twomodel == "step") {
          phase2_effect <- -log((1 / p1) - 1) - base_log_odds
          log_odds <- with(data_sim, base_log_odds + phase2_effect * phase2 + random_effect)
        } else if (twomodel == "slope") {
          slope_effect <- (-log((1 / p1) - 1) - base_log_odds) / (n_time_points - firstphase)
          log_odds <- with(data_sim, base_log_odds + slope_effect * time_slope + random_effect)
        }
        
        data_sim$prob_success <- 1 / (1 + exp(-log_odds))
        data_sim$outcome <- rbinom(n = n_obs, size = 1, prob = data_sim$prob_success)
        
        model_formula <- if (design == "three") {
          if (threemodel == "step") {
            outcome ~ phase2 + phase3 + (1 | cluster)
          }
          else {
            outcome ~ time_slope + (1 | cluster)
          }
        } else if (design == "two") { if (twomodel == "step") {
          outcome ~ phase2 + (1 | cluster)}
          else {outcome ~ time_slope + (1 | cluster)}
        }
        
        model <- glmmTMB(model_formula, family = binomial(link = "logit"), data = data_sim)
        res <- summary(model)
        
        z_score <- if (design == "three") {
          if (threemodel == "step") {
            coefs <- fixef(model)$cond
            vcov_mat <- vcov(model)$cond
            
            coef_sum <- coefs["phase2"] + coefs["phase3"]
            se_sum <- sqrt(
              vcov_mat["phase2", "phase2"] +
                vcov_mat["phase3", "phase3"] +
                2 * vcov_mat["phase2", "phase3"]
            )
            coef_sum / se_sum
          } else {
            #res$coefficients$cond["time_slope", 1] / res$coefficients$cond["time_slope", 2]
            coefs <- fixef(model)$cond
            se <- sqrt(diag(vcov(model)$cond))["time_slope"]
            coefs["time_slope"] / se
          }
        } else if (twomodel == "step") {
          #res$coefficients$cond["phase2", 1] / res$coefficients$cond["phase2", 2]
          coefs <- fixef(model)$cond
          se <- sqrt(diag(vcov(model)$cond))["phase2"]
          coefs["phase2"] / se
        } else {
          #res$coefficients$cond["time_slope", 1] / res$coefficients$cond["time_slope", 2]
          coefs <- fixef(model)$cond
          se <- sqrt(diag(vcov(model)$cond))["time_slope"]
          coefs["time_slope"] / se
        }
      } else if (outcome_type== "continuous"){
        sigma_u <- sqrt((icc * sigma_e^2) / (1 - icc))
        random_intercepts <- rnorm(n_clusters, mean = 0, sd = sigma_u)
        data_sim$random_effect <- random_intercepts[data_sim$cluster]
        
        base_mean <- p0  # passed in as y0
        target_mean <- p1  # passed in as y1
        
        if (design == "three") {
          if (threemodel == "step") {
            delta2 <- (base_mean + target_mean) / 2 - base_mean
            delta3 <- target_mean - (base_mean + target_mean) / 2
            mu <- with(data_sim, base_mean + delta2 * phase2 + delta3 * phase3 + random_effect)
          } else {
            slope <- (target_mean - base_mean) / (secondphase - firstphase)
            mu <- with(data_sim, base_mean + slope * time_slope + random_effect)
          }
        } else {
          if (twomodel == "step") {
            delta <- target_mean - base_mean
            mu <- with(data_sim, base_mean + delta * phase2 + random_effect)
          } else {
            slope <- (target_mean - base_mean) / (n_time_points - firstphase)
            mu <- with(data_sim, base_mean + slope * time_slope + random_effect)
          }
        }
        
        data_sim$outcome <- rnorm(n_obs, mean = mu, sd = sigma_e)
        
        model_formula <- if (design == "three") {
          if (threemodel == "step") {
            outcome ~ phase2 + phase3 + (1 | cluster)} else {
              outcome ~ time_slope + (1 | cluster)}
        } else if (design == "two") { 
          if (twomodel == "step") {
            outcome ~ phase2 + (1 | cluster)
          } else {
            outcome ~ time_slope + (1 | cluster) 
          }
        }
        
        model <- lmer(model_formula, data = data_sim)
        
        coefs <- fixef(model)
        vcov_mat <- vcov(model)
        
        if (design == "three" && threemodel == "step") {
          coef_sum <- coefs["phase2"] + coefs["phase3"]
          se_sum <- sqrt(
            vcov_mat["phase2", "phase2"] +
              vcov_mat["phase3", "phase3"] +
              2 * vcov_mat["phase2", "phase3"]
          )
          z_score <- coef_sum / se_sum
        } else {
          key <- ifelse(design == "three" && threemodel == "slopeplateau", "time_slope",
                        ifelse(design == "two" && twomodel == "slope", "time_slope", "phase2"))
          z_score <- coefs[key] / sqrt(vcov_mat[key, key])
        }
        
      }
      
      pval <- 2 * pnorm(q = z_score, lower.tail = FALSE)
      return(pval < 0.05)
    }
    
    
    start_time <- Sys.time()
    
    n_clusters <- input$n_clusters
    n_time_points <- input$n_time_points
    design <- input$phases
    # try to add dynamic title for power table display
    twomodel1 <- paste(input$phases, "-phase", input$twomodel)   
    threemodel1 <- paste(input$phases, "-phase", input$threemodel)
    twomodel2 <- input$twomodel
    threemodel2 <- input$threemodel
    icc_list <- as.numeric(unlist(strsplit(input$icc_range, ",")))
    icc_list <- icc_list[!is.na(icc_list) & icc_list > 0]
    icc_list <- unique(icc_list)
    #[1:min(4, length(unique(icc_list)))]
    if (input$outcome_type == "binary") {
      p0 <- input$p0
      p1 <- input$p1
    } else {
      p0 <- input$y0
      p1 <- input$y1
    }
    n_sim <- input$n_sim
    sample_sizes <- as.numeric(unlist(strsplit(input$sample_range, ",")))
    
    outcome_type= input$outcome_type
    sigma_y=input$sigma_y
    
    firstphase_list <- if (design == "two") as.numeric(unlist(strsplit(as.character(input$firstphase), ","))) else input$firstphase
    if (design == "three") {
      secondphase <- input$secondphase
    }
    
    cl <- makeCluster(parallel::detectCores() - 1)
    registerDoParallel(cl)
    
    long_table <- data.frame()
    
    withProgress(message = 'Simulation in progress...', value = 0, {
      for (icc in icc_list) {
        for (firstphase in firstphase_list) {
          for (patients_per_timepoint in sample_sizes) {
            sim_results <- foreach(i = 1:n_sim, .combine = c, .packages = c("glmmTMB","lme4")) %dopar% {
              powerf(n_clusters, n_time_points, patients_per_timepoint,
                     firstphase, if (design == "three") secondphase else NULL, 
                     icc, p0, p1, design,twomodel2,threemodel2, outcome_type, sigma_y)
            }
            power <- formatC(mean(sim_results, na.rm = TRUE), digits = 2, format = "f")
            
            total_sample <- n_clusters * n_time_points * patients_per_timepoint
            long_table <- rbind(long_table, data.frame(ICC = icc,
                                                       Totaltime=n_time_points,
                                                       Phase1End = firstphase,
                                                       SampleSize_per_time = patients_per_timepoint, 
                                                       TotalSampleSize = total_sample,
                                                       Power = power,
                                                       Model = ifelse(design == "two", 
                                                                      twomodel1, threemodel1)))
            incProgress(1 / (length(sample_sizes) * length(icc_list) * length(firstphase_list)),
                        detail = paste("ICC:", icc, "Phase1:", firstphase, "Sample size:", patients_per_timepoint))
          }
        }
      }
    })
    
    stopCluster(cl)
    end_time <- Sys.time()
    
    output$runtime <- renderPrint({
      paste("Total Runtime:", round(difftime(end_time, start_time, units = "secs"), 2), "seconds")
    })
    
    output$combinedPlot <- renderPlot({
      if (input$phases == "two") {
        if (length(firstphase_list) == 1 && length(icc_list) == 1) {
          df <- subset(long_table, ICC == icc_list[1] & Phase1End == firstphase_list[1])
          ggplot(df, aes(x = SampleSize_per_time, y = as.numeric(Power))) +
            geom_line(color = "blue") +
            geom_point(color = "blue") +
            geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
            labs(title = paste("Power vs Sample Size (ICC =", icc_list[1], ", Phase1 End =", firstphase_list[1], ")"),
                 x = "Sample size per time",
                 y = "Estimated Power") +
            theme_minimal()
        } else if (length(firstphase_list) > 1 && length(icc_list) == 1) {
          df <- subset(long_table, ICC == icc_list[1])
          ggplot(df, aes(x = SampleSize_per_time, y = as.numeric(Power), color = as.factor(Phase1End))) +
            geom_line() +
            geom_point() +
            geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
            labs(title = paste("Power vs Sample Size (ICC =", icc_list[1], ")"),
                 x = "Sample size per time",
                 y = "Estimated Power",
                 color = "Phase1 End") +
            theme_minimal()
        } else if (length(firstphase_list) == 1 && length(icc_list) > 1) {
          df <- subset(long_table, Phase1End == firstphase_list[1])
          ggplot(df, aes(x = SampleSize_per_time, y = as.numeric(Power), color = as.factor(ICC))) +
            geom_line() +
            geom_point() +
            geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
            labs(title = paste("Power vs Sample Size (Phase1 End =", firstphase_list[1], ")"),
                 x = "Sample size per time",
                 y = "Estimated Power",
                 color = "ICC") +
            theme_minimal()
        } else if (length(firstphase_list) > 1 && length(icc_list) > 1) {
          plot_list <- lapply(seq_along(firstphase_list), function(i) {
            df <- subset(long_table, Phase1End == firstphase_list[i])
            ggplot(df, aes(x = SampleSize_per_time, y = as.numeric(Power), color = as.factor(ICC))) +
              geom_line() +
              geom_point() +
              geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
              labs(title = paste("Phase1 End =", firstphase_list[i]),
                   x = "Sample size per time",
                   y = "Estimated Power",
                   color = "ICC") +
              theme_minimal()
          })
          wrap_plots(plotlist = plot_list, ncol = 2)
        }
      } else {
        df <- long_table
        ggplot(df, aes(x = SampleSize_per_time, y = as.numeric(Power), color = as.factor(ICC))) +
          geom_line() +
          geom_point() +
          geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
          labs(title = "Power vs Sample Size",
               x = "Sample size per time",
               y = "Estimated Power",
               color = "ICC") +
          theme_minimal()
      }
    })
    
    output$tableTitle <- renderText({
      paste("Power Table for", ifelse(design == "two", twomodel1, threemodel1), "model")
    })
    
    output$powerTable <- renderDT({
      dcast(long_table, Totaltime+Phase1End+SampleSize_per_time + TotalSampleSize  ~ ICC, value.var = "Power")
    })
    
    output$downloadPower <- downloadHandler(
      filename = function() {
        paste("power_table_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        write.csv(long_table, file, row.names = FALSE)
      }
    )
    
  })
}
#, options = list(launch.browser = TRUE)
shinyApp(ui = ui, server = server)
