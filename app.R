library(shiny)
library(ggplot2)
library(distributional)
library(tibble)
library(dplyr)
library(tidyr)
library(viridis)
library(ggridges)
library(ggforce)
library(patchwork)
library(purrr)
library(scales)

ui <- fluidPage(

    # Application title
    titlePanel("Bacteria/Virus Plate Dynamics"),

    # Sidebar with slider inputs 
    sidebarLayout(
        sidebarPanel(
            sliderInput("num_virus_per_colony",
                        "Number of Virus Particles per Colony:",
                        min = 1,
                        max = 100,
                        value = 20),
            checkboxInput("any_virus",
                          "Show all virus results"),
            sliderInput("num_bacteria",
                        "Number of Colony-Forming Units:",
                        min = 1,
                        max = 100,
                        value = 20),
            checkboxInput("any_bacteria",
                          "Show all colony-forming unit results"),
            sliderInput("size_plate", 
                        "Diameter of Plate (cm)",
                        min = 5,
                        max = 30,
                        value = 10),
            sliderInput("virus_within",
                        "Colony Diameter (mm)",
                        min = 1, 
                        max = 20,
                        step = 0.5,
                        value = 5)
        ),
        
        # Show a plot of the generated distribution of the plate plots
        mainPanel(
          plotOutput("distPlot", height = "600px")
        )
    )
)

server <- function(input, output) {
  
  pi <- 3.14159
  
  draw_plate <- function(num_bacteria, num_virus, size_plate, virus_within, title) {
    radius_plate <- size_plate / 2
    radius_colony <- (virus_within / 10) / 2  # Convert virus_within from mm to cm, then get radius
    effective_radius <- radius_plate - radius_colony  # Adjusted radius for colony center generation
    
    # Create data frames for bacteria colonies and virus particles
    bacteria_df <- data.frame(r = effective_radius * sqrt(runif(num_bacteria)),
                              theta = runif(num_bacteria) * 2 * pi) |> 
      transmute(x = r * cos(theta), y = r * sin(theta))
    
    virus_df <- data.frame(r = radius_plate * sqrt(runif(num_virus)),
                           theta = runif(num_virus) * 2 * pi) |> 
      transmute(x = r * cos(theta), y = r * sin(theta))
    
    # Initialize a column to track overlap
    bacteria_df$overlap <- FALSE
    
    # Check for overlap
    for (i in 1:nrow(bacteria_df)) {
      for (j in 1:nrow(virus_df)) {
        distance <- sqrt((bacteria_df$x[i] - virus_df$x[j])^2 + (bacteria_df$y[i] - virus_df$y[j])^2)
        if (distance < radius_colony) {
          bacteria_df$overlap[i] <- TRUE
          break  # Stop checking once an overlap is found
        }
      }
    }
    
    count_overlap = sum(bacteria_df$overlap)
    
    # Calculate the size for the colonies based on the specified diameter
    colony_size <- virus_within * 1  # Adjust this factor as needed for visual representation
    
    # Create the plot
    p <- ggplot() +
      geom_circle(aes(x0 = x, y0 = y, r = radius_colony, fill = overlap), data = bacteria_df, alpha = 0.5) +
      scale_fill_manual(values = c("FALSE" = "blue", "TRUE" = "red")) +
      geom_point(aes(x = x, y = y), data = virus_df, color = "red", size = 1) +
      geom_circle(aes(x0 = 0, y0 = 0, r = radius_plate)) + 
      coord_fixed(ratio = 1) +
      labs(color = "Colony-Virus Overlap")+
      xlim(-radius_plate, radius_plate) +
      ylim(-radius_plate, radius_plate) +
      theme_void() + 
      labs(title = paste0(title, ": ", count_overlap, " infected.")) 
    
    return(p)
  }
  
  n_bac <- reactive(input$num_bacteria)
  n_virus_per_colony <- reactive(input$num_virus_per_colony)
  n_virus <- reactive(input$num_virus_per_colony * input$num_bacteria)
  a_plate <- reactive(pi * (0.5 * 10 * input$size_plate)^2)
  a_virus <- reactive(pi * (0.5 * input$virus_within)^2)
  prob_one_miss <- reactive((a_plate() - a_virus())/a_plate())
  
  # for a single pair of virus, bact
  
  prob_x_miss <- reactive(prob_one_miss()^n_virus())
  prob_x_hit <- reactive(1 - prob_x_miss())
  expected_hits <- reactive(n_bac() * prob_x_hit())
  
  prob_table <- reactive(tibble(within = 0:n_bac(), 
                                prob = density(dist_binomial(n_bac(), prob_x_hit()), 0:n_bac())[[1]]))
  
  # for a specific number of colonies + any number of virus

  prob_x_hit_free_v <- reactive({
    data.frame(n_virus = n_bac() * seq(10, 100, 10)) |>
      mutate(prob_miss = prob_one_miss()^n_virus) |>
      mutate(prob_hit = 1 - prob_miss,
             expected = n_bac() * prob_hit,
             dens = density(dist_binomial(n_bac(), prob_hit), 0:n_bac())) |>
      tidyr::unnest_longer(col = dens, indices_include = TRUE)
    })
  
  
  # for a specific number of virus particles + any number of colonies
    
  prob_x_hit_free_b <- reactive({
    data.frame(n_bact = seq(10, 100, 10)) |>
      mutate(prob_miss = prob_one_miss()^n_virus()) |>
      mutate(prob_hit = 1 - prob_miss,
             expected = n_bact * prob_hit,
             dens = density(dist_binomial(n_bact, prob_hit), 0:100)) |>
      tidyr::unnest_longer(col = dens, indices_include = TRUE)
  })
  
  # for any/any
  
  prob_x_hit_free_both <- reactive({
    crossing(n_bact = seq(10, 100, 10), n_virus = n_bac() * seq(10, 100, 10)) |> 
      mutate(prob_miss = prob_one_miss()^n_virus) |>
      mutate(prob_hit = 1 - prob_miss,
           expected = n_bact * prob_hit)
  })
  
  print_as_percent <- function(number, total, digits = 0){
    return(sprintf(paste0("%1.", digits, "f%%"), 100*as.numeric(number)/as.numeric(total)))
  }
  
  # make the plot
  
  plot <- reactive({
    if(input$any_virus & !input$any_bacteria) {
      ggplot(data = prob_x_hit_free_v()) +
        geom_ridgeline(aes(y = as.factor(n_virus), x = dens_id - 1, height = dens),
                       scale = 5, fill = "#346eeb", alpha = 0.75) +
        geom_text(aes(y = as.factor(n_virus), x = expected, label = round(expected, 1)), colour = "white", vjust = -.5) + 
        theme_minimal() +
        xlim(0, 100) + 
        labs(y = "Number of Viral Particles", x = "Number of Colonies in Range of at Least One Viral Particle")
    }else if (!input$any_virus & input$any_bacteria) {
      ggplot(data = prob_x_hit_free_b()) +
        geom_ridgeline(aes(y = as.factor(n_bact), x = dens_id - 1, height = dens),
                       scale = 5, fill = "#346eeb", alpha = 0.75) +
        geom_text(aes(y = as.factor(n_bact), x = expected, label = round(expected, 1)), colour = "white", vjust = -.5) + 
        theme_minimal() +
        xlim(0, 100) +
        labs(x = "Number of Colonies in Range of at Least One Viral Particle", y = "Number of Colonies")
    }else if (input$any_virus & input$any_bacteria) {
      ggplot(data = prob_x_hit_free_both()) + 
        geom_tile(aes(x = n_virus, y = n_bact, fill = expected), alpha = 0.8) + 
        geom_text(aes(x = n_virus, y = n_bact, label = round(expected, 1))) +  
        theme_minimal() + 
        scale_fill_viridis() + 
        labs(x = "Viral Particles", y = "Number of Colonies")
    }else{
      set.seed(123)
      plate_plots <- map(c("Example 1", "Example 2", "Example 3"), 
                         .f = \(x){draw_plate(input$num_bacteria, 
                                input$num_virus_per_colony * input$num_bacteria, 
                                input$size_plate, 
                                input$virus_within, x)}) 
      
      prob_plot <- ggplot(data = prob_table()) + 
        geom_col(aes(x = as.factor(within), y = prob), fill = "#346eeb", colour = "white", alpha = .75) + 
        geom_vline(xintercept = 1 + expected_hits()) + 
        geom_text(label = expected_hits() |> round(2),
                x = 1 + expected_hits(),
                y = max(pull(prob_table(), prob))) + 
        theme_minimal() + 
        scale_x_discrete(labels = \(x){paste0(x, "\n(", print_as_percent(x, input$num_bacteria), ")")}) + 
        scale_y_continuous(labels = scales::percent) + 
        labs(x = "Number (Percentage) of Colonies in Range of a Viral Particle", y = "Probability of Outcome")
      
      prob_plot/(plate_plots[[1]] + plate_plots[[2]] + plate_plots[[3]])
    }
  })
  
  output$distPlot <- renderPlot({
      plot()
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
