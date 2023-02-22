
 ui <- fluidPage(
   pageWithSidebar(

   # Application title
   headerPanel("Bayesian analysis of a normal mean - peregrine falcons"),

   sidebarPanel(
    sliderInput("n", 
                "Sample size (n)", 
                min = 1,
                max = 100, 
                value = 12,
                step = 1),
    sliderInput("xbar", 
                "Sample mean (x-bar)", 
                min = 500,
                max = 700, 
                value = 607,
                step = 1),
    sliderInput("sigma", 
                "Population sd (sigma)", 
                min = 0.01,
                max = 50, 
                value = 27,
                step = 1),
    sliderInput("theta", 
                "Prior mean: theta", 
                min = 0,
                max = 1000, 
                value = 607,
                step = 1),
    sliderInput("tau", 
                "Prior sd: tau", 
                min = 0.01,
                max = 50, 
                value = 25,
                step = 1)
   ),

   # Show a plot of the generated land use classifications
   mainPanel(
     plotOutput("Bayes_Plot")
   )
  )
 )
