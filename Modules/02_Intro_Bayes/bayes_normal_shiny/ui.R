
 ui <- fluidPage(
   pageWithSidebar(

   # Application title
   headerPanel("Bayesian analysis of a normal mean - eucalypt trees"),

   sidebarPanel(
    sliderInput("n", 
                "Sample size (n)", 
                min = 1,
                max = 100, 
                value = 10,
                step = 5),
    sliderInput("xbar", 
                "Sample mean (x-bar)", 
                min = 35,
                max = 70, 
                value = 50,
                step = 1),
    sliderInput("sigma", 
                "Population sd (sigma)", 
                min = 0.01,
                max = 10, 
                value = 5,
                step = .01),
    sliderInput("theta", 
                "Prior mean: theta", 
                min = 0,
                max = 100, 
                value = 50,
                step = 1),
    sliderInput("tau", 
                "Prior sd: tau", 
                min = 0.01,
                max = 10, 
                value = 10,
                step = .01)
   ),

   # Show a plot of the generated land use classifications
   mainPanel(
     plotOutput("Bayes_Plot")
   )
  )
 )
