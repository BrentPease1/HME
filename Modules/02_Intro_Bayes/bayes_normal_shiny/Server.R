
 server<-function(input, output, session) {

   output$Bayes_Plot <- renderPlot({
 
     mu <- seq(-25,100,0.01)
 
     prior <- dnorm(mu,input$theta,input$tau)
     prior <- prior/sum(prior)

     VVV   <- input$n/input$sigma^2 + 1/input$tau^2
     MMM   <- input$n*input$xbar/input$sigma^2 +
              input$theta/input$tau^2 
     posterior <- dnorm(mu,MMM/VVV,1/sqrt(VVV))
     posterior <- posterior/sum(posterior)

     mx <- max(c(max(prior),max(posterior)))

     plot(mu,prior,
          type="l",lwd=2,ylim=c(0,mx),
          xlab=expression(mu),ylab="Density",cex.lab=1.5)
     abline(v=input$xbar,col=2,lwd=2)
     lines(mu,posterior,col=4,lwd=2)

     legend("topright",c("Prior","Posterior","Sample mean"),
            lwd=2,col=c(1,4,2),inset=0.0,cex=1.2,bty="n")
  })

 }
