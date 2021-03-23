output_plot <-
function(output, initial){
            if(missing(initial)) initial <- numeric(nrow(output))
            da <- format(Sys.time(), "%b %d %H:%M ")
            plot.name1   <- sprintf("Plots MCEM %s.pdf",da)
            th <- initial
            
            pdf(plot.name1, width=6, height=4)
            plot.ts(output[,1], type="l", xlab="Iterations", ylim=c(-1, 1), ylab="phi_1",  main=""); abline(h=th[1], lty=2, col="CadetBlue"); 
            plot.ts(output[,2], type="l", xlab="Iterations", ylim=c(-1, 1), ylab="psi_1",  main=""); abline(h=th[2], lty=2, col="CadetBlue"); 
            dev.off()
}
