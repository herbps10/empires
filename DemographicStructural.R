# install.packages("deSolve")

## Load in the DiffEq solver

library(deSolve)

#  First, the parameters for our model:

S_not = 10
c_param = 3
params <- c(
  r = 0.02
  beta = 0.25
)
			
#	Then, the initial values of each state variable:

init.values <- c(N = 0.1, S = 0.1)

#	The times we want to see
times <- seq(0, 2000, by = 1)
#	Now we can define the differential equation model:

k <- function(S) {
  result = 1 + c_param * (S / (S_not + S))
  
  return(result)
}

Model <- function(time, y.values, parameters) {  # Don't change this
  if(y.values[2] < 0) y.values[2] = 0
  
	with(as.list(c(y.values, parameters)), {  # or this

    
    dN.dt = r * N * (1 - (N/k(S))
    dS.dt = N * (1 - (N/k(S))) - beta * N
  
    N = N + dN.dt
    S = S + dS.dt
  
    S = S + rnorm(mean = 0, sd= 0.1 * S_not, n=1)
  
    if(N < 0) N = 0
    if(S < 0) S = 0
        
    #if(S + dS.dt < 0) dS.dt = S
		return(list(c(N, S)))		 
		})
	}
	
#	Having defined everything, now we ask the program ode
#	to actually solve the system:

Graph <- function() {

	out <- as.data.frame(ode(func = Model, y = init.values, 
							parms = params, times = times, method="iteration"))
							
	#	And now we can make a nice plot of our results:

	matplot(out$time, out[ ,2:3], type = "l", xlab = "Time", 
		ylab = "Percent of Population", main = "Demographic/Structural Model", lwd = 2)
	
	legend("right", c("Population", "State Power"),
		col = 1:3, lty = 1:2)
}
   
Graph()