# install.packages("deSolve")

## Load in the DiffEq solver

library(deSolve)

#	First, the parameters for our model:

params <- c(
  r = 0.1,     # growth rate
  K = 1,       # carrying capacity
  a = 0.1,     # hostility initiation rate
  b = 0.03,     # hostility fall off rate
  c = 0.1,     # impact war has on carrying capacity
  alpha = 0.1, # state power's effect on internal war
  rho = 0.2,   # taxation rate
  beta = 0.1,  # expenditure rate
  delta = 0.3  # effect of war on population
)
			
#	Then, the initial values of each state variable:

init.values <- c(N = 0.1, S = 0.0, W = 0.0)

#	The times we want to see
times <- seq(0, 500, by = 1)
#	Now we can define the differential equation model:

War <- function(time, y.values, parameters) {  # Don't change this
	with(as.list(c(y.values, parameters)), {  # or this
		if(N < 0) N = 0
    if(S < 0) S = 0
    if(W < 0) W = 0
    
    dN.dt = r * N * (1 - N/(K - c * W)) - delta * W * N
    dS.dt = rho * N * (1 - N/(K - c * W)) - beta * N
		dW.dt = (a * (N ^ 2)) - b * W - alpha * S
        
		return(list(c(dN.dt, dS.dt, dW.dt)))		 
		})
	}
	
#	Having defined everything, now we ask the program ode
#	to actually solve the system:

Graph <- function() {

	out <- as.data.frame(ode(func = War, y = init.values, 
							parms = params, times = times))
							
	#	And now we can make a nice plot of our results:

	matplot(out$time, out[ ,2:4], type = "l", xlab = "Time", 
		ylab = "Percent of Population", main = "Population/War Model", lwd = 2, ylim=c(0, 1))
	
	legend("right", c("Population", "State Power", "Internal War Intensity"),
		col = 1:3, lty = 1:2)
}
   
Graph()