#
# Model of dynamic civilizations over time.
#
# This is a discrete time, spatial model. The model is run on a square grid of cells. 
# Empires are made up of cells of the grid. A border cell is a cell that is on the border of an empire.
# 
# An important concept in this model is the idea of "asabiya", or the capacity of a given group for collective action.
# Turchin theorizes that asabiya is highest on imperial borders. Towards the center of an empire, asabiya decreases.
#
# The power of a given cell in an empire is related to both the average asabiya of the whole civilization
# and the distance to the capitol. Thus, while a border cell might have high asabiya, if it is far from the
# capitol its power will decrease.
#
# Cells can invade neighboring cells if they have more power.
#
# Based off of a model described by Peter Turchin in HISTORICAL DYNAMICS.
#


library(igraph)
library(animation)
library(matrix)

###
### Constants
###

len = 21          # Dimension of the grid
deltaPower = 0.05 # Difference in power necessary to invade
r0 = 0.1          # How fast asibiya grows at border
delta = 0.1       # How fast asibiya decreases in heartland
h = 0.8           # How fast power falls off as a relation of distance (1: no loss, 0: all lost)
minAsa = 0.05     # Lowest asabiya allowed before a civilization is dissolved
initialAsa = 0.5  # Initial asabiya value for cells in new empire
numInitialEmpires = 10 # Number of empires to start with
timesteps = 100    # How long to run the model

###
### Globals
###

highestEmpire = 0
capitols = numeric(0)
frontierCells = numeric(0)
empireSize <- array(0, c(timesteps,6))

###
### FUNCTIONS
###

#
# Calculates a cell's id given x, y coordinates
#
getID <- function(x, y) {
  return(y * len + x)
}

  
#
# Check if a cell is on an imperial border
#
isBorder <- function(id, g) {
  neighbors <- neighborhood(grid, 1, nodes = id)[[1]]
  
  return(sum(V(g)[id]$empire == V(g)[neighbors]$empire) != length(neighbors))
}

  
#
# Returns the smallest number of vertices between a vertex and the capitol, including the capitol and the vertex itself.
#
distToCap <- function(id, g) {   
  
  if(V(g)[id]$empire==0) {
    return(0)
  }
  
  # Look up this empire's capitol
  capitol = V(g)[capitols[V(g)[id]$empire]]
  
  # Get the shortest path from this cell to its capitol
  return(length(get.shortest.paths(g, from = id, to = capitol)[[1]]))
    
}

    
#
# Calculates the average asabiya of a given empire
#
avgAsa <- function(empireID, g) {
  sum(V(g)[V(g)$empire==empireID]$asa) / sum(V(g)$empire == empireID)
}

  
#
# Get the number of cells in an empire
#
getSize <- function(empireID, g) {
   return(length(V(g)[V(g)$empire == empireID]))
}

  
#
# Calculates a given cell's power
# 
# Power is related to the average asabiya of the empire and how far
# away the cell is from the capitol. Power falls off exponentially
# as you move farther away from the capitol.
#
# Power = (Empire Area) * (Average Asabiya) * 1/e^((Distance to Capitol) * (Constant value h))
#
# The constant h is a measure of how fast power falls off over distance
#
getPower <- function(id, g) {
  empire = V(g)[id]$empire
  
  # If the cell is in the hinterland, it has no power
  if(empire == 0) {
    return(0)
  }
  
  area = getSize(V(g)[id]$empire, g)
  dist = distToCap(id, g)
  
  power = area * avgAsa(empire, g) * exp(-1*dist/h)
    
  return(power)
}

  
#
# Randomly place an empire on the grid
# 
# Size is the number of cells in each direction to
# add to the empire
#
randomPlaceEmpire <- function(size, g) {
  # Pick a point for the capitol
  # All vertices within two jumps of the capitol are part of the empire
  capitol = 0
  repeat {
    id = sample(0:21^2, 1)
    if(V(g)[id]$empire == 0) {
      break
    }
  }
  
  return(placeEmpire(id, size, g))
}

    
#
# Places an empire with a capitol at the cell 'id'
# 
# Adds 'size' cells to the empire in all directions
# from the capitol
#
placeEmpire <- function(id, size, g) {
  highestEmpire <<- highestEmpire + 1
  
  V(g)[neighborhood(g, size, id)[[1]]]$empire = highestEmpire
  V(g)[neighborhood(g, size, id)[[1]]]$asa = initialAsa
  
  V(g)[id]$capitol = TRUE
  capitols[highestEmpire] <<- id
  
  return(g)
}

###
### Initialization
###
  
grid <- graph.lattice(c(len,len))
V(grid)$empire = 0
V(grid)$asa = 0

# Stat with two randomly placed empires
for(i in 1:numInitialEmpires) {
  grid = randomPlaceEmpire(1, grid)
  grid = randomPlaceEmpire(1, grid)
}

  
# Keep track of size and asabiya data for civilizations
sizeData = matrix(0, nrow = timesteps, ncol = 100)
asaData = matrix(0, nrow = timesteps, ncol = 100)


  
###
### Main Loop
###
  
#
# 1. Update Asabiya for each cell that is a member of an empire.
# 2. For each empire cell on a border, attempt to invade neighbors.
# 3. Dissolve any empires with average asabiya less than a threshold
#
  
#ani.start()

for(timestep in 0:timesteps) {
  
  for(empire in unique(V(grid)$empire)) {
    sizeData[timestep,empire] = getSize(empire, grid)
    asaData[timestep, empire] = avgAsa(empire, grid)
  }
  
  frontierCells = numeric(0)
  
  # Loops through and updates asibaya value for each cell that is part of an empire
  for(i in V(grid)[V(grid)$empire != 0]) {
    currAsa = V(grid)[i]$asa
    
    # Asibaya is increased for border cells
    if (isBorder(i, grid)) {
      currAsa = currAsa + r0 * currAsa * (1 - currAsa) 
      
      # Keep track of border cells for easy access later
      frontierCells = union(frontierCells, c(i))
    }
    # For non border vertices, the asabiya decreases
    else {
      currAsa = currAsa - delta * currAsa
    }
    
    V(grid)[i]$asa = currAsa
    
  }
  
  #
  # Invasion Routine
  #
  # Loop through border cells in random order. 
  # For each border cell, check neighbors in a random order. Invade the first one that has a lower power value.
  #
    
  for(i in sample(frontierCells, length(frontierCells))) {
    neighbors = neighborhood(grid, 1, i)[[1]]
    
    # Randomly shuffle so we don't bias our results by going through the neighbors in the same order every time
    neighbors = sample(neighbors, length(neighbors))
    
    for (j in neighbors) {
      
      # Invade if:
      # The other cell is in a different empire
      # We have more power then them, over the deltaPower threshold
      
      if (V(grid)[i]$empire != V(grid)[j]$empire && getPower(i, grid)-getPower(j, grid) > deltaPower) {
        V(grid)[j]$empire = V(grid)[i]$empire
        break
      }
        
    }
  }      
  
  #
  # If the average asabiya of any empire is less than minAsa,
  # dissolve the empire and return all cells to the hinterland
  # 
    
  for(e in unique(V(grid)$empire)) {
    if(e == 0) {
      next
    }
    
    if(avgAsa(e, grid) < minAsa) {
      V(grid)[V(grid)$empire == e]$empire = 0
      V(grid)[V(grid)$empire == e]$asa = 0.5
    }
  }
  
  #
  # Each timestep, have a random chance that new empire will start.
  #
  if(sample(0:100, 1) < 5) {
    #grid = randomPlaceEmpire(1, grid)
  }
  
  
  # Assign colors to each empire
  # We only have six colors, so if we have more empires than that
  # just loop around back to the beginning
  V(grid)[V(grid)$empire %% 6 == 0]$color = "red"
  V(grid)[V(grid)$empire %% 6 == 1]$color = "blue"
  V(grid)[V(grid)$empire %% 6 == 2]$color = "green"
  V(grid)[V(grid)$empire %% 6 == 3]$color = "purple"
  V(grid)[V(grid)$empire %% 6 == 4]$color = "black"
  V(grid)[V(grid)$empire %% 6 == 5]$color = "pink"
  V(grid)[V(grid)$empire==0]$color = NA

  #empireSize[timestep,1] = length(V(grid)[V(grid)$empire %% 6 == 0])
  #empireSize[timestep,2] = length(V(grid)[V(grid)$empire %% 6 == 1])
  #empireSize[timestep,3] = length(V(grid)[V(grid)$empire %% 6 == 2])
  #empireSize[timestep,4] = length(V(grid)[V(grid)$empire %% 6 == 3])
  #empireSize[timestep,5] = length(V(grid)[V(grid)$empire %% 6 == 4])
  #empireSize[timestep,6] = length(V(grid)[V(grid)$empire %% 6 == 5])
  
  #
  # Manually place each vertex for plotting
  #
  # Nick -- you'll probably be interested in this code.
  # You can pass the plot function a custom layout which is a
  # N x 2 matrix where N is the number of vertices in the graph.
  # Each vertex has a row in the matrix, with x and y coordinates.
  #
  # The code below generates such a matrix and places each vertex
  # on a grid layout.
  #
  #coords = matrix(nrow=len^2, ncol=2)
  #for(i in 1:len^2) {
  #  coords[i,] = c((i-1) %% len, floor((i-1)/len))
  #}
  
  #plot(grid, vertex.size=10, layout=coords, main = paste(avgAsa(1, grid), avgAsa(2, grid), avgAsa(3, grid)), vertex.label="", vertex.frame.color="gray", edge.color="white", vertex.shape="square")
  
  print(timestep)
}

#ani.stop()

empire_lengths = c()

colors = c("blue", "red", "green", "purple", "black")
plot(sizeData[,1], col = "blue", type="l", ylim=c(0, 30), main="Empire Size vs. Time", ylab="Empire Size", xlab="Time")
for(i in 2:100) {
  if(sum(sizeData[,i]) > 0) {
    #lines(sizeData[,i], ylim=c(0, 30), col=colors[i %% length(colors)])
    lines(sizeData[,i], ylim=c(0, 30), col="black")
    
    empire_lengths = cbind(empire_lengths, sum(sizeData[,i] != 0))
  }
}

plot(asaData[,1], type='l', ylim=c(0, 1), main="Empire Assabiya vs. Time", xlab="Time", ylab="Assabiya")
for(i in 2:100) {
  if(sum(asaData[,i]) > 0) {
    lines(asaData[,i], ylim=c(0, 1))
  }
}