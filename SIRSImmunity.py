# Import the numpy and matplotlib libraries for creating 2D-array
# and for making animations of the time-dependent lattice. Also
# use pandas library for storing the avaialble data into a .csv 
# file and seaborn as additional plotting library. Sys library is
# imported to implement values from the terminal execution line 
import numpy as np
rng = np.random.default_rng()
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
from matplotlib.colors import ListedColormap
import seaborn as sns
import sys

class SIRSImmunity(object):

	# By using the constructor, initialize a 2D grid of N rows
	# with N columns. The cells are randomly initialized to 1 
	# or 0 with equal a priori probability	
	# 1 = Infected neighbor
	# 0 = Susceptible neighbor
	def __init__(self, N, f):
		p0 = (1-f)/3
		self.grid = np.random.choice(np.arange(-1,3), size=(N,N), p = [f, *np.repeat(p0,repeats=3)])
	
	# Create function to check if there are any infected neighbors. 
	# Only the nieghbors on vertical and horizontal axis are checked,
	# the diagonals are not of interest
	def InfectedNeighbors(self, N, i, j):
		if((self.grid[i][(j+1)%N] == 1) or (self.grid[i][(j-1)%N] == 1) or (self.grid[(i-1)%N][j] == 1) or (self.grid[(i+1)%N][j] 			== 1)):
			return 1
			
		else:
			return 0
	
	# Update the cell as function of whether it is susceptible, infected
	# or recovered		
	def UpdateCell(self, N, i, j, p1, p2, p3):
	
		# Take the case where the cell is susceptible. Check the neighbors
		# and perform the MC Algorithm in order to yield an infected cell
		if(self.grid[i][j] == 0):
			if(self.InfectedNeighbors(N, i, j) == 1):
				r = np.random.random()
				if (r < p1):
					self.grid[i][j] = 1
				else:
					self.grid[i][j] = 0
			else:
				self.grid[i][j] = 0
					
		# Take the case where the cell is infected. Perform the MC algorithm
		# in order to yield a recovered cell
		elif(self.grid[i][j] == 1):
			r = np.random.random()
			if (r < p2):
				self.grid[i][j] = 2
			else:
				self.grid[i][j] = 1
				
		# Take the case of recovered cell. Perform the MC algorithm in order
		# to retrieve a susceptible cell and repeat the cycle
		elif(self.grid[i][j] == 2):
			r = np.random.random()
			if(r < p3):
				self.grid[i][j] = 0
			else:
				self.grid[i][j] = 2
	
	
	# Go through the community of suscepitble, infected and recovered patients and
	# perform the overall Monte Carlo Simulation. The simulation is done only once!!
	def UpdateGrid(self, N, p1, p2, p3):
		count = 0
		while(count < N**2):
			i = np.random.randint(N)
			j = np.random.randint(N)
			self.UpdateCell(N, i, j, p1, p2, p3)
			count += 1
			
	# Go again through the community obeying the SIRS Model. However, do it this time 
	# indefinitely and plot everything with animations included. Patterns have been found
	# for absorbing state, dynamic equilibrium and cyclic wave condition			
	def UpdateCommunity(self, N, p1, p2, p3):
		counter = 0
		i = 0
		cMap = ListedColormap(["white", "black", "red", "blue"])
		while(counter < 1):
			# Update the whole community
			self.UpdateGrid(N, p1, p2, p3)
			plt.cla()
			image = plt.imshow(self.grid, animated=True, interpolation='none', cmap=cMap, vmin=-1, vmax=2)
			plt.draw()
			plt.title("Disease spread in the community")
			plt.pause(0.001)
			if(i == 0):
				plt.colorbar()
			i += 1
			
	# Go through the whole community as follows: Do 100 sweeps (complete runs over the N*N
	# cells) for achieving equilibrium. Afterwards, do another 1000 sweeps in order to retrieve
	# the average infection rate. Return the mean value of the infection rate for the p2 value of
	# 0.5 and fixed p1 and p3		
	def ShowTheInfected(self, N, p1, p2, p3):
		p2 = 0.5
		
		# Variable used for the loop simulation
		counter = 0
		
		# Variable used to check if we are in the first 100 sweeps during the community
		# initialization
		counter_initial = 0
		
		# Variable used to count the number of infected people
		counter_infection = 0
		
		# List for the number of infected people after each sweep
		I = []
		
		# Now run through the loop until the limit is reached (or until the user wants 
		# to apply the Ctrl+C command)
		while(counter < 1):
			
			# Set the condition that for no infection rate pass directly to the 
			# next step to avoid long simulations
			if(len(self.grid[self.grid == 1]) == 0):
				return np.zeros((1,1000))
				
			# The counter_initial variable is used to measure the moment we pass
			# the intial equilibrium dynamics given by the first 100 sweeps
			counter_initial += 1
			
			# Update the grid with the susceptible, infected and recovered patients
			self.UpdateGrid(N, p1, p2, p3)
			
			# Now analyse what happens when we reach past the first 100 equiblibrium 
			# sweeps - when we actually start storing the infection rates
			if(counter_initial > 100):
			
				# Variable used to check if we are in the limit of the 1000 sweeps
				counter_infection += 1
				
				# Take the number of infected people in the community and append the
				# number to the list of infected people in each simulation
				infected = len(self.grid[self.grid == 1]) 
				I.append(infected)
				
				# Get out of the loop and finish the exercise when we are past the point
				# of 1000 sweeps. Do this by setting the initial counter variable to 1
				if(counter_infection == 1000):
					print(len(I))
					return np.array(I)
					counter = 1

###
### PART A: GENERATE A COMMUNITY OBEYING THE SIRS MODEL
### WITH A FRACTION OF POPULATION IMMUNE ALL THE TIME
###
			
def GenerateImmuneCommunity(N, p1, p2, p3):
	community = SIRSImmunity(N, 0.1)
	community.UpdateCommunity(N, p1, p2, p3)
	
###
### PART B: WORK OUT THE REQUIRED FRACTION OF THE COMMUNITY TO GET IMMUNE
### IN ORDER TO STOP THE DISEASE SPREAD FOR p1 = p2 = p3 = 0.5
###

def FindTheRequiredImmunity(N, p1, p2, p3):
	
	# Create a numpy array with the fraction required for
	# immunity - an upper limit of 0.4 suffices
	F = np.linspace(0, 0.4, num=10)
	
	# Create now lists for the average infection rate and the errors
	Y = []
	errs = []
	
	# Go through each fraction and work out, for p1 = p2 = p3 = 0.5
	p1 = 0.5
	p2 = 0.5
	p3 = 0.5
	
	# The average infection rate
	for f in F:
		
		# For each fraction, five measurements are condcuted, yielding five
		# different mean values. The errors in measurements are taken as the
		# standard of the means values
		i = 0
		arr = []
		while(i < 5):
			i += 1
			community = SIRSImmunity(N, f)
			mean_val = np.mean(community.ShowTheInfected(N, p1, p2, p3)) / N**2
			print("Simulation: " + str(i) + "fraction f: " + str(f) + " have proportion I: " + str(mean_val))
			arr.append(mean_val)
		print("Error: " + str(np.var(arr)))
		
		# Append now the values to the lists
		Y.append(np.mean(arr))
		errs.append(np.std(arr))
		
	df = pd.DataFrame({"Immunity fraction " : F, "Infection rate" : Y, "Errors on infection rate" : errs})
	df.to_csv("Immunity.csv")
		
def main():
	
	N = int(sys.argv[1])
	p1 = float(sys.argv[2])
	p2 = float(sys.argv[3])
	p3 = float(sys.argv[4])

	ok = int(input("Introduce the following command: " + "\n" + "(0) Plot the SIRS Model on (2,2) Grid" + "\n" 
		+ "(1) Average infection rate I as function of the immunity fraction: \n Place your command here: "))
		
	if(ok == 0):
		GenerateImmuneCommunity(N, p1, p2, p3)
		
	elif(ok == 1):
		FindTheRequiredImmunity(N, p1, p2, p3)
		
	else:
		raise Exception("Introduce 0 or 1!!")
main()
