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
import seaborn as sns
import sys

class SIRSModel(object):

	# By using the constructor, initialize a 2D grid of N rows
	# with N columns. The cells are randomly initialized to 1 
	# or 0 with equal a priori probability	
	# 1 = Infected neighbor
	# 0 = Susceptible neighbor
	def __init__(self, N):		
		self.grid = np.random.randint(2, size=(N,N))
		
	### Notations:
	### p1 -> probability of getting infected
	### p2 -> probability of getting recovered
	### p3 -> probability of getting susceptible again
	###
	
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
	
	###
	### PART A - IMPLEMENT THIS FUNCTION FOR GETTING AN ANIMATION OF THE DISEASE 
	### EVOLUTION IN THE COMMUNITY 
	###
	# Go again through the community obeying the SIRS Model. However, do it this time 
	# indefinitely and plot everything with animations included. Patterns have been found
	# for absorbing state, dynamic equilibrium and cyclic wave condition			
	def UpdateCommunity(self, N, p1, p2, p3):
		counter = 0
		i = 0
		while(counter < 1):
			# Update the whole community
			self.UpdateGrid(N, p1, p2, p3)
			plt.cla()
			image = plt.imshow(self.grid, animated=True, interpolation='none')
			if (i == 0):
				plt.colorbar()
			plt.draw()
			plt.title("Disease spread in the community")
			plt.pause(0.001) 
			i += 1
			
	### 
	### PROBABILITY PARAMETERS FOR (N, p1, p2, p3)
	### ABSORBING STATE: (50, 0.5, 0.6, 0.19)
	### DYNAMIC EQUILIBRIUM: (50, 0.5, 0.5, 0.5)
	### CYCLIC EQUILIBRIUM: (100, 0.8, 0.1, 0.01)
	###
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
				return np.zeros((1,10000))
				
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
				# of 1000 sweeps. Do this by setting the initial counter variable to 1.
				# NOTE: USE 1000 SWEEPS FOR MEAN AND VARIANCE CONTOURS AND USE 10000 SWEEPS
				# FOR THE ERROR BARS ON THE VARIANCE
				if(counter_infection == 10000):
					print(len(I))
					return np.array(I)
					counter = 1

###
### PART B - IMPLEMENT THIS METHOD FOR PART B
### 
					
def PlotContours():

	# Read the values from the terminal first of all
	N = int(sys.argv[1])
	p1 = float(sys.argv[2])
	p2 = float(sys.argv[3])
	p3 = float(sys.argv[4])
	
	# Define the resolution of the grid and the number of the
	# cells per one line - hence the grid dimensionality. After,
	# create the numpy array
	r = 0.05
	n = int(1/r) + 1
	isotherms = np.zeros((n,n))
	
	# Create a list for appending the grid of infected people at 
	# each step. This grid will be saved in the end and stored in 
	# a .csv file for further analysis
	InfectionAnalysis = []
	Probabilities = []
	Mean_I = []
		
	# Set the p2 parameter to a constant - 0.5 is chosen. Afterwards,
	# work through the p1-p3 plane with the given resolution in order
	# to work out the behaviour of the average infection rate
	p2 = 0.5
	for i in range(n):
		for j in range(n):
		
			# Initialize the community obeying the SIRS Model at each step
			community = SIRSModel(N)
			
			# Initialize the probabilities
			p1 = i * r
			p3 = j * r
			p = np.array([p1, p3])
			
			# Return the array of infected people at each
			infected_people = community.ShowTheInfected(N, p1, p2, p3)
			
			# The grid element takes the value of the variance in the infection rate
			# divided by the community population. Print the value in the terminal.
			# NOTE: use np.mean() or np.var() depending on the required contour plot
			isotherms[i][j] = np.mean(infected_people) / (N**2)
			print("Variance on the infection rate: " + str(isotherms[i][j]))
			
			# Update the array required for storing the infection rates from the site
			InfectionAnalysis.append(infected_people / N**2)
			Probabilities.append(p)
			Mean_I.append(isotherms[i][j])
	
	# Save the data with the infected people in a .csv file using the pandas library
	df = pd.DataFrame({"Infections " : InfectionAnalysis, "Probabilities p1-p3 " : Probabilities, "Mean" : Mean_I})
	df.to_csv("InfectedPeopleMean_copy.csv")
				
	# Plot the resulted isotherms
	plt.imshow(isotherms, extent=[0, 1, 0, 1], cmap='hot', origin='lower')
	plt.colorbar()
	plt.title("Fluctuations on the rate of infection")
	plt.show()
	
def BootstrapSampling(arr, N):

	# Bootstrap method for heat capacity in the Glauber Dynamics model
	# Declare firstly the number of samplings we would like to take.
	# Also declare the initial empty array of k samplings which will
	# be updated at each step
	# Declare the number of iterations for the bootstrap sampling
	k = 100

	# Declare the number of elements n the user wants to use for
	# each iteration
	n = len(arr)

	# Declare the number of particles
	no_particles = N**2

	# Initialize empty list of k length
	var_I_list = []

	# Iterate for k times in order to retrieve the proper bootstrap sampling
	for _ in range(k):
		# Now select n random numbers from the array. Do this using
		# the numpy random function generator
		testArr = rng.choice(arr, size=n)
			
		# Calculate now the heat capacity from the data on internal energy
		var_I_list.append(np.var(testArr) / (no_particles))
	
	# Now calculate the error for the heat capacity from the given data on the
	# retrieved energies	
	return np.std(var_I_list)

### 
### PART C - DO THIS FOR PART C
###
	
def Do_10000_Sweeps():

	# Read the values from the terminal first of all
	N = int(sys.argv[1])
	
	I = []
	p2 = 0.5
	p3 = 0.5
	P1 = np.linspace(0.2, 0.5, 31)
	Y = np.zeros(31)
	yErr = np.zeros(31)
	
	for p1 in P1:
		
		community = SIRSModel(N)
		my_arr = community.ShowTheInfected(N, p1, p2, p3)
		I.append(my_arr)
		Y[np.argwhere(P1==p1)] = np.var(my_arr) / N**2
		yerr = BootstrapSampling(my_arr, N)
		yErr[np.argwhere(P1==p1)] = yerr
		
		print("p1: " + str(p1))
		print("variance: " + str(Y[np.argwhere(P1==p1)]))
		print("variance error: " + str(yerr))
	
	df = pd.DataFrame({"Infections " : I, "Probabilities p1 " : P1, "Variances" : Y, "Error on variance" : yErr})
	df.to_csv("InfectedPeople_10000sweeps_copy.csv")	
	
	plt.errorbar(P1, Y, yerr=yErr, fmt='o')
	plt.show()	
						
def main():

	# Read the values from the terminal first of all
	N = int(sys.argv[1])
	p1 = float(sys.argv[2])
	p2 = float(sys.argv[3])
	p3 = float(sys.argv[4])
	
	ok = int(input("Introduce the following command: " + "\n" + "(0) Plot the SIRS Model on (2,2) Grid" + "\n" 
		+ "(1) Contour plot for mean and variance values: " + "\n" + "(2) Error bars on variance for p2 = p3 = 0.5 " + "\n"
		+ "Place here your command: "))
	
	# Take each case
	if(ok == 0):
		
		community = SIRSModel(N)
		community.UpdateCommunity(N, p1, p2, p3)
	
	elif(ok == 1):
		PlotContours()
		
	elif(ok == 2):
		Do_10000_Sweeps()
		
	else:
		raise Exception("Introduce an integer between 0 and 2!!")
		
main()			
