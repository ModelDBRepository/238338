"""

2DPoincareSimulation.py

Simulation of a 2D Poincaré oscillator network as a model of the Choroid Plexus

Reference: Myung J, Schmal C et al. (2018) The Choroid Plexus is an Important
Circadian Clock Component. Nat Commun, in press.

Written by Christoph Schmal, Institute for Theoretical Biology, Humboldt Universität
Correspondence: Christoph Schmal (christoph.schmal@charite.de)

February 9, 2018

"""


from numpy import *
from numpy.random import seed, normal
from scipy.integrate import odeint
from matplotlib.pyplot import figure, imshow, xticks, yticks, tight_layout, show

def phase_diff_on_circle(x, y):
     return arctan2(sin(x-y), cos(x-y))

def mean_circular_variable(x):
    return arctan2(sum(sin(x)), sum(cos(x)))

def make_adjacency_vNn(c_mask_res, image_height, image_width, grid_edge_length):
    y_grid = int(image_height/grid_edge_length)
    x_grid = int(image_width/grid_edge_length)
    coordinates_array = [[]]*y_grid*x_grid
    for i in xrange(y_grid):
        for j in xrange(x_grid):
            coordinates_array[i*x_grid + j] = [i, j]
    masked_coordinates = array(coordinates_array)[c_mask_res]
    adj_array = zeros(len(masked_coordinates)*len(masked_coordinates)).reshape(len(masked_coordinates), len(masked_coordinates))
    for i in xrange(len(masked_coordinates)):
        adj_array[i][argwhere((sum(abs(masked_coordinates-masked_coordinates[i]), axis=1)) <= 1)] = 1.
    return adj_array

# Properties of the original image
image_width       = 415
image_height      = 198
grid_element_size = 5

# Define an adjacency matrix, based on a nearest neighbor van-Neumann neighborhood and an experimentally observed geometry.
c_mask = genfromtxt("ChoroidPlexus_Mask.txt", dtype=bool)
N = sum(c_mask)
adj_array = make_adjacency_vNn(c_mask, image_height, image_width, grid_element_size)*(~identity(N, dtype=int)+2)

seed(9546) # Uncomment this line if you want to model an intrinsic period distribution different from the one used in our paper.

# Random sample intrinsic free-running periods from a normal distribution.
N = sum(c_mask)
per_mu = 25.5
per_std = 1.
per = normal(loc=per_mu, scale=per_std, size=N)

# Load initial conditions for the phases from the experimental data.
c_phases = loadtxt("Phase_InitialConditions.txt")
# Choose radial initial conditions as r=1 for all oscillators and convert the Polar coordinates into Cartesian coordinates.
X0 = array(list(cos(c_phases[:])) + list(sin(c_phases[:])) )

# Single Cell Oscillator Properties
gamma =  0.05
eps   = -0.01
A     = 1.
K     = 0.1
KArray = adj_array*K

# Define the right-hand side of the dynamical system. Note that we define a 2xN dimensional Dynamical System in numpy-array format.
def RHS(X, t):
    X1  = X[0:N]
    X2  = X[N:]
    dX1 = ( A - sqrt(X1**2 + X2**2) ) * ( gamma*X1 - eps*X2 ) - 2.*pi*X2/per + 0.5*sum(KArray*X1, axis=1)
    dX2 = ( A - sqrt(X1**2 + X2**2) ) * ( gamma*X2 + eps*X2 ) + 2.*pi*X1/per + 0.5*sum(KArray*X2, axis=1)
    return array([dX1, dX2]).flatten()

# Solve the system of differential equations numerically.
dt = 0.1
t = arange(0, 24*10, dt)
sol = odeint(RHS, X0, t)

# Easy calculation of phase and amplitude (without Hilbert transformation)
Phase = ( arctan2(sol.T[N:], sol.T[:N])     ).T
Amp   = ( sqrt( sol.T[:N]**2+sol.T[N:]**2 ) ).T

# We aim to plot the mean centered 2D phase-distribution after the decay of transient dynamics, here at t=24x6 h.
pic_num = int(24*6/dt)
PhaseDiffFromMean =  phase_diff_on_circle(Phase[pic_num], mean_circular_variable(Phase[pic_num]))

# Plotting commands that reproduce Figure S11 C of the Supplementary Text
template = array(c_mask[:], dtype=float)
template[argwhere(c_mask == 0)] = inf
template[argwhere(c_mask != 0)] = PhaseDiffFromMean.reshape(len(PhaseDiffFromMean), 1)

fig = figure(figsize=(6, 2.5))
imgplot = imshow(template.reshape(image_height/grid_element_size, image_width/grid_element_size), interpolation="nearest", vmin=-pi/2, vmax=pi/2)
imgplot.set_cmap('hsv')
xticks([], [])
yticks([], [])
tight_layout()
show()