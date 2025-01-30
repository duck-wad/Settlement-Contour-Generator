import numpy as np 
import math
<<<<<<< HEAD
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.interpolate import griddata

# Meshing Parameters
MeshLong = 0.5          # Mesh length along tunnel axis
MeshTrans = 0.5         # Mesh length transverse to tunnel axis

# Hard-coded variables
R = 2.014/2             # Tunnel radius (m) 
z0 = 7.5                # Tunnel axis depth (m)
# For now, i is a hardcoded term but it can be calculated based on the tunnel depth
i = 3.9                 # Settlement trough width parameter (m)
wmax = 0.00786          # Maximum vertical settlement (m) 
Vs = 0.0124*math.pi*R*R
# What is n?
n = 9.8

# Function to compute the x,y,z ground movements due to tunnel excavation
# Assumption that the tunnel face has advanced significantly far from the start of the tunnel (xf > 3z0) so the transverse profile behind the face develops fully
# x-xi is assumed to be infinite
def GroundMovement():

    Mesh = DefineMesh(i)

    # Initialize output variables. u is x, v is y, w is z displacement
    u = np.array(len(Mesh))
    v = np.array(len(Mesh))
    w = np.array(len(Mesh))

    # Extract x and y columns
    x = Mesh[:, 0]  
    y = Mesh[:, 1]  

    # Compute w, v, u
    w = (Vs / np.sqrt(2 * np.pi * i)) * np.exp(-np.power(y, 2) / (2 * i * i)) * (1 - norm.cdf(x / i))
    v = -n / z0 * y * w
    u = (n * Vs / (2 * np.pi * z0)) * np.exp(-np.power(y, 2) / (2 * i * i)) * (-np.exp(-np.power(x, 2) / (2 * i * i)))

    PlotContour(x, y, w, u, v)

# Take the i parameter as input and discretize a rectangular mesh around the origin
def DefineMesh(i):

   # Define x and y ranges based on i
    x_range_neg = np.arange(0, -4*i, -MeshLong)  # Left side
    x_range_pos = np.arange(MeshLong, 4*i, MeshLong)  # Right side
    y_range_neg = np.arange(0, -3*i, -MeshTrans)  # Bottom
    y_range_pos = np.arange(MeshTrans, 3*i, MeshTrans)  # Top

    # Generate grid points for each quadrant
    x1, y1 = np.meshgrid(x_range_neg, y_range_neg)  # Lower-left
    x2, y2 = np.meshgrid(x_range_neg, y_range_pos)  # Upper-left
    x3, y3 = np.meshgrid(x_range_pos, y_range_neg)  # Lower-right
    x4, y4 = np.meshgrid(x_range_pos, y_range_pos)  # Upper-right

    # Flatten and stack to create coordinate list
    Mesh = np.vstack([
        np.column_stack((x1.ravel(), y1.ravel())),
        np.column_stack((x2.ravel(), y2.ravel())),
        np.column_stack((x3.ravel(), y3.ravel())),
        np.column_stack((x4.ravel(), y4.ravel()))
    ])

    return Mesh 

def PlotContour(x, y, w, u, v):
    # Create a grid for contouring
    xi = np.linspace(min(x), max(x), 100)
    yi = np.linspace(min(y), max(y), 100)
    X, Y = np.meshgrid(xi, yi)

    # Interpolate data to fit grid for w, u, and v
    Z_w = griddata((x, y), w, (X, Y), method='cubic')
    Z_u = griddata((x, y), u, (X, Y), method='cubic')
    Z_v = griddata((x, y), v, (X, Y), method='cubic')

    fig, ax = plt.subplots(1, 3, figsize=(18, 6))

    contour_w = ax[0].contourf(Y, X, Z_w, levels=20, cmap='viridis')
    fig.colorbar(contour_w, ax=ax[0], label="Settlement (w)")
    ax[0].set_xlabel("Y Coordinate")
    ax[0].set_ylabel("X Coordinate")
    ax[0].set_title("Contour Plot of Surface Settlement")
    ax[0].invert_yaxis()

    contour_u = ax[1].contourf(Y, X, Z_u, levels=20, cmap='viridis')
    fig.colorbar(contour_u, ax=ax[1], label="Displacement (u)")
    ax[1].set_xlabel("Y Coordinate")
    ax[1].set_ylabel("X Coordinate")
    ax[1].set_title("Contour Plot of X-Displacement")
    ax[1].invert_yaxis()

    contour_v = ax[2].contourf(Y, X, Z_v, levels=20, cmap='viridis')
    fig.colorbar(contour_v, ax=ax[2], label="Displacement (v)")
    ax[2].set_xlabel("Y Coordinate")
    ax[2].set_ylabel("X Coordinate")
    ax[2].set_title("Contour Plot of Y-Displacement")
    ax[2].invert_yaxis()

    plt.tight_layout()  
    plt.show()

=======

# Meshing Parameters
MeshLong = 1         # Mesh length along tunnel axis
MeshTrans = 1         # Mesh length transverse to tunnel axis

# Hard-coded variables
R = 4                   # Tunnel radius (m) 
z0 = 7.5                # Tunnel axis depth (m)
# For now, i is a hardcoded term but it can be calculated based on the tunnel depth
i = 1                # Settlement trough width parameter (m)
wmax = 0.00786          # Maximum vertical settlement (m) 
Vs = 0.025*math.pi*R*R

# Function to compute the x,y,z ground movements due to tunnel excavation
# Assumption that the tunnel face has advanced significantly far from the start of the tunnel (xf > 3z0) so the transverse profile behind the face develops fully
def GroundMovement():

    ## Generate Mesh
    Mesh = DefineMesh(i)

    # Initialize output variables. u is x, v is y, w is z displacement
    u = []
    v = []
    w = []

    for r in range(len(Mesh)):
        G = NormalDistribution(Mesh[r][0] / i)

        w_temp = Vs / (math.sqrt(2*math.pi*i)) * math.exp(-1*math.pow(Mesh[r][1], 2))
        





# Take the i parameter as input and discretize a rectangular mesh around the origin
# Returns an array of tuples, each tuple is a (x,y) surface coordinate
def DefineMesh(i):

    Mesh = []

    yPos = 0
    xPos = 0

    # Populate the upper left quadrant
    while yPos > -3*i:
        while xPos > -4*i:
            Mesh.append((xPos, yPos))
            xPos -= MeshLong
        xPos = 0    # Re-initialize x-pos to the x-axis
        yPos -= MeshTrans
    # Populate the upper right quadrant
    xPos = 0
    yPos = MeshTrans
    while yPos < 3*i:
        while xPos > -4*i:
            Mesh.append((xPos, yPos))
            xPos -= MeshLong
        xPos = 0
        yPos += MeshTrans
    # Populate the lower left quadrant
    xPos = MeshLong
    yPos = 0
    while yPos > -3*i:
        while xPos < 4*i:
            Mesh.append((xPos, yPos))
            xPos += MeshLong
        xPos = MeshLong
        yPos -= MeshTrans
    # Populate the lower right quadrant
    xPos = MeshLong
    yPos = MeshTrans
    while yPos < 3*i:
        while xPos < 4*i:
            Mesh.append((xPos, yPos))
            xPos += MeshLong
        xPos = MeshLong
        yPos += MeshTrans

    return Mesh

def NormalDistribution(x):
    # Define the G function
    return 1
>>>>>>> 764506c0f9d0a55064de3a0e3636c11ab5fae022

GroundMovement()