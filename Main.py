import numpy as np 
import math
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
# Vs = 0.0124*math.pi*R*R
Vs = 0.0391113074
# What is n?
n = 9.867

# Function to compute the x,y,z ground movements due to tunnel excavation
# Assumption that the tunnel face has advanced significantly far from the start of the tunnel (xf > 3z0) so the transverse profile behind the face develops fully
# x-xi is assumed to be infinite
def GroundMovement():

    Mesh = DefineMesh(i)

    # Extract x and y columns
    x = Mesh[:, 0]  
    y = Mesh[:, 1]  

    # Validation using x=4, y=1.5. Expecting w=0.00113m, v=-0.00223m, u=-0.0009m
    w_temp = (Vs / np.sqrt(2.0 * np.pi * i)) * np.exp(-(1.5**2) / (2.0 * i**2)) * (1 - 0.846)
    v_temp = -n / z0 * 1.5 * w_temp
    u_temp = (n * Vs / (2.0 * np.pi * z0)) * np.exp(-(1.5**2) / (2.0 * i**2)) * (-np.exp(-(4**2) / (2.0 * i**2)))
    print(w_temp)
    print(v_temp)
    print(u_temp)
    
    # Compute w, v, u
    w = (Vs / np.sqrt(2.0 * np.pi * i)) * np.exp(-(y**2) / (2.0 * i**2)) * (1.0 - norm.cdf(x / i))
    v = -n / z0 * y * w
    u = (n * Vs / (2.0 * np.pi * z0)) * np.exp(-(y**2) / (2.0 * i**2)) * (-np.exp(-(x**2) / (2.0 * i**2)))

    x = x/i
    y = y/i
  
    PlotContour(x, y, w, u, v)
    # PlotSurfaceSettlement(x,y,w)

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

    contour_w = ax[0].contour(Y, X, Z_w, levels=10, cmap='viridis')
    fig.colorbar(contour_w, ax=ax[0], label="Settlement (w)")
    ax[0].set_xlabel("Y Coordinate")
    ax[0].set_ylabel("X Coordinate")
    ax[0].set_title("Contour Plot of Surface Settlement")
    ax[0].invert_yaxis()

    contour_u = ax[1].contour(Y, X, Z_u, levels=10, cmap='viridis')
    fig.colorbar(contour_u, ax=ax[1], label="Displacement (u)")
    ax[1].set_xlabel("Y Coordinate")
    ax[1].set_ylabel("X Coordinate")
    ax[1].set_title("Contour Plot of X-Displacement")
    ax[1].invert_yaxis()

    contour_v = ax[2].contour(Y, X, Z_v, levels=10, cmap='viridis')
    fig.colorbar(contour_v, ax=ax[2], label="Displacement (v)")
    ax[2].set_xlabel("Y Coordinate")
    ax[2].set_ylabel("X Coordinate")
    ax[2].set_title("Contour Plot of Y-Displacement")
    ax[2].invert_yaxis()

    plt.grid()
    plt.tight_layout()  
    plt.show()

def PlotSurfaceSettlement(x, y, w):
  # Create a grid for contouring
    xi = np.linspace(min(x), max(x), 100)
    yi = np.linspace(min(y), max(y), 100)
    X, Y = np.meshgrid(xi, yi)

    # Interpolate the data to fit the grid using cubic interpolation
    Z = griddata((x, y), w, (X, Y), method='cubic')

    # Create a figure and axis for 3D plotting
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Plot the surface
    surf = ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='none')

    # Add color bar
    fig.colorbar(surf, ax=ax, label="Settlement (w)")

    # Labels and title
    ax.set_xlabel("X Coordinate")
    ax.set_ylabel("Y Coordinate")
    ax.set_zlabel("Settlement (w)")
    ax.set_title("3D Contour Plot of Settlement")
    ax.invert_zaxis()

    # Show the plot
    plt.show()

GroundMovement()