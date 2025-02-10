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
z = 0                   # Sampling surface points
# For now, i is a hardcoded term but it can be calculated based on the tunnel depth
i = 3.9                 # Settlement trough width parameter (m)
n = 1
Wmax = 0.00786          # Maximum vertical settlement (m) 
Vs = Wmax * i * (2*np.pi)**0.5
# Alternatively calculate Vs as a fraction of tunnel face area
# Vs = 0.025 * np.pi * R**2

# Set position of tunnel start and face relative to coordinate system. Origin is at tunnel face
xi = -1000
xf = 0

def ValidateEquations():
    ## Validation of equations
    x_xf_temp = 4
    x_xi_temp = 1000000
    y_temp  = 1.5

    w_temp = (Vs / np.sqrt(2*np.pi) / i) * np.exp(-y_temp**2 / (2*i**2)) * (norm.cdf(x_xi_temp / i) - norm.cdf(x_xf_temp / i)) 
    v_temp = - y_temp / (z0 - z) * w_temp
    u_temp = (n * Vs / (2.0 * np.pi * (z0 - z))) * np.exp(-y_temp**2 / (2*i**2)) * (np.exp(-x_xi_temp**2 / (2.0 * i **2)) - np.exp(-x_xf_temp**2 / (2.0 * i **2)))

    print("W Expected: 1.13mm  ---  Calculated:"+  str(round(w_temp*1000,3)))
    print("V Expected:-2.23mm  ---  Calculated:"+  str(round(v_temp*1000,3)))
    print("U Expected:-0.90mm  ---  Calculated:"+  str(round(u_temp*1000,3)))

    eps_z_temp = -n * Vs / (np.sqrt(2.0*np.pi)*i*(z0-z)) * np.exp(-y_temp**2 / (2.0*i**2)) * (-1.0/np.sqrt(2.0*np.pi) * ((x_xi_temp / i) * np.exp(-x_xi_temp**2 / (2.0*i**2)) - (x_xf_temp / i)*np.exp(-x_xf_temp**2/(2.0*i**2))) + (y_temp**2 / i**2 - 1.0) * (norm.cdf(x_xi_temp / i) - norm.cdf(x_xf_temp / i))) 
    eps_y_temp = n / (z0-z) * w_temp * (y_temp**2 / i**2 - 1.0)
    eps_x_temp = -n * Vs / (2.0*np.pi * i * (z0-z)) * np.exp(-y_temp**2 / (2.0*i**2)) * ((x_xi_temp / i) * np.exp(-x_xi_temp**2 / (2.0*i**2)) - (x_xf_temp / i) * np.exp(-x_xf_temp**2 / (2.0*i**2)))

    print("Eps_z Expected: -110u  ---  Calculated:"+  str(round(eps_z_temp*1000000)))
    print("Eps_y Expected: -130u  ---  Calculated:"+  str(round(eps_y_temp*1000000)))
    print("Eps_x Expected: +240u  ---  Calculated:"+  str(round(eps_x_temp*1000000)))

# Function to compute the x,y,z ground movements due to tunnel excavation
# Origin is set at the tunnel face
# u is positive in the forward direction. w is positive downward
def GroundMovement():

    Mesh = DefineMesh(i)

    # Extract x and y columns
    x = Mesh[:, 0]  
    y = Mesh[:, 1]  

    x_xi = x - xi
    x_xf = x - xf
    
    # Compute w, v, u in m
    w = (Vs / np.sqrt(2*np.pi) / i) * np.exp(-y**2 / (2.0*i**2)) * (norm.cdf(x_xi / i) - norm.cdf(x_xf / i))
    v = -y / (z0 - z) * w
    u = (n * Vs / (2.0 * np.pi * (z0 - z))) * np.exp(-y**2 / (2*i**2)) * (np.exp(-x_xi**2 / (2.0 * i **2)) - np.exp(-x_xf**2 / (2.0 * i **2)))

    # Compute the strains in z, y, x
    eps_z = -n * Vs / (np.sqrt(2.0*np.pi)*i*(z0-z)) * np.exp(-y**2 / (2.0*i**2)) * (-1.0/np.sqrt(2.0*np.pi) * ((x_xi / i) * np.exp(-x_xi**2 / (2.0*i**2)) - (x_xf / i)*np.exp(-x_xf**2/(2.0*i**2))) + (y**2 / i**2 - 1.0) * (norm.cdf(x_xi / i) - norm.cdf(x_xf / i))) 
    eps_y = n / (z0-z) * w * (y**2 / i**2 - 1.0)
    eps_x = -n * Vs / (2.0*np.pi * i * (z0-z)) * np.exp(-y**2 / (2.0*i**2)) * ((x_xi / i) * np.exp(-x_xi**2 / (2.0*i**2)) - (x_xf / i) * np.exp(-x_xf**2 / (2.0*i**2)))

    # Scale x and y for graphing
    x = x/i
    y = y/i

    # Scale displacements to mm, scale strains up by 10e6
    w *= 1000
    v *= 1000
    u *= 1000
    eps_z *= 10e6
    eps_y *= 10e6
    eps_x *= 10e6
  
    # PlotDisplacements(x, y, w, u, v)
    PlotStrains(x, y, eps_z, eps_x, eps_y)
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

def PlotDisplacements(x, y, w, u, v):
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
    ax[0].set_title("Contour Plot of Surface Settlement (W)")
    ax[0].invert_yaxis()
    ax[0].clabel(contour_w,fontsize=6,inline=1)
    ax[0].clabel(contour_w,fontsize=6,inline=1)
    ax[0].grid()
    contour_u = ax[1].contour(Y, X, Z_u, levels=10, cmap='viridis')
    fig.colorbar(contour_u, ax=ax[1], label="Displacement (u)")
    ax[1].set_xlabel("Y Coordinate")
    ax[1].set_ylabel("X Coordinate")
    ax[1].set_title("Contour Plot of X-Displacement (U)")
    ax[1].invert_yaxis()
    ax[1].clabel(contour_u,fontsize=6,inline=1)
    ax[1].clabel(contour_u,fontsize=6,inline=1)
    ax[1].grid()
    contour_v = ax[2].contour(Y, X, Z_v, levels=10, cmap='viridis')
    fig.colorbar(contour_v, ax=ax[2], label="Displacement (V)")
    ax[2].set_xlabel("Y Coordinate")
    ax[2].set_ylabel("X Coordinate")
    ax[2].set_title("Contour Plot of Y-Displacement (V)")
    ax[2].invert_yaxis()
    ax[2].clabel(contour_v,fontsize=6,inline=1)
    ax[2].clabel(contour_v,fontsize=6,inline=1)
    ax[2].grid()
    plt.tight_layout()  
    plt.show()

def PlotStrains(x, y, eps_z, eps_x, eps_y):
    # Create a grid for contouring
    xi = np.linspace(min(x), max(x), 100)
    yi = np.linspace(min(y), max(y), 100)
    X, Y = np.meshgrid(xi, yi)

    # Interpolate data to fit grid for eps_z, eps_x, and eps_y
    Z_eps_z = griddata((x, y), eps_z, (X, Y), method='cubic')
    Z_eps_x = griddata((x, y), eps_x, (X, Y), method='cubic')
    Z_eps_y = griddata((x, y), eps_y, (X, Y), method='cubic')

    fig, ax = plt.subplots(1, 3, figsize=(18, 6))

    contour_eps_z = ax[0].contour(Y, X, Z_eps_z, levels=10, cmap='viridis')
    fig.colorbar(contour_eps_z, ax=ax[0], label="Vertical Strain (eps_z)")
    ax[0].set_xlabel("Y Coordinate")
    ax[0].set_ylabel("X Coordinate")
    ax[0].set_title("Contour Plot of Vertical Strain (eps_z)")
    ax[0].invert_yaxis()
    ax[0].clabel(contour_eps_z,fontsize=6,inline=1)
    ax[0].clabel(contour_eps_z,fontsize=6,inline=1)
    ax[0].grid()
    contour_eps_x = ax[1].contour(Y, X, Z_eps_x, levels=10, cmap='viridis')
    fig.colorbar(contour_eps_x, ax=ax[1], label="X-Horizontal Strain (eps_x)")
    ax[1].set_xlabel("Y Coordinate")
    ax[1].set_ylabel("X Coordinate")
    ax[1].set_title("Contour Plot of X-Horizontal Strain (eps_x)")
    ax[1].invert_yaxis()
    ax[1].clabel(contour_eps_x,fontsize=6,inline=1)
    ax[1].clabel(contour_eps_x,fontsize=6,inline=1)
    ax[1].grid()
    contour_eps_y = ax[2].contour(Y, X, Z_eps_y, levels=10, cmap='viridis')
    fig.colorbar(contour_eps_y, ax=ax[2], label="Y-Horizontal Strain (eps_y)")
    ax[2].set_xlabel("Y Coordinate")
    ax[2].set_ylabel("X Coordinate")
    ax[2].set_title("Contour Plot of Y-Horizontal Strain (eps_y)")
    ax[2].invert_yaxis()
    ax[2].clabel(contour_eps_y,fontsize=6,inline=1)
    ax[2].clabel(contour_eps_y,fontsize=6,inline=1)
    ax[2].grid()
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

ValidateEquations()
GroundMovement()