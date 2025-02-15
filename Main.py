import numpy as np 
import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.stats import norm
from scipy.interpolate import griddata

# Meshing Parameters
MeshLong = 0.5          # Mesh length along tunnel axis
MeshTrans = 0.5         # Mesh length transverse to tunnel axis

# Hard-coded variables
R = 1.105           # Tunnel radius (m) 
z0 = 10.5                # Tunnel axis depth (m)
z = 0.5                   # Sampling point elevation (m)
# For now, i is a hardcoded term but it can be calculated based on the tunnel depth
i = 2.7                 # Settlement trough width parameter (m)
n = 1
# Wmax = 0.00786          # Maximum vertical settlement (m) 
# Vs = Wmax * i * (2*np.pi)**0.5
# Alternatively calculate Vs as a fraction of tunnel face area
Vs = 0.05 * np.pi * R**2

# Set position of tunnel start and face relative to coordinate system. Origin is at tunnel face
xi = -10000
xf = 0

# Building foundation to be plotted on contour map
foundation_corners = [[-0.5/i, -4.0/i], [-6.0/i, -9.6/i], [-9.5/i, -6.1/i], [-4.0/i, -0.5/i]]

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

    # Plot displacement contour in mm, strain in μ
    #PlotContours(x, y, w*1e3, u*1e3, v*1e3, eps_z*1e6, eps_x*1e6, eps_y*1e6)
    PlotSurfaceSettlement(x,y,w*10e3, u*10e3, v*10e3)

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

def PlotContours(x, y, w, u, v, eps_z, eps_x, eps_y):
    # Create a grid for contouring
    xi = np.linspace(min(x), max(x), 100)
    yi = np.linspace(min(y), max(y), 100)
    X, Y = np.meshgrid(xi, yi)

    # Interpolate data to fit grid for w, u, and v
    Z_w = griddata((x, y), w, (X, Y), method='cubic')
    Z_u = griddata((x, y), u, (X, Y), method='cubic')
    Z_v = griddata((x, y), v, (X, Y), method='cubic')
    Z_eps_z = griddata((x, y), eps_z, (X, Y), method='cubic')
    Z_eps_x = griddata((x, y), eps_x, (X, Y), method='cubic')
    Z_eps_y = griddata((x, y), eps_y, (X, Y), method='cubic')

    fig, ax = plt.subplots(2, 3, figsize=(14, 9))

    contour_w = ax[0, 0].contour(Y, X, Z_w, levels=10, cmap='viridis')
    fig.colorbar(contour_w, ax=ax[0, 0], label="Settlement (mm)")
    ax[0, 0].set_xlabel("Y Coordinate (y/i)")
    ax[0, 0].set_ylabel("X Coordinate (x/i)")
    ax[0, 0].set_title("Contour Plot of Surface Settlement (W)")
    ax[0, 0].invert_yaxis()
    ax[0, 0].clabel(contour_w,fontsize=6,inline=1)
    ax[0, 0].clabel(contour_w,fontsize=6,inline=1)
    ax[0, 0].grid()

    contour_u = ax[0, 1].contour(Y, X, Z_u, levels=10, cmap='viridis')
    fig.colorbar(contour_u, ax=ax[0, 1], label="Displacement (mm)")
    ax[0, 1].set_xlabel("Y Coordinate (y/i)")
    ax[0, 1].set_ylabel("X Coordinate (x/i)")
    ax[0, 1].set_title("Contour Plot of X-Displacement (U)")
    ax[0, 1].invert_yaxis()
    ax[0, 1].clabel(contour_u,fontsize=6,inline=1)
    ax[0, 1].clabel(contour_u,fontsize=6,inline=1)
    ax[0, 1].grid()

    contour_v = ax[0, 2].contour(Y, X, Z_v, levels=10, cmap='viridis')
    fig.colorbar(contour_v, ax=ax[0, 2], label="Displacement (mm)")
    ax[0, 2].set_xlabel("Y Coordinate (y/i)")
    ax[0, 2].set_ylabel("X Coordinate (x/i)")
    ax[0, 2].set_title("Contour Plot of Y-Displacement (V)")
    ax[0, 2].invert_yaxis()
    ax[0, 2].clabel(contour_v,fontsize=6,inline=1)
    ax[0, 2].clabel(contour_v,fontsize=6,inline=1)
    ax[0, 2].grid()

    contour_eps_z = ax[1, 0].contour(Y, X, Z_eps_z, levels=10, cmap='viridis')
    fig.colorbar(contour_eps_z, ax=ax[1, 0], label="Vertical Strain (με)")
    ax[1, 0].set_xlabel("Y Coordinate (y/i)")
    ax[1, 0].set_ylabel("X Coordinate (x/i)")
    ax[1, 0].set_title("Contour Plot of Vertical Strain (εz)")
    ax[1, 0].invert_yaxis()
    ax[1, 0].clabel(contour_eps_z,fontsize=6,inline=1)
    ax[1, 0].clabel(contour_eps_z,fontsize=6,inline=1)
    ax[1, 0].grid()

    contour_eps_x = ax[1, 1].contour(Y, X, Z_eps_x, levels=10, cmap='viridis')
    fig.colorbar(contour_eps_x, ax=ax[1, 1], label="X-Horizontal Strain (με)")
    ax[1, 1].set_xlabel("Y Coordinate (y/i)")
    ax[1, 1].set_ylabel("X Coordinate (x/i)")
    ax[1, 1].set_title("Contour Plot of X-Horizontal Strain (εx)")
    ax[1, 1].invert_yaxis()
    ax[1, 1].clabel(contour_eps_x,fontsize=6,inline=1)
    ax[1, 1].clabel(contour_eps_x,fontsize=6,inline=1)
    ax[1, 1].grid()

    contour_eps_y = ax[1, 2].contour(Y, X, Z_eps_y, levels=10, cmap='viridis')
    fig.colorbar(contour_eps_y, ax=ax[1, 2], label="Y-Horizontal Strain (με)")
    ax[1, 2].set_xlabel("Y Coordinate (y/i)")
    ax[1, 2].set_ylabel("X Coordinate (x/i)")
    ax[1, 2].set_title("Contour Plot of Y-Horizontal Strain (εy)")
    ax[1, 2].invert_yaxis()
    ax[1, 2].clabel(contour_eps_y,fontsize=6,inline=1)
    ax[1, 2].clabel(contour_eps_y,fontsize=6,inline=1)
    ax[1, 2].grid()
    '''
    for row in ax:
        for subplot in row:
            subplot.add_patch(patches.Polygon(foundation_corners, closed=True, edgecolor='b', facecolor='none'))
    '''
    plt.tight_layout()  
    plt.show()

def PlotSurfaceSettlement(x, y, w, u, v):
    # Create a grid for contouring
    xi = np.linspace(min(x), max(x), 100)
    yi = np.linspace(min(y), max(y), 100)
    X, Y = np.meshgrid(xi, yi)

    # Interpolate the data to fit the grid using cubic interpolation
    Zw = griddata((x, y), w, (X, Y), method='cubic')
    Zu = griddata((x, y), u, (X, Y), method='cubic')
    Zv = griddata((x, y), v, (X, Y), method='cubic')

    fig = plt.figure(figsize=(18, 6))

    # Plot for settlement
    ax1 = fig.add_subplot(131, projection='3d')
    surf1 = ax1.plot_surface(X, Y, Zw, cmap='viridis', edgecolor='none')
    fig.colorbar(surf1, ax=ax1, label="Settlement (mm)")
    ax1.set_xlabel("X Coordinate (x/i)")
    ax1.set_ylabel("Y Coordinate (y/i)")
    ax1.set_zlabel("Settlement (mm)")
    ax1.set_title("3D Contour Plot of Settlement (W)")
    ax1.invert_zaxis()

    # Plot for U displacement
    ax2 = fig.add_subplot(132, projection='3d')
    surf2 = ax2.plot_surface(X, Y, Zu, cmap='viridis', edgecolor='none')
    fig.colorbar(surf2, ax=ax2, label="Displacement (mm)")
    ax2.set_xlabel("X Coordinate (x/i)")
    ax2.set_ylabel("Y Coordinate (y/i)")
    ax2.set_zlabel("Displacement (mm)")
    ax2.set_title("3D Contour Plot of Displacement (U)")
    ax2.invert_zaxis()

    # Plot for V displacement
    ax3 = fig.add_subplot(133, projection='3d')
    surf3 = ax3.plot_surface(X, Y, Zv, cmap='viridis', edgecolor='none')
    fig.colorbar(surf3, ax=ax3, label="Displacement (mm)")
    ax3.set_xlabel("X Coordinate (x/i)")
    ax3.set_ylabel("Y Coordinate (y/i)")
    ax3.set_zlabel("Displacement (mm)")
    ax3.set_title("3D Contour Plot of Displacement (V)")
    ax3.invert_zaxis()

    plt.tight_layout()
    plt.show()

# ValidateEquations()
GroundMovement()