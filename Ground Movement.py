import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from scipy.stats import norm
from scipy.interpolate import griddata
import sympy as sp
from sympy.stats import Normal, cdf
from sympy import Symbol, diff, lambdify

# Meshing Parameters
MeshLong = 0.5          # Mesh length along tunnel axis
MeshTrans = 0.5         # Mesh length transverse to tunnel axis

# Hard-coded variables
R = 1.105           # Tunnel radius (m) 
z0 = 10.5             # Tunnel axis depth (m)
z = 0.5                 # Sampling point elevation (m)
# For now, i is a hardcoded term but it can be calculated based on the tunnel depth
i = 2.7               # Settlement trough width parameter (m)
n = 1
#Wmax = 0.00786          # Maximum vertical settlement (m) 
#Vs = Wmax * i * (2*np.pi)**0.5
# Alternatively calculate Vs as a fraction of tunnel face area
Vs = 0.05 * np.pi * R**2

# Set position of tunnel start and face relative to coordinate system. Origin is at tunnel face
xi = -1000
xf = 0

# Wall to be plotted on contour map (unscaled)
wall_corners = np.array([[-0.5, -4.0], [-6.0, -9.6]])

# Wall angle increments
wall_angles = np.arange(0, 185, 5)
print(wall_angles)

def ValidateEquations():
    ## Validation of equations
    x_xf_temp = 4
    x_xi_temp = 1000000
    y_temp  = 1.5
    i_temp = 3.9
    Wmax_temp = 0.00786
    Vs_temp = Wmax_temp * i_temp * (2*np.pi)**0.5
    n_temp = 1
    z0_z_temp = 7.5

    w_temp = (Vs_temp / np.sqrt(2*np.pi) / i_temp) * np.exp(-y_temp**2 / (2*i_temp**2)) * (norm.cdf(x_xi_temp / i_temp) - norm.cdf(x_xf_temp / i_temp)) 
    v_temp = - y_temp / (z0_z_temp) * w_temp
    u_temp = (n_temp * Vs_temp / (2.0 * np.pi * (z0_z_temp))) * np.exp(-y_temp**2 / (2*i_temp**2)) * (np.exp(-x_xi_temp**2 / (2.0 * i_temp **2)) - np.exp(-x_xf_temp**2 / (2.0 * i_temp **2)))

    print("W Expected: 1.13mm  ---  Calculated:"+  str(round(w_temp*1000,3)))
    print("V Expected:-2.23mm  ---  Calculated:"+  str(round(v_temp*1000,3)))
    print("U Expected:-0.90mm  ---  Calculated:"+  str(round(u_temp*1000,3)))

    eps_z_temp = -n_temp * Vs_temp / (np.sqrt(2.0*np.pi)*i_temp*(z0_z_temp)) * np.exp(-y_temp**2 / (2.0*i_temp**2)) * (-1.0/np.sqrt(2.0*np.pi) * ((x_xi_temp / i_temp) * np.exp(-x_xi_temp**2 / (2.0*i_temp**2)) - (x_xf_temp / i_temp)*np.exp(-x_xf_temp**2/(2.0*i_temp**2))) + (y_temp**2 / i_temp**2 - 1.0) * (norm.cdf(x_xi_temp / i_temp) - norm.cdf(x_xf_temp / i_temp))) 
    eps_y_temp = n_temp / (z0_z_temp) * w_temp * (y_temp**2 / i_temp**2 - 1.0)
    eps_x_temp = -n_temp * Vs_temp / (2.0*np.pi * i_temp * (z0_z_temp)) * np.exp(-y_temp**2 / (2.0*i_temp**2)) * ((x_xi_temp / i_temp) * np.exp(-x_xi_temp**2 / (2.0*i_temp**2)) - (x_xf_temp / i_temp) * np.exp(-x_xf_temp**2 / (2.0*i_temp**2)))

    print("Eps_z Expected: -110u  ---  Calculated:"+  str(round(eps_z_temp*1000000)))
    print("Eps_y Expected: -130u  ---  Calculated:"+  str(round(eps_y_temp*1000000)))
    print("Eps_x Expected: +240u  ---  Calculated:"+  str(round(eps_x_temp*1000000)))

# Function to compute the x,y,z ground movements due to tunnel excavation
# Origin is set at the tunnel face
# u is positive in the forward direction. w is positive downward
def Main():

    Mesh = DefineMesh(i)

    # Extract x and y columns
    x = Mesh[:, 0]  
    y = Mesh[:, 1]  
    
    w, v, u = ComputeDisplacements(x, y)
    eps_z, eps_x, eps_y = ComputeStrains(x, y, w)
   
    # Plot displacement contour in mm, strain in μ
    # x and y are scaled for graphing
    #Plot3D(x/i, y/i, w*1e3, u*1e3, v*1e3, eps_z*1e6, eps_x*1e6, eps_y*1e6, xi, xf)
    PlotContours(x/i, y/i, w*1e3, u*1e3, v*1e3, eps_z*1e6, eps_x*1e6, eps_y*1e6)

    x_wall = wall_corners[:,0]
    y_wall = wall_corners[:,1]
    wall_length = np.sqrt((x_wall[0] - x_wall[1])**2 + (y_wall[0] - y_wall[1])**2)
    
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

def ComputeDisplacements(x, y):

    x_xi = x - xi
    x_xf = x - xf

    # Compute w, v, u in m
    w = (Vs / np.sqrt(2*np.pi) / i) * np.exp(-y**2 / (2.0*i**2)) * (norm.cdf(x_xi / i) - norm.cdf(x_xf / i))
    v = -y / (z0 - z) * w
    u = (n * Vs / (2.0 * np.pi * (z0 - z))) * np.exp(-y**2 / (2*i**2)) * (np.exp(-x_xi**2 / (2.0 * i **2)) - np.exp(-x_xf**2 / (2.0 * i **2)))

    return w, v, u

def ComputeStrains(x, y, w):

    x_xi = x - xi
    x_xf = x - xf

    # Compute the strains in z, y, x
    eps_z = -n * Vs / (np.sqrt(2.0*np.pi)*i*(z0-z)) * np.exp(-y**2 / (2.0*i**2)) * (-1.0/np.sqrt(2.0*np.pi) * ((x_xi / i) * np.exp(-x_xi**2 / (2.0*i**2)) - (x_xf / i)*np.exp(-x_xf**2/(2.0*i**2))) + (y**2 / i**2 - 1.0) * (norm.cdf(x_xi / i) - norm.cdf(x_xf / i))) 
    eps_x = -n * Vs / (2.0*np.pi * i * (z0-z)) * np.exp(-y**2 / (2.0*i**2)) * ((x_xi / i) * np.exp(-x_xi**2 / (2.0*i**2)) - (x_xf / i) * np.exp(-x_xf**2 / (2.0*i**2)))
    eps_y = n / (z0-z) * w * (y**2 / i**2 - 1.0)

    return eps_z, eps_x, eps_y

def StrainAlongTheta():
    x = Symbol

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

    for row in ax:
        for subplot in row:
            subplot.add_patch(patches.Polygon(wall_corners / i, closed=True, edgecolor='b', facecolor='none'))

    plt.tight_layout()  

def Plot3D(x, y, w, u, v, eps_z, eps_x, eps_y, xStart, xFinish):
    # Create a grid for contouring
    xi = np.linspace(min(x), max(x), 100)
    yi = np.linspace(min(y), max(y), 100)
    X, Y = np.meshgrid(xi, yi)

    # Interpolate the data to fit the grid using cubic interpolation
    Zw = griddata((x, y), w, (X, Y), method='cubic')
    Zu = griddata((x, y), u, (X, Y), method='cubic')
    Zv = griddata((x, y), v, (X, Y), method='cubic')
    Zeps_z = griddata((x, y), eps_z, (X, Y), method='cubic')
    Zeps_x = griddata((x, y), eps_x, (X, Y), method='cubic')
    Zeps_y = griddata((x, y), eps_y, (X, Y), method='cubic')

    ## Plot the displacements

    fig1 = plt.figure(figsize=(18, 6))

    # Plot for settlement
    ax1 = fig1.add_subplot(131, projection='3d')
    surf1 = ax1.plot_surface(X, Y, Zw, cmap='viridis', edgecolor='none', zorder=100)
    fig1.colorbar(surf1, ax=ax1, label="Settlement (mm)")
    ax1.set_xlabel("X Coordinate (x/i)")
    ax1.set_ylabel("Y Coordinate (y/i)")
    ax1.set_title("3D Contour Plot of Settlement (W)")
    ax1.invert_zaxis()
    ax1.plot([max(xStart/i, -4), xFinish/i], [0,0], [0,0], color='red', linewidth=2, zorder=500)

    # Plot for U displacement
    ax2 = fig1.add_subplot(132, projection='3d')
    surf2 = ax2.plot_surface(X, Y, Zu, cmap='viridis', edgecolor='none', zorder=100)
    fig1.colorbar(surf2, ax=ax2, label="Displacement (mm)")
    ax2.set_xlabel("X Coordinate (x/i)")
    ax2.set_ylabel("Y Coordinate (y/i)")
    ax2.set_title("3D Plot of Displacement Along Tunnel Axis (U)")
    ax2.plot(0, 0, zs=0, zdir='z', color='red', linewidth=2)
    ax2.plot([max(xStart/i, -4), xFinish/i], [0,0], [0,0], color='red', linewidth=2, zorder=500)

    # Plot for V displacement
    ax3 = fig1.add_subplot(133, projection='3d')
    surf3 = ax3.plot_surface(X, Y, Zv, cmap='viridis', edgecolor='none', zorder=100)
    fig1.colorbar(surf3, ax=ax3, label="Displacement (mm)")
    ax3.set_xlabel("X Coordinate (x/i)")
    ax3.set_ylabel("Y Coordinate (y/i)")
    ax3.set_title("3D Plot of Displacement Perpendicular Tunnel Axis (V)")
    ax3.plot([max(xStart/i, -4), xFinish/i], [0,0], [0,0], color='red', linewidth=2, zorder=500)

    plt.tight_layout()

    ## Plot the strains

    fig2 = plt.figure(figsize=(18, 6))

    # Plot for Z-strain
    ax4 = fig2.add_subplot(131, projection='3d')
    surf4 = ax4.plot_surface(X, Y, Zeps_z, cmap='viridis', edgecolor='none', zorder=100)
    fig2.colorbar(surf4, ax=ax4, label="Strain (με)")
    ax4.set_xlabel("X Coordinate (x/i)")
    ax4.set_ylabel("Y Coordinate (y/i)")
    ax4.set_title("3D Contour Plot of Vertical Strain")
    ax4.invert_zaxis()
    ax4.plot([max(xStart/i, -4), xFinish/i], [0,0], [0,0], color='red', linewidth=2, zorder=500)

    # Plot for X-strain
    ax5 = fig2.add_subplot(132, projection='3d')
    surf5 = ax5.plot_surface(X, Y, Zeps_x, cmap='viridis', edgecolor='none', zorder=100)
    fig2.colorbar(surf5, ax=ax5, label="Strain (με)")
    ax5.set_xlabel("X Coordinate (x/i)")
    ax5.set_ylabel("Y Coordinate (y/i)")
    ax5.set_title("3D Plot of Strain along Tunnel Axis")
    ax5.plot(0, 0, zs=0, zdir='z', color='red', linewidth=2)
    ax5.plot([max(xStart/i, -4), xFinish/i], [0,0], [0,0], color='red', linewidth=2, zorder=500)

    # Plot for Y-strain
    ax6 = fig2.add_subplot(133, projection='3d')
    surf6 = ax6.plot_surface(X, Y, Zeps_y, cmap='viridis', edgecolor='none', zorder=100)
    fig2.colorbar(surf6, ax=ax6, label="Strain (με)")
    ax6.set_xlabel("X Coordinate (x/i)")
    ax6.set_ylabel("Y Coordinate (y/i)")
    ax6.set_title("3D Plot of Strain Perpendicular Tunnel Axis")
    ax6.plot([max(xStart/i, -4), xFinish/i], [0,0], [0,0], color='red', linewidth=2, zorder=500)

    plt.tight_layout()

Main()
plt.show()