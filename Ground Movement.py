import numpy as np 
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.interpolate import griddata
import sympy as sp
from sympy.stats import Normal, cdf
from sympy import Symbol, diff, lambdify
from matplotlib.backends.backend_pdf import PdfPages

### IMPORTANT: x-direction is along the tunnel axis, runs up-down page. y-direction is perpendicular
# Mesh is defined in the typical x-y coordinates, but the contour plots show the flipped axes

# Meshing Parameters. Long is along tunnel axis, Trans is perpendicular
MeshLong = 0.5         
MeshTrans = 0.5        

# Hard-coded variables
# Radius
R = 1.105         
# Springline depth 
z0 = 7.5     
# Sampling point elevation (m)   
z = 0                
# Settlement trough width parameter (m)
# For now, i is a hardcoded term but it can be calculated based on the tunnel depth
i = 3.9             
n = 1
# Maximum vertical settlement (m) 
Wmax = 0.00786          
Vs = Wmax * i * (2*np.pi)**0.5
# Alternatively calculate Vs as a fraction of tunnel face area
#Vs = 0.05 * np.pi * R**2

# Set position of tunnel start and face relative to coordinate system. Origin is at tunnel face
xi = -1000
xf = 0

TOL = 1e-10

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

# Compute the x,y,z ground movements due to tunnel excavation
# Origin is set at the tunnel face
# u is positive in the forward direction. w is positive downward
def Main():

    Mesh = DefineMesh(i)

    x = Mesh[:, 0]  
    y = Mesh[:, 1]  

    x_sym = Symbol('x_sym')
    y_sym = Symbol('y_sym')
    
    w_eq = np.vectorize(lambdify((x_sym, y_sym), DisplacementEquations(x_sym, y_sym)[0], "numpy"))
    u_eq = np.vectorize(lambdify((x_sym, y_sym), DisplacementEquations(x_sym, y_sym)[1], "numpy"))
    v_eq = np.vectorize(lambdify((x_sym, y_sym), DisplacementEquations(x_sym, y_sym)[2], "numpy"))
    w = w_eq(x, y)
    u = u_eq(x, y)
    v = v_eq(x, y)

    eps_z_eq = np.vectorize(lambdify((x_sym, y_sym), StrainEquations(x_sym, y_sym)[0], "numpy"))
    eps_x_eq = np.vectorize(lambdify((x_sym, y_sym), StrainEquations(x_sym, y_sym)[1], "numpy"))
    eps_y_eq = np.vectorize(lambdify((x_sym, y_sym), StrainEquations(x_sym, y_sym)[2], "numpy"))
    eps_z = eps_z_eq(x, y)
    eps_x = eps_x_eq(x, y)
    eps_y = eps_y_eq(x, y)
   
    # Plot displacement contour in mm, strain in μ
    # x and y are scaled for graphing
    #Plot3D(x/i, y/i, w*1e3, u*1e3, v*1e3, eps_z*1e6, eps_x*1e6, eps_y*1e6, xi/i, xf/i)
    PlotContours(x/i, y/i, w*1e3, u*1e3, v*1e3, eps_z*1e6, eps_x*1e6, eps_y*1e6)

    StrainAlongTheta(x, y)
    DisplacementAlongTheta(x, y, u, v)

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

# Return symbolic math equations for the w, v, u
def DisplacementEquations(x, y):

    # Mean of 0, std dev of 1
    X = Normal('X', 0, 1)

    # Compute w, v, u in m
    w = (Vs / sp.sqrt(2*sp.pi) / i) * sp.exp(-y**2 / (2.0*i**2)) * (cdf(X)((x-xi)/i) - cdf(X)((x-xf)/i))
    u = (n * Vs / (2.0 * sp.pi * (z0 - z))) * sp.exp(-y**2 / (2*i**2)) * (sp.exp(-(x-xi)**2 / (2.0 * i **2)) - sp.exp(-(x-xf)**2 / (2.0 * i **2)))
    v = -y / (z0 - z) * w


    return w, u, v

# Return the symbolic math equations for the strains
def StrainEquations(x, y):

    w = DisplacementEquations(x, y)[0]
    X = Normal('X', 0, 1)

    # Compute the strains in z, y, x
    eps_z = -n * Vs / (sp.sqrt(2.0*sp.pi)*i*(z0-z)) * sp.exp(-y**2 / (2.0*i**2)) * (-1.0/sp.sqrt(2.0*sp.pi) * (((x-xi) / i) * sp.exp(-(x-xi)**2 / (2.0*i**2)) - ((x-xf) / i)*sp.exp(-(x-xf)**2/(2.0*i**2))) + (y**2 / i**2 - 1.0) * (cdf(X)((x-xi)/i) - cdf(X)((x-xf)/i))) 
    eps_x = -n * Vs / (2.0*sp.pi * i * (z0-z)) * sp.exp(-y**2 / (2.0*i**2)) * ((x-xi / i) * sp.exp(-(x-xi)**2 / (2.0*i**2)) - ((x-xf) / i) * sp.exp(-(x-xf)**2 / (2.0*i**2)))
    eps_y = n / (z0-z) * w * (y**2 / i**2 - 1.0)
    return eps_z, eps_x, eps_y

# Determine the strains in direction angle theta from x-axis at 5 degree increments
# Compute the strain tensor by doing the partial derivatives of u and v wrt x and y
# Then perform the strain transformation to get strain in the theta direction
def StrainAlongTheta(x, y):

    # Wall angle increments
    wall_angle = np.arange(0, 185, 5)
    #wall_angle = np.array([0,90])
    wall_angle_rad = np.deg2rad(wall_angle)

    x_sym = Symbol('x_sym')
    y_sym = Symbol('y_sym')

    u_sym = DisplacementEquations(x_sym, y_sym)[1]
    v_sym = DisplacementEquations(x_sym, y_sym)[2]

    eps_x_sym = StrainEquations(x_sym, y_sym)[1]
    eps_y_sym = StrainEquations(x_sym, y_sym)[2]
    eps_xy_sym = 0.5 * (diff(u_sym, y_sym) + diff(v_sym, x_sym))

    # Strain transformation equation to direction angle theta from x-axis
    theta_sym = Symbol('theta_sym')
    eps_theta_sym = ((eps_x_sym + eps_y_sym) / 2.0) + ((eps_x_sym - eps_y_sym) / 2.0 * sp.cos(2.0 * theta_sym)) + (eps_xy_sym * sp.sin(2.0 * theta_sym))
    eps_theta = lambdify((x_sym, y_sym, theta_sym), eps_theta_sym, "numpy")

    # Loop over every angle, and calculate strain at each point in mesh
    strain_results = np.zeros((len(x), len(wall_angle_rad)))
    for it1 in range(len(x)):
        for it2 in range(len(wall_angle_rad)):
            strain_results[it1, it2] = eps_theta(x[it1], y[it1], wall_angle_rad[it2])

    xi = np.linspace(min(x/i), max(x/i), 100)
    yi = np.linspace(min(y/i), max(y/i), 100)
    X, Y = np.meshgrid(xi, yi)
    
    with PdfPages('Strain Plots.pdf') as pdf:
        for it1 in range(len(wall_angle)):
            Z_eps = griddata((x/i, y/i), (strain_results[:,it1]*1e6), (X, Y), method='cubic')
            plt.figure(figsize = (10, 10))
            contour = plt.contour(Y, X, Z_eps, levels=10, cmap='viridis')
            plt.colorbar(contour, label="Horizontal Strain (με)", ax=plt.gca())
            plt.xlabel("Y Coordinate (y/i)")
            plt.ylabel("X Coordinate (x/i)")
            plt.title("Contour Plot of Horizontal Strain: " + str(wall_angle[it1]) + " Degrees from Tunnel Axis")
            plt.gca().invert_yaxis()
            plt.clabel(contour,fontsize=6,inline=1)
            plt.clabel(contour,fontsize=6,inline=1)
            plt.grid()
            pdf.savefig()
            plt.close()
    
    # Find the angle producing maximum strain at each offset (y-distance) from tunnel
    offsets = np.unique(y)
    indices = np.empty(len(offsets), dtype=object)
    for it1 in range(len(offsets)):
        indices[it1] = np.where(abs(y - offsets[it1]) <= TOL)[0]

    # Array to store the angles producing highest strain at each cross-section
    max_angles = np.zeros(len(offsets))
    # Corresponding max strain for each cross-section
    max_strain_for_offset = np.zeros(len(offsets))

    for it1 in range(len(offsets)):
        
        temp_max_strain = np.zeros(len(wall_angle))

        for it2 in range(len(wall_angle)):
            temp_strains = strain_results[:,it2]
            temp_indices = indices[it1]
            # strain at every point at tunnel offset cross-section for the specified angle
            strains_at_angle = temp_strains[temp_indices]
            # are we considering absolute or just positive strain?
            temp_max_strain[it2] = np.max(np.absolute(strains_at_angle))
        
        # Get the angle which causes the highest strain at the specified cross-section
        max_angles[it1] = wall_angle[temp_max_strain.argmax()]
        max_strain_for_offset[it1] = temp_max_strain[temp_max_strain.argmax()]

    fig, ax = plt.subplots(2, 1, figsize = (14, 9))
    ax[0].scatter(offsets, max_angles)
    ax[0].set_title("Angles Producing Maximum Strain vs. Tunnel Offset")
    ax[0].set_xlabel("Offset from Tunnel Axis (m)")
    ax[0].set_ylabel("Angle for Max Strain (degree)")

    ax[1].scatter(offsets, max_strain_for_offset)
    ax[1].set_title("Maximum Strain at Each Tunnel Offset")
    ax[1].set_xlabel("Offset from Tunnel Axis (m)")
    ax[1].set_ylabel("Maximum Horizontal Strain at Offset (με)")

    plt.tight_layout(pad=2.0)
            
def DisplacementAlongTheta(x, y, u, v):
    #wall_angle = np.array([0, 90])
    wall_angle = np.arange(0, 185, 5)
    wall_angle_rad = np.deg2rad(wall_angle)
    displacement_results = np.zeros((len(u),len(wall_angle_rad)))
    for it1 in range(len(u)):
        for it2 in range(len(wall_angle_rad)):
            displacement_results[it1, it2] = u[it1] * np.cos(wall_angle_rad[it2]) + v[it1] * np.cos(np.pi/2.0 - wall_angle_rad[it2])

    xi = np.linspace(min(x/i), max(x/i), 100)
    yi = np.linspace(min(y/i), max(y/i), 100)
    X, Y = np.meshgrid(xi, yi)
    
    with PdfPages('Displacement Plots.pdf') as pdf:
        for it1 in range(len(wall_angle)):
            Z_disp = griddata((x/i, y/i), (displacement_results[:,it1]*1e3), (X, Y), method='cubic')
            plt.figure(figsize = (10, 10))
            contour = plt.contour(Y, X, Z_disp, levels=10, cmap='viridis')
            plt.colorbar(contour, label="Displacement (mm)", ax=plt.gca())
            plt.xlabel("Y Coordinate (y/i)")
            plt.ylabel("X Coordinate (x/i)")
            plt.title("Contour Plot of Displacement: " + str(wall_angle[it1]) + " Degrees from Tunnel Axis")
            plt.gca().invert_yaxis()
            plt.clabel(contour,fontsize=6,inline=1)
            plt.clabel(contour,fontsize=6,inline=1)
            plt.grid()
            pdf.savefig()
            plt.close()  

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
    ax1.plot([max(xStart, -4), xFinish], [0,0], [0,0], color='red', linewidth=2, zorder=500)

    # Plot for U displacement
    ax2 = fig1.add_subplot(132, projection='3d')
    surf2 = ax2.plot_surface(X, Y, Zu, cmap='viridis', edgecolor='none', zorder=100)
    fig1.colorbar(surf2, ax=ax2, label="Displacement (mm)")
    ax2.set_xlabel("X Coordinate (x/i)")
    ax2.set_ylabel("Y Coordinate (y/i)")
    ax2.set_title("3D Plot of Displacement Along Tunnel Axis (U)")
    ax2.plot(0, 0, zs=0, zdir='z', color='red', linewidth=2)
    ax2.plot([max(xStart, -4), xFinish], [0,0], [0,0], color='red', linewidth=2, zorder=500)

    # Plot for V displacement
    ax3 = fig1.add_subplot(133, projection='3d')
    surf3 = ax3.plot_surface(X, Y, Zv, cmap='viridis', edgecolor='none', zorder=100)
    fig1.colorbar(surf3, ax=ax3, label="Displacement (mm)")
    ax3.set_xlabel("X Coordinate (x/i)")
    ax3.set_ylabel("Y Coordinate (y/i)")
    ax3.set_title("3D Plot of Displacement Perpendicular Tunnel Axis (V)")
    ax3.plot([max(xStart, -4), xFinish], [0,0], [0,0], color='red', linewidth=2, zorder=500)

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
    ax4.plot([max(xStart, -4), xFinish], [0,0], [0,0], color='red', linewidth=2, zorder=500)

    # Plot for X-strain
    ax5 = fig2.add_subplot(132, projection='3d')
    surf5 = ax5.plot_surface(X, Y, Zeps_x, cmap='viridis', edgecolor='none', zorder=100)
    fig2.colorbar(surf5, ax=ax5, label="Strain (με)")
    ax5.set_xlabel("X Coordinate (x/i)")
    ax5.set_ylabel("Y Coordinate (y/i)")
    ax5.set_title("3D Plot of Strain along Tunnel Axis")
    ax5.plot(0, 0, zs=0, zdir='z', color='red', linewidth=2)
    ax5.plot([max(xStart, -4), xFinish], [0,0], [0,0], color='red', linewidth=2, zorder=500)

    # Plot for Y-strain
    ax6 = fig2.add_subplot(133, projection='3d')
    surf6 = ax6.plot_surface(X, Y, Zeps_y, cmap='viridis', edgecolor='none', zorder=100)
    fig2.colorbar(surf6, ax=ax6, label="Strain (με)")
    ax6.set_xlabel("X Coordinate (x/i)")
    ax6.set_ylabel("Y Coordinate (y/i)")
    ax6.set_title("3D Plot of Strain Perpendicular Tunnel Axis")
    ax6.plot([max(xStart, -4), xFinish], [0,0], [0,0], color='red', linewidth=2, zorder=500)

    plt.tight_layout()

Main()
plt.show()