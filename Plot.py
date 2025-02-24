import numpy as np 
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

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