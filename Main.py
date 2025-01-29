import numpy as np 
import math

# Meshing Parameters
MeshLong = 1         # Mesh length along tunnel axis
MeshTrans = 1         # Mesh length transverse to tunnel axis

# Hard-coded variables
R = 4                   # Tunnel radius (m) 
z0 = 7.5                # Tunnel axis depth (m)
i = 1                # Settlement trough width parameter (m)
wmax = 0.00786          # Maximum vertical settlement (m) 

# Function to compute the x,y,z ground movements due to tunnel excavation
# Assumption that the tunnel face has advanced significantly far from the start of the tunnel (xf > 3z0) so the transverse profile behind the face develops fully
def GroundMovement():

    # For now, i is a hardcoded term but it can be calculated based on the tunnel depth

    ## Generate Mesh
    
    Mesh = DefineMesh(i)
    print(Mesh)
    

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

GroundMovement()