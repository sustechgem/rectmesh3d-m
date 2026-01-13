# rectmesh3d-m (v1.0)
A collection of utilities for 3D rectilinear mesh and model manipulation in geology and geophysics.

### Features
- Fully functional programming in MATLAB
- Minimum dependency on advanced/external toolboxes
- For the purpose of algorithm archiving and prototyping
- Designed for convenience, not for speed

### Definitions and Concepts
- A rectilinear mesh discretizes 3D space into many cuboidal **cells**;
the interfaces between adjacent cells are **faces**; each face is bounded by **edges**;
edges connect **nodes**.
- A **grid** is a lattice structure of points. Nodes of a mesh can form a grid, so do cell centers, face centers and edge centers.
- A right-handed coordinate system is adopted; A recommended geographical correspondence is x positive for east, y positive for north and z positive for upward.
- A mesh can be fully specified by three vectors **nodeX**, **nodeY** and **nodeZ**. Such a vector lists all possible x, y or z-coordinates of grid points, without recurrence of the same number. The numbers are ascending in nodeX (from west to east) and nodeY (from south to north), and descending in nodeZ (from top to bottom).
- By default, a **model** refers to a vector of numbers, each of which is assigned to a mesh cell. The convention of reshaping a 3D array to a vector is always that the z-coordinate of cell centers varies the most frequently, then the x-coordinate, and finally the y-coordinate. It is like painting walls from top to bottom, from left to right and from front to back when you are facing the southernmost wall of the mesh.

### Functions

#### addBlock
- Add cuboidal blocks to an existing model.

#### addCylinder
- Add round cylinders to an existing model.

#### addEllipsoid
- Add an ellipsoid to an existing model.

#### addSlab
- Add a slab to an existing model.

#### addSurface
- Add an arbitrary surface to an existing model. 

#### addDisk
- Add round disks to an existing model.

#### addPolyhedron
- Add a polyhedron to an existing model.

#### addSphere
- Add a sphere to an existing model.

#### CellIndex2PointXYZ
- Get the x-y-z locations of cell centers.

#### cellVolume
- Get cell volumes.

#### cropModel
- Cut out a portion of a model.

#### DirectionalIndex2GlobalIndex
- Get global indices of grid points with their directional indices in x, y, z given.

#### flatTopoModel
- Make a flat-topography model.

#### formLatticeTrilinearInterpMatrix
- Form an interpolation matrix that tri-linearly interpolates values on a lattice grid to
arbitrary points.

#### formMeshConversionMatrixCellCenterInterp
- Form a matrix that converts model1 defined on mesh1 to model2 defined on mesh2 using cell-center interpolation.

#### formMeshConversionMatrixVolumeWeighted
- Form a matrix that converts model1 defined on mesh1 to model2 defined on mesh2 using volume-weighting averaging.

#### getMeshPara
- Get basic mesh parameters from nodeX, nodeY and nodeZ.

#### GlobalIndex2DirectionalIndex
- Get directional indices of grid points in x, y, z with their global indices given.

#### makeRectMesh3D
- Make a 3D rectilinear mesh; the mesh has fine cells in the center and coarser cells toward the boundaries.

#### node2center
- Get locations of cell centers along x/y/z direction when nodes are given.

#### node2size
- Get cell sizes along x/y/z direction when nodes are given.

#### PointTopo2MeshTopo
- Convert topography represented by scattered points to topography model (0=air, 1=earth) in a 3D mesh.

#### PointXYZ2CellIndex
- Find the global indices of cells that encloses given points.

#### readMeshFileUBC
- Read 3D mesh file in UBC-GIF format from disk.

#### size2node
- Get nodes in x/y/z direction when cell sizes and the first node are given.

#### viewModelBlock
- View entire or slices or subdomains of 3D mesh/model.

#### writeMeshFileUBC
- Write 3D mesh file in UBC-GIF format to disk.
