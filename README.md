<h1 align="center">
   Project for CSE528
  <br>
</h1>

<h2 align="center">Project On Surface Parameterization</h4>

![screenshot](https://github.com/xuan-li/GraphicsProject/blob/master/images/gui.png)

## Project Overview

In this project, I explore two kinds of parameterization algorithms. The first kind is from the view of Tutte embedding called Orbifold Tutte embedding (OTE). The other kind is from the view of differential geometry called boundary first flattening (BFF).

The first part follows my [proposal](https://github.com/xuan-li/GraphicsProject/blob/master/Documents/Project-Proposal/Proposal.pdf). The other part is for the bonus points.

[Full report on this project.](https://github.com/xuan-li/GraphicsProject/blob/master/Documents/Final-Report/main.pdf)

[Software Demo](https://github.com/xuan-li/GraphicsProject/blob/master/SoftwareDemo.mp4)


## References

- Aigerman, N. and Lipman, Y. (2015). Orbifold tutte embeddings. ACM Trans. Graph., 34(6):190:1–190:12.

- Aigerman, N. and Lipman, Y. (2016). Hyperbolic orbifold tutte embeddings. ACM Trans. Graph., 35(6):217:1–217:14.

- Sawhney, R. and Crane, K. (2017). Boundary first flattening. https://arxiv.org/abs/1704.06873. 

## Building Platform 

* Windows 10 x64
  
* Visual Studio 2017 

## Needed Library

* [OpenMesh](https://www.openmesh.org): Data structure used to represent polygonal meshes.

* [Libigl](http://libigl.github.io/libigl/): Generate GUI and solve quadratic programming.

* [Eigen](http://eigen.tuxfamily.org): Handle matrix operations.

* [LBFGS++](https://github.com/yixuan/LBFGSpp): Solve first order optimization.

Needed libraries are all included in Code/external/ and compiled into static libraries (*.lib).


## How To Use

### Folder structure

Released .exe file is in experiment/bin/ 

Test data are included in experiment/. 

There are seven folders in experiment/: 

- bin: released exe file

- Texture: some textures to show texture mappings.

- Euclidean Orbifold / Hyperbolic Orbifold: test data for OTE

- BoundaryFree / Polygon / ConeParameterization: test data for BFF


### Viewer Options

Show Option:

- Original: show loaded model.

- Sliced: show model after cut, which is a disk. 

- Embedding: show results of algorithms

- Covering Space: show tiled plane for orbifolds. (Only enabled in orbifold embeddings).

Show Slices and Cones: whether to show added cones and slices.



### Vertex Selection and Cutting System

Both algorithm needs to cut mesh and set angle sum for some vertices.

The loaded mesh is shown in "Original" mode.  In this mode, you can select vertices while pressing down 'S'. Input cone angle sum in the unit of <img src="https://latex.codecogs.com/gif.latex?\pi" />. Select one vertex to add a cone and select two vertices to add a slice.

Current marker can be saved and loaded. There's a .mark file associated with each model file.


### Orbifold Tutte Embedding

Only sphere-type Euclidean orbifolds with three cones supported: 
<img src="https://latex.codecogs.com/gif.latex?(\frac{\pi}{2},\pi,\frac{\pi}{2})" />,
<img src="https://latex.codecogs.com/gif.latex?(\frac{2\pi}{3},\frac{2\pi}{3},\frac{2\pi}{3})" />,
<img src="https://latex.codecogs.com/gif.latex?(\frac{\pi}{3},\frac{2\pi}{3},\pi)" />


Only sphere-type hyperbolic orbifolds with cone angles all <img src="https://latex.codecogs.com/gif.latex?\pi" /> supported. The number of cones for this kind should be larger than 4.

Needed .mark file is included for test.

Use Covering Space flag to see tiled plane: part of Euclidean plane and part of Poincare disk.


### Boundary First Flattening

Three settings: free boundary (Harmonic/Hilbert BFF with free B), polygonal boundary (Harmonic/Hilbert BFF with K), cone parameterization (Harmonic/Hilbert BFF with Cones).

Needed .mark file is included for test.


### Texture Mapping

After computation, load a texture. Choose "Show texture" flag.  See it in "Sliced" mode. 

## Build From Source

All codes are included in Code/.

Open .sln file in VS2017 and build. Only x64 mode supported. Needed libraries' path are set in relative mode, you needn't modify them.

#### License

MIT
