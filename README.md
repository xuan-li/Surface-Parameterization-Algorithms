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


### Vertex Selection System

Both algorithm needs to cut mesh and set angle sum for some vertices.

The loaded mesh is shown in "Original" mode.  In this mode, you can select vertices while pressing down 'S'. Input cone angle sum in the unit of $\pi$  


You can [download](https://github.com/amitmerchant1990/electron-markdownify/releases/tag/v1.2.0) latest installable version of Markdownify for Windows, macOS and Linux.

## Credits

This software uses code from several open source packages.

- [Electron](http://electron.atom.io/)
- [Node.js](https://nodejs.org/)
- [Marked - a markdown parser](https://github.com/chjj/marked)
- [showdown](http://showdownjs.github.io/showdown/)
- [CodeMirror](http://codemirror.net/)
- Emojis are taken from [here](https://github.com/arvida/emoji-cheat-sheet.com)
- [highlight.js](https://highlightjs.org/)

## Related

[markdownify-web](https://github.com/amitmerchant1990/markdownify-web) - Web version of Markdownify

## You may also like...

- [Pomolectron](https://github.com/amitmerchant1990/pomolectron) - A pomodoro app
- [Correo](https://github.com/amitmerchant1990/correo) - A menubar/taskbar Gmail App for Windows and macOS

#### License

MIT

---

> [amitmerchant.com](https://www.amitmerchant.com) &nbsp;&middot;&nbsp;
> GitHub [@amitmerchant1990](https://github.com/amitmerchant1990) &nbsp;&middot;&nbsp;
> Twitter [@amit_merchant](https://twitter.com/amit_merchant)
