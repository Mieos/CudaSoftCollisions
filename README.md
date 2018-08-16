# CudaSoftCollisions
CudaSoftCollisions is a library that allow you to detect self collisions in a tetrahedrical mesh that is deformed over time.

## Dependencies
* [OpenCV](https://github.com/opencv/opencv)
* [CUDA](https://developer.nvidia.com/cuda-downloads)
* [VTK](https://www.vtk.org/)

## Build the project
* Get all the dependencies you need
* Clone the project :
```bash
git clone https://github.com/Mieos/CudaSoftCollisions.git
```
* Build the project :
```bash
mkdir build 
cd build
cmake ..
make
```

## How does this work ?
This library is doing fast model collision detection using two main tricks:
* First, we reduce the search space dimension by only considering tetrahedron which have a face on the surface of the polyhedra
* Second for each tetrahedron independantly (on GPU), we first check the intersection of circumsphere before checking intersection between tetrahedron (much faster) using [the method of separating axes](https://www.geometrictools.com/Documentation/MethodOfSeparatingAxes.pdf)

## How to test the library
This library provides different test in bin/tests/unit to check if the collision is working properly
* uTest_extractSurfaceVolume : extract the surface of a polyhedra and save it as data/meshes/test.ply (in the exemple it is a torus)
![Torus model](https://github.com/Mieos/CudaSoftCollisions/data/img/tor_model.png)
* uTest_extractModelFromFile : extract all the tetrahedron of a polyhedra that have a face on the surface and save it as data/meshes/test3Dmodel.ply
* uTest_simpleCollision : simple collision test
* uTest_collisionTwoTet : collision between two tetrahedra (with and without collision)
* uTest_trickyCollision : tricky configuration without collision 
* uTest_trickyCollision2 : tricky configuration 2 without collision 
* uTest_sphereCollision : check the collisions in a sphere (before and after adding collision) and save the result as data/meshes/intersectingSphere.vtk
* uTest_bunnyCollision : check the collisions in a bunny (before and after adding collision) and save the result as data/meshes/intersectingBunny.vtk

## Future Improvements
* Compute the epsilon used in the intersection test instead of using a default one
* Faster collision function but that use more gpu memory (not the priority for now)
