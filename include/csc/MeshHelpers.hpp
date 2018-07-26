#ifndef _HPP_MESHHELPERS_
#define _HPP_MESHHELPERS_

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

class MeshHelpers {

   public:

      static bool readVolumeMeshVTK(const std::string & fileName, vtkSmartPointer<vtkUnstructuredGrid> & mesh3d);
      static bool getSurfaceOfVolMesh(const vtkSmartPointer<vtkUnstructuredGrid> & mesh3d);

};

#endif
