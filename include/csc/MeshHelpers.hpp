#ifndef _HPP_MESHHELPERS_
#define _HPP_MESHHELPERS_

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>

class MeshHelpers {

   public:

      static bool readVolumeMeshVTK(const std::string & fileName, vtkSmartPointer<vtkUnstructuredGrid> & mesh3d);
      static bool getSurfaceOfVolMesh(const vtkSmartPointer<vtkUnstructuredGrid> & mesh3d, vtkSmartPointer<vtkPolyData> & resultPolyData);
      static bool checkDirectVectorOrientation(const double* v1, const double* v2, const double* v3);
};

#endif
