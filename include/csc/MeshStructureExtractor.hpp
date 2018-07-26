#ifndef _HPP_MESHSTRUCTUREEXTRACTOR_
#define _HPP_MESHSTRUCTUREEXTRACTOR_

#include <string>

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

class MeshStructureExtractor {

   public:

      static bool extractModelFromFile(std::string fileName);
      static bool extractModelFromMesh(const vtkSmartPointer<vtkUnstructuredGrid> & mesh3d);
};

#endif
