#ifndef _HPP_MESHSTRUCTUREEXTRACTOR_
#define _HPP_MESHSTRUCTUREEXTRACTOR_

#include <string>

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include <opencv2/core/core.hpp>

class MeshStructureExtractor {

   public:

      static bool extractModelFromFile(std::string fileName, cv::Mat & resultStructMatTet, std::vector<std::vector<size_t> > & tetIdVector, std::vector<size_t> & associationVectorResult, std::vector<size_t> & tetSelected);
      static bool extractModelFromMesh(const vtkSmartPointer<vtkUnstructuredGrid> & mesh3d, cv::Mat & resultStructMatTet, std::vector<std::vector<size_t>> & tetIdVector, std::vector<size_t> & associationVectorResult, std::vector<size_t> & tetSelected);

};

#endif
