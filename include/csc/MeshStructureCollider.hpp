#ifndef _HPP_MESHSTRUCTURECOLLIDER_
#define _HPP_MESHSTRUCTURECOLLIDER_

#include <string>
#include <opencv2/core/core.hpp>

class MeshStructureCollider {

   public:

      MeshStructureCollider(const cv::Mat & dataMesh, const std::vector<std::vector<size_t> > & tetIdVector, const std::vector<size_t> & associationVectorU);
      ~MeshStructureCollider();
      bool isProperlyInitialized();

   private:
      bool initialized;
      bool verbose;
      size_t numPoints;
      size_t numTets;
      float* data_d;
      float* tetId_d;
      std::vector<size_t> associationVector;

};

#endif
