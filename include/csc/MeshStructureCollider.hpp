#ifndef _HPP_MESHSTRUCTURECOLLIDER_
#define _HPP_MESHSTRUCTURECOLLIDER_

#include <string>
#include <opencv2/core/core.hpp>

class MeshStructureCollider {

   public:

      MeshStructureCollider(const cv::Mat & dataMesh, const std::vector<std::vector<size_t> > & tetIdVector, const std::vector<size_t> & associationVectorU);
      ~MeshStructureCollider();
      bool isProperlyInitialized();
      bool updatePointsPositions(const cv::Mat & newPositions);
      bool collide(std::vector<bool> & collisionList);
   private:
      bool initialized;
      bool verbose;
      size_t numPoints;
      size_t numTets;
      float* dataArrayBuff;
      float* data_d;
      size_t* tetId_d;
      float* sphereBuf_d;
      bool* collideVectorArray;
      bool* collideVectorArray_d;
      std::vector<size_t> associationVector;

};

#endif
