#ifndef _HPP_MESHSTRUCTURECOLLIDER_
#define _HPP_MESHSTRUCTURECOLLIDER_

#include <string>
#include <opencv2/core/core.hpp>

class MeshStructureCollider {

   public:

      MeshStructureCollider(const cv::Mat & dataMesh, const std::vector<std::vector<size_t> > & tetIdVector, const std::vector<size_t> & associationVectorU, const std::vector<bool> & ignoredInversionTet);

      ~MeshStructureCollider();
      bool isProperlyInitialized();
      bool updatePointsPositions(const cv::Mat & newPositions);
      bool collide(std::vector<bool> & collisionList);
      bool collideAndGetMovements(std::vector<bool> & collisionList, std::vector<std::vector<float> > & movVect);
      
      //new TODO
      bool collideAndGetMovementsAndInversions(std::vector<bool> & collisionList, std::vector<std::vector<float> > & movVect, std::vector<bool> & invertionList, std::vector<std::vector<float> > & movVectInvertion);

      //bool collideAndGetInversions(std::vector<bool> & collisionList, std::vector<bool> & inversionList);
   
   private:
      /////////////
      //Functions//
      /////////////
      bool collide();
      static bool cudaAllocation(void ** pointer_d, size_t sizeUsed, std::string errorName);
      static bool copyTodevice(void* pointer_gpu, void* pointer_cpu, size_t sizeUsed , std::string errorName);
      static bool copyTohost(void* pointer_cpu, void* pointer_gpu, size_t sizeUsed , std::string errorName);
      static bool checkGPUerrors(std::string errorName);
      //////////////
      //Parameters//
      //////////////
      bool initialized;
      bool verbose;
      size_t numPoints;
      size_t numTets;
      float* dataArrayBuff;
      float* data_d;
      size_t* tetId_d;
      float* sphereBuf_d;
      float* normalBuf_d;
      //Collisions vectors (results)
      bool* collideVectorArray;
      bool* collideVectorArray_d;
      //Inversions vectors (results)
      bool * inversionTetVectorArray;
      bool * inversionTetVectorArray_d;
      //Subdivision(loop)
      size_t numberSub;
      bool useSubdivision;
      bool* subdividedCollisionVector_d;
      std::vector<size_t> associationVector;
      
      //Spacial subdivision //Not working FIXME
      bool useSpatialSubdivision;
      size_t numberOfSpatialSubdivision;
      float* subdivisionXYZArray;
      float* subdivisionXYZArray_d;
      bool* spatialSub_d; //Maybe it has to be re-thought
      
      //Movements arrays 
      float* movementsArray;
      float* movementsArray_d;

      //Ignored inversions
      std::vector<bool> inversionsIgnored;
      float* movementsInversions;
      float* movementsInversions_d;

};

#endif
