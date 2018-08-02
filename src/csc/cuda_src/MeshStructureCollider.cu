#include "csc/MeshStructureCollider.hpp"
#include <iostream>

MeshStructureCollider::MeshStructureCollider(const cv::Mat & dataMesh, const std::vector<std::vector<size_t> > & tetIdVector, const std::vector<size_t> & associationVectorU) : initialized(false), numPoints(0), verbose(true){

   std::vector<float> array;
   float * arrayP = &array[0];

   if(dataMesh.type()==CV_32FC1){

      std::cout << "Type detected : CV_32FC1 : OK" << std::endl;
      this->numPoints=dataMesh.rows;

      //Copy the points
      if(dataMesh.isContinuous()){
         std::cout << "Copying data to buffer : should be fast" << std::endl;
         array.assign((float*)dataMesh.datastart, (float*)dataMesh.dataend);

      } else { 
         std::cout << "Copying data to buffer : can be slow, data not continuous" << std::endl;
         for (int i = 0; i < dataMesh.rows; ++i) {
            array.insert(array.end(), dataMesh.ptr<float>(i), dataMesh.ptr<float>(i)+dataMesh.cols);
         } 
      }

      size_t size_data = 3*this->numPoints*sizeof(float);
      cudaMalloc((void **) &(this->data_d), size_data);
      cudaMemcpy(this->data_d, arrayP, size_data, cudaMemcpyHostToDevice);

      //Copy the index of the tetrahedrons
      std::vector<size_t> tetVectorArray;
      size_t* tetVectorArrayP = &tetVectorArray[0];
      this->numTets=tetIdVector.size();
      for(size_t k=0; k<tetIdVector.size(); k++){
         
         if(tetIdVector.at(k).size()!=4){
            std::cout << "Issue with association vector in collider initialisation, abort.."<< std::endl;
            break;
         }
         
         tetVectorArray.push_back(tetIdVector.at(k).at(0));
         tetVectorArray.push_back(tetIdVector.at(k).at(1));
         tetVectorArray.push_back(tetIdVector.at(k).at(2));
         tetVectorArray.push_back(tetIdVector.at(k).at(3));
         
      }

      size_t size_tetVector = tetVectorArray.size()*sizeof(float);
      cudaMalloc((void **) &(this->tetId_d), size_tetVector);
      cudaMemcpy(this->tetId_d, tetVectorArrayP, size_tetVector, cudaMemcpyHostToDevice);

      //Copy the association vector
      this->associationVector = associationVectorU;

      //Update state
      this->initialized=true;

   } else { 
      std::cout << "Type not taken into account in MeshStructureCollider, please convert to CV_32FC1" << std::endl;
   }

}

MeshStructureCollider::~MeshStructureCollider(){

   cudaFree(this->data_d);
   cudaFree(this->tetId_d);

}

bool MeshStructureCollider::isProperlyInitialized(){
   
   return this->initialized;

}
