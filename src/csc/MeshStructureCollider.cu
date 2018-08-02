#include "csc/MeshStructureCollider.hpp"
#include <iostream>

MeshStructureCollider::MeshStructureCollider(const cv::Mat & dataMesh, const std::vector<std::vector<size_t> > & tetIdVector, const std::vector<size_t> & associationVector) : initialized(false), numPoints(0){

   std::vector<float> array;
   float * arrayP = &array[0];

   if(dataMesh.type()==CV_32FC1){

      std::cout << "Type detected : CV_32FC1 : OK" << std::endl;
      this->numPoints=dataMesh.rows;

      if(dataMesh.isContinuous()){
         std::cout << "Copying data to buffer : should be fast" << std::endl;
         array.assign((float*)dataMesh.datastart, (float*)dataMesh.dataend);

      } else { 
         std::cout << "Copying data to buffer : can be slow, data not continuous" << std::endl;
         for (int i = 0; i < dataMesh.rows; ++i) {
            array.insert(array.end(), dataMesh.ptr<float>(i), dataMesh.ptr<float>(i)+dataMesh.cols);
         } 
      }

      //std::cout << "debug : " << 3*numPoints << " == " << array.size() << std::endl;
      //std::cout << "debug : " << array.at(3) << " == " << dataMesh.at<float>(1,0) << std::endl;

      size_t size_data = 3*this->numPoints*sizeof(float);
      cudaMalloc((void **) &(this->data_d), size_data);
      cudaMemcpy(this->data_d, arrayP, sizeNM, cudaMemcpyHostToDevice) ;

      std::vector<size_t> associationVectorArray;
      size_t* associationVectorArrayP = &associationVectorArray[0];
      for(size_t k=0; k<associationVector.size()){
         if(associationVector.at(k).size()!=4){
            std::cout << "Issue with association vector in collider initialisation, abort.."<< std::endl;
            break;
         }
         associationVectorArray.push_back(associationVector.at(k).at(0));
         associationVectorArray.push_back(associationVector.at(k).at(1));
         associationVectorArray.push_back(associationVector.at(k).at(2));
         associationVectorArray.push_back(associationVector.at(k).at(3));
      }




   } else { 
      std::cout << "Type not taken into account in MeshStructureCollider, please convert to CV_32FC1" << std::endl;
   }

}

MeshStructureCollider::~MeshStructureCollider(){

   cudaFree(this->data_d);

}
