#include "csc/MeshStructureCollider.hpp"
#include <iostream>
//#include <string>

MeshStructureCollider::MeshStructureCollider(const cv::Mat & dataMesh, const std::vector<std::vector<size_t> > & tetIdVector, const std::vector<size_t> & associationVectorU) : initialized(false), numPoints(0), verbose(true){

   std::vector<float> array;

   if(dataMesh.type()==CV_32FC1){

      //std::cout << "Type detected : CV_32FC1 : OK" << std::endl;

      this->numPoints=dataMesh.rows;

      if(dataMesh.cols==3){

         this->dataArrayBuff = new float[3*numPoints];
         
         //Copy the points
         for(size_t i=0; i<dataMesh.rows; i++){
            for(size_t j=0; j<3; j++){
               this->dataArrayBuff[3*i+j]=dataMesh.at<float>(i,j);   
            }
         }

         //Copy points to gpu
         size_t size_data = 3*this->numPoints*sizeof(float);
         cudaMalloc((void **) &(this->data_d), size_data);
         cudaMemcpy(this->data_d, this->dataArrayBuff, size_data, cudaMemcpyHostToDevice);

         //Copy the index of the tetrahedrons
         this->numTets=tetIdVector.size();
         float* tetVectorPointer = new float[4*tetIdVector.size()];

         for(size_t k=0; k<tetIdVector.size(); k++){

            if(tetIdVector.at(k).size()!=4){
               std::cout << "Issue with association vector in collider initialisation, abort.."<< std::endl;
               break;
            }

            tetVectorPointer[4*k] = tetIdVector.at(k).at(0);
            tetVectorPointer[4*k+1] = tetIdVector.at(k).at(1);
            tetVectorPointer[4*k+2] = tetIdVector.at(k).at(2);
            tetVectorPointer[4*k+3] = tetIdVector.at(k).at(3);

         }

         size_t size_tetVector = 4*this->numTets*sizeof(float);
         cudaMalloc((void **) &(this->tetId_d), size_tetVector);
         cudaMemcpy(this->tetId_d, tetVectorPointer, size_tetVector, cudaMemcpyHostToDevice);

         //Copy the association vector
         this->associationVector = associationVectorU;

         //Update state
         this->initialized=true;

         delete[] tetVectorPointer;

      } else {

         std::cout << "Wrong number of collumns" << std::endl;

      }

   } else { 
      std::cout << "Type not taken into account in MeshStructureCollider, please convert to CV_32FC1" << std::endl;
   }

}

MeshStructureCollider::~MeshStructureCollider(){

   delete[] this->dataArrayBuff;
   cudaFree(this->data_d);
   cudaFree(this->tetId_d);

}

bool MeshStructureCollider::isProperlyInitialized(){

   return this->initialized;

}

bool MeshStructureCollider::updatePointsPositions(const cv::Mat & newPositions){

   //Not init
   if(!this->isProperlyInitialized()){
      return false;
   }

   //Bad number of cols
   if(newPositions.cols!=3){
      return false;
   }

   //Test type
   if(newPositions.type()==CV_32FC1){
   
      for(size_t k=0; k<this->associationVector.size(); k++){
         this->dataArrayBuff[3*k]=newPositions.at<float>(associationVector.at(k),0); 
         this->dataArrayBuff[3*k+1]=newPositions.at<float>(associationVector.at(k),1); 
         this->dataArrayBuff[3*k+2]=newPositions.at<float>(associationVector.at(k),2); 
      }

   } else if(newPositions.type()==CV_64FC1){
   
      for(size_t k=0; k<this->associationVector.size(); k++){
         this->dataArrayBuff[3*k]=float(newPositions.at<double>(associationVector.at(k),0)); 
         this->dataArrayBuff[3*k+1]=float(newPositions.at<double>(associationVector.at(k),1)); 
         this->dataArrayBuff[3*k+2]=float(newPositions.at<double>(associationVector.at(k),2)); 
      }

   } 

   //Copy to buffer
   size_t size_data = 3*this->numPoints*sizeof(float);
   cudaMemcpy(this->data_d, this->dataArrayBuff, size_data, cudaMemcpyHostToDevice);

   return true;

}
