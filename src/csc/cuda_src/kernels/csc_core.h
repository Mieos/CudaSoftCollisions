#ifndef _H_KERNELS_CSC_CORE_
#define _H_KERNELS_CSC_CORE_

#include "csc_tetahedron_handling.h"

//Check the intersection between one tetrahedron and all the others (kernels)
__global__ void checkForIntersection(float*  dataPointsD, size_t* idArrayD, float* centerSphereB, float* normalsB, size_t numberTets, bool* intersectionVector){

   size_t numTet = blockIdx.x*blockDim.x*blockDim.y +blockDim.x*threadIdx.y+threadIdx.x;

   if(numTet<numberTets){

      intersectionVector[numTet]=false;

      if(numTet>0){ //We have a "small" issue if numTet==0 in the next loop

         //First part
         for(size_t k=0; k<numTet; k++){

            if(!checkSameTet(idArrayD,numTet,k)){

               if(checkSphereIntersection(centerSphereB,numTet,k)){
                  if(checkTetraIntersection(dataPointsD,idArrayD, normalsB, numTet,k)){
                     intersectionVector[numTet]=true;
                     //break;
                  }

               }

            }
         }

      }

      //Second part
      if(!intersectionVector[numTet]){
         for(size_t k=numTet+1; k<numberTets; k++){ 


            if(!checkSameTet(idArrayD,numTet,k)){
               if(checkSphereIntersection(centerSphereB,numTet,k)){ 
                  if(checkTetraIntersection(dataPointsD,idArrayD,normalsB, numTet,k)){
                     intersectionVector[numTet]=true;
                     //break;
                  }
               } 
            }

         }
      }
      

   }

}

#endif
