#ifndef _H_KERNELS_CSC_CORE_
#define _H_KERNELS_CSC_CORE_

#include "csc_tetahedron_handling.h"

//Check the intersection between one tetrahedron and all the others (kernels)
__global__ void checkForIntersectionV0(float*  dataPointsD, size_t* idArrayD, float* centerSphereB, float* normalsB, size_t numberTets, bool* intersectionVector){

   size_t numTet = blockIdx.x*blockDim.x*blockDim.y +blockDim.x*threadIdx.y+threadIdx.x;
   //printf("%u\n",subD);

   if(numTet<numberTets){

      intersectionVector[numTet]=false;

      for(size_t k=0; k<numberTets; k++){

         if(!checkSameTet(idArrayD,numTet,k)){

            if(checkSphereIntersection(centerSphereB,numTet,k)){

               if(checkTetraIntersection(dataPointsD,idArrayD, normalsB, numTet,k)){
                  intersectionVector[numTet]=true;
               }

            }

         }

      }

   }

}

//Check the intersection between one tetrahedron and a part of all the others (kernels)
__global__ void checkForIntersectionV1(float*  dataPointsD, size_t* idArrayD, float* centerSphereB, float* normalsB, size_t numberTets, bool* subIntersectionVector){

   size_t numTet = blockIdx.x*blockDim.x*blockDim.y +blockDim.x*threadIdx.y+threadIdx.x;
   size_t subD = gridDim.y;
   size_t numS = blockIdx.y;


   if(numTet<numberTets){

      subIntersectionVector[subD*numTet+numS]=false;

      for(size_t k=numS*subD; k<(numS+1)*subD; k++){

         if(k<numberTets){

            if(!checkSameTet(idArrayD,numTet,k)){

               if(checkSphereIntersection(centerSphereB,numTet,k)){

                  if(checkTetraIntersection(dataPointsD,idArrayD, normalsB, numTet,k)){
                     subIntersectionVector[subD*numTet+numS]=true;
                  }

               }

            }
         }

      }

   }

}

__global__ void reduceIntersectionVector(size_t numberTets, size_t subStep, bool* subIntersectionVector, bool* intersectionVector){

   size_t numTet = blockIdx.x*blockDim.x*blockDim.y +blockDim.x*threadIdx.y+threadIdx.x;

   if(numTet<numberTets){

      subIntersectionVector[numTet]=false;

      for(size_t k=0; k<subStep; k++){

         if(subIntersectionVector[subStep*numTet+k]){
            intersectionVector[numTet]=true;
         }

      }

   }

}

#endif