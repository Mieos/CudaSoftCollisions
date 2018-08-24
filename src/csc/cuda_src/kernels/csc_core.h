#ifndef _H_KERNELS_CSC_CORE_
#define _H_KERNELS_CSC_CORE_

#include "csc_tetahedron_handling.h"

//Check the intersection between one tetrahedron and all the others (kernels)
__global__ void checkForIntersectionV0_withmovements(float*  dataPointsD, size_t* idArrayD, float* centerSphereB, float* normalsB, size_t numberTets, bool* intersectionVector, bool* inversionVector, float* movementArray){

   size_t numTet = blockIdx.x*blockDim.x*blockDim.y +blockDim.x*threadIdx.y+threadIdx.x;
   //printf("%u\n",subD);

   if(numTet<numberTets){

      intersectionVector[numTet]=false;
      /*
      float dx,dy,dz;
      float normV;
      */

      //Inversions
      bool invert1=inversionVector[numTet];
      bool invert2;

      //Update movement
      movementArray[3*numTet]=0;
      movementArray[3*numTet+1]=0;
      movementArray[3*numTet+2]=0;

      size_t numberColisions = 0;

      for(size_t k=0; k<numberTets; k++){

         if(!checkSameTet(idArrayD,&numTet,&k)){

            if(checkSphereIntersection(centerSphereB,&numTet,&k)){

               numberColisions++;

               if(inversionVector[k]){
                  invert2=true;
               } else {
                  invert2=false;
               }

               //intersectionVector[numTet]=true;
               if(checkTetraIntersection(dataPointsD,idArrayD, normalsB, numTet, k, invert1, invert2)){
                  intersectionVector[numTet]=true;

                  /*
                  dx = centerSphereB[4*numTet] - centerSphereB[4*k];
                  dy = centerSphereB[4*numTet+1] - centerSphereB[4*k+1];
                  dz = centerSphereB[4*numTet+2] - centerSphereB[4*k+2];

                  normV=dx*dx + dy*dy + dz*dz;
                  normV=sqrt(normV);

                  dx=dx/normV;
                  dy=dy/normV;
                  dz=dz/normV;

                  normV = 0.5*(centerSphereB[4*numTet+3] + centerSphereB[4*k+3] - normV);

                  movementArray[3*numTet]=movementArray[3*numTet] + normV*dx;
                  movementArray[3*numTet+1]=movementArray[3*numTet] + normV*dy;
                  movementArray[3*numTet+2]=movementArray[3*numTet] + normV*dz;
                  */

               }

            }

         }

      }

      if(intersectionVector[numTet]){

         float norm = 0.5*centerSphereB[4*numTet+3] ;

         if (invert1){
            movementArray[3*numTet] = norm * normalsB[4*6*numTet];
            movementArray[3*numTet+1] = norm * normalsB[4*6*numTet+1];
            movementArray[3*numTet+2] = norm * normalsB[4*6*numTet+2];
         } else {
            movementArray[3*numTet] = - norm * normalsB[4*6*numTet];
            movementArray[3*numTet+1] = -norm * normalsB[4*6*numTet+1];
            movementArray[3*numTet+2] = -norm * normalsB[4*6*numTet+2];
         }
         /*
         movementArray[3*numTet]=movementArray[3*numTet]/numberColisions;
         movementArray[3*numTet+1]=movementArray[3*numTet+1]/numberColisions;
         movementArray[3*numTet+2]=movementArray[3*numTet+2]/numberColisions;
         */
      }

   }

}

//Check the intersection between one tetrahedron and all the others (kernels)
__global__ void checkForIntersectionV0(float*  dataPointsD, size_t* idArrayD, float* centerSphereB, float* normalsB, size_t numberTets, bool* intersectionVector){

   size_t numTet = blockIdx.x*blockDim.x*blockDim.y +blockDim.x*threadIdx.y+threadIdx.x;
   //printf("%u\n",subD);

   if(numTet<numberTets){

      intersectionVector[numTet]=false;

      for(size_t k=0; k<numberTets; k++){

         if(!checkSameTet(idArrayD,&numTet,&k)){

            if(checkSphereIntersection(centerSphereB,&numTet,&k)){

               //intersectionVector[numTet]=true;
               if(checkTetraIntersection(dataPointsD,idArrayD, normalsB, numTet,k,false, false)){
                  intersectionVector[numTet]=true;
               }

            }

         }

      }

   }

}

//Check the intersection between one tetrahedron and a part of all the others (kernels)
__global__ void checkForIntersectionV1(float*  dataPointsD, size_t* idArrayD, float* centerSphereB, float* normalsB, size_t numberTets, size_t size_loop, bool* subIntersectionVector){

   size_t numTet = blockIdx.x*blockDim.x*blockDim.y +blockDim.x*threadIdx.y+threadIdx.x;
   size_t subD = gridDim.y;
   size_t numS = blockIdx.y;

   //size_t size_loop = size_t(float(numberTets)/float(subD));

   bool foundInter=false;

   if(numTet<numberTets){

      subIntersectionVector[subD*numTet+numS]=false;

      //V1
      size_t beginLoopK=numS*size_loop;
      size_t endLoopK;
      if(numS==subD){
         endLoopK=numberTets;
      } else {
         endLoopK=(numS+1)*size_loop;
      }

      /*
         size_t debugValue = endLoopK-beginLoopK;
         if(debugValue>696){
         printf("%lu\n",debugValue);
         }
         */

      for(size_t k=beginLoopK; k<endLoopK; k++){

         if(!checkSameTet(idArrayD,&numTet,&k)){

            /*

               if(!foundInter){
               subIntersectionVector[subD*numTet+numS]=true;
               foundInter=true;
               }

*/

            if(checkSphereIntersection(centerSphereB,&numTet,&k)){

               /*
                  if(!foundInter){
                  subIntersectionVector[subD*numTet+numS]=true;
                  foundInter=true;
                  }
                  */


               if(checkTetraIntersection(dataPointsD,idArrayD, normalsB, numTet,k,false,false)){
                  if(!foundInter){
                     subIntersectionVector[subD*numTet+numS]=true;
                     foundInter=true;
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

__global__ void checkTetOrientations(float*  dataPointsD, size_t* idArrayD, size_t numberTets, bool* intersectionVector){

   size_t numTet = blockIdx.x*blockDim.x*blockDim.y +blockDim.x*threadIdx.y+threadIdx.x;

   if(numTet<numberTets){

      intersectionVector[numTet]=false;

      size_t idTetConsidered[4];
      idTetConsidered[0] = idArrayD[4*numTet];
      idTetConsidered[1] = idArrayD[4*numTet+1];
      idTetConsidered[2] = idArrayD[4*numTet+2];
      idTetConsidered[3] = idArrayD[4*numTet+3];

      //Points of oriented surface + last point
      size_t idOrientedTetFaces[] = {
         0,1,2,3,
         0,2,3,1,
         0,3,1,2,
         2,1,3,0
      };

      for(size_t k=0; k<4; k++){

         //Get id faces
         size_t id1,id2,id3, idOut;
         id1=idTetConsidered[idOrientedTetFaces[4*k]];
         id2=idTetConsidered[idOrientedTetFaces[4*k+1]];
         id3=idTetConsidered[idOrientedTetFaces[4*k+2]];
         idOut=idTetConsidered[idOrientedTetFaces[4*k+3]];

         //The points
         float P1[3];
         float P2[3];
         float P3[3];
         float POut[3];

         for(size_t i=0; i<3; i++){
            P1[i] = dataPointsD[3*id1+i];
            P2[i] = dataPointsD[3*id2+i];
            P3[i] = dataPointsD[3*id3+i];
            POut[i] = dataPointsD[3*idOut+i];
         }

         //The vector
         float v1[3], v2[3], v3[3];
         for(size_t i=0; i<3; i++){
            v1[i]=P2[i]-P1[i];
            v2[i]=P3[i]-P1[i];
            v3[i]=POut[i]-P1[i];
         }

         float v1Vv2[3];
         v1Vv2[0] = v1[1]*v2[2] - v1[2]*v2[1];
         v1Vv2[1] = v1[2]*v2[0] - v1[0]*v2[2];
         v1Vv2[2] = v1[0]*v2[1] - v1[1]*v2[0];

         float dotP = v1Vv2[0]*v3[0] + v1Vv2[1]*v3[1] + v1Vv2[2]*v3[2];

         if(dotP>0){
            intersectionVector[numTet]=true;
         }

      }

   }

}

__global__ void updateSpatialSub(float*  dataPointsD, size_t* idArrayD, size_t numberTets, float* espaceDistribustion, size_t numberOfSpatialSubdivision, bool* spatialSub){

   size_t numTet = blockIdx.x*blockDim.x*blockDim.y +blockDim.x*threadIdx.y+threadIdx.x;

   if(numTet<numberTets){

      //Get the ids
      size_t id1 = idArrayD[4*numTet];
      size_t id2 = idArrayD[4*numTet+1];
      size_t id3 = idArrayD[4*numTet+2];
      size_t id4 = idArrayD[4*numTet+3];

      //Get the points
      float p1[3], p2[3], p3[3], p4[3];
      //p1
      p1[0] = dataPointsD[3*id1];
      p1[1] = dataPointsD[3*id1+1];
      p1[2] = dataPointsD[3*id1+2];
      //p2
      p2[0] = dataPointsD[3*id2];
      p2[1] = dataPointsD[3*id2+1];
      p2[2] = dataPointsD[3*id2+2];
      //p3
      p3[0] = dataPointsD[3*id3];
      p3[1] = dataPointsD[3*id3+1];
      p3[2] = dataPointsD[3*id3+2];
      //p4
      p4[0] = dataPointsD[3*id4];
      p4[1] = dataPointsD[3*id4+1];
      p4[2] = dataPointsD[3*id4+2];
      
      bool foundInitXYZ[3];
      foundInitXYZ[0]=false;
      foundInitXYZ[1]=false;
      foundInitXYZ[2]=false;
      
      size_t beginXYZ[3];
      size_t endXYZ[3];

      for(size_t k=0; k<numberOfSpatialSubdivision; k++){
      
         for(size_t i=0; i<3; i++){
        
            //P1
            if((espaceDistribustion[i*numberOfSpatialSubdivision+k] <= p1[i]) &&
               (espaceDistribustion[i*numberOfSpatialSubdivision+k+1] >= p1[i])){
               
               if(!foundInitXYZ[i]){
                 beginXYZ[i]=k;
                 foundInitXYZ[i]=true;
               }
               endXYZ[i]=k;

            }
   
            //P2
            if((espaceDistribustion[i*numberOfSpatialSubdivision+k] <= p2[i]) &&
               (espaceDistribustion[i*numberOfSpatialSubdivision+k+1] >= p2[i])){
               
               if(!foundInitXYZ[i]){
                 beginXYZ[i]=k;
                 foundInitXYZ[i]=true;
               }
               endXYZ[i]=k;

            }


            //P3
            if((espaceDistribustion[i*numberOfSpatialSubdivision+k] <= p3[i]) &&
               (espaceDistribustion[i*numberOfSpatialSubdivision+k+1] >= p3[i])){
               
               if(!foundInitXYZ[i]){
                 beginXYZ[i]=k;
                 foundInitXYZ[i]=true;
               }
               endXYZ[i]=k;

            }


            //P4
            if((espaceDistribustion[i*numberOfSpatialSubdivision+k] <= p4[i]) &&
               (espaceDistribustion[i*numberOfSpatialSubdivision+k+1] >= p4[i])){
               
               if(!foundInitXYZ[i]){
                 beginXYZ[i]=k;
                 foundInitXYZ[i]=true;
               }
               endXYZ[i]=k;

            }

         }

      }

      size_t ySkip = numberOfSpatialSubdivision;
      size_t zSkip = numberOfSpatialSubdivision*numberOfSpatialSubdivision;
      size_t tSkip = numberOfSpatialSubdivision*ySkip;

      for(size_t i=beginXYZ[0]; i<endXYZ[0]; i++){
      
         for(size_t j=beginXYZ[1]; j<endXYZ[1]; j++){
         
            for(size_t k=beginXYZ[2]; k<endXYZ[2]; k++){
            
               spatialSub[numTet*tSkip+i*zSkip+j*ySkip+k]=true;

            }

         }

      }
      

   }

}
#endif
