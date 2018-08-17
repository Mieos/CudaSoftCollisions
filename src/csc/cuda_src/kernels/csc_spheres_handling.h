#ifndef _H_KERNELS_CSC_SPHERES_HANDLING_
#define _H_KERNELS_CSC_SPHERES_HANDLING_

//Check the intersection between two sphere in the buffers

__device__ bool checkSphereIntersection(float* centerSphereID, size_t s1, size_t s2){

   float distCenter = 
      (centerSphereID[4*s1] - centerSphereID[4*s2])*
      (centerSphereID[4*s1] - centerSphereID[4*s2])+
      (centerSphereID[4*s1+1] - centerSphereID[4*s2+1])*
      (centerSphereID[4*s1+1] - centerSphereID[4*s2+1])+
      (centerSphereID[4*s1+2] - centerSphereID[4*s2+2])*
      (centerSphereID[4*s1+2] - centerSphereID[4*s2+2]);


   distCenter=sqrt(distCenter);
   float addRadius = centerSphereID[4*s1+3] + centerSphereID[4*s2+3]; 

   if(distCenter<addRadius){
      return true;
   } else {
      return false;
   }

} 

__global__ void updateCenterSphere(float*  dataPointsD, size_t* idArrayD, float* sphereBuffer, size_t numberTets){

   size_t numTet = blockIdx.x*blockDim.x*blockDim.y +blockDim.x*threadIdx.y+threadIdx.x;

   if(numTet<numberTets){

      size_t id1 = idArrayD[4*numTet];
      size_t id2 = idArrayD[4*numTet+1];
      size_t id3 = idArrayD[4*numTet+2];
      size_t id4 = idArrayD[4*numTet+3];


      //Coordinate
      sphereBuffer[4*numTet] = (dataPointsD[3*id1] + dataPointsD[3*id2] + dataPointsD[3*id3] + dataPointsD[3*id4])/4.0;
      sphereBuffer[4*numTet+1] = (dataPointsD[3*id1+1] + dataPointsD[3*id2+1] + dataPointsD[3*id3+1] + dataPointsD[3*id4+1])/4.0;
      sphereBuffer[4*numTet+2] = (dataPointsD[3*id1+2] + dataPointsD[3*id2+2] + dataPointsD[3*id3+2] + dataPointsD[3*id4+2])/4.0;

      //Radius
      float radius = 
         (sphereBuffer[4*numTet] - dataPointsD[3*id1])*(sphereBuffer[4*numTet] - dataPointsD[3*id1])+
         (sphereBuffer[4*numTet+1] - dataPointsD[3*id1+1])*(sphereBuffer[4*numTet+1] - dataPointsD[3*id1+1])+
         (sphereBuffer[4*numTet+2] - dataPointsD[3*id1+2])*(sphereBuffer[4*numTet+2] - dataPointsD[3*id1+2]);
      sphereBuffer[4*numTet+3] = sqrt(radius);

   } 

}

#endif


