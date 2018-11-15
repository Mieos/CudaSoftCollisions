#ifndef _H_KERNELS_CSC_SPHERES_HANDLING_
#define _H_KERNELS_CSC_SPHERES_HANDLING_

#include "csc_geometry_helpers.h"

//Check the intersection between two sphere in the buffers

__device__ bool checkSphereIntersection(float* centerSphereID, size_t *s1, size_t *s2){

   float distCenter = 
      (centerSphereID[4*(*s1)] - centerSphereID[4*(*s2)])*
      (centerSphereID[4*(*s1)] - centerSphereID[4*(*s2)])+
      (centerSphereID[4*(*s1)+1] - centerSphereID[4*(*s2)+1])*
      (centerSphereID[4*(*s1)+1] - centerSphereID[4*(*s2)+1])+
      (centerSphereID[4*(*s1)+2] - centerSphereID[4*(*s2)+2])*
      (centerSphereID[4*(*s1)+2] - centerSphereID[4*(*s2)+2]);


   distCenter=sqrt(distCenter);
   float addRadius = centerSphereID[4*(*s1)+3] + centerSphereID[4*(*s2)+3]; 

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

      float ax = dataPointsD[3*id1];
      float ay = dataPointsD[3*id1+1];
      float az = dataPointsD[3*id1+2];

      float bx = dataPointsD[3*id2];
      float by = dataPointsD[3*id2+1];
      float bz = dataPointsD[3*id2+2];

      float cx = dataPointsD[3*id3];
      float cy = dataPointsD[3*id3+1];
      float cz = dataPointsD[3*id3+2];

      float dx = dataPointsD[3*id4];
      float dy = dataPointsD[3*id4+1];
      float dz = dataPointsD[3*id4+2];

      float Mal[16];
      Mal[0] = ax; Mal[1] = ay; Mal[2] = az; Mal[3] = 1;
      Mal[4] = bx; Mal[5] = by; Mal[6] = bz; Mal[7] = 1;
      Mal[8] = cx; Mal[9] = cy; Mal[10] = cz; Mal[11] = 1;
      Mal[12] = dx; Mal[13] = dy; Mal[14] = dz; Mal[15] = 1;
      float Dal = cDet4(Mal);

      float Mx[16]; 
      Mx[0] = ax*ax + ay*ay + az*az; Mx[1] = ay; Mx[2] = az; Mx[3] = 1;
      Mx[4] = bx*bx + by*by + bz*bz; Mx[5] = by; Mx[6] = bz; Mx[7] = 1;
      Mx[8] = cx*cx + cy*cy + cz*cz; Mx[9] = cy; Mx[10] = cz; Mx[11] = 1;
      Mx[12] = dx*dx + dy*dy + dz*dz; Mx[13] = dy; Mx[14] = dz; Mx[15] = 1;
      float Dx = cDet4(Mx);

      float My[16]; 
      My[0] = ax*ax + ay*ay + az*az; My[1] = ax; My[2] = az; My[3] = 1;
      My[4] = bx*bx + by*by + bz*bz; My[5] = bx; My[6] = bz; My[7] = 1;
      My[8] = cx*cx + cy*cy + cz*cz; My[9] = cx; My[10] = cz; My[11] = 1;
      My[12] = dx*dx + dy*dy + dz*dz; My[13] = dx; My[14] = dz; My[15] = 1;
      float Dy = cDet4(My);

      float Mz[16]; 
      Mz[0] = ax*ax + ay*ay + az*az; Mz[1] = ax; Mz[2] = ay; Mz[3] = 1;
      Mz[4] = bx*bx + by*by + bz*bz; Mz[5] = bx; Mz[6] = by; Mz[7] = 1;
      Mz[8] = cx*cx + cy*cy + cz*cz; Mz[9] = cx; Mz[10] = cy; Mz[11] = 1;
      Mz[12] = dx*dx + dy*dy + dz*dz; Mz[13] = dx; Mz[14] = dy; Mz[15] = 1;
      float Dz = cDet4(Mz);

      //Coordinate
      //sphereBuffer[4*numTet] = (dataPointsD[3*id1] + dataPointsD[3*id2] + dataPointsD[3*id3] + dataPointsD[3*id4])/4.0;
      //sphereBuffer[4*numTet+1] = (dataPointsD[3*id1+1] + dataPointsD[3*id2+1] + dataPointsD[3*id3+1] + dataPointsD[3*id4+1])/4.0;
      //sphereBuffer[4*numTet+2] = (dataPointsD[3*id1+2] + dataPointsD[3*id2+2] + dataPointsD[3*id3+2] + dataPointsD[3*id4+2])/4.0;
      sphereBuffer[4*numTet] = Dx/2*Dal;
      sphereBuffer[4*numTet+1] = Dx/2*Dal;
      sphereBuffer[4*numTet+2] = Dx/2*Dal;

      //Radius
      float radius = 
         (sphereBuffer[4*numTet] - dataPointsD[3*id1])*(sphereBuffer[4*numTet] - dataPointsD[3*id1])+
         (sphereBuffer[4*numTet+1] - dataPointsD[3*id1+1])*(sphereBuffer[4*numTet+1] - dataPointsD[3*id1+1])+
         (sphereBuffer[4*numTet+2] - dataPointsD[3*id1+2])*(sphereBuffer[4*numTet+2] - dataPointsD[3*id1+2]);
      sphereBuffer[4*numTet+3] = sqrt(radius);

   } 

}

#endif


