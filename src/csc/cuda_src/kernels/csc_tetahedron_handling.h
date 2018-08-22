#ifndef _H_KERNELS_CSC_TETRAHEDRON_HANDLING_
#define _H_KERNELS_CSC_TETRAHEDRON_HANDLING_

#include "csc_geometry_helpers.h"
#include "csc_spheres_handling.h"

//__constant__ size_t edgesID_d[12];

//Check if two tetrahedron are the same one (if they have the same indexed points)
__device__ bool checkSameTet(size_t* idArrayD, size_t* s1, size_t* s2){

   bool isSameInter;

   //Avoid too many reading (not a good idea)
   /*
   size_t id1[4], id2[4];
   id1[0] = idArrayD[4*(*s1)];
   id1[1] = idArrayD[4*(*s1)+1];
   id1[2] = idArrayD[4*(*s1)+2];
   id1[3] = idArrayD[4*(*s1)+3]; 
   id2[0] = idArrayD[4*(*s2)];
   id2[1] = idArrayD[4*(*s2)+1];
   id2[2] = idArrayD[4*(*s2)+2];
   id2[3] = idArrayD[4*(*s2)+3];
   */

   for(size_t k=0; k<4; k++){
      isSameInter=false;
      for(size_t i=0; i<4; i++){
         if(idArrayD[4*(*s1)+k]==idArrayD[4*(*s2)+i]){
         //if(id1[k]==id2[i]){
            isSameInter=true;
            //break;
         }
      }
      if(!isSameInter){
         return false;
      }
   }

   return true;
}

__device__ bool checkTetraIntersection(float*  dataPointsD, size_t* idArrayD, float* normalBuf, size_t s1, size_t s2, bool inversionTet1, bool inversionTet2){

   /*
      printf("s1 = %u\n", s1);
      printf("s2 = %u\n", s2);


      printf("DEBUG ID TET1 : %u ",idArrayD[4*s1]);
      printf(" %u ",idArrayD[4*s1+1]);
      printf(" %u ",idArrayD[4*s1+2]);
      printf(" %u\n",idArrayD[4*s1+3]);

      printf("DEBUG ID TET2 : %u ",idArrayD[4*s2]);
      printf(" %u ",idArrayD[4*s2+1]);
      printf(" %u ",idArrayD[4*s2+2]);
      printf(" %u\n",idArrayD[4*s2+3]);

      printf("TET 1 p1 = %f ",dataPointsD[3*idArrayD[4*s1]]);
      printf(" %f ",dataPointsD[3*idArrayD[4*s1]+1]);
      printf(" %f\n",dataPointsD[3*idArrayD[4*s1]+2]);

      printf("TET 1 p2 = %f ",dataPointsD[3*idArrayD[4*s1+1]]);
      printf(" %f ",dataPointsD[3*idArrayD[4*s1+1]+1]);
      printf(" %f\n",dataPointsD[3*idArrayD[4*s1+1]+2]);

      printf("TET 1 p3 = %f ",dataPointsD[3*idArrayD[4*s1+2]]);
      printf(" %f ",dataPointsD[3*idArrayD[4*s1+2]+1]);
      printf(" %f\n",dataPointsD[3*idArrayD[4*s1+2]+2]);

      printf("TET 1 p4 = %f ",dataPointsD[3*idArrayD[4*s1+3]]);
      printf(" %f ",dataPointsD[3*idArrayD[4*s1+3]+1]);
      printf(" %f\n",dataPointsD[3*idArrayD[4*s1+3]+2]);

      printf("TET 2 p1 = %f ",dataPointsD[3*idArrayD[4*s2]]);
      printf(" %f ",dataPointsD[3*idArrayD[4*s2]+1]);
      printf(" %f\n",dataPointsD[3*idArrayD[4*s2]+2]);

      printf("TET 2 p2 = %f ",dataPointsD[3*idArrayD[4*s2+1]]);
      printf(" %f ",dataPointsD[3*idArrayD[4*s2+1]+1]);
      printf(" %f\n",dataPointsD[3*idArrayD[4*s2+1]+2]);

      printf("TET 2 p3 = %f ",dataPointsD[3*idArrayD[4*s2+2]]);
      printf(" %f ",dataPointsD[3*idArrayD[4*s2+2]+1]);
      printf(" %f\n",dataPointsD[3*idArrayD[4*s2+2]+2]);

      printf("TET 2 p4 = %f ",dataPointsD[3*idArrayD[4*s2+3]]);
      printf(" %f ",dataPointsD[3*idArrayD[4*s2+3]+1]);
      printf(" %f\n",dataPointsD[3*idArrayD[4*s2+3]+2]);

    */

   //Real test
   float normalU[3];
   float pU[3];
   float pointTested[3];

   bool outTestSuccess=false;

   //Test all points of tet2 with faces of tet1
   size_t numPos=0;
   for(size_t k=0; k<4; k++){

      if(inversionTet1){
         normalU[0] = -normalBuf[4*6*s1+k*6];
         normalU[1] = -normalBuf[4*6*s1+k*6+1];
         normalU[2] = -normalBuf[4*6*s1+k*6+2]; 
      } else {
         normalU[0] = normalBuf[4*6*s1+k*6];
         normalU[1] = normalBuf[4*6*s1+k*6+1];
         normalU[2] = normalBuf[4*6*s1+k*6+2];
      }
      pU[0] = normalBuf[4*6*s1+k*6+3];
      pU[1] = normalBuf[4*6*s1+k*6+4];
      pU[2] = normalBuf[4*6*s1+k*6+5];

      size_t numNegInter = 0;

      for(size_t i=0; i<4; i++){

         size_t idP = idArrayD[4*s2+i];
         pointTested[0]=dataPointsD[3*idP];
         pointTested[1]=dataPointsD[3*idP+1];
         pointTested[2]=dataPointsD[3*idP+2];

         int res_test = checkOrientation(pointTested,normalU,pU);
         if(res_test<0){
            numNegInter++;
         } else if(res_test>0){
            numPos++;
         }

      }

      if(numNegInter==0){ //All points are at least in the positive side of one face 
         //printf("OUSIDE s1 %u\n",k);
         outTestSuccess=true;
      }


   }


   if(numPos == 0){ //The tet is inside the other
      return true;      
   }

   //Test all points of tet1 with faces of tet2
   numPos=0;
   for(size_t k=0; k<4;k++){
      if(inversionTet2){
         normalU[0] = -normalBuf[4*6*s2+k*6];
         normalU[1] = -normalBuf[4*6*s2+k*6+1];
         normalU[2] = -normalBuf[4*6*s2+k*6+2];
      } else {
         normalU[0] = normalBuf[4*6*s2+k*6];
         normalU[1] = normalBuf[4*6*s2+k*6+1];
         normalU[2] = normalBuf[4*6*s2+k*6+2];
      }
      pU[0] = normalBuf[4*6*s2+k*6+3];
      pU[1] = normalBuf[4*6*s2+k*6+4];
      pU[2] = normalBuf[4*6*s2+k*6+5];

      size_t numNegInter = 0;

      for(size_t i=0; i<4; i++){

         size_t idP = idArrayD[4*s1+i];
         pointTested[0]=dataPointsD[3*idP];
         pointTested[1]=dataPointsD[3*idP+1];
         pointTested[2]=dataPointsD[3*idP+2];

         int res_test = checkOrientation(pointTested,normalU,pU);
         if(res_test<0){
            numNegInter++;   
         } else if (res_test>0){
            numPos++;
         }

      }

      if(numNegInter==0){ //All points are at least in the positive side of one face 
         outTestSuccess=true;
      }

   }

   if(numPos == 0){ //The tet is inside the other
      return true;      
   }

   //This test has to be done only after the test concerning the inside test is done
   if(outTestSuccess){
      return false;
   }

   size_t edgesID_d[] = {
      0,1,
      1,2,
      2,0,
      0,3,
      1,3,
      2,3
   };

   //Cross product between edges
   float p1_v[3];
   float p2_v[3];
   float crossP_v[3];
   for(size_t m=0; m<6; m++){

      size_t idP1_1 = idArrayD[4*s1+edgesID_d[2*m]];
      size_t idP1_2 = idArrayD[4*s1+edgesID_d[2*m+1]];

      p1_v[0]=dataPointsD[3*idP1_2] - dataPointsD[3*idP1_1];
      p1_v[1]=dataPointsD[3*idP1_2+1] - dataPointsD[3*idP1_1+1];
      p1_v[2]=dataPointsD[3*idP1_2+2] - dataPointsD[3*idP1_1+2];

      pU[0] = dataPointsD[3*idP1_2];
      pU[1] = dataPointsD[3*idP1_2+1];
      pU[2] = dataPointsD[3*idP1_2+2];;

      for(size_t n=0; n<6; n++){

         size_t idP2_1 = idArrayD[4*s2+edgesID_d[2*n]];
         size_t idP2_2 = idArrayD[4*s2+edgesID_d[2*n+1]];

         p2_v[0]=dataPointsD[3*idP2_2] - dataPointsD[3*idP2_1];
         p2_v[1]=dataPointsD[3*idP2_2+1] - dataPointsD[3*idP2_1+1];
         p2_v[2]=dataPointsD[3*idP2_2+2] - dataPointsD[3*idP2_1+2];

         //Cross product
         crossP(p1_v,p2_v,crossP_v);

         float normCross = crossP_v[0]*crossP_v[0] + crossP_v[1]*crossP_v[1] + crossP_v[2]*crossP_v[2];

         if(normCross<EPSILON_d_2){
            //Skip
            continue;
         } else {
            //Normalize
            crossP_v[0] = crossP_v[0]/normCross;
            crossP_v[1] = crossP_v[1]/normCross;
            crossP_v[2] = crossP_v[2]/normCross;
         }

         //Test 1
         int results_inter1 = -2;

         for(size_t i=0; i<4; i++){

            size_t idP = idArrayD[4*s1+i];
            pointTested[0]=dataPointsD[3*idP];
            pointTested[1]=dataPointsD[3*idP+1];
            pointTested[2]=dataPointsD[3*idP+2];

            //Update tests
            int res_check = checkOrientation(pointTested,crossP_v,pU);
            if(res_check>0){
               if(results_inter1==-2){
                  results_inter1=1;
               } else if(results_inter1==-1){
                  results_inter1=0;
                  //break;
               } else if(results_inter1==1){
                  //OK nothing
               }
            } else if(res_check<0){
               if(results_inter1==-2){
                  results_inter1=-1;
               } else if(results_inter1==-1){
                  //OK nothing
               } else if(results_inter1==1){
                  results_inter1=0;
                  //break;
               }
            }

         }

         //If points are on one sides
         if(results_inter1==0){
            continue;
         }

         //Test2
         int results_inter2 = -2;

         for(size_t i=0; i<4; i++){

            size_t idP = idArrayD[4*s2+i];
            pointTested[0]=dataPointsD[3*idP];
            pointTested[1]=dataPointsD[3*idP+1];
            pointTested[2]=dataPointsD[3*idP+2];

            //Update tests 
            int res_check = checkOrientation(pointTested,crossP_v,pU);
            if(res_check>0){
               if(results_inter2==-2){
                  results_inter2=1;
               } else if(results_inter2==-1){
                  results_inter2=0;
                  //break;
               } else if(results_inter2==1){
                  //OK nothing
               }
            } else if(res_check<0){
               if(results_inter2==-2){
                  results_inter2=-1;
               } else if(results_inter2==-1){
                  //OK nothing
               } else if(results_inter2==1){
                  results_inter2=0;
                  //break;
               }
            }

         }

         //If points are on one sides
         if(results_inter2==0){
            continue;
         }

         //We do not want printf in the kernel
         /*
         if((results_inter1==-2)&&(results_inter2==-2)){ //All points are on the plane
            printf("BUG (important), please check the kernel : %f \n", normCross);
         }
         */

         if(results_inter1*results_inter2<0){
            return false;
         }

      }

   }

   return true;

}

__global__ void updatePlanesTet(float*  dataPointsD, size_t* idArrayD, float* normalBuf, size_t numberTets){

   size_t numTet = blockIdx.x*blockDim.x*blockDim.y +blockDim.x*threadIdx.y+threadIdx.x;

   if(numTet<numberTets){

      size_t idTetConsidered[4];
      idTetConsidered[0] = idArrayD[4*numTet];
      idTetConsidered[1] = idArrayD[4*numTet+1];
      idTetConsidered[2] = idArrayD[4*numTet+2];
      idTetConsidered[3] = idArrayD[4*numTet+3];

      size_t idOrientedTetFaces[] = {
         0,1,2,
         0,2,3,
         0,3,1,
         2,1,3
      };


      for(size_t k=0; k<4; k++){

         size_t id1,id2,id3;
         id1=idTetConsidered[idOrientedTetFaces[3*k]];
         id2=idTetConsidered[idOrientedTetFaces[3*k+1]];
         id3=idTetConsidered[idOrientedTetFaces[3*k+2]];

         //Point of the face (oriented)
         float P1[3];
         float P2[3];
         float P3[3];
         for(size_t i=0; i<3; i++){
            P1[i] = dataPointsD[3*id1+i];
            P2[i] = dataPointsD[3*id2+i];
            P3[i] = dataPointsD[3*id3+i];
         }

         //Vectors
         float v1[3], v2[3];
         for(size_t i=0; i<3; i++){
            v1[i]=P2[i]-P1[i];
            v2[i]=P3[i]-P1[i];
         }

         //Vect product
         float v1Vv2[3];
         v1Vv2[0] = v1[1]*v2[2] - v1[2]*v2[1];
         v1Vv2[1] = v1[2]*v2[0] - v1[0]*v2[2];
         v1Vv2[2] = v1[0]*v2[1] - v1[1]*v2[0];

         //Update buffer
         float normVP = v1Vv2[0]*v1Vv2[0] + v1Vv2[1]*v1Vv2[1] + v1Vv2[2]*v1Vv2[2];
         normVP = sqrt(normVP);

         //Update buffer
         normalBuf[4*6*numTet+k*6] = v1Vv2[0]/normVP; //Normal
         normalBuf[4*6*numTet+k*6+1] = v1Vv2[1]/normVP; 
         normalBuf[4*6*numTet+k*6+2] = v1Vv2[2]/normVP;
         normalBuf[4*6*numTet+k*6+3] = dataPointsD[3*id1]; //Point
         normalBuf[4*6*numTet+k*6+4] = dataPointsD[3*id1+1];
         normalBuf[4*6*numTet+k*6+5] = dataPointsD[3*id1+2];

      }

   }

}

#endif

