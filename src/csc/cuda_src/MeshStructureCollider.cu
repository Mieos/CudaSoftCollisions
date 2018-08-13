#include "csc/MeshStructureCollider.hpp"
#include <iostream>
#include <math.h> 

#define EPSILON_d 0.001f

__constant__ size_t edgesID_d[12];

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

__device__ int checkOrientation(float* Ptest, float* Nm, float* Pu){

   float vt2[3];
   for(size_t k=0; k<3; k++){
      vt2[k]=Ptest[k]-Pu[k];
   }


   float dotP = Nm[0]*vt2[0] + Nm[1]*vt2[1] + Nm[2]*vt2[2];


   if(abs(dotP)<EPSILON_d){ 
      return 0;
   } else {
      if(dotP>0){
         return 1;
      } else {
         return -1;
      }
   }

}

__device__ void crossP(float* v1, float* v2, float* res){

   res[0] = v1[1]*v2[2] - v1[2]*v2[1];
   res[1] = v1[2]*v2[0] - v1[0]*v2[2];
   res[2] = v1[0]*v2[1] - v1[1]*v2[0];

}

__device__ bool checkSameTet(size_t* idArrayD, size_t s1, size_t s2){

   bool isSame=true;
   for(size_t k=0; k<4; k++){
      bool isSameInter=false;
      for(size_t i=0; i<4; i++){
         if(idArrayD[4*s1+k]==idArrayD[4*s2+i]){
            isSameInter=true;
            break;
         }
      }
      if(!isSameInter){
         isSame=false;
         break;
      }
   }

   return isSame;
}

__device__ bool checkTetraIntersection(float*  dataPointsD, size_t* idArrayD, float* normalBuf, size_t s1, size_t s2){

   //printf("OK 0 !\n");

   //printf("s1 = %u\n", s1);
   //printf("s2 = %u\n", s2);

   /*
   printf ("TET 1 || %u , ", idArrayD[4*s1]);
   printf (" %u , ", idArrayD[4*s1+1]);
   printf (" %u , ", idArrayD[4*s1+2]);
   printf (" %u \n", idArrayD[4*s1+3]);
 
   */
   /*
   if(s2==1109){
      printf ("TET 2 || %u , ", idArrayD[4*s2]);
      printf (" %u , ", idArrayD[4*s2+1]);
      printf (" %u , ", idArrayD[4*s2+2]);
      printf (" %u \n", idArrayD[4*s2+3]);   
   }
   */

   /*
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

   //printf("Test init tetra\n");
   //printf("Tet1 : %u , %u , %u, %u \n",idArrayD[4*s1], idArrayD[4*s1+1], idArrayD[4*s1+2], idArrayD[4*s1+3]);
   //printf("Tet2 : %u , %u , %u, %u \n",idArrayD[4*s2], idArrayD[4*s2+1], idArrayD[4*s2+2], idArrayD[4*s2+3]);

   //Real test
   float normalU[3];
   float pU[3];
   float pointTested[3];

   bool outTestSuccess=false;

   //Test all points of tet2 with faces of tet1
   size_t numPos=0;
   for(size_t k=0; k<4;k++){
      normalU[0] = normalBuf[4*6*s1+k*6];
      normalU[1] = normalBuf[4*6*s1+k*6+1];
      normalU[2] = normalBuf[4*6*s1+k*6+2];
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

      /*
      printf("DEBUG nX = %f",normalU[0]);
      printf(" , nY = %f",normalU[1]);
      printf(" , nZ = %f\n",normalU[2]);
      */

      if(numNegInter==0){ //All points are at least in the positive side of one face 
         outTestSuccess=true;
      }

   }

   if(numPos == 0){ //The tet is inside the other
      //printf("INSIDE 0 \n");
      return true;      
   }

   //Test all points of tet1 with faces of tet2
   numPos=0;
   for(size_t k=0; k<4;k++){
      normalU[0] = normalBuf[4*6*s2+k*6];
      normalU[1] = normalBuf[4*6*s2+k*6+1];
      normalU[2] = normalBuf[4*6*s2+k*6+2];
      pU[0] = normalBuf[4*6*s2+k*6+3];
      pU[1] = normalBuf[4*6*s2+k*6+4];
      pU[2] = normalBuf[4*6*s2+k*6+5];

      size_t numNegInter = 0;
     
      //printf("Face = %u\n", k);
      for(size_t i=0; i<4; i++){

         size_t idP = idArrayD[4*s1+i];
         pointTested[0]=dataPointsD[3*idP];
         pointTested[1]=dataPointsD[3*idP+1];
         pointTested[2]=dataPointsD[3*idP+2];

         int res_test = checkOrientation(pointTested,normalU,pU);
         if(res_test<0){
            //printf("Num point Neg = %u\n", i);
            numNegInter++;   
         } else if (res_test>0){
            //printf("Num point Pos = %u\n", i);
            numPos++;
         }

      }

      if(numNegInter==0){ //All points are at least in the positive side of one face 
         //printf("numPos = %u\n",numPos);
         //printf("face = %u\n",k);
         outTestSuccess=true;
      }

   }

   if(numPos == 0){ //The tet is inside the other
      //printf("INSIDE\n");
      return true;      
   }

   //This test has to be done only after the test concerning the inside test is done
   if(outTestSuccess){
      //printf("OUSIDE \n");
      return false;
   }

   //printf("WTF! %u\n",s2);

   size_t edgesID_d[] = {
         0,1,
         1,2,
         2,0,
         0,3,
         1,3,
         2,3
   };

   //TODO test that part 
   //Cross product

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
                  break;
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
                  break;
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
                  break;
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
                  break;
               }
            }

         }

         //If points are on one sides
         if(results_inter2==0){
            continue;
         }

         if((results_inter1==-2)||(results_inter2==-2)){
            //printf("BUG (important)\n");//FIXME All the points are ON the plane WTF
            continue;
         }

         //printf("At least, there is a check\n");

         if(results_inter1*results_inter2<0){
            return false;
         }

      }

   }

   //printf("WTFn2\n");

   return true;

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

      //printf("numTet = %u\n",numTet);

      for(size_t k=0; k<4; k++){
      //for(size_t k=3; k<4; k++){
      
         /*
         printf("ID1 = %u\n",idTetConsidered[idOrientedTetFaces[3*k]]);
         printf("ID2 = %u\n",idTetConsidered[idOrientedTetFaces[3*k+1]]);
         printf("ID3 = %u\n",idTetConsidered[idOrientedTetFaces[3*k+2]]);
         */

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

__global__ void checkForIntersection(float*  dataPointsD, size_t* idArrayD, float* centerSphereB, float* normalsB, size_t numberTets, bool* intersectionVector){

   size_t numTet = blockIdx.x*blockDim.x*blockDim.y +blockDim.x*threadIdx.y+threadIdx.x;


   if(numTet<numberTets){

      intersectionVector[numTet]=false;

      //if(numTet ==0){ //FIXME
      //if(numTet ==1066){ //FIXME

      
      if(numTet>0){ //We have a "small" issue if numTet==0 in the next loop

         //First part
         for(size_t k=0; k<numTet; k++){

            if(!checkSameTet(idArrayD,numTet,k)){

               if(checkSphereIntersection(centerSphereB,numTet,k)){
                  if(checkTetraIntersection(dataPointsD,idArrayD, normalsB, numTet,k)){
                     intersectionVector[numTet]=true;
                     //printf("DEBUG = %u\n",k);
                     break;
                  }

               }

            }

         }

      }
       //FIXME

      //Second part
      if(!intersectionVector[numTet]){
         for(size_t k=numTet+1; k<numberTets; k++){ 

            //if(k==1109){ //FIXME

            if(!checkSameTet(idArrayD,numTet,k)){

               if(checkSphereIntersection(centerSphereB,numTet,k)){ 
                  if(checkTetraIntersection(dataPointsD,idArrayD,normalsB, numTet,k)){
                     intersectionVector[numTet]=true;
                     //printf("DEBUG = %u\n",k);
                     break;
                  }
               } 
            }

            //} //END FIXME
         }
      }

      //} //FIXME

      /*
      if(intersectionVector[numTet]){
         printf("Intersecting : %u\n", numTet);
      } else {
         printf("NOT intersecting : %u\n", numTet);
      }
      */
      

   }

}

__global__ void debugTestKernel(size_t* tetIndex, size_t numberTets){

   size_t numTet = blockIdx.x*blockDim.x*blockDim.y +blockDim.x*threadIdx.y+threadIdx.x;

   if(numTet<numberTets){

      size_t i1, i2, i3, i4;
      i1 = tetIndex[4*numTet];
      i2 = tetIndex[4*numTet+1];
      i3 = tetIndex[4*numTet+2];
      i4 = tetIndex[4*numTet+3];
      //printf("Numtet = %u || %u , %u , %u , %u \n", numTet, i1, i2, i3, i4);
      printf(" i1 = %u\n", i1);
      printf(" i2 = %u\n", i2);
      printf(" i3 = %u\n", i3);
      printf(" i4 = %u\n", i4);

      /*
         boolV[numTet]=false;

         if(boolV[numTet]){
         printf("DEBUG\n");
         }
       */

   }

}

MeshStructureCollider::MeshStructureCollider(const cv::Mat & dataMesh, const std::vector<std::vector<size_t> > & tetIdVector, const std::vector<size_t> & associationVectorU) : initialized(false), numPoints(0), verbose(true){

   //Cuda error
   cudaError_t error;

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
         error=cudaMalloc((void **) &(this->data_d), size_data);
         if (error != cudaSuccess) {
            fprintf(stderr, "Issue while initializing data buffer\n");
            fprintf(stderr, "cudaMalloc failed: %s\n", cudaGetErrorString(error));
         }
         error=cudaMemcpy(this->data_d, this->dataArrayBuff, size_data, cudaMemcpyHostToDevice);
         if (error != cudaSuccess) {
            fprintf(stderr, "Issue while loading data on GPU\n");
            fprintf(stderr, "cudaMemcpy failed: %s\n", cudaGetErrorString(error));
         }


         //Copy the index of the tetrahedrons
         this->numTets=tetIdVector.size();
         size_t* tetVectorPointer = new size_t[4*tetIdVector.size()];

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

         //DEBUG
         /*
            for(size_t k=0; k<this->numTets; k++){ 
            std::cout << "NUM = " << k << " || Id = " << tetVectorPointer[4*k] << " , " << tetVectorPointer[4*k+1] << " , " << tetVectorPointer[4*k+2] << " , " << tetVectorPointer[4*k+3] << std::endl;
            }
          */

         size_t size_tetVector = 4*this->numTets*sizeof(size_t);
         error=cudaMalloc((void **) &(this->tetId_d), size_tetVector);
         if (error != cudaSuccess) {
            fprintf(stderr, "Issue while initializing data buffer (structure)\n");
            fprintf(stderr, "cudaMalloc failed: %s\n", cudaGetErrorString(error));
         }

         error=cudaMemcpy(this->tetId_d, tetVectorPointer, size_tetVector, cudaMemcpyHostToDevice);
         if (error != cudaSuccess) {
            fprintf(stderr, "Issue while loading data (structure) on GPU\n");
            fprintf(stderr, "cudaMemcpy failed: %s\n", cudaGetErrorString(error));
         }

         //Copy the association vector
         this->associationVector = associationVectorU;

         //Create sphere buffer: (x y z radius)
         size_t size_sphereBuffer = 4*this->numTets*sizeof(float);
         error=cudaMalloc((void **) &(this->sphereBuf_d), size_sphereBuffer);
         if (error != cudaSuccess) {
            fprintf(stderr, "Issue while initializing buffer of spheres\n");
            fprintf(stderr, "cudaMalloc failed: %s\n", cudaGetErrorString(error));
         }


         //Create normal buffer for tetrahedron (n1 n2 n3 n4 P0)
         size_t size_normalBuffer = 4*6*this->numTets*sizeof(float);
         error=cudaMalloc((void **) &(this->normalBuf_d), size_normalBuffer);
         if (error != cudaSuccess) {
            fprintf(stderr, "Issue while initializing normal buffer\n");
            fprintf(stderr, "cudaMalloc failed: %s\n", cudaGetErrorString(error));
         }

         //Create collision vector (CPU and GPU)
         this->collideVectorArray = new bool[this->numTets];
         size_t size_collideVectorArray = this->numTets*sizeof(bool); 
         error=cudaMalloc((void **) &(this->collideVectorArray_d), size_collideVectorArray);
         if (error != cudaSuccess) {
            fprintf(stderr, "Issue while creating output vector\n");
            fprintf(stderr, "cudaMalloc failed: %s\n", cudaGetErrorString(error));
         }

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
   delete[] this->collideVectorArray;
   cudaFree(this->data_d);
   cudaFree(this->tetId_d);
   cudaFree(this->sphereBuf_d);
   cudaFree(this->normalBuf_d);

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

bool MeshStructureCollider::collide(std::vector<bool> & collisionList){

   cudaError_t error;

   cudaDeviceProp prop;
   cudaGetDeviceProperties(&prop,0);
   size_t sizeBlockToUse = sqrt (prop.maxThreadsDim[0]);
   size_t sizeGridx = ceil(float(this->numTets)/float(sizeBlockToUse*sizeBlockToUse));
   size_t sizeGridy = 1;

   //std::cout << "Tets = " << this->numTets << std::endl;
   //std::cout << "Block size = " << sizeBlockToUse << std::endl;
   //std::cout << "Grid size = " << sizeGridx << std::endl;

   dim3 dimGrid(sizeGridx,sizeGridy);
   dim3 dimBlock(sizeBlockToUse, sizeBlockToUse);

   //Update sphere
   updateCenterSphere<<<dimGrid, dimBlock>>>(this->data_d, this->tetId_d, this->sphereBuf_d, this->numTets);
   if ( cudaSuccess != cudaGetLastError() ){ 
      fprintf(stderr, "Error updating spheres\n");
   }
   cudaDeviceSynchronize();
   if ( cudaSuccess != cudaGetLastError() ){ 
      fprintf(stderr, "Error updating spheres (synchro)\n");
   }

   //Update Normals Buff
   updatePlanesTet<<<dimGrid, dimBlock>>>(this->data_d, this->tetId_d, this->normalBuf_d, this->numTets);
   if ( cudaSuccess != cudaGetLastError() ){ 
      fprintf(stderr, "Error updating face buffers\n");
   }
   cudaDeviceSynchronize();
   if ( cudaSuccess != cudaGetLastError() ){ 
      fprintf(stderr, "Error updating face buffers(synchro)\n");
   }

   //Check collision
   checkForIntersection<<<dimGrid, dimBlock>>>(this->data_d, this->tetId_d, this->sphereBuf_d, this->normalBuf_d, this->numTets, this->collideVectorArray_d); 
   if ( cudaSuccess != cudaGetLastError() ){ 
      fprintf(stderr, "Error collision\n");
   }
   error = cudaDeviceSynchronize();
   if ( cudaSuccess != error ){ 
      fprintf(stderr, "Error collision(synchro)\n");
      fprintf(stderr, "Synchro: %s\n", cudaGetErrorString(error));

   }


   //debugTestKernel<<<dimGrid, dimBlock>>>(this->tetId_d, this->numTets);
   //cudaDeviceSynchronize();

   //Get results
   size_t sizeResult = this->numTets * sizeof(bool);

   error = cudaMemcpy(collideVectorArray, this->collideVectorArray_d, sizeResult, cudaMemcpyDeviceToHost);
   if (error != cudaSuccess) {
      fprintf(stderr, "Issue while saving the collision vector\n");
      fprintf(stderr, "cudaMemcpy failed: %s\n", cudaGetErrorString(error));
   }

   for(size_t k=0; k<this->numTets;k++){
      collisionList.at(k)=collideVectorArray[k];
   }

   return false; 

}
