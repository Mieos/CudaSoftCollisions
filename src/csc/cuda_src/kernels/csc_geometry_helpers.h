#ifndef _H_KERNELS_CSC_GEOMETRY_HELPERS_
#define _H_KERNELS_CSC_GEOMETRY_HELPERS_

#define EPSILON_d 0.00001f
#define EPSILON_d_2 0.00000001f

//Cross product
__device__ void crossP(float* v1, float* v2, float* res){

   res[0] = v1[1]*v2[2] - v1[2]*v2[1];
   res[1] = v1[2]*v2[0] - v1[0]*v2[2];
   res[2] = v1[0]*v2[1] - v1[1]*v2[0];

}

//Check position of a point toward a plane (one side, the other side or on the plane)
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


__device__ float cDet4(float* M){

   float result = 
      M[0]*(
            M[5]*(
               M[10]*M[15]-M[14]*M[11]
               )
            -M[9]*(
               M[6]*M[15]-M[14]*M[7]
               )
            +M[13]*(
               M[6]*M[11]-M[10]*M[7]
               )
            )
      - M[4]*(
            M[1]*(
               M[10]*M[15]-M[14]*M[11]
               )
            -M[9]*(
               M[2]*M[15]-M[14]*M[3]
               )
            +M[13]*(
               M[2]*M[11]-M[10]*M[3]
               )
            )
      + M[8]*(
            M[1]*(
               M[6]*M[15]-M[14]*M[7]
               )
            -M[5]*(
               M[2]*M[15]-M[14]*M[3]
               )
            +M[13]*(
               M[2]*M[7]-M[6]*M[3]
               )
            )
      - M[12]*(
            M[1]*(
               M[6]*M[11]-M[10]*M[7]
               )
            -M[5]*(
               M[2]*M[11]-M[10]*M[3]
               )
            +M[9]*(
               M[2]*M[7]-M[6]*M[3]
               )
            );

               

   return result;

}

#endif


