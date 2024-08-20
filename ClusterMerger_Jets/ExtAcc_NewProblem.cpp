#include "CUPOT.h"
#ifdef __CUDACC__
#include "CUDA_CheckError.h"
#endif

#ifdef GRAVITY

/********************************************************
1. Point-mass external acceleration
   --> It can be regarded as a template for implementing
       other external acceleration

2. This file is shared by both CPU and GPU

   GPU_Gravity/CUPOT_ExtAcc_NFW.cu -> CPU_Gravity/CPU_ExtAcc_NFW.cpp

3. Three steps are required to implement external acceleration

   I.   Set an auxiliary array
        --> SetExtAccAuxArray_NFW()

   II.  Specify external acceleration
        --> ExtAcc_NFW()

   III. Set initialization functions
        --> SetGPUExtAcc_NFW()
            SetCPUExtAcc_NFW()
            Init_ExtAcc_NFW()

4. The external acceleration major routine, ExtAcc_NFW(),
   must be thread-safe and not use any global variable

5. Reference: https://github.com/gamer-project/gamer/wiki/Gravity#external-accelerationpotential
********************************************************/



// soften length implementation
#  define SOFTEN_PLUMMER
//#  define SOFTEN_RUFFERT



// =================================
// I. Set an auxiliary array
// =================================

#ifndef __CUDACC__
//-------------------------------------------------------------------------------------------------------
// Function    :  SetExtAccAuxArray_NFW
// Description :  Set the auxiliary array ExtAcc_AuxArray[] used by ExtAcc_NFW()
//
// Note        :  1. Invoked by Init_ExtAcc_NFW()
//                2. AuxArray[] has the size of EXT_ACC_NAUX_MAX defined in Macro.h (default = 20)
//                3. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  AuxArray : Array to be filled up
//                Time     : Target physical time
//
// Return      :  AuxArray[]
//-------------------------------------------------------------------------------------------------------
void SetExtAccAuxArray_NFW( double AuxArray[], const double Time )
{

// example parameters
   // the mass in only for Perseus cluster. virial radius in function below is also for that (xianshu)
   const double coef= 1.0; // a ratio between total mass by integration and by analytical function (not used, finally)
   const double M   = 8.5E14 * Const_Msun / UNIT_M * coef; // 1.9885e47 as UNIT_M
   const double GM  = NEWTON_G*M;
   const double Eps = 0.0;

   const double r200     = 1.9355819525244697 * Const_Mpc / UNIT_L; // in unit Mpc
   const double c        = 6.81;
   const double rs       = r200 / c; // in unit Mpc
   const double a        = 1.5*rs;
   const double mp       = 1 - (1-(2+3*c)/2/POW(1+c, 1.5)); // ratio between M200 and M_total

   AuxArray[0] = 0.5*amr->BoxSize[0];  // x coordinate of the external acceleration center
   AuxArray[1] = 0.5*amr->BoxSize[1];  // y ...
   AuxArray[2] = 0.5*amr->BoxSize[2];  // z ...
   AuxArray[3] = GM;                   // gravitational_constant*point_source_mass
   AuxArray[4] = Eps;                  // soften_length (<=0.0 --> disable)
   AuxArray[5] = rs;                   // specific radius
   AuxArray[6] = c;                    // concentration coefficient
   AuxArray[7] = a;                    // 1.5 rs defined in sNFW acceleration
   AuxArray[8] = mp;                   // mass ratio between virial mass and total mass

} // FUNCTION : SetExtAccAuxArray_NFW
#endif // #ifndef __CUDACC__



// =================================
// II. Specify external acceleration
// =================================

//-----------------------------------------------------------------------------------------
// Function    :  ExtAcc_NFW
// Description :  Calculate the external acceleration at the given coordinates and time
//
// Note        :  1. This function is shared by CPU and GPU
//                2. Auxiliary array UserArray[] is set by SetExtAccAuxArray_NFW(), where
//                      UserArray[0] = x coordinate of the external acceleration center
//                      UserArray[1] = y ...
//                      UserArray[2] = z ..
//                      UserArray[3] = gravitational_constant*point_source_mass
//                      UserArray[4] = soften_length (<=0.0 --> disable)
//                3. Two different soften length implementations are supported
//                   --> SOFTEN_PLUMMER & SOFTEN_RUFFERT
//
// Parameter   :  Acc       : Array to store the output external acceleration
//                x/y/z     : Target spatial coordinates
//                Time      : Target physical time
//                UserArray : User-provided auxiliary array
//
// Return      :  External acceleration Acc[] at (x,y,z,Time)
//-----------------------------------------------------------------------------------------
GPU_DEVICE_NOINLINE
static void ExtAcc_NFW( real Acc[], const double x, const double y, const double z, const double Time,
                              const double UserArray[] )
{

   const real Cen[3]   = { (real)UserArray[0], (real)UserArray[1], (real)UserArray[2] };
   const real GM       = (real)UserArray[3];
   const real eps      = (real)UserArray[4];
   const real dx       = (real)(x - Cen[0]);
   const real dy       = (real)(y - Cen[1]);
   const real dz       = (real)(z - Cen[2]);
   const real r        = SQRT( dx*dx + dy*dy + dz*dz );
   const real rs       = (real)UserArray[5];
   const real c        = (real)UserArray[6];
   const real a        = (real)UserArray[7];
   const real mp       = (real)UserArray[8];
   // const real r200     = 2.440; // in unit Mpc
   // const real r200     = 1.9355819525244697; // in unit Mpc
   // const real c        = 6.81;
   // const real rs       = r200 / c; // in unit Mpc
   // const real a        = 1.5*rs;
   // const real mp       = 1 - (1-(2+3*c)/2/POW(1+c, (real)1.5)); // ratio between M200 and M_total

// Plummer
#  if   ( defined SOFTEN_PLUMMER )
   const real _r3 = ( eps <= (real)0.0 ) ? (real)1.0/CUBE(r) : POW( SQR(r)+SQR(eps), (real)-1.5 );

// Ruffert 1994
#  elif ( defined SOFTEN_RUFFERT )
   const real tmp = EXP( -SQR(r)/SQR(eps) );
   const real _r3 = ( eps <= (real)0.0 ) ? (real)1.0/CUBE(r) : POW( SQR(r)+SQR(eps)*tmp, (real)-1.5 )*( (real)1.0 - tmp );

#  else
   const real _r3 = (real)1.0/CUBE(r);
#  endif
   // This coefficient follows the mass profile in cluster generator (sNFW)
   // reference: Lilley et al. 2018, MNRAS 476, 2086-2091
   // Acc[0] = -GM/r*dx * (1/POW(r,2) - POW(a,0.5)/POW(r,2)/POW(a+r,0.5) - POW(a,0.5)/2/r/POW(a+r,1.5)) / mp;
   // Acc[1] = -GM/r*dy * (1/POW(r,2) - POW(a,0.5)/POW(r,2)/POW(a+r,0.5) - POW(a,0.5)/2/r/POW(a+r,1.5)) / mp;
   // Acc[2] = -GM/r*dz * (1/POW(r,2) - POW(a,0.5)/POW(r,2)/POW(a+r,0.5) - POW(a,0.5)/2/r/POW(a+r,1.5)) / mp;
   // This coefficient follows the mass profile in cluster generator (sNFW version)
   Acc[0] = -GM/mp*_r3*dx * ((real)1 - ((real)2+(real)3*r/rs)/(real)2/POW((real)1+r/rs, (real)1.5));
   Acc[1] = -GM/mp*_r3*dy * ((real)1 - ((real)2+(real)3*r/rs)/(real)2/POW((real)1+r/rs, (real)1.5));
   Acc[2] = -GM/mp*_r3*dz * ((real)1 - ((real)2+(real)3*r/rs)/(real)2/POW((real)1+r/rs, (real)1.5));
   // This coefficient follows the mass profile (NFW version)
   // Acc[0] = GM*_r3*dx / (log(1+c)-c/(1+c)) * (r/(r+rs)-log(1+r/rs));
   // Acc[1] = GM*_r3*dy / (log(1+c)-c/(1+c)) * (r/(r+rs)-log(1+r/rs));
   // Acc[2] = GM*_r3*dz / (log(1+c)-c/(1+c)) * (r/(r+rs)-log(1+r/rs));

} // FUNCTION : ExtAcc_NFW



// =================================
// III. Set initialization functions
// =================================

#ifdef __CUDACC__
#  define FUNC_SPACE __device__ static
#else
#  define FUNC_SPACE            static
#endif

FUNC_SPACE ExtAcc_t ExtAcc_Ptr = ExtAcc_NFW;

//-----------------------------------------------------------------------------------------
// Function    :  SetCPU/GPUExtAcc_NFW
// Description :  Return the function pointers of the CPU/GPU external acceleration routines
//
// Note        :  1. Invoked by Init_ExtAcc_NFW()
//                2. Must obtain the CPU and GPU function pointers by **separate** routines
//                   since CPU and GPU functions are compiled completely separately in GAMER
//                   --> In other words, a unified routine like the following won't work
//
//                      SetExtAcc_NFW( ExtAcc_t &CPUExtAcc_Ptr, ExtAcc_t &GPUExtAcc_Ptr )
//
// Parameter   :  CPU/GPUExtAcc_Ptr (call-by-reference)
//
// Return      :  CPU/GPUExtAcc_Ptr
//-----------------------------------------------------------------------------------------
#ifdef __CUDACC__
__host__
void SetGPUExtAcc_NFW( ExtAcc_t &GPUExtAcc_Ptr )
{
   CUDA_CHECK_ERROR(  cudaMemcpyFromSymbol( &GPUExtAcc_Ptr, ExtAcc_Ptr, sizeof(ExtAcc_t) )  );
}

#else // #ifdef __CUDACC__

void SetCPUExtAcc_NFW( ExtAcc_t &CPUExtAcc_Ptr )
{
   CPUExtAcc_Ptr = ExtAcc_Ptr;
}

#endif // #ifdef __CUDACC__ ... else ...



#ifndef __CUDACC__

// local function prototypes
void SetExtAccAuxArray_NFW( double [], const double );
void SetCPUExtAcc_NFW( ExtAcc_t & );
#ifdef GPU
void SetGPUExtAcc_NFW( ExtAcc_t & );
#endif

//-----------------------------------------------------------------------------------------
// Function    :  Init_ExtAcc_NFW
// Description :  Initialize external acceleration
//
// Note        :  1. Set an auxiliary array by invoking SetExtAccAuxArray_*()
//                   --> It will be copied to GPU automatically in CUAPI_SetConstMemory()
//                2. Set the CPU/GPU external acceleration major routines by invoking SetCPU/GPUExtAcc_*()
//                3. Invoked by Init_ExtAccPot()
//                   --> Enable it by linking to the function pointer "Init_ExtAcc_Ptr"
//                4. Add "#ifndef __CUDACC__" since this routine is only useful on CPU
//
// Parameter   :  None
//
// Return      :  None
//-----------------------------------------------------------------------------------------
void Init_ExtAcc_NFW()
{

   SetExtAccAuxArray_NFW( ExtAcc_AuxArray, Time[0] );

   SetCPUExtAcc_NFW( CPUExtAcc_Ptr );
#  ifdef GPU
   SetGPUExtAcc_NFW( GPUExtAcc_Ptr );
#  endif

} // FUNCTION : Init_ExtAcc_NFW

#endif // #ifndef __CUDACC__



#endif // #ifdef GRAVITY
