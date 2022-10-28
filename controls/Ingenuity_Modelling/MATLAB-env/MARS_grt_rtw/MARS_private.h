/*
 * MARS_private.h
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "MARS".
 *
 * Model version              : 1.6
 * Simulink Coder version : 9.7 (R2022a) 13-Nov-2021
 * C++ source code generated on : Thu Oct 20 23:51:39 2022
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objective: Execution efficiency
 * Validation result: Not run
 */

#ifndef RTW_HEADER_MARS_private_h_
#define RTW_HEADER_MARS_private_h_
#include "rtwtypes.h"
#include "multiword_types.h"

/* Private macros used by the generated code to access rtModel */
#ifndef rtmSetFirstInitCond
#define rtmSetFirstInitCond(rtm, val)  ((rtm)->Timing.firstInitCondFlag = (val))
#endif

#ifndef rtmIsFirstInitCond
#define rtmIsFirstInitCond(rtm)        ((rtm)->Timing.firstInitCondFlag)
#endif

#ifndef rtmIsMajorTimeStep
#define rtmIsMajorTimeStep(rtm)        (((rtm)->Timing.simTimeStep) == MAJOR_TIME_STEP)
#endif

#ifndef rtmIsMinorTimeStep
#define rtmIsMinorTimeStep(rtm)        (((rtm)->Timing.simTimeStep) == MINOR_TIME_STEP)
#endif

#ifndef rtmSetTPtr
#define rtmSetTPtr(rtm, val)           ((rtm)->Timing.t = (val))
#endif

extern real_T rt_atan2d_snf(real_T u0, real_T u1);
extern real_T rt_modd_snf(real_T u0, real_T u1);
extern real_T rt_remd_snf(real_T u0, real_T u1);
extern real_T rt_powd_snf(real_T u0, real_T u1);
extern void rt_mrdivide_U1d1x3_U2d_9vOrDY9Z(const real_T u0[3], const real_T u1
  [9], real_T y[3]);

/* private model entry point functions */
extern void MARS_derivatives();

#endif                                 /* RTW_HEADER_MARS_private_h_ */
