/*
 * MARS.h
 *
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * Code generation for model "MARS".
 *
 * Model version              : 1.28
 * Simulink Coder version : 9.7 (R2022a) 13-Nov-2021
 * C source code generated on : Fri Oct 28 08:09:17 2022
 *
 * Target selection: grt.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Intel->x86-64 (Windows64)
 * Code generation objective: Execution efficiency
 * Validation result: Not run
 */

#ifndef RTW_HEADER_MARS_h_
#define RTW_HEADER_MARS_h_
#ifndef MARS_COMMON_INCLUDES_
#define MARS_COMMON_INCLUDES_
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#endif                                 /* MARS_COMMON_INCLUDES_ */

#include "MARS_types.h"
#include <string.h>
#include "rt_nonfinite.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetContStateDisabled
#define rtmGetContStateDisabled(rtm)   ((rtm)->contStateDisabled)
#endif

#ifndef rtmSetContStateDisabled
#define rtmSetContStateDisabled(rtm, val) ((rtm)->contStateDisabled = (val))
#endif

#ifndef rtmGetContStates
#define rtmGetContStates(rtm)          ((rtm)->contStates)
#endif

#ifndef rtmSetContStates
#define rtmSetContStates(rtm, val)     ((rtm)->contStates = (val))
#endif

#ifndef rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm) ((rtm)->CTOutputIncnstWithState)
#endif

#ifndef rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag
#define rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm, val) ((rtm)->CTOutputIncnstWithState = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
#define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
#define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetIntgData
#define rtmGetIntgData(rtm)            ((rtm)->intgData)
#endif

#ifndef rtmSetIntgData
#define rtmSetIntgData(rtm, val)       ((rtm)->intgData = (val))
#endif

#ifndef rtmGetOdeF
#define rtmGetOdeF(rtm)                ((rtm)->odeF)
#endif

#ifndef rtmSetOdeF
#define rtmSetOdeF(rtm, val)           ((rtm)->odeF = (val))
#endif

#ifndef rtmGetOdeY
#define rtmGetOdeY(rtm)                ((rtm)->odeY)
#endif

#ifndef rtmSetOdeY
#define rtmSetOdeY(rtm, val)           ((rtm)->odeY = (val))
#endif

#ifndef rtmGetPeriodicContStateIndices
#define rtmGetPeriodicContStateIndices(rtm) ((rtm)->periodicContStateIndices)
#endif

#ifndef rtmSetPeriodicContStateIndices
#define rtmSetPeriodicContStateIndices(rtm, val) ((rtm)->periodicContStateIndices = (val))
#endif

#ifndef rtmGetPeriodicContStateRanges
#define rtmGetPeriodicContStateRanges(rtm) ((rtm)->periodicContStateRanges)
#endif

#ifndef rtmSetPeriodicContStateRanges
#define rtmSetPeriodicContStateRanges(rtm, val) ((rtm)->periodicContStateRanges = (val))
#endif

#ifndef rtmGetZCCacheNeedsReset
#define rtmGetZCCacheNeedsReset(rtm)   ((rtm)->zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
#define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->zCCacheNeedsReset = (val))
#endif

#ifndef rtmGetdX
#define rtmGetdX(rtm)                  ((rtm)->derivs)
#endif

#ifndef rtmSetdX
#define rtmSetdX(rtm, val)             ((rtm)->derivs = (val))
#endif

#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
#define rtmGetStopRequested(rtm)       ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
#define rtmSetStopRequested(rtm, val)  ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
#define rtmGetStopRequestedPtr(rtm)    (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
#define rtmGetT(rtm)                   (rtmGetTPtr((rtm))[0])
#endif

#ifndef rtmGetTPtr
#define rtmGetTPtr(rtm)                ((rtm)->Timing.t)
#endif

/* Block signals (default storage) */
typedef struct {
  real_T Merge[4];                     /* '<S74>/Merge' */
  real_T Product4[9];                  /* '<S61>/Product4' */
  real_T Output;                       /* '<S5>/Output' */
  real_T IntegralGain;                 /* '<S35>/Integral Gain' */
  real_T FilterCoefficient;            /* '<S41>/Filter Coefficient' */
  real_T Sum2[3];                      /* '<S57>/Sum2' */
  real_T TmpSignalConversionAtqInport1[4];/* '<S61>/qdot' */
  real_T Product1[3];                  /* '<S62>/Product1' */
  real_T Product5[3];                  /* '<S62>/Product5' */
  real_T Sum2_e[3];                    /* '<S62>/Sum2' */
  real_T DeadZone[3];                  /* '<S63>/Dead Zone' */
  real_T Product2[3];                  /* '<S64>/Product2' */
  real_T Merge_j;                      /* '<S104>/Merge' */
  real_T Merge_j2;                     /* '<S108>/Merge' */
} B_MARS_T;

/* Block states (default storage) for system '<Root>' */
typedef struct {
  real_T NextOutput;                   /* '<S5>/White Noise' */
  real_T Product2_DWORK4[9];           /* '<S64>/Product2' */
  uint32_T RandSeed;                   /* '<S5>/White Noise' */
  int_T p1_IWORK;                      /* '<S62>/p1' */
  int_T q_IWORK;                       /* '<S61>/q' */
  int_T pqr_IWORK;                     /* '<S57>/p,q,r ' */
  int_T p_IWORK;                       /* '<S62>/p' */
  int8_T If_ActiveSubsystem;           /* '<S74>/If' */
  int8_T If_ActiveSubsystem_h;         /* '<S73>/If' */
  int8_T If1_ActiveSubsystem;          /* '<S103>/If1' */
  int8_T If1_ActiveSubsystem_l;        /* '<S126>/If1' */
  int8_T If_ActiveSubsystem_c;         /* '<S104>/If' */
  int8_T If_ActiveSubsystem_c0;        /* '<S108>/If' */
} DW_MARS_T;

/* Continuous states (default storage) */
typedef struct {
  real_T p1_CSTATE[3];                 /* '<S62>/p1' */
  real_T q_CSTATE[4];                  /* '<S61>/q' */
  real_T Filter_CSTATE;                /* '<S33>/Filter' */
  real_T Integrator_CSTATE;            /* '<S38>/Integrator' */
  real_T Integrator_CSTATE_f;          /* '<S65>/Integrator' */
  real_T ubvbwb_CSTATE[3];             /* '<S63>/ub,vb,wb' */
  real_T pqr_CSTATE[3];                /* '<S57>/p,q,r ' */
  real_T p_CSTATE[3];                  /* '<S62>/p' */
} X_MARS_T;

/* State derivatives (default storage) */
typedef struct {
  real_T p1_CSTATE[3];                 /* '<S62>/p1' */
  real_T q_CSTATE[4];                  /* '<S61>/q' */
  real_T Filter_CSTATE;                /* '<S33>/Filter' */
  real_T Integrator_CSTATE;            /* '<S38>/Integrator' */
  real_T Integrator_CSTATE_f;          /* '<S65>/Integrator' */
  real_T ubvbwb_CSTATE[3];             /* '<S63>/ub,vb,wb' */
  real_T pqr_CSTATE[3];                /* '<S57>/p,q,r ' */
  real_T p_CSTATE[3];                  /* '<S62>/p' */
} XDot_MARS_T;

/* State disabled  */
typedef struct {
  boolean_T p1_CSTATE[3];              /* '<S62>/p1' */
  boolean_T q_CSTATE[4];               /* '<S61>/q' */
  boolean_T Filter_CSTATE;             /* '<S33>/Filter' */
  boolean_T Integrator_CSTATE;         /* '<S38>/Integrator' */
  boolean_T Integrator_CSTATE_f;       /* '<S65>/Integrator' */
  boolean_T ubvbwb_CSTATE[3];          /* '<S63>/ub,vb,wb' */
  boolean_T pqr_CSTATE[3];             /* '<S57>/p,q,r ' */
  boolean_T p_CSTATE[3];               /* '<S62>/p' */
} XDis_MARS_T;

/* Invariant block signals (default storage) */
typedef struct {
  const real_T Product2[9];            /* '<S61>/Product2' */
  const real_T Selector[9];            /* '<S64>/Selector' */
  const real_T Selector2[9];           /* '<S64>/Selector2' */
  const real_T Abs2[9];                /* '<S153>/Abs2' */
  const real_T Product_n[3];           /* '<S128>/Product' */
  const real_T Product_l[3];           /* '<S130>/Product' */
  const real_T Product_d[3];           /* '<S129>/Product' */
  const real_T Add_o;                  /* '<S147>/Add' */
  const real_T Add_jx;                 /* '<S148>/Add' */
  const real_T Add_n;                  /* '<S149>/Add' */
} ConstB_MARS_T;

#ifndef ODE3_INTG
#define ODE3_INTG

/* ODE3 Integration Data */
typedef struct {
  real_T *y;                           /* output */
  real_T *f[3];                        /* derivatives */
} ODE3_IntgData;

#endif

/* Constant parameters (default storage) */
typedef struct {
  /* Pooled Parameter (Expression: -eye(3))
   * Referenced by:
   *   '<S115>/Bias1'
   *   '<S153>/Bias1'
   */
  real_T pooled5[9];
} ConstP_MARS_T;

/* Real-time Model Data Structure */
struct tag_RTM_MARS_T {
  const char_T *errorStatus;
  RTWSolverInfo solverInfo;
  X_MARS_T *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  boolean_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T CTOutputIncnstWithState;
  real_T odeY[19];
  real_T odeF[3][19];
  ODE3_IntgData intgData;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    int_T numContStates;
    int_T numPeriodicContStates;
    int_T numSampTimes;
  } Sizes;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    uint32_T clockTick0;
    time_T stepSize0;
    uint32_T clockTick1;
    boolean_T firstInitCondFlag;
    struct {
      uint8_T TID[3];
    } TaskCounters;

    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *t;
    time_T tArray[3];
  } Timing;
};

/* Block signals (default storage) */
extern B_MARS_T MARS_B;

/* Continuous states (default storage) */
extern X_MARS_T MARS_X;

/* Block states (default storage) */
extern DW_MARS_T MARS_DW;
extern const ConstB_MARS_T MARS_ConstB;/* constant block i/o */

/* Constant parameters (default storage) */
extern const ConstP_MARS_T MARS_ConstP;

/* Model entry point functions */
extern void MARS_initialize(void);
extern void MARS_step(void);
extern void MARS_terminate(void);

/* Real-time Model object */
extern RT_MODEL_MARS_T *const MARS_M;

/*-
 * These blocks were eliminated from the model due to optimizations:
 *
 * Block '<S1>/Constant' : Unused code path elimination
 * Block '<S2>/Display' : Unused code path elimination
 * Block '<S2>/Scope' : Unused code path elimination
 * Block '<S2>/Scope1' : Unused code path elimination
 * Block '<S57>/Math Function1' : Unused code path elimination
 * Block '<S57>/Product1' : Unused code path elimination
 * Block '<S57>/Sum1' : Unused code path elimination
 * Block '<S68>/Unit Conversion' : Unused code path elimination
 * Block '<S69>/Unit Conversion' : Unused code path elimination
 * Block '<S4>/Height' : Unused code path elimination
 * Block '<S4>/Index' : Unused code path elimination
 * Block '<S58>/Divide' : Unused code path elimination
 * Block '<S58>/Exp' : Unused code path elimination
 * Block '<S58>/Gamma*R' : Unused code path elimination
 * Block '<S58>/Height Limiter' : Unused code path elimination
 * Block '<S58>/Multiply' : Unused code path elimination
 * Block '<S58>/Pres Lapse' : Unused code path elimination
 * Block '<S58>/R' : Unused code path elimination
 * Block '<S58>/Square Root' : Unused code path elimination
 * Block '<S58>/Sum' : Unused code path elimination
 * Block '<S58>/Surface Temp' : Unused code path elimination
 * Block '<S58>/TempLapse' : Unused code path elimination
 * Block '<S3>/Rate Transition' : Eliminated since input and output rates are identical
 * Block '<S88>/Reshape (9) to [3x3] column-major' : Reshape block reduction
 * Block '<S99>/Reshape (9) to [3x3] column-major' : Reshape block reduction
 * Block '<S102>/Reshape' : Reshape block reduction
 * Block '<S102>/Reshape1' : Reshape block reduction
 * Block '<S102>/Reshape2' : Reshape block reduction
 * Block '<S115>/Reshape' : Reshape block reduction
 * Block '<S74>/Reshape 3x3 -> 9' : Reshape block reduction
 * Block '<S153>/Reshape' : Reshape block reduction
 * Block '<S160>/Reshape' : Reshape block reduction
 * Block '<S171>/Reshape (9) to [3x3] column-major' : Reshape block reduction
 * Block '<S175>/Reshape (9) to [3x3] column-major' : Reshape block reduction
 * Block '<S192>/Reshape1' : Reshape block reduction
 * Block '<S192>/Reshape2' : Reshape block reduction
 * Block '<S193>/Reshape1' : Reshape block reduction
 * Block '<S193>/Reshape2' : Reshape block reduction
 * Block '<S64>/Reshape' : Reshape block reduction
 * Block '<S64>/Reshape1' : Reshape block reduction
 * Block '<S196>/Reshape (9) to [3x3] column-major' : Reshape block reduction
 */

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'MARS'
 * '<S1>'   : 'MARS/Altimeter'
 * '<S2>'   : 'MARS/Altitude Hold Control'
 * '<S3>'   : 'MARS/FlightGear Preconfigured 6DoF Animation'
 * '<S4>'   : 'MARS/M.A.R.S.'
 * '<S5>'   : 'MARS/Altimeter/Band-Limited White Noise'
 * '<S6>'   : 'MARS/Altitude Hold Control/PID Controller'
 * '<S7>'   : 'MARS/Altitude Hold Control/PID Controller/Anti-windup'
 * '<S8>'   : 'MARS/Altitude Hold Control/PID Controller/D Gain'
 * '<S9>'   : 'MARS/Altitude Hold Control/PID Controller/Filter'
 * '<S10>'  : 'MARS/Altitude Hold Control/PID Controller/Filter ICs'
 * '<S11>'  : 'MARS/Altitude Hold Control/PID Controller/I Gain'
 * '<S12>'  : 'MARS/Altitude Hold Control/PID Controller/Ideal P Gain'
 * '<S13>'  : 'MARS/Altitude Hold Control/PID Controller/Ideal P Gain Fdbk'
 * '<S14>'  : 'MARS/Altitude Hold Control/PID Controller/Integrator'
 * '<S15>'  : 'MARS/Altitude Hold Control/PID Controller/Integrator ICs'
 * '<S16>'  : 'MARS/Altitude Hold Control/PID Controller/N Copy'
 * '<S17>'  : 'MARS/Altitude Hold Control/PID Controller/N Gain'
 * '<S18>'  : 'MARS/Altitude Hold Control/PID Controller/P Copy'
 * '<S19>'  : 'MARS/Altitude Hold Control/PID Controller/Parallel P Gain'
 * '<S20>'  : 'MARS/Altitude Hold Control/PID Controller/Reset Signal'
 * '<S21>'  : 'MARS/Altitude Hold Control/PID Controller/Saturation'
 * '<S22>'  : 'MARS/Altitude Hold Control/PID Controller/Saturation Fdbk'
 * '<S23>'  : 'MARS/Altitude Hold Control/PID Controller/Sum'
 * '<S24>'  : 'MARS/Altitude Hold Control/PID Controller/Sum Fdbk'
 * '<S25>'  : 'MARS/Altitude Hold Control/PID Controller/Tracking Mode'
 * '<S26>'  : 'MARS/Altitude Hold Control/PID Controller/Tracking Mode Sum'
 * '<S27>'  : 'MARS/Altitude Hold Control/PID Controller/Tsamp - Integral'
 * '<S28>'  : 'MARS/Altitude Hold Control/PID Controller/Tsamp - Ngain'
 * '<S29>'  : 'MARS/Altitude Hold Control/PID Controller/postSat Signal'
 * '<S30>'  : 'MARS/Altitude Hold Control/PID Controller/preSat Signal'
 * '<S31>'  : 'MARS/Altitude Hold Control/PID Controller/Anti-windup/Passthrough'
 * '<S32>'  : 'MARS/Altitude Hold Control/PID Controller/D Gain/Internal Parameters'
 * '<S33>'  : 'MARS/Altitude Hold Control/PID Controller/Filter/Cont. Filter'
 * '<S34>'  : 'MARS/Altitude Hold Control/PID Controller/Filter ICs/Internal IC - Filter'
 * '<S35>'  : 'MARS/Altitude Hold Control/PID Controller/I Gain/Internal Parameters'
 * '<S36>'  : 'MARS/Altitude Hold Control/PID Controller/Ideal P Gain/Passthrough'
 * '<S37>'  : 'MARS/Altitude Hold Control/PID Controller/Ideal P Gain Fdbk/Disabled'
 * '<S38>'  : 'MARS/Altitude Hold Control/PID Controller/Integrator/Continuous'
 * '<S39>'  : 'MARS/Altitude Hold Control/PID Controller/Integrator ICs/Internal IC'
 * '<S40>'  : 'MARS/Altitude Hold Control/PID Controller/N Copy/Disabled'
 * '<S41>'  : 'MARS/Altitude Hold Control/PID Controller/N Gain/Internal Parameters'
 * '<S42>'  : 'MARS/Altitude Hold Control/PID Controller/P Copy/Disabled'
 * '<S43>'  : 'MARS/Altitude Hold Control/PID Controller/Parallel P Gain/Internal Parameters'
 * '<S44>'  : 'MARS/Altitude Hold Control/PID Controller/Reset Signal/Disabled'
 * '<S45>'  : 'MARS/Altitude Hold Control/PID Controller/Saturation/Enabled'
 * '<S46>'  : 'MARS/Altitude Hold Control/PID Controller/Saturation Fdbk/Disabled'
 * '<S47>'  : 'MARS/Altitude Hold Control/PID Controller/Sum/Sum_PID'
 * '<S48>'  : 'MARS/Altitude Hold Control/PID Controller/Sum Fdbk/Disabled'
 * '<S49>'  : 'MARS/Altitude Hold Control/PID Controller/Tracking Mode/Disabled'
 * '<S50>'  : 'MARS/Altitude Hold Control/PID Controller/Tracking Mode Sum/Passthrough'
 * '<S51>'  : 'MARS/Altitude Hold Control/PID Controller/Tsamp - Integral/Passthrough'
 * '<S52>'  : 'MARS/Altitude Hold Control/PID Controller/Tsamp - Ngain/Passthrough'
 * '<S53>'  : 'MARS/Altitude Hold Control/PID Controller/postSat Signal/Forward_Path'
 * '<S54>'  : 'MARS/Altitude Hold Control/PID Controller/preSat Signal/Forward_Path'
 * '<S55>'  : 'MARS/FlightGear Preconfigured 6DoF Animation/Angle Conversion'
 * '<S56>'  : 'MARS/FlightGear Preconfigured 6DoF Animation/Angle Conversion1'
 * '<S57>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)'
 * '<S58>'  : 'MARS/M.A.R.S./Mars Atmospheric Model'
 * '<S59>'  : 'MARS/M.A.R.S./Subsystem1'
 * '<S60>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Angle Conversion'
 * '<S61>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles'
 * '<S62>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate Position in EI '
 * '<S63>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate Velocity in Body Axes'
 * '<S64>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate omega_dot'
 * '<S65>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Celestial Longitude of Greenwich'
 * '<S66>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Determine Force,  Mass & Inertia'
 * '<S67>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/ECEF to Inertial'
 * '<S68>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Velocity Conversion'
 * '<S69>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Velocity Conversion2'
 * '<S70>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/w_ned'
 * '<S71>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED'
 * '<S72>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED1'
 * '<S73>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles'
 * '<S74>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions'
 * '<S75>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix'
 * '<S76>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Rotation Angles to Direction Cosine Matrix'
 * '<S77>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/qdot'
 * '<S78>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED/A11'
 * '<S79>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED/A12'
 * '<S80>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED/A13'
 * '<S81>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED/A21'
 * '<S82>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED/A22'
 * '<S83>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED/A23'
 * '<S84>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED/A31'
 * '<S85>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED/A32'
 * '<S86>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED/A33'
 * '<S87>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED/Angle Conversion'
 * '<S88>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED/Create Transformation Matrix'
 * '<S89>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED1/A11'
 * '<S90>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED1/A12'
 * '<S91>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED1/A13'
 * '<S92>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED1/A21'
 * '<S93>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED1/A22'
 * '<S94>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED1/A23'
 * '<S95>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED1/A31'
 * '<S96>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED1/A32'
 * '<S97>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED1/A33'
 * '<S98>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED1/Angle Conversion'
 * '<S99>'  : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED1/Create Transformation Matrix'
 * '<S100>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/AxisRotDefault'
 * '<S101>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/AxisRotZeroR3'
 * '<S102>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/Get DCM Values'
 * '<S103>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/Validate DCM'
 * '<S104>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/AxisRotDefault/Protect asincos input'
 * '<S105>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/AxisRotDefault/Protect asincos input/If Action Subsystem'
 * '<S106>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/AxisRotDefault/Protect asincos input/If Action Subsystem1'
 * '<S107>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/AxisRotDefault/Protect asincos input/If Action Subsystem2'
 * '<S108>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/AxisRotZeroR3/Protect asincos input'
 * '<S109>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/AxisRotZeroR3/Protect asincos input/If Action Subsystem'
 * '<S110>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/AxisRotZeroR3/Protect asincos input/If Action Subsystem1'
 * '<S111>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/AxisRotZeroR3/Protect asincos input/If Action Subsystem2'
 * '<S112>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/Validate DCM/If Warning//Error'
 * '<S113>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/Validate DCM/If Warning//Error/Else If Not Orthogonal'
 * '<S114>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/Validate DCM/If Warning//Error/If Not Proper'
 * '<S115>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/Validate DCM/If Warning//Error/isNotOrthogonal'
 * '<S116>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/Validate DCM/If Warning//Error/isNotProper'
 * '<S117>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/Validate DCM/If Warning//Error/Else If Not Orthogonal/Error'
 * '<S118>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/Validate DCM/If Warning//Error/Else If Not Orthogonal/Warning'
 * '<S119>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/Validate DCM/If Warning//Error/If Not Proper/Error'
 * '<S120>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/Validate DCM/If Warning//Error/If Not Proper/Warning'
 * '<S121>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/Validate DCM/If Warning//Error/isNotOrthogonal/transpose*dcm ~= eye(3)'
 * '<S122>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/Validate DCM/If Warning//Error/isNotProper/Determinant of 3x3 Matrix'
 * '<S123>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/Validate DCM/If Warning//Error/isNotProper/determinant does not equal 1'
 * '<S124>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace'
 * '<S125>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Positive Trace'
 * '<S126>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Validate DCM'
 * '<S127>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/trace(DCM)'
 * '<S128>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(1,1)'
 * '<S129>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(2,2)'
 * '<S130>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(3,3)'
 * '<S131>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/diag(DCM)'
 * '<S132>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(1,1)/cos(phi)sin(theta)cos(psi) + sin(phi)sin(psi) -sin(theta)'
 * '<S133>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(1,1)/cos(theta)sin(phi) - (cos(phi)sin(theta)sin(psi) - sin(phi)cos(psi))'
 * '<S134>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(1,1)/cos(theta)sin(psi) + (sin(phi)sin(theta)cos(psi) - cos(phi)sin(psi))'
 * '<S135>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(1,1)/if s~=0; s=0.5//s'
 * '<S136>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(1,1)/u(1) -(u(5)+u(9)) +1'
 * '<S137>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(2,2)/cos(phi)sin(theta)cos(psi) + sin(phi)sin(psi) +sin(theta)'
 * '<S138>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(2,2)/cos(theta)sin(phi) + (cos(phi)sin(theta)sin(psi) - sin(phi)cos(psi))'
 * '<S139>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(2,2)/cos(theta)sin(psi) + (sin(phi)sin(theta)cos(psi) - cos(phi)sin(psi))'
 * '<S140>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(2,2)/if s~=0; s=0.5//s'
 * '<S141>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(2,2)/u(5) -(u(1)+u(9)) +1'
 * '<S142>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(3,3)/cos(phi)sin(theta)cos(psi) + sin(phi)sin(psi) -sin(theta)'
 * '<S143>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(3,3)/cos(theta)sin(phi) + (cos(phi)sin(theta)sin(psi) - sin(phi)cos(psi))'
 * '<S144>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(3,3)/cos(theta)sin(psi) - (sin(phi)sin(theta)cos(psi) - cos(phi)sin(psi))'
 * '<S145>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(3,3)/if s~=0; s=0.5//s'
 * '<S146>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(3,3)/u(9) -(u(1)+u(5)) +1'
 * '<S147>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Positive Trace/cos(phi)sin(theta)cos(psi) + sin(phi)sin(psi) +sin(theta)'
 * '<S148>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Positive Trace/cos(theta)sin(phi) - (cos(phi)sin(theta)sin(psi) - sin(phi)cos(psi))'
 * '<S149>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Positive Trace/cos(theta)sin(psi) - (sin(phi)sin(theta)cos(psi) - cos(phi)sin(psi))'
 * '<S150>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error'
 * '<S151>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/Else If Not Orthogonal'
 * '<S152>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/If Not Proper'
 * '<S153>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/isNotOrthogonal'
 * '<S154>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/isNotProper'
 * '<S155>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/Else If Not Orthogonal/Error'
 * '<S156>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/Else If Not Orthogonal/Warning'
 * '<S157>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/If Not Proper/Error'
 * '<S158>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/If Not Proper/Warning'
 * '<S159>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/isNotOrthogonal/transpose*dcm ~= eye(3)'
 * '<S160>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/isNotProper/Determinant of 3x3 Matrix'
 * '<S161>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/isNotProper/determinant does not equal 1'
 * '<S162>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/A11'
 * '<S163>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/A12'
 * '<S164>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/A13'
 * '<S165>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/A21'
 * '<S166>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/A22'
 * '<S167>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/A23'
 * '<S168>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/A31'
 * '<S169>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/A32'
 * '<S170>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/A33'
 * '<S171>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S172>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/Quaternion Normalize'
 * '<S173>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/Quaternion Normalize/Quaternion Modulus'
 * '<S174>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/Quaternion Normalize/Quaternion Modulus/Quaternion Norm'
 * '<S175>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S176>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/qdot/Quaternion Normalize'
 * '<S177>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/qdot/Quaternion Normalize/Quaternion Modulus'
 * '<S178>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/qdot/Quaternion Normalize/Quaternion Modulus/Quaternion Norm'
 * '<S179>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate Position in EI /pxwe'
 * '<S180>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate Position in EI /pxwe/Subsystem'
 * '<S181>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate Position in EI /pxwe/Subsystem1'
 * '<S182>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate Velocity in Body Axes/Vbxwb'
 * '<S183>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate Velocity in Body Axes/wex(wexp)'
 * '<S184>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate Velocity in Body Axes/wexp'
 * '<S185>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate Velocity in Body Axes/Vbxwb/Subsystem'
 * '<S186>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate Velocity in Body Axes/Vbxwb/Subsystem1'
 * '<S187>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate Velocity in Body Axes/wex(wexp)/Subsystem'
 * '<S188>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate Velocity in Body Axes/wex(wexp)/Subsystem1'
 * '<S189>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate Velocity in Body Axes/wexp/Subsystem'
 * '<S190>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate Velocity in Body Axes/wexp/Subsystem1'
 * '<S191>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate omega_dot/3x3 Cross Product'
 * '<S192>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate omega_dot/I x w'
 * '<S193>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate omega_dot/I x w1'
 * '<S194>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate omega_dot/3x3 Cross Product/Subsystem'
 * '<S195>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/Calculate omega_dot/3x3 Cross Product/Subsystem1'
 * '<S196>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/ECEF to Inertial/Create 3x3 Matrix'
 * '<S197>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/w_ned/Angle Conversion'
 * '<S198>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/w_ned/M+h'
 * '<S199>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/w_ned/N+h'
 * '<S200>' : 'MARS/M.A.R.S./6DOF ECEF (Quaternion)/w_ned/e2'
 * '<S201>' : 'MARS/M.A.R.S./Subsystem1/Components'
 * '<S202>' : 'MARS/M.A.R.S./Subsystem1/Subsystem'
 */
#endif                                 /* RTW_HEADER_MARS_h_ */
