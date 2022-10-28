/*
 * MARS.h
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

#ifndef RTW_HEADER_MARS_h_
#define RTW_HEADER_MARS_h_
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#include "MARS_types.h"
#include <cstring>

extern "C" {

#include "rt_nonfinite.h"

}
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
struct B_MARS_T {
  real_T Merge[4];                     /* '<S17>/Merge' */
  real_T Product4[9];                  /* '<S4>/Product4' */
  real_T Sum2[3];                      /* '<S1>/Sum2' */
  real_T TmpSignalConversionAtqInport1[4];/* '<S4>/qdot' */
  real_T Product1[3];                  /* '<S5>/Product1' */
  real_T Product5[3];                  /* '<S5>/Product5' */
  real_T Sum2_e[3];                    /* '<S5>/Sum2' */
  real_T DeadZone[3];                  /* '<S6>/Dead Zone' */
  real_T Product2[3];                  /* '<S7>/Product2' */
  real_T Merge_j;                      /* '<S47>/Merge' */
  real_T Merge_j2;                     /* '<S51>/Merge' */
};

/* Block states (default storage) for system '<Root>' */
struct DW_MARS_T {
  real_T Product2_DWORK4[9];           /* '<S7>/Product2' */
  int_T p1_IWORK;                      /* '<S5>/p1' */
  int_T q_IWORK;                       /* '<S4>/q' */
  int_T pqr_IWORK;                     /* '<S1>/p,q,r ' */
  int_T p_IWORK;                       /* '<S5>/p' */
  int8_T If_ActiveSubsystem;           /* '<S17>/If' */
  int8_T If_ActiveSubsystem_h;         /* '<S16>/If' */
  int8_T If1_ActiveSubsystem;          /* '<S46>/If1' */
  int8_T If1_ActiveSubsystem_l;        /* '<S69>/If1' */
  int8_T If_ActiveSubsystem_c;         /* '<S47>/If' */
  int8_T If_ActiveSubsystem_c0;        /* '<S51>/If' */
};

/* Continuous states (default storage) */
struct X_MARS_T {
  real_T p1_CSTATE[3];                 /* '<S5>/p1' */
  real_T q_CSTATE[4];                  /* '<S4>/q' */
  real_T Integrator_CSTATE;            /* '<S8>/Integrator' */
  real_T ubvbwb_CSTATE[3];             /* '<S6>/ub,vb,wb' */
  real_T pqr_CSTATE[3];                /* '<S1>/p,q,r ' */
  real_T p_CSTATE[3];                  /* '<S5>/p' */
};

/* State derivatives (default storage) */
struct XDot_MARS_T {
  real_T p1_CSTATE[3];                 /* '<S5>/p1' */
  real_T q_CSTATE[4];                  /* '<S4>/q' */
  real_T Integrator_CSTATE;            /* '<S8>/Integrator' */
  real_T ubvbwb_CSTATE[3];             /* '<S6>/ub,vb,wb' */
  real_T pqr_CSTATE[3];                /* '<S1>/p,q,r ' */
  real_T p_CSTATE[3];                  /* '<S5>/p' */
};

/* State disabled  */
struct XDis_MARS_T {
  boolean_T p1_CSTATE[3];              /* '<S5>/p1' */
  boolean_T q_CSTATE[4];               /* '<S4>/q' */
  boolean_T Integrator_CSTATE;         /* '<S8>/Integrator' */
  boolean_T ubvbwb_CSTATE[3];          /* '<S6>/ub,vb,wb' */
  boolean_T pqr_CSTATE[3];             /* '<S1>/p,q,r ' */
  boolean_T p_CSTATE[3];               /* '<S5>/p' */
};

/* Invariant block signals (default storage) */
struct ConstB_MARS_T {
  real_T Product2[9];                  /* '<S4>/Product2' */
  real_T Selector[9];                  /* '<S7>/Selector' */
  real_T Selector2[9];                 /* '<S7>/Selector2' */
  real_T Product[3];                   /* '<S1>/Product' */
  real_T Abs2[9];                      /* '<S96>/Abs2' */
  real_T Product_n[3];                 /* '<S71>/Product' */
  real_T Product_l[3];                 /* '<S73>/Product' */
  real_T Product_d[3];                 /* '<S72>/Product' */
  real_T Add_o;                        /* '<S90>/Add' */
  real_T Add_jx;                       /* '<S91>/Add' */
  real_T Add_n;                        /* '<S92>/Add' */
};

#ifndef ODE3_INTG
#define ODE3_INTG

/* ODE3 Integration Data */
struct ODE3_IntgData {
  real_T *y;                           /* output */
  real_T *f[3];                        /* derivatives */
};

#endif

/* Constant parameters (default storage) */
struct ConstP_MARS_T {
  /* Pooled Parameter (Expression: -eye(3))
   * Referenced by:
   *   '<S58>/Bias1'
   *   '<S96>/Bias1'
   */
  real_T pooled5[9];
};

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
  real_T odeY[17];
  real_T odeF[3][17];
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
    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *t;
    time_T tArray[2];
  } Timing;
};

extern const ConstB_MARS_T MARS_ConstB;/* constant block i/o */

/* Constant parameters (default storage) */
extern const ConstP_MARS_T MARS_ConstP;

/* Class declaration for model MARS */
class MARS final
{
  /* public data and function members */
 public:
  /* Copy Constructor */
  MARS(MARS const&) = delete;

  /* Assignment Operator */
  MARS& operator= (MARS const&) & = delete;

  /* Move Constructor */
  MARS(MARS &&) = delete;

  /* Move Assignment Operator */
  MARS& operator= (MARS &&) = delete;

  /* Real-Time Model get method */
  RT_MODEL_MARS_T * getRTM();

  /* model start function */
  void start();

  /* Initial conditions function */
  void initialize();

  /* model step function */
  void step();

  /* model terminate function */
  static void terminate();

  /* Constructor */
  MARS();

  /* Destructor */
  ~MARS();

  /* private data and function members */
 private:
  /* Block signals */
  B_MARS_T MARS_B;

  /* Block states */
  DW_MARS_T MARS_DW;

  /* Block continuous states */
  X_MARS_T MARS_X;

  /* private member function(s) for subsystem '<S51>/If Action Subsystem2'*/
  static void MARS_IfActionSubsystem2(real_T rtu_In, real_T *rty_OutOrig);

  /* private member function(s) for subsystem '<S55>/If Not Proper'*/
  static void MARS_IfNotProper(real_T rtp_action);

  /* private member function(s) for subsystem '<S55>/Else If Not Orthogonal'*/
  static void MARS_ElseIfNotOrthogonal(real_T rtp_action);

  /* Continuous states update member function*/
  void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si );

  /* Derivatives member function */
  void MARS_derivatives();

  /* Real-Time Model */
  RT_MODEL_MARS_T MARS_M;
};

/*-
 * These blocks were eliminated from the model due to optimizations:
 *
 * Block '<S1>/Math Function1' : Unused code path elimination
 * Block '<S1>/Product1' : Unused code path elimination
 * Block '<S1>/Sum1' : Unused code path elimination
 * Block '<S11>/Unit Conversion' : Unused code path elimination
 * Block '<S12>/Unit Conversion' : Unused code path elimination
 * Block '<S31>/Reshape (9) to [3x3] column-major' : Reshape block reduction
 * Block '<S42>/Reshape (9) to [3x3] column-major' : Reshape block reduction
 * Block '<S45>/Reshape' : Reshape block reduction
 * Block '<S45>/Reshape1' : Reshape block reduction
 * Block '<S45>/Reshape2' : Reshape block reduction
 * Block '<S58>/Reshape' : Reshape block reduction
 * Block '<S17>/Reshape 3x3 -> 9' : Reshape block reduction
 * Block '<S96>/Reshape' : Reshape block reduction
 * Block '<S103>/Reshape' : Reshape block reduction
 * Block '<S114>/Reshape (9) to [3x3] column-major' : Reshape block reduction
 * Block '<S118>/Reshape (9) to [3x3] column-major' : Reshape block reduction
 * Block '<S135>/Reshape1' : Reshape block reduction
 * Block '<S135>/Reshape2' : Reshape block reduction
 * Block '<S136>/Reshape1' : Reshape block reduction
 * Block '<S136>/Reshape2' : Reshape block reduction
 * Block '<S7>/Reshape' : Reshape block reduction
 * Block '<S7>/Reshape1' : Reshape block reduction
 * Block '<S139>/Reshape (9) to [3x3] column-major' : Reshape block reduction
 * Block '<S2>/Rate Transition' : Eliminated since input and output rates are identical
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
 * '<S1>'   : 'MARS/6DOF ECEF (Quaternion)'
 * '<S2>'   : 'MARS/FlightGear Preconfigured 6DoF Animation'
 * '<S3>'   : 'MARS/6DOF ECEF (Quaternion)/Angle Conversion'
 * '<S4>'   : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles'
 * '<S5>'   : 'MARS/6DOF ECEF (Quaternion)/Calculate Position in EI '
 * '<S6>'   : 'MARS/6DOF ECEF (Quaternion)/Calculate Velocity in Body Axes'
 * '<S7>'   : 'MARS/6DOF ECEF (Quaternion)/Calculate omega_dot'
 * '<S8>'   : 'MARS/6DOF ECEF (Quaternion)/Celestial Longitude of Greenwich'
 * '<S9>'   : 'MARS/6DOF ECEF (Quaternion)/Determine Force,  Mass & Inertia'
 * '<S10>'  : 'MARS/6DOF ECEF (Quaternion)/ECEF to Inertial'
 * '<S11>'  : 'MARS/6DOF ECEF (Quaternion)/Velocity Conversion'
 * '<S12>'  : 'MARS/6DOF ECEF (Quaternion)/Velocity Conversion2'
 * '<S13>'  : 'MARS/6DOF ECEF (Quaternion)/w_ned'
 * '<S14>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED'
 * '<S15>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED1'
 * '<S16>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles'
 * '<S17>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions'
 * '<S18>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix'
 * '<S19>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Rotation Angles to Direction Cosine Matrix'
 * '<S20>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/qdot'
 * '<S21>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED/A11'
 * '<S22>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED/A12'
 * '<S23>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED/A13'
 * '<S24>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED/A21'
 * '<S25>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED/A22'
 * '<S26>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED/A23'
 * '<S27>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED/A31'
 * '<S28>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED/A32'
 * '<S29>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED/A33'
 * '<S30>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED/Angle Conversion'
 * '<S31>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED/Create Transformation Matrix'
 * '<S32>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED1/A11'
 * '<S33>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED1/A12'
 * '<S34>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED1/A13'
 * '<S35>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED1/A21'
 * '<S36>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED1/A22'
 * '<S37>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED1/A23'
 * '<S38>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED1/A31'
 * '<S39>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED1/A32'
 * '<S40>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED1/A33'
 * '<S41>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED1/Angle Conversion'
 * '<S42>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix ECEF to NED1/Create Transformation Matrix'
 * '<S43>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/AxisRotDefault'
 * '<S44>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/AxisRotZeroR3'
 * '<S45>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/Get DCM Values'
 * '<S46>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/Validate DCM'
 * '<S47>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/AxisRotDefault/Protect asincos input'
 * '<S48>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/AxisRotDefault/Protect asincos input/If Action Subsystem'
 * '<S49>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/AxisRotDefault/Protect asincos input/If Action Subsystem1'
 * '<S50>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/AxisRotDefault/Protect asincos input/If Action Subsystem2'
 * '<S51>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/AxisRotZeroR3/Protect asincos input'
 * '<S52>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/AxisRotZeroR3/Protect asincos input/If Action Subsystem'
 * '<S53>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/AxisRotZeroR3/Protect asincos input/If Action Subsystem1'
 * '<S54>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/AxisRotZeroR3/Protect asincos input/If Action Subsystem2'
 * '<S55>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/Validate DCM/If Warning//Error'
 * '<S56>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/Validate DCM/If Warning//Error/Else If Not Orthogonal'
 * '<S57>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/Validate DCM/If Warning//Error/If Not Proper'
 * '<S58>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/Validate DCM/If Warning//Error/isNotOrthogonal'
 * '<S59>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/Validate DCM/If Warning//Error/isNotProper'
 * '<S60>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/Validate DCM/If Warning//Error/Else If Not Orthogonal/Error'
 * '<S61>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/Validate DCM/If Warning//Error/Else If Not Orthogonal/Warning'
 * '<S62>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/Validate DCM/If Warning//Error/If Not Proper/Error'
 * '<S63>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/Validate DCM/If Warning//Error/If Not Proper/Warning'
 * '<S64>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/Validate DCM/If Warning//Error/isNotOrthogonal/transpose*dcm ~= eye(3)'
 * '<S65>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/Validate DCM/If Warning//Error/isNotProper/Determinant of 3x3 Matrix'
 * '<S66>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix to Rotation Angles/Validate DCM/If Warning//Error/isNotProper/determinant does not equal 1'
 * '<S67>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace'
 * '<S68>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Positive Trace'
 * '<S69>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Validate DCM'
 * '<S70>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/trace(DCM)'
 * '<S71>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(1,1)'
 * '<S72>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(2,2)'
 * '<S73>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(3,3)'
 * '<S74>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/diag(DCM)'
 * '<S75>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(1,1)/cos(phi)sin(theta)cos(psi) + sin(phi)sin(psi) -sin(theta)'
 * '<S76>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(1,1)/cos(theta)sin(phi) - (cos(phi)sin(theta)sin(psi) - sin(phi)cos(psi))'
 * '<S77>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(1,1)/cos(theta)sin(psi) + (sin(phi)sin(theta)cos(psi) - cos(phi)sin(psi))'
 * '<S78>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(1,1)/if s~=0; s=0.5//s'
 * '<S79>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(1,1)/u(1) -(u(5)+u(9)) +1'
 * '<S80>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(2,2)/cos(phi)sin(theta)cos(psi) + sin(phi)sin(psi) +sin(theta)'
 * '<S81>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(2,2)/cos(theta)sin(phi) + (cos(phi)sin(theta)sin(psi) - sin(phi)cos(psi))'
 * '<S82>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(2,2)/cos(theta)sin(psi) + (sin(phi)sin(theta)cos(psi) - cos(phi)sin(psi))'
 * '<S83>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(2,2)/if s~=0; s=0.5//s'
 * '<S84>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(2,2)/u(5) -(u(1)+u(9)) +1'
 * '<S85>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(3,3)/cos(phi)sin(theta)cos(psi) + sin(phi)sin(psi) -sin(theta)'
 * '<S86>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(3,3)/cos(theta)sin(phi) + (cos(phi)sin(theta)sin(psi) - sin(phi)cos(psi))'
 * '<S87>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(3,3)/cos(theta)sin(psi) - (sin(phi)sin(theta)cos(psi) - cos(phi)sin(psi))'
 * '<S88>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(3,3)/if s~=0; s=0.5//s'
 * '<S89>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Negative Trace/Maximum Value at DCM(3,3)/u(9) -(u(1)+u(5)) +1'
 * '<S90>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Positive Trace/cos(phi)sin(theta)cos(psi) + sin(phi)sin(psi) +sin(theta)'
 * '<S91>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Positive Trace/cos(theta)sin(phi) - (cos(phi)sin(theta)sin(psi) - sin(phi)cos(psi))'
 * '<S92>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Positive Trace/cos(theta)sin(psi) - (sin(phi)sin(theta)cos(psi) - cos(phi)sin(psi))'
 * '<S93>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error'
 * '<S94>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/Else If Not Orthogonal'
 * '<S95>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/If Not Proper'
 * '<S96>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/isNotOrthogonal'
 * '<S97>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/isNotProper'
 * '<S98>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/Else If Not Orthogonal/Error'
 * '<S99>'  : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/Else If Not Orthogonal/Warning'
 * '<S100>' : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/If Not Proper/Error'
 * '<S101>' : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/If Not Proper/Warning'
 * '<S102>' : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/isNotOrthogonal/transpose*dcm ~= eye(3)'
 * '<S103>' : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/isNotProper/Determinant of 3x3 Matrix'
 * '<S104>' : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Direction Cosine Matrix  to Quaternions/Validate DCM/If Warning//Error/isNotProper/determinant does not equal 1'
 * '<S105>' : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/A11'
 * '<S106>' : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/A12'
 * '<S107>' : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/A13'
 * '<S108>' : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/A21'
 * '<S109>' : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/A22'
 * '<S110>' : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/A23'
 * '<S111>' : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/A31'
 * '<S112>' : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/A32'
 * '<S113>' : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/A33'
 * '<S114>' : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S115>' : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/Quaternion Normalize'
 * '<S116>' : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/Quaternion Normalize/Quaternion Modulus'
 * '<S117>' : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Quaternions to  Direction Cosine Matrix/Quaternion Normalize/Quaternion Modulus/Quaternion Norm'
 * '<S118>' : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/Rotation Angles to Direction Cosine Matrix/Create 3x3 Matrix'
 * '<S119>' : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/qdot/Quaternion Normalize'
 * '<S120>' : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/qdot/Quaternion Normalize/Quaternion Modulus'
 * '<S121>' : 'MARS/6DOF ECEF (Quaternion)/Calculate DCM & Euler Angles/qdot/Quaternion Normalize/Quaternion Modulus/Quaternion Norm'
 * '<S122>' : 'MARS/6DOF ECEF (Quaternion)/Calculate Position in EI /pxwe'
 * '<S123>' : 'MARS/6DOF ECEF (Quaternion)/Calculate Position in EI /pxwe/Subsystem'
 * '<S124>' : 'MARS/6DOF ECEF (Quaternion)/Calculate Position in EI /pxwe/Subsystem1'
 * '<S125>' : 'MARS/6DOF ECEF (Quaternion)/Calculate Velocity in Body Axes/Vbxwb'
 * '<S126>' : 'MARS/6DOF ECEF (Quaternion)/Calculate Velocity in Body Axes/wex(wexp)'
 * '<S127>' : 'MARS/6DOF ECEF (Quaternion)/Calculate Velocity in Body Axes/wexp'
 * '<S128>' : 'MARS/6DOF ECEF (Quaternion)/Calculate Velocity in Body Axes/Vbxwb/Subsystem'
 * '<S129>' : 'MARS/6DOF ECEF (Quaternion)/Calculate Velocity in Body Axes/Vbxwb/Subsystem1'
 * '<S130>' : 'MARS/6DOF ECEF (Quaternion)/Calculate Velocity in Body Axes/wex(wexp)/Subsystem'
 * '<S131>' : 'MARS/6DOF ECEF (Quaternion)/Calculate Velocity in Body Axes/wex(wexp)/Subsystem1'
 * '<S132>' : 'MARS/6DOF ECEF (Quaternion)/Calculate Velocity in Body Axes/wexp/Subsystem'
 * '<S133>' : 'MARS/6DOF ECEF (Quaternion)/Calculate Velocity in Body Axes/wexp/Subsystem1'
 * '<S134>' : 'MARS/6DOF ECEF (Quaternion)/Calculate omega_dot/3x3 Cross Product'
 * '<S135>' : 'MARS/6DOF ECEF (Quaternion)/Calculate omega_dot/I x w'
 * '<S136>' : 'MARS/6DOF ECEF (Quaternion)/Calculate omega_dot/I x w1'
 * '<S137>' : 'MARS/6DOF ECEF (Quaternion)/Calculate omega_dot/3x3 Cross Product/Subsystem'
 * '<S138>' : 'MARS/6DOF ECEF (Quaternion)/Calculate omega_dot/3x3 Cross Product/Subsystem1'
 * '<S139>' : 'MARS/6DOF ECEF (Quaternion)/ECEF to Inertial/Create 3x3 Matrix'
 * '<S140>' : 'MARS/6DOF ECEF (Quaternion)/w_ned/Angle Conversion'
 * '<S141>' : 'MARS/6DOF ECEF (Quaternion)/w_ned/M+h'
 * '<S142>' : 'MARS/6DOF ECEF (Quaternion)/w_ned/N+h'
 * '<S143>' : 'MARS/6DOF ECEF (Quaternion)/w_ned/e2'
 * '<S144>' : 'MARS/FlightGear Preconfigured 6DoF Animation/Angle Conversion'
 * '<S145>' : 'MARS/FlightGear Preconfigured 6DoF Animation/Angle Conversion1'
 */
#endif                                 /* RTW_HEADER_MARS_h_ */
