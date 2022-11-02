/*
 * MARS.c
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

#include "MARS.h"
#include "rtwtypes.h"
#include "MARS_private.h"
#include "rt_assert.h"
#include <math.h>
#include "rt_nonfinite.h"
#include <emmintrin.h>
#include <string.h>
#include "rt_defines.h"
#include <float.h>

/* Block signals (default storage) */
B_MARS_T MARS_B;

/* Continuous states */
X_MARS_T MARS_X;

/* Block states (default storage) */
DW_MARS_T MARS_DW;

/* Real-time model */
static RT_MODEL_MARS_T MARS_M_;
RT_MODEL_MARS_T *const MARS_M = &MARS_M_;
static void rate_scheduler(void);

/*
 *         This function updates active task flag for each subrate.
 *         The function is called at model base rate, hence the
 *         generated code self-manages all its subrates.
 */
static void rate_scheduler(void)
{
  /* Compute which subrates run during the next base time step.  Subrates
   * are an integer multiple of the base rate counter.  Therefore, the subtask
   * counter is reset when it reaches its limit (zero means run).
   */
  (MARS_M->Timing.TaskCounters.TID[2])++;
  if ((MARS_M->Timing.TaskCounters.TID[2]) > 5) {/* Sample time: [0.1s, 0.0s] */
    MARS_M->Timing.TaskCounters.TID[2] = 0;
  }
}

/*
 * This function updates continuous states using the ODE3 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  /* Solver Matrices */
  static const real_T rt_ODE3_A[3] = {
    1.0/2.0, 3.0/4.0, 1.0
  };

  static const real_T rt_ODE3_B[3][3] = {
    { 1.0/2.0, 0.0, 0.0 },

    { 0.0, 3.0/4.0, 0.0 },

    { 2.0/9.0, 1.0/3.0, 4.0/9.0 }
  };

  time_T t = rtsiGetT(si);
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE3_IntgData *id = (ODE3_IntgData *)rtsiGetSolverData(si);
  real_T *y = id->y;
  real_T *f0 = id->f[0];
  real_T *f1 = id->f[1];
  real_T *f2 = id->f[2];
  real_T hB[3];
  int_T i;
  int_T nXc = 19;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) memcpy(y, x,
                (uint_T)nXc*sizeof(real_T));

  /* Assumes that rtsiSetT and ModelOutputs are up-to-date */
  /* f0 = f(t,y) */
  rtsiSetdX(si, f0);
  MARS_derivatives();

  /* f(:,2) = feval(odefile, t + hA(1), y + f*hB(:,1), args(:)(*)); */
  hB[0] = h * rt_ODE3_B[0][0];
  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[0]);
  rtsiSetdX(si, f1);
  MARS_step();
  MARS_derivatives();

  /* f(:,3) = feval(odefile, t + hA(2), y + f*hB(:,2), args(:)(*)); */
  for (i = 0; i <= 1; i++) {
    hB[i] = h * rt_ODE3_B[1][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1]);
  }

  rtsiSetT(si, t + h*rt_ODE3_A[1]);
  rtsiSetdX(si, f2);
  MARS_step();
  MARS_derivatives();

  /* tnew = t + hA(3);
     ynew = y + f*hB(:,3); */
  for (i = 0; i <= 2; i++) {
    hB[i] = h * rt_ODE3_B[2][i];
  }

  for (i = 0; i < nXc; i++) {
    x[i] = y[i] + (f0[i]*hB[0] + f1[i]*hB[1] + f2[i]*hB[2]);
  }

  rtsiSetT(si, tnew);
  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/*
 * Output and update for action system:
 *    '<S108>/If Action Subsystem2'
 *    '<S104>/If Action Subsystem2'
 */
void MARS_IfActionSubsystem2(real_T rtu_In, real_T *rty_OutOrig)
{
  /* SignalConversion generated from: '<S111>/In' */
  *rty_OutOrig = rtu_In;
}

/*
 * Output and update for action system:
 *    '<S112>/If Not Proper'
 *    '<S150>/If Not Proper'
 */
void MARS_IfNotProper(real_T rtp_action)
{
  /* If: '<S114>/If' incorporates:
   *  Constant: '<S114>/Constant'
   */
  if (rtp_action == 2.0) {
    /* Outputs for IfAction SubSystem: '<S114>/Warning' incorporates:
     *  ActionPort: '<S120>/Action Port'
     */
    /* Assertion: '<S120>/Assertion' incorporates:
     *  Constant: '<S114>/Constant1'
     */
    utAssert(false);

    /* End of Outputs for SubSystem: '<S114>/Warning' */
  } else if (rtp_action == 3.0) {
    /* Outputs for IfAction SubSystem: '<S114>/Error' incorporates:
     *  ActionPort: '<S119>/Action Port'
     */
    /* Assertion: '<S119>/Assertion' incorporates:
     *  Constant: '<S114>/Constant1'
     */
    utAssert(false);

    /* End of Outputs for SubSystem: '<S114>/Error' */
  }

  /* End of If: '<S114>/If' */
}

/*
 * Output and update for action system:
 *    '<S112>/Else If Not Orthogonal'
 *    '<S150>/Else If Not Orthogonal'
 */
void MARS_ElseIfNotOrthogonal(real_T rtp_action)
{
  /* If: '<S113>/If' incorporates:
   *  Constant: '<S113>/Constant'
   */
  if (rtp_action == 2.0) {
    /* Outputs for IfAction SubSystem: '<S113>/Warning' incorporates:
     *  ActionPort: '<S118>/Action Port'
     */
    /* Assertion: '<S118>/Assertion' incorporates:
     *  Constant: '<S113>/Constant1'
     */
    utAssert(false);

    /* End of Outputs for SubSystem: '<S113>/Warning' */
  } else if (rtp_action == 3.0) {
    /* Outputs for IfAction SubSystem: '<S113>/Error' incorporates:
     *  ActionPort: '<S117>/Action Port'
     */
    /* Assertion: '<S117>/Assertion' incorporates:
     *  Constant: '<S113>/Constant1'
     */
    utAssert(false);

    /* End of Outputs for SubSystem: '<S113>/Error' */
  }

  /* End of If: '<S113>/If' */
}

real_T rt_urand_Upu32_Yd_f_pw_snf(uint32_T *u)
{
  uint32_T hi;
  uint32_T lo;

  /* Uniform random number generator (random number between 0 and 1)

     #define IA      16807                      magic multiplier = 7^5
     #define IM      2147483647                 modulus = 2^31-1
     #define IQ      127773                     IM div IA
     #define IR      2836                       IM modulo IA
     #define S       4.656612875245797e-10      reciprocal of 2^31-1
     test = IA * (seed % IQ) - IR * (seed/IQ)
     seed = test < 0 ? (test + IM) : test
     return (seed*S)
   */
  lo = *u % 127773U * 16807U;
  hi = *u / 127773U * 2836U;
  if (lo < hi) {
    *u = 2147483647U - (hi - lo);
  } else {
    *u = lo - hi;
  }

  return (real_T)*u * 4.6566128752457969E-10;
}

real_T rt_nrand_Upu32_Yd_f_pw_snf(uint32_T *u)
{
  real_T si;
  real_T sr;
  real_T y;

  /* Normal (Gaussian) random number generator */
  do {
    sr = 2.0 * rt_urand_Upu32_Yd_f_pw_snf(u) - 1.0;
    si = 2.0 * rt_urand_Upu32_Yd_f_pw_snf(u) - 1.0;
    si = sr * sr + si * si;
  } while (si > 1.0);

  y = sqrt(-2.0 * log(si) / si) * sr;
  return y;
}

real_T rt_atan2d_snf(real_T u0, real_T u1)
{
  real_T y;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else if (rtIsInf(u0) && rtIsInf(u1)) {
    int32_T u0_0;
    int32_T u1_0;
    if (u0 > 0.0) {
      u0_0 = 1;
    } else {
      u0_0 = -1;
    }

    if (u1 > 0.0) {
      u1_0 = 1;
    } else {
      u1_0 = -1;
    }

    y = atan2(u0_0, u1_0);
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = atan2(u0, u1);
  }

  return y;
}

real_T rt_modd_snf(real_T u0, real_T u1)
{
  real_T y;
  y = u0;
  if (u1 == 0.0) {
    if (u0 == 0.0) {
      y = u1;
    }
  } else if (rtIsNaN(u0) || rtIsNaN(u1) || rtIsInf(u0)) {
    y = (rtNaN);
  } else if (u0 == 0.0) {
    y = 0.0 / u1;
  } else if (rtIsInf(u1)) {
    if ((u1 < 0.0) != (u0 < 0.0)) {
      y = u1;
    }
  } else {
    boolean_T yEq;
    y = fmod(u0, u1);
    yEq = (y == 0.0);
    if ((!yEq) && (u1 > floor(u1))) {
      real_T q;
      q = fabs(u0 / u1);
      yEq = !(fabs(q - floor(q + 0.5)) > DBL_EPSILON * q);
    }

    if (yEq) {
      y = u1 * 0.0;
    } else if ((u0 < 0.0) != (u1 < 0.0)) {
      y += u1;
    }
  }

  return y;
}

real_T rt_remd_snf(real_T u0, real_T u1)
{
  real_T y;
  if (rtIsNaN(u0) || rtIsNaN(u1) || rtIsInf(u0)) {
    y = (rtNaN);
  } else if (rtIsInf(u1)) {
    y = u0;
  } else if ((u1 != 0.0) && (u1 != trunc(u1))) {
    real_T q;
    q = fabs(u0 / u1);
    if (!(fabs(q - floor(q + 0.5)) > DBL_EPSILON * q)) {
      y = 0.0 * u0;
    } else {
      y = fmod(u0, u1);
    }
  } else {
    y = fmod(u0, u1);
  }

  return y;
}

real_T rt_powd_snf(real_T u0, real_T u1)
{
  real_T y;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else {
    real_T tmp;
    real_T tmp_0;
    tmp = fabs(u0);
    tmp_0 = fabs(u1);
    if (rtIsInf(u1)) {
      if (tmp == 1.0) {
        y = 1.0;
      } else if (tmp > 1.0) {
        if (u1 > 0.0) {
          y = (rtInf);
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = (rtInf);
      }
    } else if (tmp_0 == 0.0) {
      y = 1.0;
    } else if (tmp_0 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = (rtNaN);
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

void rt_mrdivide_U1d1x3_U2d_9vOrDY9Z(const real_T u0[3], const real_T u1[9],
  real_T y[3])
{
  real_T A[9];
  real_T a21;
  real_T maxval;
  int32_T r1;
  int32_T r2;
  int32_T r3;
  memcpy(&A[0], &u1[0], 9U * sizeof(real_T));
  r1 = 0;
  r2 = 1;
  r3 = 2;
  maxval = fabs(u1[0]);
  a21 = fabs(u1[1]);
  if (a21 > maxval) {
    maxval = a21;
    r1 = 1;
    r2 = 0;
  }

  if (fabs(u1[2]) > maxval) {
    r1 = 2;
    r2 = 1;
    r3 = 0;
  }

  A[r2] = u1[r2] / u1[r1];
  A[r3] /= A[r1];
  A[r2 + 3] -= A[r1 + 3] * A[r2];
  A[r3 + 3] -= A[r1 + 3] * A[r3];
  A[r2 + 6] -= A[r1 + 6] * A[r2];
  A[r3 + 6] -= A[r1 + 6] * A[r3];
  if (fabs(A[r3 + 3]) > fabs(A[r2 + 3])) {
    int32_T rtemp;
    rtemp = r2 + 1;
    r2 = r3;
    r3 = rtemp - 1;
  }

  A[r3 + 3] /= A[r2 + 3];
  A[r3 + 6] -= A[r3 + 3] * A[r2 + 6];
  y[r1] = u0[0] / A[r1];
  y[r2] = u0[1] - A[r1 + 3] * y[r1];
  y[r3] = u0[2] - A[r1 + 6] * y[r1];
  y[r2] /= A[r2 + 3];
  y[r3] -= A[r2 + 6] * y[r2];
  y[r3] /= A[r3 + 6];
  y[r2] -= A[r3 + 3] * y[r3];
  y[r1] -= y[r3] * A[r3];
  y[r1] -= y[r2] * A[r2];
}

/* Model step function */
void MARS_step(void)
{
  /* local block i/o variables */
  real_T rtb_VectorConcatenate[7];
  __m128d tmp_1;
  __m128d tmp_2;
  __m128d tmp_3;
  __m128d tmp_4;
  real_T rtb_Product4_g_tmp[9];
  real_T rtb_VectorConcatenate_b[9];
  real_T rtb_VectorConcatenate_f[9];
  real_T rtb_Gain_a[3];
  real_T rtb_Sum2[3];
  real_T rtb_Sum2_0[3];
  real_T rtb_Sum_hp[3];
  real_T rtb_ECEFPositiontoLLA_o1_idx_0;
  real_T rtb_ECEFPositiontoLLA_o2;
  real_T rtb_Gain_f;
  real_T rtb_VectorConcatenate_f_tmp;
  real_T rtb_VectorConcatenate_f_tmp_0;
  real_T rtb_VectorConcatenate_f_tmp_1;
  real_T rtb_VectorConcatenate_f_tmp_2;
  real_T rtb_VectorConcatenate_f_tmp_3;
  real_T rtb_VectorConcatenate_f_tmp_4;
  real_T rtb_VectorConcatenate_f_tmp_5;
  real_T rtb_ixk;
  real_T rtb_jxi;
  real_T rtb_kxj;
  real_T rtb_sincos_o2_idx_1;
  real_T uTmp_idx_0;
  real_T uTmp_idx_1;
  int32_T Product4_tmp;
  int32_T i;
  int32_T i_0;
  int8_T rtAction;
  boolean_T rtb_Compare_a[9];
  if (rtmIsMajorTimeStep(MARS_M)) {
    /* set solver stop time */
    rtsiSetSolverStopTime(&MARS_M->solverInfo,((MARS_M->Timing.clockTick0+1)*
      MARS_M->Timing.stepSize0));
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(MARS_M)) {
    MARS_M->Timing.t[0] = rtsiGetT(&MARS_M->solverInfo);
  }

  /* Integrator: '<S62>/p1' */
  if (MARS_DW.p1_IWORK != 0) {
    MARS_X.p1_CSTATE[0] = 2.5914814724705652E+6;
    MARS_X.p1_CSTATE[1] = -1.0790161949738523E+6;
    MARS_X.p1_CSTATE[2] = 5.7079879931331929E+6;
  }

  /* ECEF2LLA: '<S57>/ECEF Position to LLA' incorporates:
   *  Integrator: '<S62>/p1'
   */
  uTmp_idx_0 = MARS_X.p1_CSTATE[0];
  uTmp_idx_1 = MARS_X.p1_CSTATE[1];
  rtb_jxi = sqrt(MARS_X.p1_CSTATE[0] * MARS_X.p1_CSTATE[0] + MARS_X.p1_CSTATE[1]
                 * MARS_X.p1_CSTATE[1]);
  rtb_ixk = rt_atan2d_snf(MARS_X.p1_CSTATE[2], 0.99664718933525254 * rtb_jxi);
  rtb_Gain_f = sin(rtb_ixk);
  rtb_ECEFPositiontoLLA_o2 = cos(rtb_ixk);
  rtb_Gain_f = rt_atan2d_snf(42841.311513313573 * rtb_Gain_f * rtb_Gain_f *
    rtb_Gain_f + MARS_X.p1_CSTATE[2], rtb_jxi - 42697.672707179969 *
    rtb_ECEFPositiontoLLA_o2 * rtb_ECEFPositiontoLLA_o2 *
    rtb_ECEFPositiontoLLA_o2);
  rtb_ECEFPositiontoLLA_o2 = rt_atan2d_snf(0.99664718933525254 * sin(rtb_Gain_f),
    cos(rtb_Gain_f));
  i = 0;
  while ((rtb_ixk != rtb_ECEFPositiontoLLA_o2) && (i < 5)) {
    rtb_ixk = rtb_ECEFPositiontoLLA_o2;
    rtb_Gain_f = sin(rtb_ECEFPositiontoLLA_o2);
    rtb_ECEFPositiontoLLA_o2 = cos(rtb_ECEFPositiontoLLA_o2);
    rtb_Gain_f = rt_atan2d_snf(42841.311513313573 * rtb_Gain_f * rtb_Gain_f *
      rtb_Gain_f + MARS_X.p1_CSTATE[2], sqrt(uTmp_idx_0 * uTmp_idx_0 +
      uTmp_idx_1 * uTmp_idx_1) - 42697.672707179969 * rtb_ECEFPositiontoLLA_o2 *
      rtb_ECEFPositiontoLLA_o2 * rtb_ECEFPositiontoLLA_o2);
    rtb_ECEFPositiontoLLA_o2 = rt_atan2d_snf(0.99664718933525254 * sin
      (rtb_Gain_f), cos(rtb_Gain_f));
    i++;
  }

  rtb_kxj = fabs(rtb_Gain_f);
  rtb_ixk = rtb_Gain_f;
  rtb_ECEFPositiontoLLA_o2 = rt_atan2d_snf(MARS_X.p1_CSTATE[1],
    MARS_X.p1_CSTATE[0]);
  if (rtb_kxj > 3.1415926535897931) {
    rtb_ixk = rt_modd_snf(rtb_Gain_f + 3.1415926535897931, 6.2831853071795862) -
      3.1415926535897931;
    rtb_kxj = fabs(rtb_ixk);
  }

  if (rtb_kxj > 1.5707963267948966) {
    rtb_ECEFPositiontoLLA_o2 += 3.1415926535897931;
    if (!rtIsNaN(rtb_ixk)) {
      if (rtb_ixk < 0.0) {
        rtb_ixk = -1.0;
      } else {
        rtb_ixk = (rtb_ixk > 0.0);
      }
    }

    rtb_ixk *= 1.5707963267948966 - (rtb_kxj - 1.5707963267948966);
  }

  if (fabs(rtb_ECEFPositiontoLLA_o2) > 3.1415926535897931) {
    rtb_ECEFPositiontoLLA_o2 = rt_remd_snf(rtb_ECEFPositiontoLLA_o2,
      6.2831853071795862);
    rtb_ECEFPositiontoLLA_o2 -= trunc(rtb_ECEFPositiontoLLA_o2 /
      3.1415926535897931) * 6.2831853071795862;
  }

  rtb_ECEFPositiontoLLA_o1_idx_0 = rtb_ixk * 180.0 / 3.1415926535897931;
  rtb_kxj = rtb_ECEFPositiontoLLA_o2 * 180.0 / 3.1415926535897931;
  rtb_ixk = sin(rtb_Gain_f);
  rtb_ECEFPositiontoLLA_o2 = ((6.378137E+6 / sqrt(1.0 - rtb_ixk * rtb_ixk *
    0.0066943799901413165) * 0.0066943799901413165 * rtb_ixk + MARS_X.p1_CSTATE
    [2]) * rtb_ixk + rtb_jxi * cos(rtb_Gain_f)) - 6.378137E+6 / sqrt(1.0 - sin
    (rtb_Gain_f) * sin(rtb_Gain_f) * 0.0066943799901413165);

  /* End of ECEF2LLA: '<S57>/ECEF Position to LLA' */
  if (rtmIsMajorTimeStep(MARS_M) &&
      MARS_M->Timing.TaskCounters.TID[1] == 0) {
    /* If: '<S74>/If' */
    if (rtsiIsModeUpdateTimeStep(&MARS_M->solverInfo)) {
      rtAction = 1;
      MARS_DW.If_ActiveSubsystem = 1;
    } else {
      rtAction = MARS_DW.If_ActiveSubsystem;
    }

    switch (rtAction) {
     case 0:
      /* Outputs for IfAction SubSystem: '<S74>/Positive Trace' incorporates:
       *  ActionPort: '<S125>/Action Port'
       */
      /* Gain: '<S125>/Gain' incorporates:
       *  Merge: '<S74>/Merge'
       */
      MARS_B.Merge[0] = 0.22088421489473928;

      /* Product: '<S125>/Product' incorporates:
       *  Merge: '<S74>/Merge'
       */
      MARS_B.Merge[1] = MARS_ConstB.Add_jx / 0.88353685957895711;
      MARS_B.Merge[2] = MARS_ConstB.Add_o / 0.88353685957895711;
      MARS_B.Merge[3] = MARS_ConstB.Add_n / 0.88353685957895711;

      /* End of Outputs for SubSystem: '<S74>/Positive Trace' */
      break;

     case 1:
      /* Outputs for IfAction SubSystem: '<S74>/Negative Trace' incorporates:
       *  ActionPort: '<S124>/Action Port'
       */
      /* If: '<S124>/Find Maximum Diagonal Value' */
      if ((MARS_ConstB.Product2[4] > MARS_ConstB.Product2[0]) &&
          (MARS_ConstB.Product2[4] > MARS_ConstB.Product2[8])) {
        /* Outputs for IfAction SubSystem: '<S124>/Maximum Value at DCM(2,2)' incorporates:
         *  ActionPort: '<S129>/Action Port'
         */
        /* Gain: '<S129>/Gain' incorporates:
         *  Merge: '<S74>/Merge'
         */
        MARS_B.Merge[2] = 0.955404187046654;

        /* Gain: '<S129>/Gain1' incorporates:
         *  Merge: '<S74>/Merge'
         */
        MARS_B.Merge[1] = MARS_ConstB.Product_d[0];

        /* Gain: '<S129>/Gain3' incorporates:
         *  Merge: '<S74>/Merge'
         */
        MARS_B.Merge[3] = MARS_ConstB.Product_d[1];

        /* Gain: '<S129>/Gain4' incorporates:
         *  Merge: '<S74>/Merge'
         */
        MARS_B.Merge[0] = MARS_ConstB.Product_d[2];

        /* End of Outputs for SubSystem: '<S124>/Maximum Value at DCM(2,2)' */
      } else if (MARS_ConstB.Product2[8] > MARS_ConstB.Product2[0]) {
        /* Outputs for IfAction SubSystem: '<S124>/Maximum Value at DCM(3,3)' incorporates:
         *  ActionPort: '<S130>/Action Port'
         */
        /* Gain: '<S130>/Gain' incorporates:
         *  Merge: '<S74>/Merge'
         */
        MARS_B.Merge[3] = 0.044147846096478843;

        /* Gain: '<S130>/Gain1' incorporates:
         *  Merge: '<S74>/Merge'
         */
        MARS_B.Merge[1] = MARS_ConstB.Product_l[0];

        /* Gain: '<S130>/Gain2' incorporates:
         *  Merge: '<S74>/Merge'
         */
        MARS_B.Merge[2] = MARS_ConstB.Product_l[1];

        /* Gain: '<S130>/Gain3' incorporates:
         *  Merge: '<S74>/Merge'
         */
        MARS_B.Merge[0] = MARS_ConstB.Product_l[2];

        /* End of Outputs for SubSystem: '<S124>/Maximum Value at DCM(3,3)' */
      } else {
        /* Outputs for IfAction SubSystem: '<S124>/Maximum Value at DCM(1,1)' incorporates:
         *  ActionPort: '<S128>/Action Port'
         */
        /* Gain: '<S128>/Gain' incorporates:
         *  Merge: '<S74>/Merge'
         */
        MARS_B.Merge[1] = 0.19095541539610353;

        /* Gain: '<S128>/Gain1' incorporates:
         *  Merge: '<S74>/Merge'
         */
        MARS_B.Merge[2] = MARS_ConstB.Product_n[0];

        /* Gain: '<S128>/Gain2' incorporates:
         *  Merge: '<S74>/Merge'
         */
        MARS_B.Merge[3] = MARS_ConstB.Product_n[1];

        /* Gain: '<S128>/Gain3' incorporates:
         *  Merge: '<S74>/Merge'
         */
        MARS_B.Merge[0] = MARS_ConstB.Product_n[2];

        /* End of Outputs for SubSystem: '<S124>/Maximum Value at DCM(1,1)' */
      }

      /* End of If: '<S124>/Find Maximum Diagonal Value' */
      /* End of Outputs for SubSystem: '<S74>/Negative Trace' */
      break;
    }

    /* End of If: '<S74>/If' */
  }

  /* Integrator: '<S61>/q' */
  if (MARS_DW.q_IWORK != 0) {
    MARS_X.q_CSTATE[0] = MARS_B.Merge[0];
    MARS_X.q_CSTATE[1] = MARS_B.Merge[1];
    MARS_X.q_CSTATE[2] = MARS_B.Merge[2];
    MARS_X.q_CSTATE[3] = MARS_B.Merge[3];
  }

  /* Sqrt: '<S173>/sqrt' incorporates:
   *  Integrator: '<S61>/q'
   *  Product: '<S174>/Product'
   *  Product: '<S174>/Product1'
   *  Product: '<S174>/Product2'
   *  Product: '<S174>/Product3'
   *  Sqrt: '<S177>/sqrt'
   *  Sum: '<S174>/Sum'
   */
  rtb_Gain_f = sqrt(((MARS_X.q_CSTATE[0] * MARS_X.q_CSTATE[0] + MARS_X.q_CSTATE
                      [1] * MARS_X.q_CSTATE[1]) + MARS_X.q_CSTATE[2] *
                     MARS_X.q_CSTATE[2]) + MARS_X.q_CSTATE[3] * MARS_X.q_CSTATE
                    [3]);

  /* Product: '<S172>/Product' incorporates:
   *  Integrator: '<S61>/q'
   *  Sqrt: '<S173>/sqrt'
   */
  rtb_ixk = MARS_X.q_CSTATE[0] / rtb_Gain_f;

  /* Product: '<S172>/Product1' incorporates:
   *  Integrator: '<S61>/q'
   *  Sqrt: '<S173>/sqrt'
   */
  rtb_jxi = MARS_X.q_CSTATE[1] / rtb_Gain_f;

  /* Product: '<S172>/Product2' incorporates:
   *  Integrator: '<S61>/q'
   *  Product: '<S176>/Product2'
   *  Sqrt: '<S173>/sqrt'
   */
  uTmp_idx_0 = MARS_X.q_CSTATE[2] / rtb_Gain_f;

  /* Product: '<S172>/Product3' incorporates:
   *  Integrator: '<S61>/q'
   *  Product: '<S176>/Product3'
   *  Sqrt: '<S173>/sqrt'
   */
  uTmp_idx_1 = MARS_X.q_CSTATE[3] / rtb_Gain_f;

  /* Product: '<S162>/Product3' incorporates:
   *  Product: '<S166>/Product3'
   */
  rtb_sincos_o2_idx_1 = rtb_ixk * rtb_ixk;

  /* Product: '<S162>/Product2' incorporates:
   *  Product: '<S166>/Product2'
   */
  rtb_VectorConcatenate_f_tmp_1 = rtb_jxi * rtb_jxi;

  /* Product: '<S162>/Product1' incorporates:
   *  Product: '<S166>/Product1'
   *  Product: '<S170>/Product1'
   *  Product: '<S172>/Product2'
   */
  rtb_VectorConcatenate_f_tmp_2 = uTmp_idx_0 * uTmp_idx_0;

  /* Product: '<S162>/Product' incorporates:
   *  Product: '<S166>/Product'
   *  Product: '<S170>/Product'
   *  Product: '<S172>/Product3'
   */
  rtb_VectorConcatenate_f_tmp_3 = uTmp_idx_1 * uTmp_idx_1;

  /* Sum: '<S162>/Sum' incorporates:
   *  Product: '<S162>/Product'
   *  Product: '<S162>/Product1'
   *  Product: '<S162>/Product2'
   *  Product: '<S162>/Product3'
   */
  rtb_VectorConcatenate_f[0] = ((rtb_sincos_o2_idx_1 +
    rtb_VectorConcatenate_f_tmp_1) - rtb_VectorConcatenate_f_tmp_2) -
    rtb_VectorConcatenate_f_tmp_3;

  /* Product: '<S165>/Product3' incorporates:
   *  Product: '<S163>/Product3'
   *  Product: '<S172>/Product3'
   */
  rtb_VectorConcatenate_f_tmp = uTmp_idx_1 * rtb_ixk;

  /* Product: '<S165>/Product2' incorporates:
   *  Product: '<S163>/Product2'
   *  Product: '<S172>/Product2'
   */
  rtb_VectorConcatenate_f_tmp_0 = rtb_jxi * uTmp_idx_0;

  /* Gain: '<S165>/Gain' incorporates:
   *  Product: '<S165>/Product2'
   *  Product: '<S165>/Product3'
   *  Sum: '<S165>/Sum'
   */
  rtb_VectorConcatenate_f[1] = (rtb_VectorConcatenate_f_tmp_0 -
    rtb_VectorConcatenate_f_tmp) * 2.0;

  /* Product: '<S168>/Product2' incorporates:
   *  Product: '<S164>/Product2'
   *  Product: '<S172>/Product3'
   */
  rtb_VectorConcatenate_f_tmp_4 = rtb_jxi * uTmp_idx_1;

  /* Product: '<S168>/Product1' incorporates:
   *  Product: '<S164>/Product1'
   *  Product: '<S172>/Product2'
   */
  rtb_VectorConcatenate_f_tmp_5 = rtb_ixk * uTmp_idx_0;

  /* Gain: '<S168>/Gain' incorporates:
   *  Product: '<S168>/Product1'
   *  Product: '<S168>/Product2'
   *  Sum: '<S168>/Sum'
   */
  rtb_VectorConcatenate_f[2] = (rtb_VectorConcatenate_f_tmp_5 +
    rtb_VectorConcatenate_f_tmp_4) * 2.0;

  /* Gain: '<S163>/Gain' incorporates:
   *  Sum: '<S163>/Sum'
   */
  rtb_VectorConcatenate_f[3] = (rtb_VectorConcatenate_f_tmp +
    rtb_VectorConcatenate_f_tmp_0) * 2.0;

  /* Sum: '<S166>/Sum' incorporates:
   *  Sum: '<S170>/Sum'
   */
  rtb_sincos_o2_idx_1 -= rtb_VectorConcatenate_f_tmp_1;
  rtb_VectorConcatenate_f[4] = (rtb_sincos_o2_idx_1 +
    rtb_VectorConcatenate_f_tmp_2) - rtb_VectorConcatenate_f_tmp_3;

  /* Product: '<S169>/Product1' incorporates:
   *  Product: '<S167>/Product1'
   */
  rtb_VectorConcatenate_f_tmp_1 = rtb_ixk * rtb_jxi;

  /* Product: '<S169>/Product2' incorporates:
   *  Product: '<S167>/Product2'
   *  Product: '<S172>/Product2'
   *  Product: '<S172>/Product3'
   */
  rtb_VectorConcatenate_f_tmp = uTmp_idx_0 * uTmp_idx_1;

  /* Gain: '<S169>/Gain' incorporates:
   *  Product: '<S169>/Product1'
   *  Product: '<S169>/Product2'
   *  Sum: '<S169>/Sum'
   */
  rtb_VectorConcatenate_f[5] = (rtb_VectorConcatenate_f_tmp -
    rtb_VectorConcatenate_f_tmp_1) * 2.0;

  /* Gain: '<S164>/Gain' incorporates:
   *  Sum: '<S164>/Sum'
   */
  rtb_VectorConcatenate_f[6] = (rtb_VectorConcatenate_f_tmp_4 -
    rtb_VectorConcatenate_f_tmp_5) * 2.0;

  /* Gain: '<S167>/Gain' incorporates:
   *  Sum: '<S167>/Sum'
   */
  rtb_VectorConcatenate_f[7] = (rtb_VectorConcatenate_f_tmp_1 +
    rtb_VectorConcatenate_f_tmp) * 2.0;

  /* Sum: '<S170>/Sum' */
  rtb_VectorConcatenate_f[8] = (rtb_sincos_o2_idx_1 -
    rtb_VectorConcatenate_f_tmp_2) + rtb_VectorConcatenate_f_tmp_3;

  /* UnitConversion: '<S87>/Unit Conversion' incorporates:
   *  UnitConversion: '<S197>/Unit Conversion'
   */
  /* Unit Conversion - from: deg to: rad
     Expression: output = (0.0174533*input) + (0) */
  rtb_ECEFPositiontoLLA_o1_idx_0 *= 0.017453292519943295;

  /* Trigonometry: '<S71>/sincos' incorporates:
   *  Trigonometry: '<S70>/sine'
   *  UnitConversion: '<S87>/Unit Conversion'
   */
  rtb_jxi = cos(rtb_ECEFPositiontoLLA_o1_idx_0);
  rtb_VectorConcatenate_f_tmp_1 = sin(rtb_ECEFPositiontoLLA_o1_idx_0);

  /* UnitConversion: '<S87>/Unit Conversion' */
  rtb_kxj *= 0.017453292519943295;

  /* Trigonometry: '<S71>/sincos' */
  rtb_sincos_o2_idx_1 = cos(rtb_kxj);
  rtb_kxj = sin(rtb_kxj);

  /* UnaryMinus: '<S78>/Unary Minus' incorporates:
   *  Product: '<S78>/u(1)*u(4)'
   *  Trigonometry: '<S71>/sincos'
   */
  rtb_VectorConcatenate_b[0] = -(rtb_VectorConcatenate_f_tmp_1 *
    rtb_sincos_o2_idx_1);

  /* UnaryMinus: '<S81>/Unary Minus' */
  rtb_VectorConcatenate_b[1] = -rtb_kxj;

  /* UnaryMinus: '<S84>/Unary Minus' incorporates:
   *  Product: '<S84>/u(3)*u(4)'
   */
  rtb_VectorConcatenate_b[2] = -(rtb_jxi * rtb_sincos_o2_idx_1);

  /* UnaryMinus: '<S79>/Unary Minus' incorporates:
   *  Product: '<S79>/u(1)*u(2)'
   *  Trigonometry: '<S71>/sincos'
   */
  rtb_VectorConcatenate_b[3] = -(rtb_VectorConcatenate_f_tmp_1 * rtb_kxj);

  /* SignalConversion generated from: '<S88>/Vector Concatenate' */
  rtb_VectorConcatenate_b[4] = rtb_sincos_o2_idx_1;

  /* UnaryMinus: '<S85>/Unary Minus' incorporates:
   *  Product: '<S85>/u(2)*u(3)'
   */
  rtb_VectorConcatenate_b[5] = -(rtb_jxi * rtb_kxj);

  /* SignalConversion generated from: '<S88>/Vector Concatenate' */
  rtb_VectorConcatenate_b[6] = rtb_jxi;

  /* SignalConversion generated from: '<S88>/Vector Concatenate' incorporates:
   *  Constant: '<S83>/Constant'
   */
  rtb_VectorConcatenate_b[7] = 0.0;

  /* UnaryMinus: '<S86>/Unary Minus' incorporates:
   *  Trigonometry: '<S71>/sincos'
   */
  rtb_VectorConcatenate_b[8] = -rtb_VectorConcatenate_f_tmp_1;
  for (i = 0; i < 3; i++) {
    for (i_0 = 0; i_0 <= 0; i_0 += 2) {
      /* Product: '<S61>/Product4' incorporates:
       *  Math: '<S61>/Math Function2'
       */
      Product4_tmp = 3 * i + i_0;
      _mm_storeu_pd(&MARS_B.Product4[Product4_tmp], _mm_set1_pd(0.0));

      /* Product: '<S61>/Product4' incorporates:
       *  Concatenate: '<S171>/Vector Concatenate'
       *  Concatenate: '<S196>/Vector Concatenate'
       *  Math: '<S61>/Math Function2'
       */
      tmp_3 = _mm_loadu_pd(&rtb_VectorConcatenate_f[i_0]);
      tmp_4 = _mm_loadu_pd(&MARS_B.Product4[Product4_tmp]);
      _mm_storeu_pd(&MARS_B.Product4[Product4_tmp], _mm_add_pd(tmp_4, _mm_mul_pd
        (tmp_3, _mm_set1_pd(rtb_VectorConcatenate_b[i]))));
      tmp_3 = _mm_loadu_pd(&rtb_VectorConcatenate_f[i_0 + 3]);
      tmp_4 = _mm_loadu_pd(&MARS_B.Product4[Product4_tmp]);
      _mm_storeu_pd(&MARS_B.Product4[Product4_tmp], _mm_add_pd(tmp_4, _mm_mul_pd
        (tmp_3, _mm_set1_pd(rtb_VectorConcatenate_b[i + 3]))));
      tmp_3 = _mm_loadu_pd(&rtb_VectorConcatenate_f[i_0 + 6]);
      tmp_4 = _mm_loadu_pd(&MARS_B.Product4[Product4_tmp]);
      _mm_storeu_pd(&MARS_B.Product4[Product4_tmp], _mm_add_pd(tmp_4, _mm_mul_pd
        (tmp_3, _mm_set1_pd(rtb_VectorConcatenate_b[i + 6]))));
    }

    for (i_0 = 2; i_0 < 3; i_0++) {
      /* Product: '<S61>/Product4' incorporates:
       *  Concatenate: '<S171>/Vector Concatenate'
       *  Concatenate: '<S196>/Vector Concatenate'
       *  Math: '<S61>/Math Function2'
       */
      Product4_tmp = 3 * i + i_0;
      MARS_B.Product4[Product4_tmp] = 0.0;
      MARS_B.Product4[Product4_tmp] += rtb_VectorConcatenate_f[i_0] *
        rtb_VectorConcatenate_b[i];
      MARS_B.Product4[Product4_tmp] += rtb_VectorConcatenate_f[i_0 + 3] *
        rtb_VectorConcatenate_b[i + 3];
      MARS_B.Product4[Product4_tmp] += rtb_VectorConcatenate_f[i_0 + 6] *
        rtb_VectorConcatenate_b[i + 6];
    }
  }

  /* Gain: '<S102>/Gain1' incorporates:
   *  Concatenate: '<S102>/Vector Concatenate'
   *  Product: '<S61>/Product4'
   *  Selector: '<S102>/Selector1'
   */
  rtb_VectorConcatenate[0] = MARS_B.Product4[3];
  rtb_VectorConcatenate[1] = MARS_B.Product4[0];
  rtb_VectorConcatenate[2] = -MARS_B.Product4[6];

  /* Gain: '<S102>/Gain2' incorporates:
   *  Concatenate: '<S102>/Vector Concatenate'
   *  Product: '<S61>/Product4'
   *  Selector: '<S102>/Selector2'
   */
  rtb_VectorConcatenate[3] = MARS_B.Product4[7];

  /* Gain: '<S102>/Gain3' incorporates:
   *  Concatenate: '<S102>/Vector Concatenate'
   *  Product: '<S61>/Product4'
   *  Selector: '<S102>/Selector3'
   */
  rtb_VectorConcatenate[5] = -MARS_B.Product4[1];

  /* Gain: '<S102>/Gain2' incorporates:
   *  Concatenate: '<S102>/Vector Concatenate'
   *  Product: '<S61>/Product4'
   *  Selector: '<S102>/Selector2'
   */
  rtb_VectorConcatenate[4] = MARS_B.Product4[8];

  /* Gain: '<S102>/Gain3' incorporates:
   *  Concatenate: '<S102>/Vector Concatenate'
   *  Product: '<S61>/Product4'
   *  Selector: '<S102>/Selector3'
   */
  rtb_VectorConcatenate[6] = MARS_B.Product4[4];

  /* If: '<S73>/If' */
  if (rtsiIsModeUpdateTimeStep(&MARS_M->solverInfo)) {
    rtAction = (int8_T)((!(rtb_VectorConcatenate[2] >= 1.0)) &&
                        (!(rtb_VectorConcatenate[2] <= -1.0)));
    MARS_DW.If_ActiveSubsystem_h = rtAction;
  } else {
    rtAction = MARS_DW.If_ActiveSubsystem_h;
  }

  switch (rtAction) {
   case 0:
    /* Outputs for IfAction SubSystem: '<S73>/AxisRotZeroR3' incorporates:
     *  ActionPort: '<S101>/Action Port'
     */
    /* If: '<S108>/If' */
    if (rtsiIsModeUpdateTimeStep(&MARS_M->solverInfo)) {
      if (rtb_VectorConcatenate[2] > 1.0) {
        rtAction = 0;
      } else if (rtb_VectorConcatenate[2] < -1.0) {
        rtAction = 1;
      } else {
        rtAction = 2;
      }

      MARS_DW.If_ActiveSubsystem_c0 = rtAction;
    } else {
      rtAction = MARS_DW.If_ActiveSubsystem_c0;
    }

    switch (rtAction) {
     case 0:
      /* Outputs for IfAction SubSystem: '<S108>/If Action Subsystem' incorporates:
       *  ActionPort: '<S109>/Action Port'
       */
      if (rtmIsMajorTimeStep(MARS_M) &&
          MARS_M->Timing.TaskCounters.TID[1] == 0) {
        /* Merge: '<S108>/Merge' incorporates:
         *  Constant: '<S109>/Constant'
         */
        MARS_B.Merge_j2 = 1.0;
      }

      /* End of Outputs for SubSystem: '<S108>/If Action Subsystem' */
      break;

     case 1:
      /* Outputs for IfAction SubSystem: '<S108>/If Action Subsystem1' incorporates:
       *  ActionPort: '<S110>/Action Port'
       */
      if (rtmIsMajorTimeStep(MARS_M) &&
          MARS_M->Timing.TaskCounters.TID[1] == 0) {
        /* Merge: '<S108>/Merge' incorporates:
         *  Constant: '<S110>/Constant'
         */
        MARS_B.Merge_j2 = 1.0;
      }

      /* End of Outputs for SubSystem: '<S108>/If Action Subsystem1' */
      break;

     case 2:
      /* Outputs for IfAction SubSystem: '<S108>/If Action Subsystem2' incorporates:
       *  ActionPort: '<S111>/Action Port'
       */
      MARS_IfActionSubsystem2(rtb_VectorConcatenate[2], &MARS_B.Merge_j2);

      /* End of Outputs for SubSystem: '<S108>/If Action Subsystem2' */
      break;
    }

    /* End of If: '<S108>/If' */
    /* End of Outputs for SubSystem: '<S73>/AxisRotZeroR3' */
    break;

   case 1:
    /* Outputs for IfAction SubSystem: '<S73>/AxisRotDefault' incorporates:
     *  ActionPort: '<S100>/Action Port'
     */
    /* If: '<S104>/If' */
    if (rtsiIsModeUpdateTimeStep(&MARS_M->solverInfo)) {
      if (rtb_VectorConcatenate[2] > 1.0) {
        rtAction = 0;
      } else if (rtb_VectorConcatenate[2] < -1.0) {
        rtAction = 1;
      } else {
        rtAction = 2;
      }

      MARS_DW.If_ActiveSubsystem_c = rtAction;
    } else {
      rtAction = MARS_DW.If_ActiveSubsystem_c;
    }

    switch (rtAction) {
     case 0:
      /* Outputs for IfAction SubSystem: '<S104>/If Action Subsystem' incorporates:
       *  ActionPort: '<S105>/Action Port'
       */
      if (rtmIsMajorTimeStep(MARS_M) &&
          MARS_M->Timing.TaskCounters.TID[1] == 0) {
        /* Merge: '<S104>/Merge' incorporates:
         *  Constant: '<S105>/Constant'
         */
        MARS_B.Merge_j = 1.0;
      }

      /* End of Outputs for SubSystem: '<S104>/If Action Subsystem' */
      break;

     case 1:
      /* Outputs for IfAction SubSystem: '<S104>/If Action Subsystem1' incorporates:
       *  ActionPort: '<S106>/Action Port'
       */
      if (rtmIsMajorTimeStep(MARS_M) &&
          MARS_M->Timing.TaskCounters.TID[1] == 0) {
        /* Merge: '<S104>/Merge' incorporates:
         *  Constant: '<S106>/Constant'
         */
        MARS_B.Merge_j = 1.0;
      }

      /* End of Outputs for SubSystem: '<S104>/If Action Subsystem1' */
      break;

     case 2:
      /* Outputs for IfAction SubSystem: '<S104>/If Action Subsystem2' incorporates:
       *  ActionPort: '<S107>/Action Port'
       */
      MARS_IfActionSubsystem2(rtb_VectorConcatenate[2], &MARS_B.Merge_j);

      /* End of Outputs for SubSystem: '<S104>/If Action Subsystem2' */
      break;
    }

    /* End of If: '<S104>/If' */
    /* End of Outputs for SubSystem: '<S73>/AxisRotDefault' */
    break;
  }

  /* End of If: '<S73>/If' */
  if (rtmIsMajorTimeStep(MARS_M) &&
      MARS_M->Timing.TaskCounters.TID[2] == 0) {
    /* Gain: '<S5>/Output' incorporates:
     *  RandomNumber: '<S5>/White Noise'
     */
    MARS_B.Output = 2.23606797749979 * MARS_DW.NextOutput;
  }

  /* Sum: '<S2>/Subtract' incorporates:
   *  Constant: '<S2>/Ref Altitude'
   *  MultiPortSwitch: '<S1>/Height'
   *  Sum: '<S1>/Sum'
   */
  rtb_ixk = 100.0 - (MARS_B.Output + rtb_ECEFPositiontoLLA_o2);

  /* Gain: '<S35>/Integral Gain' */
  MARS_B.IntegralGain = 0.002 * rtb_ixk;

  /* Gain: '<S41>/Filter Coefficient' incorporates:
   *  Gain: '<S32>/Derivative Gain'
   *  Integrator: '<S33>/Filter'
   *  Sum: '<S33>/SumD'
   */
  MARS_B.FilterCoefficient = (0.064118 * rtb_ixk - MARS_X.Filter_CSTATE) * 360.0;

  /* Sum: '<S47>/Sum' incorporates:
   *  Gain: '<S43>/Proportional Gain'
   *  Integrator: '<S38>/Integrator'
   */
  rtb_jxi = (0.02 * rtb_ixk + MARS_X.Integrator_CSTATE) +
    MARS_B.FilterCoefficient;
  if (rtmIsMajorTimeStep(MARS_M) &&
      MARS_M->Timing.TaskCounters.TID[1] == 0) {
    boolean_T tmp;

    /* If: '<S103>/If1' */
    rtAction = -1;
    if (rtsiIsModeUpdateTimeStep(&MARS_M->solverInfo)) {
      MARS_DW.If1_ActiveSubsystem = -1;
    } else {
      rtAction = MARS_DW.If1_ActiveSubsystem;
    }

    if (rtAction == 0) {
      /* Outputs for IfAction SubSystem: '<S103>/If Warning//Error' incorporates:
       *  ActionPort: '<S112>/if'
       */
      /* Bias: '<S115>/Bias1' incorporates:
       *  Math: '<S115>/Math Function'
       *  Product: '<S115>/Product'
       *  Product: '<S61>/Product4'
       */
      for (i = 0; i < 3; i++) {
        for (i_0 = 0; i_0 < 3; i_0++) {
          Product4_tmp = 3 * i + i_0;
          rtb_Product4_g_tmp[Product4_tmp] = ((MARS_B.Product4[3 * i_0 + 1] *
            MARS_B.Product4[3 * i + 1] + MARS_B.Product4[3 * i_0] *
            MARS_B.Product4[3 * i]) + MARS_B.Product4[3 * i_0 + 2] *
            MARS_B.Product4[3 * i + 2]) + MARS_ConstP.pooled5[Product4_tmp];
        }
      }

      /* End of Bias: '<S115>/Bias1' */

      /* RelationalOperator: '<S121>/Compare' incorporates:
       *  Abs: '<S115>/Abs2'
       *  Constant: '<S121>/Constant'
       */
      for (i = 0; i < 9; i++) {
        rtb_Compare_a[i] = (fabs(rtb_Product4_g_tmp[i]) > 4.4408920985006262E-16);
      }

      /* End of RelationalOperator: '<S121>/Compare' */

      /* Logic: '<S115>/Logical Operator1' incorporates:
       *  RelationalOperator: '<S121>/Compare'
       */
      tmp = rtb_Compare_a[0];
      for (i = 0; i < 8; i++) {
        tmp = (tmp || rtb_Compare_a[i + 1]);
      }

      /* If: '<S112>/If' incorporates:
       *  Abs: '<S116>/Abs1'
       *  Bias: '<S116>/Bias'
       *  Constant: '<S123>/Constant'
       *  Logic: '<S115>/Logical Operator1'
       *  Product: '<S122>/Product'
       *  Product: '<S122>/Product1'
       *  Product: '<S122>/Product2'
       *  Product: '<S122>/Product3'
       *  Product: '<S122>/Product4'
       *  Product: '<S122>/Product5'
       *  Product: '<S61>/Product4'
       *  RelationalOperator: '<S123>/Compare'
       *  Reshape: '<S122>/Reshape'
       *  Sum: '<S122>/Sum'
       */
      if (fabs((((((MARS_B.Product4[0] * MARS_B.Product4[4] * MARS_B.Product4[8]
                    - MARS_B.Product4[0] * MARS_B.Product4[5] * MARS_B.Product4
                    [7]) - MARS_B.Product4[1] * MARS_B.Product4[3] *
                   MARS_B.Product4[8]) + MARS_B.Product4[2] * MARS_B.Product4[3]
                  * MARS_B.Product4[7]) + MARS_B.Product4[1] * MARS_B.Product4[5]
                 * MARS_B.Product4[6]) - MARS_B.Product4[2] * MARS_B.Product4[4]
                * MARS_B.Product4[6]) + -1.0) > 4.4408920985006262E-16) {
        /* Outputs for IfAction SubSystem: '<S112>/If Not Proper' incorporates:
         *  ActionPort: '<S114>/Action Port'
         */
        MARS_IfNotProper(1.0);

        /* End of Outputs for SubSystem: '<S112>/If Not Proper' */
      } else if (tmp) {
        /* Outputs for IfAction SubSystem: '<S112>/Else If Not Orthogonal' incorporates:
         *  ActionPort: '<S113>/Action Port'
         */
        MARS_ElseIfNotOrthogonal(1.0);

        /* End of Outputs for SubSystem: '<S112>/Else If Not Orthogonal' */
      }

      /* End of If: '<S112>/If' */
      /* End of Outputs for SubSystem: '<S103>/If Warning//Error' */
    }

    /* End of If: '<S103>/If1' */

    /* If: '<S126>/If1' */
    rtAction = -1;
    if (rtsiIsModeUpdateTimeStep(&MARS_M->solverInfo)) {
      MARS_DW.If1_ActiveSubsystem_l = -1;
    } else {
      rtAction = MARS_DW.If1_ActiveSubsystem_l;
    }

    if (rtAction == 0) {
      /* Outputs for IfAction SubSystem: '<S126>/If Warning//Error' incorporates:
       *  ActionPort: '<S150>/if'
       */
      /* RelationalOperator: '<S159>/Compare' incorporates:
       *  Abs: '<S153>/Abs2'
       *  Constant: '<S159>/Constant'
       */
      for (i = 0; i < 9; i++) {
        rtb_Compare_a[i] = (MARS_ConstB.Abs2[i] > 4.4408920985006262E-16);
      }

      /* End of RelationalOperator: '<S159>/Compare' */

      /* Logic: '<S153>/Logical Operator1' incorporates:
       *  RelationalOperator: '<S159>/Compare'
       */
      tmp = rtb_Compare_a[0];
      for (i = 0; i < 8; i++) {
        tmp = (tmp || rtb_Compare_a[i + 1]);
      }

      /* If: '<S150>/If' incorporates:
       *  Logic: '<S153>/Logical Operator1'
       */
      if (tmp) {
        /* Outputs for IfAction SubSystem: '<S150>/Else If Not Orthogonal' incorporates:
         *  ActionPort: '<S151>/Action Port'
         */
        MARS_ElseIfNotOrthogonal(1.0);

        /* End of Outputs for SubSystem: '<S150>/Else If Not Orthogonal' */
      }

      /* End of If: '<S150>/If' */
      /* End of Outputs for SubSystem: '<S126>/If Warning//Error' */
    }

    /* End of If: '<S126>/If1' */
  }

  /* Trigonometry: '<S67>/sincos' incorporates:
   *  Integrator: '<S65>/Integrator'
   */
  rtb_ixk = sin(MARS_X.Integrator_CSTATE_f);
  rtb_kxj = cos(MARS_X.Integrator_CSTATE_f);

  /* SignalConversion generated from: '<S196>/Vector Concatenate' */
  rtb_VectorConcatenate_b[0] = rtb_kxj;

  /* SignalConversion generated from: '<S196>/Vector Concatenate' */
  rtb_VectorConcatenate_b[1] = rtb_ixk;

  /* SignalConversion generated from: '<S196>/Vector Concatenate' incorporates:
   *  Constant: '<S67>/Zero'
   */
  rtb_VectorConcatenate_b[2] = 0.0;

  /* UnaryMinus: '<S67>/Unary Minus' */
  rtb_VectorConcatenate_b[3] = -rtb_ixk;

  /* SignalConversion generated from: '<S196>/Vector Concatenate' */
  rtb_VectorConcatenate_b[4] = rtb_kxj;

  /* SignalConversion generated from: '<S196>/Vector Concatenate' incorporates:
   *  Constant: '<S67>/Zero'
   */
  rtb_VectorConcatenate_b[5] = 0.0;

  /* SignalConversion generated from: '<S196>/Vector Concatenate' incorporates:
   *  Constant: '<S67>/Zero'
   */
  rtb_VectorConcatenate_b[6] = 0.0;

  /* SignalConversion generated from: '<S196>/Vector Concatenate' incorporates:
   *  Constant: '<S67>/Zero'
   */
  rtb_VectorConcatenate_b[7] = 0.0;

  /* SignalConversion generated from: '<S196>/Vector Concatenate' incorporates:
   *  Constant: '<S67>/Zero1'
   */
  rtb_VectorConcatenate_b[8] = 1.0;
  for (i = 0; i < 3; i++) {
    /* Product: '<S57>/Product3' incorporates:
     *  Concatenate: '<S171>/Vector Concatenate'
     *  Product: '<S63>/Product2'
     */
    rtb_Sum_hp[i] = (rtb_VectorConcatenate_f[i + 3] * 0.0 +
                     rtb_VectorConcatenate_f[i] * 0.0) +
      rtb_VectorConcatenate_f[i + 6] * 7.292115E-5;

    /* Product: '<S70>/Product3' incorporates:
     *  Integrator: '<S63>/ub,vb,wb'
     *  Math: '<S70>/Math Function2'
     *  Product: '<S61>/Product4'
     */
    rtb_Sum2[i] = (MARS_B.Product4[3 * i + 1] * MARS_X.ubvbwb_CSTATE[1] +
                   MARS_B.Product4[3 * i] * MARS_X.ubvbwb_CSTATE[0]) +
      MARS_B.Product4[3 * i + 2] * MARS_X.ubvbwb_CSTATE[2];
  }

  /* Sum: '<S70>/Sum2' incorporates:
   *  Constant: '<S70>/f2'
   *  Product: '<S70>/Product'
   */
  /* Unit Conversion - from: deg to: rad
     Expression: output = (0.0174533*input) + (0) */
  rtb_kxj = 1.0 - rtb_VectorConcatenate_f_tmp_1 * rtb_VectorConcatenate_f_tmp_1 *
    0.00669437999014133;

  /* Product: '<S70>/w1' incorporates:
   *  Constant: '<S199>/f4'
   *  Product: '<S199>/N'
   *  Sqrt: '<S199>/sqrt'
   *  Sum: '<S199>/Sum4'
   */
  rtb_ixk = rtb_Sum2[1] / (6.378137E+6 / sqrt(rtb_kxj) +
    rtb_ECEFPositiontoLLA_o2);

  /* SignalConversion generated from: '<S57>/Product2' incorporates:
   *  Constant: '<S198>/f1'
   *  Constant: '<S198>/f3'
   *  Gain: '<S70>/Gain'
   *  Gain: '<S70>/Gain1'
   *  Math: '<S198>/Math Function'
   *  Product: '<S198>/M'
   *  Product: '<S70>/w2'
   *  Product: '<S70>/w3'
   *  Sum: '<S198>/Sum1'
   *  Trigonometry: '<S70>/tan'
   */
  rtb_ECEFPositiontoLLA_o2 = -(rtb_Sum2[0] / (6.3354393272928195E+6 /
    rt_powd_snf(rtb_kxj, 1.5) + rtb_ECEFPositiontoLLA_o2));
  rtb_kxj = -(rtb_ixk * tan(rtb_ECEFPositiontoLLA_o1_idx_0));
  for (i = 0; i <= 0; i += 2) {
    /* Sum: '<S57>/Sum2' incorporates:
     *  Product: '<S57>/Product2'
     *  Product: '<S61>/Product4'
     */
    tmp_3 = _mm_loadu_pd(&MARS_B.Product4[i + 3]);
    tmp_4 = _mm_loadu_pd(&MARS_B.Product4[i]);
    tmp_1 = _mm_loadu_pd(&MARS_B.Product4[i + 6]);

    /* Product: '<S57>/Product3' incorporates:
     *  Product: '<S57>/Product2'
     *  Sum: '<S57>/Sum2'
     */
    tmp_2 = _mm_loadu_pd(&rtb_Sum_hp[i]);

    /* Sum: '<S57>/Sum2' incorporates:
     *  Product: '<S57>/Product2'
     *  SignalConversion generated from: '<S57>/Product2'
     */
    _mm_storeu_pd(&MARS_B.Sum2[i], _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd
      (tmp_3, _mm_set1_pd(rtb_ECEFPositiontoLLA_o2)), _mm_mul_pd(tmp_4,
      _mm_set1_pd(rtb_ixk))), _mm_mul_pd(tmp_1, _mm_set1_pd(rtb_kxj))), tmp_2));
  }

  for (i = 2; i < 3; i++) {
    /* Sum: '<S57>/Sum2' incorporates:
     *  Product: '<S57>/Product2'
     *  Product: '<S57>/Product3'
     *  Product: '<S61>/Product4'
     *  SignalConversion generated from: '<S57>/Product2'
     *  Sum: '<S57>/Sum3'
     */
    MARS_B.Sum2[i] = ((MARS_B.Product4[i + 3] * rtb_ECEFPositiontoLLA_o2 +
                       MARS_B.Product4[i] * rtb_ixk) + MARS_B.Product4[i + 6] *
                      rtb_kxj) + rtb_Sum_hp[i];
  }

  /* Integrator: '<S57>/p,q,r ' */
  if (MARS_DW.pqr_IWORK != 0) {
    MARS_X.pqr_CSTATE[0] = MARS_B.Sum2[0];
    MARS_X.pqr_CSTATE[1] = MARS_B.Sum2[1];
    MARS_X.pqr_CSTATE[2] = MARS_B.Sum2[2];
  }

  /* Sum: '<S57>/Sum4' incorporates:
   *  Integrator: '<S57>/p,q,r '
   *  Product: '<S57>/Product3'
   */
  rtb_Sum2[0] = MARS_X.pqr_CSTATE[0] - rtb_Sum_hp[0];
  rtb_Sum2[1] = MARS_X.pqr_CSTATE[1] - rtb_Sum_hp[1];
  rtb_Sum2[2] = MARS_X.pqr_CSTATE[2] - rtb_Sum_hp[2];

  /* Product: '<S176>/Product' incorporates:
   *  Integrator: '<S61>/q'
   */
  rtb_ECEFPositiontoLLA_o2 = MARS_X.q_CSTATE[0] / rtb_Gain_f;

  /* Product: '<S176>/Product1' incorporates:
   *  Integrator: '<S61>/q'
   */
  rtb_kxj = MARS_X.q_CSTATE[1] / rtb_Gain_f;

  /* SignalConversion generated from: '<S61>/q' incorporates:
   *  Fcn: '<S77>/q0dot'
   *  Fcn: '<S77>/q1dot'
   *  Fcn: '<S77>/q2dot'
   *  Fcn: '<S77>/q3dot'
   */
  MARS_B.TmpSignalConversionAtqInport1[0] = ((rtb_kxj * rtb_Sum2[0] + uTmp_idx_0
    * rtb_Sum2[1]) + uTmp_idx_1 * rtb_Sum2[2]) * -0.5;
  MARS_B.TmpSignalConversionAtqInport1[1] = ((rtb_ECEFPositiontoLLA_o2 *
    rtb_Sum2[0] + uTmp_idx_0 * rtb_Sum2[2]) - uTmp_idx_1 * rtb_Sum2[1]) * 0.5;
  MARS_B.TmpSignalConversionAtqInport1[2] = ((rtb_ECEFPositiontoLLA_o2 *
    rtb_Sum2[1] + uTmp_idx_1 * rtb_Sum2[0]) - rtb_kxj * rtb_Sum2[2]) * 0.5;
  MARS_B.TmpSignalConversionAtqInport1[3] = ((rtb_ECEFPositiontoLLA_o2 *
    rtb_Sum2[2] + rtb_kxj * rtb_Sum2[1]) - uTmp_idx_0 * rtb_Sum2[0]) * 0.5;
  for (i = 0; i < 3; i++) {
    /* Product: '<S62>/Product1' */
    MARS_B.Product1[i] = 0.0;

    /* Product: '<S62>/Product1' incorporates:
     *  Concatenate: '<S196>/Vector Concatenate'
     */
    rtb_Gain_f = rtb_VectorConcatenate_b[i];

    /* Product: '<S62>/Product1' incorporates:
     *  LLA2ECEF: '<S62>/LLA to ECEF Position'
     */
    MARS_B.Product1[i] += rtb_Gain_f * 2.5914814724705652E+6;

    /* Math: '<S62>/Math Function1' incorporates:
     *  Math: '<S61>/Math Function'
     */
    rtb_Product4_g_tmp[3 * i] = rtb_Gain_f;

    /* Product: '<S62>/Product1' incorporates:
     *  Concatenate: '<S196>/Vector Concatenate'
     */
    rtb_Gain_f = rtb_VectorConcatenate_b[i + 3];

    /* Product: '<S62>/Product1' incorporates:
     *  LLA2ECEF: '<S62>/LLA to ECEF Position'
     */
    MARS_B.Product1[i] += rtb_Gain_f * -1.0790161949738523E+6;

    /* Math: '<S62>/Math Function1' incorporates:
     *  Math: '<S61>/Math Function'
     */
    rtb_Product4_g_tmp[3 * i + 1] = rtb_Gain_f;

    /* Product: '<S62>/Product1' incorporates:
     *  Concatenate: '<S196>/Vector Concatenate'
     */
    rtb_Gain_f = rtb_VectorConcatenate_b[i + 6];

    /* Product: '<S62>/Product1' incorporates:
     *  LLA2ECEF: '<S62>/LLA to ECEF Position'
     */
    MARS_B.Product1[i] += rtb_Gain_f * 5.7079879931331929E+6;

    /* Math: '<S62>/Math Function1' incorporates:
     *  Math: '<S61>/Math Function'
     */
    rtb_Product4_g_tmp[3 * i + 2] = rtb_Gain_f;
  }

  for (i = 0; i < 3; i++) {
    /* Product: '<S62>/Product5' incorporates:
     *  Concatenate: '<S171>/Vector Concatenate'
     *  Integrator: '<S63>/ub,vb,wb'
     *  Math: '<S62>/Math Function2'
     */
    MARS_B.Product5[i] = 0.0;
    MARS_B.Product5[i] += rtb_VectorConcatenate_f[3 * i] * MARS_X.ubvbwb_CSTATE
      [0];
    MARS_B.Product5[i] += rtb_VectorConcatenate_f[3 * i + 1] *
      MARS_X.ubvbwb_CSTATE[1];
    MARS_B.Product5[i] += rtb_VectorConcatenate_f[3 * i + 2] *
      MARS_X.ubvbwb_CSTATE[2];

    /* Product: '<S62>/Product4' incorporates:
     *  Math: '<S62>/Math Function1'
     */
    rtb_Sum2[i] = (rtb_Product4_g_tmp[i + 3] * 0.0 + rtb_Product4_g_tmp[i] * 0.0)
      + rtb_Product4_g_tmp[i + 6] * 7.292115E-5;
  }

  /* Integrator: '<S62>/p' */
  if (MARS_DW.p_IWORK != 0) {
    MARS_X.p_CSTATE[0] = MARS_B.Product1[0];
    MARS_X.p_CSTATE[1] = MARS_B.Product1[1];
    MARS_X.p_CSTATE[2] = MARS_B.Product1[2];
  }

  /* Sum: '<S179>/Sum' incorporates:
   *  Integrator: '<S62>/p'
   *  Product: '<S180>/i x j'
   *  Product: '<S180>/j x k'
   *  Product: '<S180>/k x i'
   *  Product: '<S181>/i x k'
   *  Product: '<S181>/j x i'
   *  Product: '<S181>/k x j'
   */
  rtb_Gain_a[0] = MARS_X.p_CSTATE[1] * rtb_Sum2[2];
  rtb_Gain_a[1] = rtb_Sum2[0] * MARS_X.p_CSTATE[2];
  rtb_Gain_a[2] = MARS_X.p_CSTATE[0] * rtb_Sum2[1];
  rtb_Sum2_0[0] = rtb_Sum2[1] * MARS_X.p_CSTATE[2];
  rtb_Sum2_0[1] = MARS_X.p_CSTATE[0] * rtb_Sum2[2];
  rtb_Sum2_0[2] = rtb_Sum2[0] * MARS_X.p_CSTATE[1];
  for (i = 0; i < 3; i++) {
    /* Sum: '<S62>/Sum2' */
    rtb_Gain_f = 0.0;
    for (i_0 = 0; i_0 < 3; i_0++) {
      /* Math: '<S62>/Math Function' incorporates:
       *  Concatenate: '<S171>/Vector Concatenate'
       *  Math: '<S61>/Math Function'
       *  Product: '<S61>/Product1'
       *  Product: '<S62>/Product2'
       */
      Product4_tmp = 3 * i_0 + i;
      rtb_VectorConcatenate_b[Product4_tmp] = 0.0;
      rtb_VectorConcatenate_b[Product4_tmp] += rtb_Product4_g_tmp[3 * i] *
        rtb_VectorConcatenate_f[i_0];
      rtb_VectorConcatenate_b[Product4_tmp] += rtb_Product4_g_tmp[3 * i + 1] *
        rtb_VectorConcatenate_f[i_0 + 3];
      rtb_VectorConcatenate_b[Product4_tmp] += rtb_Product4_g_tmp[3 * i + 2] *
        rtb_VectorConcatenate_f[i_0 + 6];

      /* Sum: '<S62>/Sum2' incorporates:
       *  Integrator: '<S63>/ub,vb,wb'
       *  Product: '<S62>/Product2'
       */
      rtb_Gain_f += rtb_VectorConcatenate_b[Product4_tmp] *
        MARS_X.ubvbwb_CSTATE[i_0];
    }

    /* Sum: '<S62>/Sum2' incorporates:
     *  Product: '<S62>/Product2'
     *  Sum: '<S179>/Sum'
     */
    MARS_B.Sum2_e[i] = rtb_Gain_f - (rtb_Gain_a[i] - rtb_Sum2_0[i]);

    /* Sum: '<S63>/Sum2' incorporates:
     *  Integrator: '<S57>/p,q,r '
     *  Product: '<S63>/Product2'
     */
    rtb_Sum_hp[i] += MARS_X.pqr_CSTATE[i];
  }

  /* Sum: '<S184>/Sum' incorporates:
   *  Constant: '<S57>/omega_earth2'
   *  Constant: '<S57>/omega_earth3'
   *  Integrator: '<S62>/p1'
   *  Product: '<S189>/i x j'
   *  Product: '<S189>/j x k'
   *  Product: '<S189>/k x i'
   *  Product: '<S190>/i x k'
   *  Product: '<S190>/j x i'
   *  Product: '<S190>/k x j'
   */
  rtb_Sum2[0] = 0.0 * MARS_X.p1_CSTATE[2] - 7.292115E-5 * MARS_X.p1_CSTATE[1];
  rtb_Sum2[1] = 7.292115E-5 * MARS_X.p1_CSTATE[0] - 0.0 * MARS_X.p1_CSTATE[2];
  rtb_Sum2[2] = 0.0 * MARS_X.p1_CSTATE[1] - 0.0 * MARS_X.p1_CSTATE[0];

  /* Saturate: '<S45>/Saturation' */
  if (rtb_jxi > 1.0) {
    rtb_jxi = 1.0;
  } else if (rtb_jxi < 0.2) {
    rtb_jxi = 0.2;
  }

  /* End of Saturate: '<S45>/Saturation' */

  /* Product: '<S59>/Product2' incorporates:
   *  Constant: '<S4>/Max Thrust'
   */
  rtb_jxi *= 40.0;

  /* Product: '<S202>/Tsin(delta)' */
  rtb_ECEFPositiontoLLA_o2 = rtb_jxi * 0.0;

  /* Product: '<S202>/Product' */
  rtb_Gain_f = rtb_ECEFPositiontoLLA_o2;

  /* Product: '<S202>/Product1' */
  rtb_ECEFPositiontoLLA_o2 *= 0.0;

  /* Sum: '<S182>/Sum' incorporates:
   *  Integrator: '<S63>/ub,vb,wb'
   *  Product: '<S185>/i x j'
   *  Product: '<S185>/j x k'
   *  Product: '<S185>/k x i'
   *  Product: '<S186>/i x k'
   *  Product: '<S186>/j x i'
   *  Product: '<S186>/k x j'
   */
  rtb_Gain_a[0] = MARS_X.ubvbwb_CSTATE[1] * rtb_Sum_hp[2];
  rtb_Gain_a[1] = rtb_Sum_hp[0] * MARS_X.ubvbwb_CSTATE[2];
  rtb_Gain_a[2] = MARS_X.ubvbwb_CSTATE[0] * rtb_Sum_hp[1];
  rtb_Sum2_0[0] = rtb_Sum_hp[1] * MARS_X.ubvbwb_CSTATE[2];
  rtb_Sum2_0[1] = MARS_X.ubvbwb_CSTATE[0] * rtb_Sum_hp[2];
  rtb_Sum2_0[2] = rtb_Sum_hp[0] * MARS_X.ubvbwb_CSTATE[1];

  /* Sum: '<S183>/Sum' incorporates:
   *  Constant: '<S57>/omega_earth2'
   *  Constant: '<S57>/omega_earth3'
   *  Product: '<S187>/i x j'
   *  Product: '<S187>/j x k'
   *  Product: '<S187>/k x i'
   *  Product: '<S188>/i x k'
   *  Product: '<S188>/j x i'
   *  Product: '<S188>/k x j'
   */
  uTmp_idx_0 = 0.0 * rtb_Sum2[2] - 7.292115E-5 * rtb_Sum2[1];
  uTmp_idx_1 = 7.292115E-5 * rtb_Sum2[0] - 0.0 * rtb_Sum2[2];
  rtb_ixk = 0.0 * rtb_Sum2[1] - 0.0 * rtb_Sum2[0];
  for (i = 0; i <= 0; i += 2) {
    __m128d tmp_0;

    /* Sum: '<S182>/Sum' incorporates:
     *  Product: '<S63>/Product1'
     *  Sum: '<S63>/Sum'
     */
    tmp_3 = _mm_loadu_pd(&rtb_Gain_a[i]);
    tmp_4 = _mm_loadu_pd(&rtb_Sum2_0[i]);

    /* Sum: '<S63>/Sum' incorporates:
     *  Concatenate: '<S171>/Vector Concatenate'
     *  Product: '<S63>/Product1'
     */
    tmp_1 = _mm_loadu_pd(&rtb_VectorConcatenate_f[i + 3]);
    tmp_2 = _mm_loadu_pd(&rtb_VectorConcatenate_f[i]);
    tmp_0 = _mm_loadu_pd(&rtb_VectorConcatenate_f[i + 6]);
    _mm_storeu_pd(&rtb_Sum2[i], _mm_sub_pd(_mm_sub_pd(tmp_3, tmp_4), _mm_add_pd
      (_mm_add_pd(_mm_mul_pd(tmp_1, _mm_set1_pd(uTmp_idx_1)), _mm_mul_pd(tmp_2,
      _mm_set1_pd(uTmp_idx_0))), _mm_mul_pd(tmp_0, _mm_set1_pd(rtb_ixk)))));
  }

  /* Sum: '<S63>/Sum' incorporates:
   *  Concatenate: '<S171>/Vector Concatenate'
   *  Product: '<S63>/Product1'
   *  Sum: '<S182>/Sum'
   */
  for (i = 2; i < 3; i++) {
    rtb_Sum2[i] = (rtb_Gain_a[i] - rtb_Sum2_0[i]) - ((rtb_VectorConcatenate_f[i
      + 3] * uTmp_idx_1 + rtb_VectorConcatenate_f[i] * uTmp_idx_0) +
      rtb_VectorConcatenate_f[i + 6] * rtb_ixk);
  }

  /* Product: '<S57>/Product' incorporates:
   *  Constant: '<S66>/Constant'
   *  Gain: '<S202>/Gain2'
   *  Gain: '<S202>/Gain4'
   *  Product: '<S202>/Tcos(delta)'
   *  Sum: '<S202>/Sum'
   */
  rtb_Gain_a[0] = rtb_Gain_f / 1.6;
  rtb_Gain_a[1] = -rtb_ECEFPositiontoLLA_o2 / 1.6;
  rtb_Gain_a[2] = -(rtb_jxi + -18.605) / 1.6;
  for (i = 0; i < 3; i++) {
    /* Sum: '<S63>/Sum' */
    rtb_jxi = rtb_Sum2[i] + rtb_Gain_a[i];

    /* DeadZone: '<S63>/Dead Zone' */
    if (rtb_jxi > 2.2204460492503131E-16) {
      /* DeadZone: '<S63>/Dead Zone' */
      MARS_B.DeadZone[i] = rtb_jxi - 2.2204460492503131E-16;
    } else if (rtb_jxi >= -2.2204460492503131E-16) {
      /* DeadZone: '<S63>/Dead Zone' */
      MARS_B.DeadZone[i] = 0.0;
    } else {
      /* DeadZone: '<S63>/Dead Zone' */
      MARS_B.DeadZone[i] = rtb_jxi - -2.2204460492503131E-16;
    }

    /* End of DeadZone: '<S63>/Dead Zone' */

    /* Sum: '<S63>/Sum' incorporates:
     *  Integrator: '<S57>/p,q,r '
     *  Product: '<S192>/Product'
     *  Selector: '<S64>/Selector'
     */
    rtb_Sum_hp[i] = (MARS_ConstB.Selector[i + 3] * MARS_X.pqr_CSTATE[1] +
                     MARS_ConstB.Selector[i] * MARS_X.pqr_CSTATE[0]) +
      MARS_ConstB.Selector[i + 6] * MARS_X.pqr_CSTATE[2];
  }

  /* Sum: '<S191>/Sum' incorporates:
   *  Integrator: '<S57>/p,q,r '
   *  Product: '<S194>/i x j'
   *  Product: '<S194>/j x k'
   *  Product: '<S194>/k x i'
   *  Product: '<S195>/i x k'
   *  Product: '<S195>/j x i'
   *  Product: '<S195>/k x j'
   */
  rtb_Gain_a[0] = MARS_X.pqr_CSTATE[1] * rtb_Sum_hp[2];
  rtb_Gain_a[1] = rtb_Sum_hp[0] * MARS_X.pqr_CSTATE[2];
  rtb_Gain_a[2] = MARS_X.pqr_CSTATE[0] * rtb_Sum_hp[1];
  rtb_Sum2_0[0] = rtb_Sum_hp[1] * MARS_X.pqr_CSTATE[2];
  rtb_Sum2_0[1] = MARS_X.pqr_CSTATE[0] * rtb_Sum_hp[2];
  rtb_Sum2_0[2] = rtb_Sum_hp[0] * MARS_X.pqr_CSTATE[1];
  for (i = 0; i <= 0; i += 2) {
    /* Sum: '<S191>/Sum' */
    tmp_3 = _mm_loadu_pd(&rtb_Gain_a[i]);
    tmp_4 = _mm_loadu_pd(&rtb_Sum2_0[i]);
    _mm_storeu_pd(&rtb_Sum_hp[i], _mm_sub_pd(tmp_3, tmp_4));

    /* Product: '<S193>/Product' incorporates:
     *  Integrator: '<S57>/p,q,r '
     *  Sum: '<S191>/Sum'
     *  Sum: '<S64>/Sum2'
     */
    _mm_storeu_pd(&rtb_Sum2[i], _mm_set1_pd(0.0 * MARS_X.pqr_CSTATE[2] + (0.0 *
      MARS_X.pqr_CSTATE[1] + 0.0 * MARS_X.pqr_CSTATE[0])));
  }

  for (i = 2; i < 3; i++) {
    /* Sum: '<S191>/Sum' */
    rtb_Sum_hp[i] = rtb_Gain_a[i] - rtb_Sum2_0[i];

    /* Product: '<S193>/Product' incorporates:
     *  Integrator: '<S57>/p,q,r '
     *  Sum: '<S64>/Sum2'
     */
    rtb_Sum2[i] = (0.0 * MARS_X.pqr_CSTATE[1] + 0.0 * MARS_X.pqr_CSTATE[0]) +
      0.0 * MARS_X.pqr_CSTATE[2];
  }

  /* Sum: '<S64>/Sum2' incorporates:
   *  Constant: '<S4>/M_x'
   *  Gain: '<S4>/Gain'
   *  Gain: '<S4>/Gain1'
   */
  rtb_Gain_a[0] = (0.1 * rtb_Gain_f - rtb_Sum2[0]) - rtb_Sum_hp[0];
  rtb_Gain_a[1] = (0.1 * rtb_ECEFPositiontoLLA_o2 - rtb_Sum2[1]) - rtb_Sum_hp[1];
  rtb_Gain_a[2] = (0.0 - rtb_Sum2[2]) - rtb_Sum_hp[2];

  /* Product: '<S64>/Product2' incorporates:
   *  Selector: '<S64>/Selector2'
   */
  rt_mrdivide_U1d1x3_U2d_9vOrDY9Z(rtb_Gain_a, MARS_ConstB.Selector2,
    MARS_B.Product2);
  if (rtmIsMajorTimeStep(MARS_M)) {
    /* Update for Integrator: '<S62>/p1' */
    MARS_DW.p1_IWORK = 0;

    /* Update for Integrator: '<S61>/q' */
    MARS_DW.q_IWORK = 0;
    if (rtmIsMajorTimeStep(MARS_M) &&
        MARS_M->Timing.TaskCounters.TID[2] == 0) {
      /* Update for RandomNumber: '<S5>/White Noise' */
      MARS_DW.NextOutput = rt_nrand_Upu32_Yd_f_pw_snf(&MARS_DW.RandSeed);
    }

    /* Update for Integrator: '<S57>/p,q,r ' */
    MARS_DW.pqr_IWORK = 0;

    /* Update for Integrator: '<S62>/p' */
    MARS_DW.p_IWORK = 0;
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep(MARS_M)) {
    rt_ertODEUpdateContinuousStates(&MARS_M->solverInfo);

    /* Update absolute time for base rate */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     */
    ++MARS_M->Timing.clockTick0;
    MARS_M->Timing.t[0] = rtsiGetSolverStopTime(&MARS_M->solverInfo);

    {
      /* Update absolute timer for sample time: [0.016666666666666666s, 0.0s] */
      /* The "clockTick1" counts the number of times the code of this task has
       * been executed. The resolution of this integer timer is 0.016666666666666666, which is the step size
       * of the task. Size of "clockTick1" ensures timer will not overflow during the
       * application lifespan selected.
       */
      MARS_M->Timing.clockTick1++;
    }

    rate_scheduler();
  }                                    /* end MajorTimeStep */
}

/* Derivatives for root system: '<Root>' */
void MARS_derivatives(void)
{
  XDot_MARS_T *_rtXdot;
  _rtXdot = ((XDot_MARS_T *) MARS_M->derivs);

  /* Derivatives for Integrator: '<S62>/p1' */
  _rtXdot->p1_CSTATE[0] = MARS_B.Product5[0];
  _rtXdot->p1_CSTATE[1] = MARS_B.Product5[1];
  _rtXdot->p1_CSTATE[2] = MARS_B.Product5[2];

  /* Derivatives for Integrator: '<S61>/q' */
  _rtXdot->q_CSTATE[0] = MARS_B.TmpSignalConversionAtqInport1[0];
  _rtXdot->q_CSTATE[1] = MARS_B.TmpSignalConversionAtqInport1[1];
  _rtXdot->q_CSTATE[2] = MARS_B.TmpSignalConversionAtqInport1[2];
  _rtXdot->q_CSTATE[3] = MARS_B.TmpSignalConversionAtqInport1[3];

  /* Derivatives for Integrator: '<S33>/Filter' */
  _rtXdot->Filter_CSTATE = MARS_B.FilterCoefficient;

  /* Derivatives for Integrator: '<S38>/Integrator' */
  _rtXdot->Integrator_CSTATE = MARS_B.IntegralGain;

  /* Derivatives for Integrator: '<S65>/Integrator' incorporates:
   *  Constant: '<S65>/Constant1'
   */
  _rtXdot->Integrator_CSTATE_f = 7.292115E-5;

  /* Derivatives for Integrator: '<S63>/ub,vb,wb' */
  _rtXdot->ubvbwb_CSTATE[0] = MARS_B.DeadZone[0];

  /* Derivatives for Integrator: '<S57>/p,q,r ' */
  _rtXdot->pqr_CSTATE[0] = MARS_B.Product2[0];

  /* Derivatives for Integrator: '<S62>/p' */
  _rtXdot->p_CSTATE[0] = MARS_B.Sum2_e[0];

  /* Derivatives for Integrator: '<S63>/ub,vb,wb' */
  _rtXdot->ubvbwb_CSTATE[1] = MARS_B.DeadZone[1];

  /* Derivatives for Integrator: '<S57>/p,q,r ' */
  _rtXdot->pqr_CSTATE[1] = MARS_B.Product2[1];

  /* Derivatives for Integrator: '<S62>/p' */
  _rtXdot->p_CSTATE[1] = MARS_B.Sum2_e[1];

  /* Derivatives for Integrator: '<S63>/ub,vb,wb' */
  _rtXdot->ubvbwb_CSTATE[2] = MARS_B.DeadZone[2];

  /* Derivatives for Integrator: '<S57>/p,q,r ' */
  _rtXdot->pqr_CSTATE[2] = MARS_B.Product2[2];

  /* Derivatives for Integrator: '<S62>/p' */
  _rtXdot->p_CSTATE[2] = MARS_B.Sum2_e[2];
}

/* Model initialize function */
void MARS_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* initialize real-time model */
  (void) memset((void *)MARS_M, 0,
                sizeof(RT_MODEL_MARS_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&MARS_M->solverInfo, &MARS_M->Timing.simTimeStep);
    rtsiSetTPtr(&MARS_M->solverInfo, &rtmGetTPtr(MARS_M));
    rtsiSetStepSizePtr(&MARS_M->solverInfo, &MARS_M->Timing.stepSize0);
    rtsiSetdXPtr(&MARS_M->solverInfo, &MARS_M->derivs);
    rtsiSetContStatesPtr(&MARS_M->solverInfo, (real_T **) &MARS_M->contStates);
    rtsiSetNumContStatesPtr(&MARS_M->solverInfo, &MARS_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&MARS_M->solverInfo,
      &MARS_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&MARS_M->solverInfo,
      &MARS_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&MARS_M->solverInfo,
      &MARS_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&MARS_M->solverInfo, (&rtmGetErrorStatus(MARS_M)));
    rtsiSetRTModelPtr(&MARS_M->solverInfo, MARS_M);
  }

  rtsiSetSimTimeStep(&MARS_M->solverInfo, MAJOR_TIME_STEP);
  MARS_M->intgData.y = MARS_M->odeY;
  MARS_M->intgData.f[0] = MARS_M->odeF[0];
  MARS_M->intgData.f[1] = MARS_M->odeF[1];
  MARS_M->intgData.f[2] = MARS_M->odeF[2];
  MARS_M->contStates = ((X_MARS_T *) &MARS_X);
  rtsiSetSolverData(&MARS_M->solverInfo, (void *)&MARS_M->intgData);
  rtsiSetIsMinorTimeStepWithModeChange(&MARS_M->solverInfo, false);
  rtsiSetSolverName(&MARS_M->solverInfo,"ode3");
  rtmSetTPtr(MARS_M, &MARS_M->Timing.tArray[0]);
  MARS_M->Timing.stepSize0 = 0.016666666666666666;
  rtmSetFirstInitCond(MARS_M, 1);

  /* block I/O */
  (void) memset(((void *) &MARS_B), 0,
                sizeof(B_MARS_T));

  /* states (continuous) */
  {
    (void) memset((void *)&MARS_X, 0,
                  sizeof(X_MARS_T));
  }

  /* states (dwork) */
  (void) memset((void *)&MARS_DW, 0,
                sizeof(DW_MARS_T));

  /* Start for If: '<S74>/If' */
  MARS_DW.If_ActiveSubsystem = -1;

  /* Start for If: '<S73>/If' */
  MARS_DW.If_ActiveSubsystem_h = -1;

  /* Start for If: '<S103>/If1' */
  MARS_DW.If1_ActiveSubsystem = -1;

  /* Start for If: '<S126>/If1' */
  MARS_DW.If1_ActiveSubsystem_l = -1;

  /* InitializeConditions for Integrator: '<S62>/p1' incorporates:
   *  Integrator: '<S61>/q'
   */
  if (rtmIsFirstInitCond(MARS_M)) {
    MARS_X.p1_CSTATE[0] = 0.0;
    MARS_X.p1_CSTATE[1] = 0.0;
    MARS_X.p1_CSTATE[2] = 0.0;
    MARS_X.q_CSTATE[0] = 0.0;
    MARS_X.q_CSTATE[1] = 0.0;
    MARS_X.q_CSTATE[2] = 0.0;
    MARS_X.q_CSTATE[3] = 0.0;
  }

  MARS_DW.p1_IWORK = 1;

  /* End of InitializeConditions for Integrator: '<S62>/p1' */

  /* InitializeConditions for Integrator: '<S61>/q' */
  MARS_DW.q_IWORK = 1;

  /* InitializeConditions for RandomNumber: '<S5>/White Noise' */
  MARS_DW.RandSeed = 1529675776U;
  MARS_DW.NextOutput = rt_nrand_Upu32_Yd_f_pw_snf(&MARS_DW.RandSeed);

  /* InitializeConditions for Integrator: '<S33>/Filter' */
  MARS_X.Filter_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S38>/Integrator' */
  MARS_X.Integrator_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S65>/Integrator' */
  MARS_X.Integrator_CSTATE_f = 0.0;

  /* InitializeConditions for Integrator: '<S63>/ub,vb,wb' */
  MARS_X.ubvbwb_CSTATE[0] = 0.0;
  MARS_X.ubvbwb_CSTATE[1] = 0.0;
  MARS_X.ubvbwb_CSTATE[2] = 0.0;

  /* InitializeConditions for Integrator: '<S57>/p,q,r ' incorporates:
   *  Integrator: '<S62>/p'
   */
  if (rtmIsFirstInitCond(MARS_M)) {
    MARS_X.pqr_CSTATE[0] = 0.0;
    MARS_X.pqr_CSTATE[1] = 0.0;
    MARS_X.pqr_CSTATE[2] = 0.0;
    MARS_X.p_CSTATE[0] = 0.0;
    MARS_X.p_CSTATE[1] = 0.0;
    MARS_X.p_CSTATE[2] = 0.0;
  }

  MARS_DW.pqr_IWORK = 1;

  /* End of InitializeConditions for Integrator: '<S57>/p,q,r ' */

  /* InitializeConditions for Integrator: '<S62>/p' */
  MARS_DW.p_IWORK = 1;

  /* SystemInitialize for Merge: '<S74>/Merge' */
  MARS_B.Merge[0] = 1.0;
  MARS_B.Merge[1] = 0.0;
  MARS_B.Merge[2] = 0.0;
  MARS_B.Merge[3] = 0.0;

  /* SystemInitialize for IfAction SubSystem: '<S73>/AxisRotZeroR3' */
  /* Start for If: '<S108>/If' */
  MARS_DW.If_ActiveSubsystem_c0 = -1;

  /* End of SystemInitialize for SubSystem: '<S73>/AxisRotZeroR3' */

  /* SystemInitialize for IfAction SubSystem: '<S73>/AxisRotDefault' */
  /* Start for If: '<S104>/If' */
  MARS_DW.If_ActiveSubsystem_c = -1;

  /* End of SystemInitialize for SubSystem: '<S73>/AxisRotDefault' */

  /* set "at time zero" to false */
  if (rtmIsFirstInitCond(MARS_M)) {
    rtmSetFirstInitCond(MARS_M, 0);
  }
}

/* Model terminate function */
void MARS_terminate(void)
{
  /* (no terminate code required) */
}
