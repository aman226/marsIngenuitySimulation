/*
 * MARS.cpp
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

#include "MARS.h"
#include "rtwtypes.h"
#include "rt_assert.h"
#include <cmath>
#include "MARS_private.h"
#include <emmintrin.h>
#include <cstring>
#include "rt_defines.h"
#include <cfloat>

extern "C" {

#include "rt_nonfinite.h"

}
/*
 * This function updates continuous states using the ODE3 fixed-step
 * solver algorithm
 */
  void MARS::rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  /* Solver Matrices */
  static const real_T rt_ODE3_A[3]{
    1.0/2.0, 3.0/4.0, 1.0
  };

  static const real_T rt_ODE3_B[3][3]{
    { 1.0/2.0, 0.0, 0.0 },

    { 0.0, 3.0/4.0, 0.0 },

    { 2.0/9.0, 1.0/3.0, 4.0/9.0 }
  };

  time_T t { rtsiGetT(si) };

  time_T tnew { rtsiGetSolverStopTime(si) };

  time_T h { rtsiGetStepSize(si) };

  real_T *x { rtsiGetContStates(si) };

  ODE3_IntgData *id { static_cast<ODE3_IntgData *>(rtsiGetSolverData(si)) };

  real_T *y { id->y };

  real_T *f0 { id->f[0] };

  real_T *f1 { id->f[1] };

  real_T *f2 { id->f[2] };

  real_T hB[3];
  int_T i;
  int_T nXc { 17 };

  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);

  /* Save the state values at time t in y, we'll use x as ynew. */
  (void) std::memcpy(y, x,
                     static_cast<uint_T>(nXc)*sizeof(real_T));

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
  this->step();
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
  this->step();
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
 *    '<S51>/If Action Subsystem2'
 *    '<S47>/If Action Subsystem2'
 */
void MARS::MARS_IfActionSubsystem2(real_T rtu_In, real_T *rty_OutOrig)
{
  /* SignalConversion generated from: '<S54>/In' */
  *rty_OutOrig = rtu_In;
}

/*
 * Output and update for action system:
 *    '<S55>/If Not Proper'
 *    '<S93>/If Not Proper'
 */
void MARS::MARS_IfNotProper(real_T rtp_action)
{
  /* If: '<S57>/If' incorporates:
   *  Constant: '<S57>/Constant'
   */
  if (rtp_action == 2.0) {
    /* Outputs for IfAction SubSystem: '<S57>/Warning' incorporates:
     *  ActionPort: '<S63>/Action Port'
     */
    /* Assertion: '<S63>/Assertion' incorporates:
     *  Constant: '<S57>/Constant1'
     */
    utAssert(false);

    /* End of Outputs for SubSystem: '<S57>/Warning' */
  } else if (rtp_action == 3.0) {
    /* Outputs for IfAction SubSystem: '<S57>/Error' incorporates:
     *  ActionPort: '<S62>/Action Port'
     */
    /* Assertion: '<S62>/Assertion' incorporates:
     *  Constant: '<S57>/Constant1'
     */
    utAssert(false);

    /* End of Outputs for SubSystem: '<S57>/Error' */
  }

  /* End of If: '<S57>/If' */
}

/*
 * Output and update for action system:
 *    '<S55>/Else If Not Orthogonal'
 *    '<S93>/Else If Not Orthogonal'
 */
void MARS::MARS_ElseIfNotOrthogonal(real_T rtp_action)
{
  /* If: '<S56>/If' incorporates:
   *  Constant: '<S56>/Constant'
   */
  if (rtp_action == 2.0) {
    /* Outputs for IfAction SubSystem: '<S56>/Warning' incorporates:
     *  ActionPort: '<S61>/Action Port'
     */
    /* Assertion: '<S61>/Assertion' incorporates:
     *  Constant: '<S56>/Constant1'
     */
    utAssert(false);

    /* End of Outputs for SubSystem: '<S56>/Warning' */
  } else if (rtp_action == 3.0) {
    /* Outputs for IfAction SubSystem: '<S56>/Error' incorporates:
     *  ActionPort: '<S60>/Action Port'
     */
    /* Assertion: '<S60>/Assertion' incorporates:
     *  Constant: '<S56>/Constant1'
     */
    utAssert(false);

    /* End of Outputs for SubSystem: '<S56>/Error' */
  }

  /* End of If: '<S56>/If' */
}

real_T rt_atan2d_snf(real_T u0, real_T u1)
{
  real_T y;
  if (std::isnan(u0) || std::isnan(u1)) {
    y = (rtNaN);
  } else if (std::isinf(u0) && std::isinf(u1)) {
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

    y = std::atan2(static_cast<real_T>(u0_0), static_cast<real_T>(u1_0));
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = std::atan2(u0, u1);
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
  } else if (std::isnan(u0) || std::isnan(u1) || std::isinf(u0)) {
    y = (rtNaN);
  } else if (u0 == 0.0) {
    y = 0.0 / u1;
  } else if (std::isinf(u1)) {
    if ((u1 < 0.0) != (u0 < 0.0)) {
      y = u1;
    }
  } else {
    boolean_T yEq;
    y = std::fmod(u0, u1);
    yEq = (y == 0.0);
    if ((!yEq) && (u1 > std::floor(u1))) {
      real_T q;
      q = std::abs(u0 / u1);
      yEq = !(std::abs(q - std::floor(q + 0.5)) > DBL_EPSILON * q);
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
  if (std::isnan(u0) || std::isnan(u1) || std::isinf(u0)) {
    y = (rtNaN);
  } else if (std::isinf(u1)) {
    y = u0;
  } else if ((u1 != 0.0) && (u1 != std::trunc(u1))) {
    real_T q;
    q = std::abs(u0 / u1);
    if (!(std::abs(q - std::floor(q + 0.5)) > DBL_EPSILON * q)) {
      y = 0.0 * u0;
    } else {
      y = std::fmod(u0, u1);
    }
  } else {
    y = std::fmod(u0, u1);
  }

  return y;
}

real_T rt_powd_snf(real_T u0, real_T u1)
{
  real_T y;
  if (std::isnan(u0) || std::isnan(u1)) {
    y = (rtNaN);
  } else {
    real_T tmp;
    real_T tmp_0;
    tmp = std::abs(u0);
    tmp_0 = std::abs(u1);
    if (std::isinf(u1)) {
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
      y = std::sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > std::floor(u1))) {
      y = (rtNaN);
    } else {
      y = std::pow(u0, u1);
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
  std::memcpy(&A[0], &u1[0], 9U * sizeof(real_T));
  r1 = 0;
  r2 = 1;
  r3 = 2;
  maxval = std::abs(u1[0]);
  a21 = std::abs(u1[1]);
  if (a21 > maxval) {
    maxval = a21;
    r1 = 1;
    r2 = 0;
  }

  if (std::abs(u1[2]) > maxval) {
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
  if (std::abs(A[r3 + 3]) > std::abs(A[r2 + 3])) {
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
void MARS::step()
{
  /* local block i/o variables */
  real_T rtb_VectorConcatenate[7];
  __m128d tmp_3;
  __m128d tmp_4;
  real_T rtb_Product4_g_tmp[9];
  real_T rtb_VectorConcatenate_b[9];
  real_T rtb_VectorConcatenate_f[9];
  real_T rtb_Sum2[3];
  real_T rtb_Sum2_0[3];
  real_T rtb_Sum_hp[3];
  real_T tmp_0[3];
  real_T phi_tmp;
  real_T rtb_ECEFPositiontoLLA_o1_idx_0;
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
  real_T rtb_q0dot;
  real_T rtb_sincos_o2_idx_1;
  real_T uTmp_idx_0;
  real_T uTmp_idx_1;
  int32_T Product4_tmp;
  int32_T i;
  int8_T rtAction;
  boolean_T rtb_Compare_a[9];
  if (rtmIsMajorTimeStep((&MARS_M))) {
    /* set solver stop time */
    rtsiSetSolverStopTime(&(&MARS_M)->solverInfo,(((&MARS_M)->Timing.clockTick0+
      1)*(&MARS_M)->Timing.stepSize0));
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep((&MARS_M))) {
    (&MARS_M)->Timing.t[0] = rtsiGetT(&(&MARS_M)->solverInfo);
  }

  /* Integrator: '<S5>/p1' */
  if (MARS_DW.p1_IWORK != 0) {
    MARS_X.p1_CSTATE[0] = 2.5914855245392346E+6;
    MARS_X.p1_CSTATE[1] = -1.0790178821353326E+6;
    MARS_X.p1_CSTATE[2] = 5.7079969783558194E+6;
  }

  /* ECEF2LLA: '<S1>/ECEF Position to LLA' incorporates:
   *  Integrator: '<S5>/p1'
   */
  uTmp_idx_0 = MARS_X.p1_CSTATE[0];
  uTmp_idx_1 = MARS_X.p1_CSTATE[1];
  rtb_q0dot = std::sqrt(MARS_X.p1_CSTATE[0] * MARS_X.p1_CSTATE[0] +
                        MARS_X.p1_CSTATE[1] * MARS_X.p1_CSTATE[1]);
  rtb_ixk = rt_atan2d_snf(MARS_X.p1_CSTATE[2], 0.99664718933525254 * rtb_q0dot);
  phi_tmp = std::sin(rtb_ixk);
  rtb_jxi = std::cos(rtb_ixk);
  phi_tmp = rt_atan2d_snf(42841.311513313573 * phi_tmp * phi_tmp * phi_tmp +
    MARS_X.p1_CSTATE[2], rtb_q0dot - 42697.672707179969 * rtb_jxi * rtb_jxi *
    rtb_jxi);
  rtb_jxi = rt_atan2d_snf(0.99664718933525254 * std::sin(phi_tmp), std::cos
    (phi_tmp));
  i = 0;
  while ((rtb_ixk != rtb_jxi) && (i < 5)) {
    rtb_ixk = rtb_jxi;
    phi_tmp = std::sin(rtb_jxi);
    rtb_jxi = std::cos(rtb_jxi);
    phi_tmp = rt_atan2d_snf(42841.311513313573 * phi_tmp * phi_tmp * phi_tmp +
      MARS_X.p1_CSTATE[2], std::sqrt(uTmp_idx_0 * uTmp_idx_0 + uTmp_idx_1 *
      uTmp_idx_1) - 42697.672707179969 * rtb_jxi * rtb_jxi * rtb_jxi);
    rtb_jxi = rt_atan2d_snf(0.99664718933525254 * std::sin(phi_tmp), std::cos
      (phi_tmp));
    i++;
  }

  rtb_kxj = std::abs(phi_tmp);
  rtb_ixk = phi_tmp;
  rtb_jxi = rt_atan2d_snf(MARS_X.p1_CSTATE[1], MARS_X.p1_CSTATE[0]);
  if (rtb_kxj > 3.1415926535897931) {
    rtb_ixk = rt_modd_snf(phi_tmp + 3.1415926535897931, 6.2831853071795862) -
      3.1415926535897931;
    rtb_kxj = std::abs(rtb_ixk);
  }

  if (rtb_kxj > 1.5707963267948966) {
    rtb_jxi += 3.1415926535897931;
    if (!std::isnan(rtb_ixk)) {
      if (rtb_ixk < 0.0) {
        rtb_ixk = -1.0;
      } else {
        rtb_ixk = (rtb_ixk > 0.0);
      }
    }

    rtb_ixk *= 1.5707963267948966 - (rtb_kxj - 1.5707963267948966);
  }

  if (std::abs(rtb_jxi) > 3.1415926535897931) {
    rtb_jxi = rt_remd_snf(rtb_jxi, 6.2831853071795862);
    rtb_jxi -= std::trunc(rtb_jxi / 3.1415926535897931) * 6.2831853071795862;
  }

  rtb_ECEFPositiontoLLA_o1_idx_0 = rtb_ixk * 180.0 / 3.1415926535897931;
  rtb_kxj = rtb_jxi * 180.0 / 3.1415926535897931;
  rtb_ixk = std::sin(phi_tmp);
  rtb_q0dot = ((6.378137E+6 / std::sqrt(1.0 - rtb_ixk * rtb_ixk *
    0.0066943799901413165) * 0.0066943799901413165 * rtb_ixk + MARS_X.p1_CSTATE
                [2]) * rtb_ixk + rtb_q0dot * std::cos(phi_tmp)) - 6.378137E+6 /
    std::sqrt(1.0 - std::sin(phi_tmp) * std::sin(phi_tmp) *
              0.0066943799901413165);

  /* End of ECEF2LLA: '<S1>/ECEF Position to LLA' */
  if (rtmIsMajorTimeStep((&MARS_M))) {
    /* If: '<S17>/If' */
    if (rtsiIsModeUpdateTimeStep(&(&MARS_M)->solverInfo)) {
      rtAction = 1;
      MARS_DW.If_ActiveSubsystem = 1;
    } else {
      rtAction = MARS_DW.If_ActiveSubsystem;
    }

    switch (rtAction) {
     case 0:
      /* Outputs for IfAction SubSystem: '<S17>/Positive Trace' incorporates:
       *  ActionPort: '<S68>/Action Port'
       */
      /* Gain: '<S68>/Gain' incorporates:
       *  Merge: '<S17>/Merge'
       */
      MARS_B.Merge[0] = 0.22088421489473928;

      /* Product: '<S68>/Product' incorporates:
       *  Merge: '<S17>/Merge'
       */
      MARS_B.Merge[1] = MARS_ConstB.Add_jx / 0.88353685957895711;
      MARS_B.Merge[2] = MARS_ConstB.Add_o / 0.88353685957895711;
      MARS_B.Merge[3] = MARS_ConstB.Add_n / 0.88353685957895711;

      /* End of Outputs for SubSystem: '<S17>/Positive Trace' */
      break;

     case 1:
      /* Outputs for IfAction SubSystem: '<S17>/Negative Trace' incorporates:
       *  ActionPort: '<S67>/Action Port'
       */
      /* If: '<S67>/Find Maximum Diagonal Value' */
      if ((MARS_ConstB.Product2[4] > MARS_ConstB.Product2[0]) &&
          (MARS_ConstB.Product2[4] > MARS_ConstB.Product2[8])) {
        /* Outputs for IfAction SubSystem: '<S67>/Maximum Value at DCM(2,2)' incorporates:
         *  ActionPort: '<S72>/Action Port'
         */
        /* Gain: '<S72>/Gain' incorporates:
         *  Merge: '<S17>/Merge'
         */
        MARS_B.Merge[2] = 0.955404187046654;

        /* Gain: '<S72>/Gain1' incorporates:
         *  Merge: '<S17>/Merge'
         */
        MARS_B.Merge[1] = MARS_ConstB.Product_d[0];

        /* Gain: '<S72>/Gain3' incorporates:
         *  Merge: '<S17>/Merge'
         */
        MARS_B.Merge[3] = MARS_ConstB.Product_d[1];

        /* Gain: '<S72>/Gain4' incorporates:
         *  Merge: '<S17>/Merge'
         */
        MARS_B.Merge[0] = MARS_ConstB.Product_d[2];

        /* End of Outputs for SubSystem: '<S67>/Maximum Value at DCM(2,2)' */
      } else if (MARS_ConstB.Product2[8] > MARS_ConstB.Product2[0]) {
        /* Outputs for IfAction SubSystem: '<S67>/Maximum Value at DCM(3,3)' incorporates:
         *  ActionPort: '<S73>/Action Port'
         */
        /* Gain: '<S73>/Gain' incorporates:
         *  Merge: '<S17>/Merge'
         */
        MARS_B.Merge[3] = 0.044147846096478843;

        /* Gain: '<S73>/Gain1' incorporates:
         *  Merge: '<S17>/Merge'
         */
        MARS_B.Merge[1] = MARS_ConstB.Product_l[0];

        /* Gain: '<S73>/Gain2' incorporates:
         *  Merge: '<S17>/Merge'
         */
        MARS_B.Merge[2] = MARS_ConstB.Product_l[1];

        /* Gain: '<S73>/Gain3' incorporates:
         *  Merge: '<S17>/Merge'
         */
        MARS_B.Merge[0] = MARS_ConstB.Product_l[2];

        /* End of Outputs for SubSystem: '<S67>/Maximum Value at DCM(3,3)' */
      } else {
        /* Outputs for IfAction SubSystem: '<S67>/Maximum Value at DCM(1,1)' incorporates:
         *  ActionPort: '<S71>/Action Port'
         */
        /* Gain: '<S71>/Gain' incorporates:
         *  Merge: '<S17>/Merge'
         */
        MARS_B.Merge[1] = 0.19095541539610353;

        /* Gain: '<S71>/Gain1' incorporates:
         *  Merge: '<S17>/Merge'
         */
        MARS_B.Merge[2] = MARS_ConstB.Product_n[0];

        /* Gain: '<S71>/Gain2' incorporates:
         *  Merge: '<S17>/Merge'
         */
        MARS_B.Merge[3] = MARS_ConstB.Product_n[1];

        /* Gain: '<S71>/Gain3' incorporates:
         *  Merge: '<S17>/Merge'
         */
        MARS_B.Merge[0] = MARS_ConstB.Product_n[2];

        /* End of Outputs for SubSystem: '<S67>/Maximum Value at DCM(1,1)' */
      }

      /* End of If: '<S67>/Find Maximum Diagonal Value' */
      /* End of Outputs for SubSystem: '<S17>/Negative Trace' */
      break;
    }

    /* End of If: '<S17>/If' */
  }

  /* Integrator: '<S4>/q' */
  if (MARS_DW.q_IWORK != 0) {
    MARS_X.q_CSTATE[0] = MARS_B.Merge[0];
    MARS_X.q_CSTATE[1] = MARS_B.Merge[1];
    MARS_X.q_CSTATE[2] = MARS_B.Merge[2];
    MARS_X.q_CSTATE[3] = MARS_B.Merge[3];
  }

  /* Sqrt: '<S116>/sqrt' incorporates:
   *  Integrator: '<S4>/q'
   *  Product: '<S117>/Product'
   *  Product: '<S117>/Product1'
   *  Product: '<S117>/Product2'
   *  Product: '<S117>/Product3'
   *  Sqrt: '<S120>/sqrt'
   *  Sum: '<S117>/Sum'
   */
  phi_tmp = std::sqrt(((MARS_X.q_CSTATE[0] * MARS_X.q_CSTATE[0] +
                        MARS_X.q_CSTATE[1] * MARS_X.q_CSTATE[1]) +
                       MARS_X.q_CSTATE[2] * MARS_X.q_CSTATE[2]) +
                      MARS_X.q_CSTATE[3] * MARS_X.q_CSTATE[3]);

  /* Product: '<S115>/Product' incorporates:
   *  Integrator: '<S4>/q'
   *  Sqrt: '<S116>/sqrt'
   */
  rtb_ixk = MARS_X.q_CSTATE[0] / phi_tmp;

  /* Product: '<S115>/Product1' incorporates:
   *  Integrator: '<S4>/q'
   *  Sqrt: '<S116>/sqrt'
   */
  rtb_jxi = MARS_X.q_CSTATE[1] / phi_tmp;

  /* Product: '<S115>/Product2' incorporates:
   *  Integrator: '<S4>/q'
   *  Product: '<S119>/Product2'
   *  Sqrt: '<S116>/sqrt'
   */
  uTmp_idx_0 = MARS_X.q_CSTATE[2] / phi_tmp;

  /* Product: '<S115>/Product3' incorporates:
   *  Integrator: '<S4>/q'
   *  Product: '<S119>/Product3'
   *  Sqrt: '<S116>/sqrt'
   */
  uTmp_idx_1 = MARS_X.q_CSTATE[3] / phi_tmp;

  /* Product: '<S105>/Product3' incorporates:
   *  Product: '<S109>/Product3'
   */
  rtb_sincos_o2_idx_1 = rtb_ixk * rtb_ixk;

  /* Product: '<S105>/Product2' incorporates:
   *  Product: '<S109>/Product2'
   */
  rtb_VectorConcatenate_f_tmp_1 = rtb_jxi * rtb_jxi;

  /* Product: '<S105>/Product1' incorporates:
   *  Product: '<S109>/Product1'
   *  Product: '<S113>/Product1'
   *  Product: '<S115>/Product2'
   */
  rtb_VectorConcatenate_f_tmp_2 = uTmp_idx_0 * uTmp_idx_0;

  /* Product: '<S105>/Product' incorporates:
   *  Product: '<S109>/Product'
   *  Product: '<S113>/Product'
   *  Product: '<S115>/Product3'
   */
  rtb_VectorConcatenate_f_tmp_3 = uTmp_idx_1 * uTmp_idx_1;

  /* Sum: '<S105>/Sum' incorporates:
   *  Product: '<S105>/Product'
   *  Product: '<S105>/Product1'
   *  Product: '<S105>/Product2'
   *  Product: '<S105>/Product3'
   */
  rtb_VectorConcatenate_f[0] = ((rtb_sincos_o2_idx_1 +
    rtb_VectorConcatenate_f_tmp_1) - rtb_VectorConcatenate_f_tmp_2) -
    rtb_VectorConcatenate_f_tmp_3;

  /* Product: '<S108>/Product3' incorporates:
   *  Product: '<S106>/Product3'
   *  Product: '<S115>/Product3'
   */
  rtb_VectorConcatenate_f_tmp = uTmp_idx_1 * rtb_ixk;

  /* Product: '<S108>/Product2' incorporates:
   *  Product: '<S106>/Product2'
   *  Product: '<S115>/Product2'
   */
  rtb_VectorConcatenate_f_tmp_0 = rtb_jxi * uTmp_idx_0;

  /* Gain: '<S108>/Gain' incorporates:
   *  Product: '<S108>/Product2'
   *  Product: '<S108>/Product3'
   *  Sum: '<S108>/Sum'
   */
  rtb_VectorConcatenate_f[1] = (rtb_VectorConcatenate_f_tmp_0 -
    rtb_VectorConcatenate_f_tmp) * 2.0;

  /* Product: '<S111>/Product2' incorporates:
   *  Product: '<S107>/Product2'
   *  Product: '<S115>/Product3'
   */
  rtb_VectorConcatenate_f_tmp_4 = rtb_jxi * uTmp_idx_1;

  /* Product: '<S111>/Product1' incorporates:
   *  Product: '<S107>/Product1'
   *  Product: '<S115>/Product2'
   */
  rtb_VectorConcatenate_f_tmp_5 = rtb_ixk * uTmp_idx_0;

  /* Gain: '<S111>/Gain' incorporates:
   *  Product: '<S111>/Product1'
   *  Product: '<S111>/Product2'
   *  Sum: '<S111>/Sum'
   */
  rtb_VectorConcatenate_f[2] = (rtb_VectorConcatenate_f_tmp_5 +
    rtb_VectorConcatenate_f_tmp_4) * 2.0;

  /* Gain: '<S106>/Gain' incorporates:
   *  Sum: '<S106>/Sum'
   */
  rtb_VectorConcatenate_f[3] = (rtb_VectorConcatenate_f_tmp +
    rtb_VectorConcatenate_f_tmp_0) * 2.0;

  /* Sum: '<S109>/Sum' incorporates:
   *  Sum: '<S113>/Sum'
   */
  rtb_sincos_o2_idx_1 -= rtb_VectorConcatenate_f_tmp_1;
  rtb_VectorConcatenate_f[4] = (rtb_sincos_o2_idx_1 +
    rtb_VectorConcatenate_f_tmp_2) - rtb_VectorConcatenate_f_tmp_3;

  /* Product: '<S112>/Product1' incorporates:
   *  Product: '<S110>/Product1'
   */
  rtb_VectorConcatenate_f_tmp_1 = rtb_ixk * rtb_jxi;

  /* Product: '<S112>/Product2' incorporates:
   *  Product: '<S110>/Product2'
   *  Product: '<S115>/Product2'
   *  Product: '<S115>/Product3'
   */
  rtb_VectorConcatenate_f_tmp = uTmp_idx_0 * uTmp_idx_1;

  /* Gain: '<S112>/Gain' incorporates:
   *  Product: '<S112>/Product1'
   *  Product: '<S112>/Product2'
   *  Sum: '<S112>/Sum'
   */
  rtb_VectorConcatenate_f[5] = (rtb_VectorConcatenate_f_tmp -
    rtb_VectorConcatenate_f_tmp_1) * 2.0;

  /* Gain: '<S107>/Gain' incorporates:
   *  Sum: '<S107>/Sum'
   */
  rtb_VectorConcatenate_f[6] = (rtb_VectorConcatenate_f_tmp_4 -
    rtb_VectorConcatenate_f_tmp_5) * 2.0;

  /* Gain: '<S110>/Gain' incorporates:
   *  Sum: '<S110>/Sum'
   */
  rtb_VectorConcatenate_f[7] = (rtb_VectorConcatenate_f_tmp_1 +
    rtb_VectorConcatenate_f_tmp) * 2.0;

  /* Sum: '<S113>/Sum' */
  rtb_VectorConcatenate_f[8] = (rtb_sincos_o2_idx_1 -
    rtb_VectorConcatenate_f_tmp_2) + rtb_VectorConcatenate_f_tmp_3;

  /* UnitConversion: '<S30>/Unit Conversion' incorporates:
   *  UnitConversion: '<S140>/Unit Conversion'
   */
  /* Unit Conversion - from: deg to: rad
     Expression: output = (0.0174533*input) + (0) */
  rtb_ECEFPositiontoLLA_o1_idx_0 *= 0.017453292519943295;

  /* Trigonometry: '<S14>/sincos' incorporates:
   *  Trigonometry: '<S13>/sine'
   *  UnitConversion: '<S30>/Unit Conversion'
   */
  rtb_ixk = std::cos(rtb_ECEFPositiontoLLA_o1_idx_0);
  rtb_VectorConcatenate_f_tmp_1 = std::sin(rtb_ECEFPositiontoLLA_o1_idx_0);

  /* UnitConversion: '<S30>/Unit Conversion' */
  rtb_kxj *= 0.017453292519943295;

  /* Trigonometry: '<S14>/sincos' */
  rtb_sincos_o2_idx_1 = std::cos(rtb_kxj);
  rtb_kxj = std::sin(rtb_kxj);

  /* UnaryMinus: '<S21>/Unary Minus' incorporates:
   *  Product: '<S21>/u(1)*u(4)'
   *  Trigonometry: '<S14>/sincos'
   */
  rtb_VectorConcatenate_b[0] = -(rtb_VectorConcatenate_f_tmp_1 *
    rtb_sincos_o2_idx_1);

  /* UnaryMinus: '<S24>/Unary Minus' */
  rtb_VectorConcatenate_b[1] = -rtb_kxj;

  /* UnaryMinus: '<S27>/Unary Minus' incorporates:
   *  Product: '<S27>/u(3)*u(4)'
   */
  rtb_VectorConcatenate_b[2] = -(rtb_ixk * rtb_sincos_o2_idx_1);

  /* UnaryMinus: '<S22>/Unary Minus' incorporates:
   *  Product: '<S22>/u(1)*u(2)'
   *  Trigonometry: '<S14>/sincos'
   */
  rtb_VectorConcatenate_b[3] = -(rtb_VectorConcatenate_f_tmp_1 * rtb_kxj);

  /* SignalConversion generated from: '<S31>/Vector Concatenate' */
  rtb_VectorConcatenate_b[4] = rtb_sincos_o2_idx_1;

  /* UnaryMinus: '<S28>/Unary Minus' incorporates:
   *  Product: '<S28>/u(2)*u(3)'
   */
  rtb_VectorConcatenate_b[5] = -(rtb_ixk * rtb_kxj);

  /* SignalConversion generated from: '<S31>/Vector Concatenate' */
  rtb_VectorConcatenate_b[6] = rtb_ixk;

  /* SignalConversion generated from: '<S31>/Vector Concatenate' incorporates:
   *  Constant: '<S26>/Constant'
   */
  rtb_VectorConcatenate_b[7] = 0.0;

  /* UnaryMinus: '<S29>/Unary Minus' incorporates:
   *  Trigonometry: '<S14>/sincos'
   */
  rtb_VectorConcatenate_b[8] = -rtb_VectorConcatenate_f_tmp_1;
  for (i = 0; i < 3; i++) {
    for (int32_T i_0{0}; i_0 <= 0; i_0 += 2) {
      /* Product: '<S4>/Product4' incorporates:
       *  Math: '<S4>/Math Function2'
       */
      Product4_tmp = 3 * i + i_0;
      _mm_storeu_pd(&MARS_B.Product4[Product4_tmp], _mm_set1_pd(0.0));

      /* Product: '<S4>/Product4' incorporates:
       *  Concatenate: '<S114>/Vector Concatenate'
       *  Concatenate: '<S139>/Vector Concatenate'
       *  Math: '<S4>/Math Function2'
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

    for (int32_T i_0{2}; i_0 < 3; i_0++) {
      /* Product: '<S4>/Product4' incorporates:
       *  Concatenate: '<S114>/Vector Concatenate'
       *  Concatenate: '<S139>/Vector Concatenate'
       *  Math: '<S4>/Math Function2'
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

  /* Gain: '<S45>/Gain1' incorporates:
   *  Concatenate: '<S45>/Vector Concatenate'
   *  Product: '<S4>/Product4'
   *  Selector: '<S45>/Selector1'
   */
  rtb_VectorConcatenate[0] = MARS_B.Product4[3];
  rtb_VectorConcatenate[1] = MARS_B.Product4[0];
  rtb_VectorConcatenate[2] = -MARS_B.Product4[6];

  /* Gain: '<S45>/Gain2' incorporates:
   *  Concatenate: '<S45>/Vector Concatenate'
   *  Product: '<S4>/Product4'
   *  Selector: '<S45>/Selector2'
   */
  rtb_VectorConcatenate[3] = MARS_B.Product4[7];

  /* Gain: '<S45>/Gain3' incorporates:
   *  Concatenate: '<S45>/Vector Concatenate'
   *  Product: '<S4>/Product4'
   *  Selector: '<S45>/Selector3'
   */
  rtb_VectorConcatenate[5] = -MARS_B.Product4[1];

  /* Gain: '<S45>/Gain2' incorporates:
   *  Concatenate: '<S45>/Vector Concatenate'
   *  Product: '<S4>/Product4'
   *  Selector: '<S45>/Selector2'
   */
  rtb_VectorConcatenate[4] = MARS_B.Product4[8];

  /* Gain: '<S45>/Gain3' incorporates:
   *  Concatenate: '<S45>/Vector Concatenate'
   *  Product: '<S4>/Product4'
   *  Selector: '<S45>/Selector3'
   */
  rtb_VectorConcatenate[6] = MARS_B.Product4[4];

  /* If: '<S16>/If' */
  if (rtsiIsModeUpdateTimeStep(&(&MARS_M)->solverInfo)) {
    rtAction = static_cast<int8_T>((!(rtb_VectorConcatenate[2] >= 1.0)) &&
      (!(rtb_VectorConcatenate[2] <= -1.0)));
    MARS_DW.If_ActiveSubsystem_h = rtAction;
  } else {
    rtAction = MARS_DW.If_ActiveSubsystem_h;
  }

  switch (rtAction) {
   case 0:
    /* Outputs for IfAction SubSystem: '<S16>/AxisRotZeroR3' incorporates:
     *  ActionPort: '<S44>/Action Port'
     */
    /* If: '<S51>/If' */
    if (rtsiIsModeUpdateTimeStep(&(&MARS_M)->solverInfo)) {
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
      /* Outputs for IfAction SubSystem: '<S51>/If Action Subsystem' incorporates:
       *  ActionPort: '<S52>/Action Port'
       */
      if (rtmIsMajorTimeStep((&MARS_M))) {
        /* Merge: '<S51>/Merge' incorporates:
         *  Constant: '<S52>/Constant'
         */
        MARS_B.Merge_j2 = 1.0;
      }

      /* End of Outputs for SubSystem: '<S51>/If Action Subsystem' */
      break;

     case 1:
      /* Outputs for IfAction SubSystem: '<S51>/If Action Subsystem1' incorporates:
       *  ActionPort: '<S53>/Action Port'
       */
      if (rtmIsMajorTimeStep((&MARS_M))) {
        /* Merge: '<S51>/Merge' incorporates:
         *  Constant: '<S53>/Constant'
         */
        MARS_B.Merge_j2 = 1.0;
      }

      /* End of Outputs for SubSystem: '<S51>/If Action Subsystem1' */
      break;

     case 2:
      /* Outputs for IfAction SubSystem: '<S51>/If Action Subsystem2' incorporates:
       *  ActionPort: '<S54>/Action Port'
       */
      MARS_IfActionSubsystem2(rtb_VectorConcatenate[2], &MARS_B.Merge_j2);

      /* End of Outputs for SubSystem: '<S51>/If Action Subsystem2' */
      break;
    }

    /* End of If: '<S51>/If' */
    /* End of Outputs for SubSystem: '<S16>/AxisRotZeroR3' */
    break;

   case 1:
    /* Outputs for IfAction SubSystem: '<S16>/AxisRotDefault' incorporates:
     *  ActionPort: '<S43>/Action Port'
     */
    /* If: '<S47>/If' */
    if (rtsiIsModeUpdateTimeStep(&(&MARS_M)->solverInfo)) {
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
      /* Outputs for IfAction SubSystem: '<S47>/If Action Subsystem' incorporates:
       *  ActionPort: '<S48>/Action Port'
       */
      if (rtmIsMajorTimeStep((&MARS_M))) {
        /* Merge: '<S47>/Merge' incorporates:
         *  Constant: '<S48>/Constant'
         */
        MARS_B.Merge_j = 1.0;
      }

      /* End of Outputs for SubSystem: '<S47>/If Action Subsystem' */
      break;

     case 1:
      /* Outputs for IfAction SubSystem: '<S47>/If Action Subsystem1' incorporates:
       *  ActionPort: '<S49>/Action Port'
       */
      if (rtmIsMajorTimeStep((&MARS_M))) {
        /* Merge: '<S47>/Merge' incorporates:
         *  Constant: '<S49>/Constant'
         */
        MARS_B.Merge_j = 1.0;
      }

      /* End of Outputs for SubSystem: '<S47>/If Action Subsystem1' */
      break;

     case 2:
      /* Outputs for IfAction SubSystem: '<S47>/If Action Subsystem2' incorporates:
       *  ActionPort: '<S50>/Action Port'
       */
      MARS_IfActionSubsystem2(rtb_VectorConcatenate[2], &MARS_B.Merge_j);

      /* End of Outputs for SubSystem: '<S47>/If Action Subsystem2' */
      break;
    }

    /* End of If: '<S47>/If' */
    /* End of Outputs for SubSystem: '<S16>/AxisRotDefault' */
    break;
  }

  /* End of If: '<S16>/If' */
  if (rtmIsMajorTimeStep((&MARS_M))) {
    boolean_T tmp;

    /* If: '<S46>/If1' */
    rtAction = -1;
    if (rtsiIsModeUpdateTimeStep(&(&MARS_M)->solverInfo)) {
      MARS_DW.If1_ActiveSubsystem = -1;
    } else {
      rtAction = MARS_DW.If1_ActiveSubsystem;
    }

    if (rtAction == 0) {
      /* Outputs for IfAction SubSystem: '<S46>/If Warning//Error' incorporates:
       *  ActionPort: '<S55>/if'
       */
      /* Bias: '<S58>/Bias1' incorporates:
       *  Math: '<S58>/Math Function'
       *  Product: '<S4>/Product4'
       *  Product: '<S58>/Product'
       */
      for (i = 0; i < 3; i++) {
        for (int32_T i_0{0}; i_0 < 3; i_0++) {
          Product4_tmp = 3 * i + i_0;
          rtb_Product4_g_tmp[Product4_tmp] = ((MARS_B.Product4[3 * i_0 + 1] *
            MARS_B.Product4[3 * i + 1] + MARS_B.Product4[3 * i_0] *
            MARS_B.Product4[3 * i]) + MARS_B.Product4[3 * i_0 + 2] *
            MARS_B.Product4[3 * i + 2]) + MARS_ConstP.pooled5[Product4_tmp];
        }
      }

      /* End of Bias: '<S58>/Bias1' */

      /* RelationalOperator: '<S64>/Compare' incorporates:
       *  Abs: '<S58>/Abs2'
       *  Constant: '<S64>/Constant'
       */
      for (i = 0; i < 9; i++) {
        rtb_Compare_a[i] = (std::abs(rtb_Product4_g_tmp[i]) >
                            4.4408920985006262E-16);
      }

      /* End of RelationalOperator: '<S64>/Compare' */

      /* Logic: '<S58>/Logical Operator1' incorporates:
       *  RelationalOperator: '<S64>/Compare'
       */
      tmp = rtb_Compare_a[0];
      for (i = 0; i < 8; i++) {
        tmp = (tmp || rtb_Compare_a[i + 1]);
      }

      /* If: '<S55>/If' incorporates:
       *  Abs: '<S59>/Abs1'
       *  Bias: '<S59>/Bias'
       *  Constant: '<S66>/Constant'
       *  Logic: '<S58>/Logical Operator1'
       *  Product: '<S4>/Product4'
       *  Product: '<S65>/Product'
       *  Product: '<S65>/Product1'
       *  Product: '<S65>/Product2'
       *  Product: '<S65>/Product3'
       *  Product: '<S65>/Product4'
       *  Product: '<S65>/Product5'
       *  RelationalOperator: '<S66>/Compare'
       *  Reshape: '<S65>/Reshape'
       *  Sum: '<S65>/Sum'
       */
      if (std::abs((((((MARS_B.Product4[0] * MARS_B.Product4[4] *
                        MARS_B.Product4[8] - MARS_B.Product4[0] *
                        MARS_B.Product4[5] * MARS_B.Product4[7]) -
                       MARS_B.Product4[1] * MARS_B.Product4[3] *
                       MARS_B.Product4[8]) + MARS_B.Product4[2] *
                      MARS_B.Product4[3] * MARS_B.Product4[7]) +
                     MARS_B.Product4[1] * MARS_B.Product4[5] * MARS_B.Product4[6])
                    - MARS_B.Product4[2] * MARS_B.Product4[4] * MARS_B.Product4
                    [6]) + -1.0) > 4.4408920985006262E-16) {
        /* Outputs for IfAction SubSystem: '<S55>/If Not Proper' incorporates:
         *  ActionPort: '<S57>/Action Port'
         */
        MARS_IfNotProper(1.0);

        /* End of Outputs for SubSystem: '<S55>/If Not Proper' */
      } else if (tmp) {
        /* Outputs for IfAction SubSystem: '<S55>/Else If Not Orthogonal' incorporates:
         *  ActionPort: '<S56>/Action Port'
         */
        MARS_ElseIfNotOrthogonal(1.0);

        /* End of Outputs for SubSystem: '<S55>/Else If Not Orthogonal' */
      }

      /* End of If: '<S55>/If' */
      /* End of Outputs for SubSystem: '<S46>/If Warning//Error' */
    }

    /* End of If: '<S46>/If1' */

    /* If: '<S69>/If1' */
    rtAction = -1;
    if (rtsiIsModeUpdateTimeStep(&(&MARS_M)->solverInfo)) {
      MARS_DW.If1_ActiveSubsystem_l = -1;
    } else {
      rtAction = MARS_DW.If1_ActiveSubsystem_l;
    }

    if (rtAction == 0) {
      /* Outputs for IfAction SubSystem: '<S69>/If Warning//Error' incorporates:
       *  ActionPort: '<S93>/if'
       */
      /* RelationalOperator: '<S102>/Compare' incorporates:
       *  Abs: '<S96>/Abs2'
       *  Constant: '<S102>/Constant'
       */
      for (i = 0; i < 9; i++) {
        rtb_Compare_a[i] = (MARS_ConstB.Abs2[i] > 4.4408920985006262E-16);
      }

      /* End of RelationalOperator: '<S102>/Compare' */

      /* Logic: '<S96>/Logical Operator1' incorporates:
       *  RelationalOperator: '<S102>/Compare'
       */
      tmp = rtb_Compare_a[0];
      for (i = 0; i < 8; i++) {
        tmp = (tmp || rtb_Compare_a[i + 1]);
      }

      /* If: '<S93>/If' incorporates:
       *  Logic: '<S96>/Logical Operator1'
       */
      if (tmp) {
        /* Outputs for IfAction SubSystem: '<S93>/Else If Not Orthogonal' incorporates:
         *  ActionPort: '<S94>/Action Port'
         */
        MARS_ElseIfNotOrthogonal(1.0);

        /* End of Outputs for SubSystem: '<S93>/Else If Not Orthogonal' */
      }

      /* End of If: '<S93>/If' */
      /* End of Outputs for SubSystem: '<S69>/If Warning//Error' */
    }

    /* End of If: '<S69>/If1' */
  }

  /* Trigonometry: '<S10>/sincos' incorporates:
   *  Integrator: '<S8>/Integrator'
   */
  rtb_ixk = std::sin(MARS_X.Integrator_CSTATE);
  rtb_jxi = std::cos(MARS_X.Integrator_CSTATE);

  /* SignalConversion generated from: '<S139>/Vector Concatenate' */
  rtb_VectorConcatenate_b[0] = rtb_jxi;

  /* SignalConversion generated from: '<S139>/Vector Concatenate' */
  rtb_VectorConcatenate_b[1] = rtb_ixk;

  /* SignalConversion generated from: '<S139>/Vector Concatenate' incorporates:
   *  Constant: '<S10>/Zero'
   */
  rtb_VectorConcatenate_b[2] = 0.0;

  /* UnaryMinus: '<S10>/Unary Minus' */
  rtb_VectorConcatenate_b[3] = -rtb_ixk;

  /* SignalConversion generated from: '<S139>/Vector Concatenate' */
  rtb_VectorConcatenate_b[4] = rtb_jxi;

  /* SignalConversion generated from: '<S139>/Vector Concatenate' incorporates:
   *  Constant: '<S10>/Zero'
   */
  rtb_VectorConcatenate_b[5] = 0.0;

  /* SignalConversion generated from: '<S139>/Vector Concatenate' incorporates:
   *  Constant: '<S10>/Zero'
   */
  rtb_VectorConcatenate_b[6] = 0.0;

  /* SignalConversion generated from: '<S139>/Vector Concatenate' incorporates:
   *  Constant: '<S10>/Zero'
   */
  rtb_VectorConcatenate_b[7] = 0.0;

  /* SignalConversion generated from: '<S139>/Vector Concatenate' incorporates:
   *  Constant: '<S10>/Zero1'
   */
  rtb_VectorConcatenate_b[8] = 1.0;
  for (i = 0; i < 3; i++) {
    /* Product: '<S1>/Product3' incorporates:
     *  Concatenate: '<S114>/Vector Concatenate'
     *  Product: '<S6>/Product2'
     */
    rtb_Sum_hp[i] = (rtb_VectorConcatenate_f[i + 3] * 0.0 +
                     rtb_VectorConcatenate_f[i] * 0.0) +
      rtb_VectorConcatenate_f[i + 6] * 7.292115E-5;

    /* Product: '<S13>/Product3' incorporates:
     *  Integrator: '<S6>/ub,vb,wb'
     *  Math: '<S13>/Math Function2'
     *  Product: '<S4>/Product4'
     */
    rtb_Sum2[i] = (MARS_B.Product4[3 * i + 1] * MARS_X.ubvbwb_CSTATE[1] +
                   MARS_B.Product4[3 * i] * MARS_X.ubvbwb_CSTATE[0]) +
      MARS_B.Product4[3 * i + 2] * MARS_X.ubvbwb_CSTATE[2];
  }

  /* Sum: '<S13>/Sum2' incorporates:
   *  Constant: '<S13>/f2'
   *  Product: '<S13>/Product'
   */
  /* Unit Conversion - from: deg to: rad
     Expression: output = (0.0174533*input) + (0) */
  rtb_jxi = 1.0 - rtb_VectorConcatenate_f_tmp_1 * rtb_VectorConcatenate_f_tmp_1 *
    0.00669437999014133;

  /* Product: '<S13>/w1' incorporates:
   *  Constant: '<S142>/f4'
   *  Product: '<S142>/N'
   *  Sqrt: '<S142>/sqrt'
   *  Sum: '<S142>/Sum4'
   */
  rtb_kxj = rtb_Sum2[1] / (6.378137E+6 / std::sqrt(rtb_jxi) + rtb_q0dot);

  /* SignalConversion generated from: '<S1>/Product2' incorporates:
   *  Constant: '<S141>/f1'
   *  Constant: '<S141>/f3'
   *  Gain: '<S13>/Gain'
   *  Gain: '<S13>/Gain1'
   *  Math: '<S141>/Math Function'
   *  Product: '<S13>/w2'
   *  Product: '<S13>/w3'
   *  Product: '<S141>/M'
   *  Sum: '<S141>/Sum1'
   *  Trigonometry: '<S13>/tan'
   */
  rtb_q0dot = -(rtb_Sum2[0] / (6.3354393272928195E+6 / rt_powd_snf(rtb_jxi, 1.5)
    + rtb_q0dot));
  rtb_ixk = -(rtb_kxj * std::tan(rtb_ECEFPositiontoLLA_o1_idx_0));
  for (i = 0; i <= 0; i += 2) {
    __m128d tmp_1;
    __m128d tmp_2;

    /* Sum: '<S1>/Sum2' incorporates:
     *  Product: '<S1>/Product2'
     *  Product: '<S4>/Product4'
     */
    tmp_3 = _mm_loadu_pd(&MARS_B.Product4[i + 3]);
    tmp_4 = _mm_loadu_pd(&MARS_B.Product4[i]);
    tmp_1 = _mm_loadu_pd(&MARS_B.Product4[i + 6]);

    /* Product: '<S1>/Product3' incorporates:
     *  Product: '<S1>/Product2'
     *  Sum: '<S1>/Sum2'
     */
    tmp_2 = _mm_loadu_pd(&rtb_Sum_hp[i]);

    /* Sum: '<S1>/Sum2' incorporates:
     *  Product: '<S1>/Product2'
     *  SignalConversion generated from: '<S1>/Product2'
     */
    _mm_storeu_pd(&MARS_B.Sum2[i], _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd
      (tmp_3, _mm_set1_pd(rtb_q0dot)), _mm_mul_pd(tmp_4, _mm_set1_pd(rtb_kxj))),
      _mm_mul_pd(tmp_1, _mm_set1_pd(rtb_ixk))), tmp_2));
  }

  for (i = 2; i < 3; i++) {
    /* Sum: '<S1>/Sum2' incorporates:
     *  Product: '<S1>/Product2'
     *  Product: '<S1>/Product3'
     *  Product: '<S4>/Product4'
     *  SignalConversion generated from: '<S1>/Product2'
     *  Sum: '<S1>/Sum3'
     */
    MARS_B.Sum2[i] = ((MARS_B.Product4[i + 3] * rtb_q0dot + MARS_B.Product4[i] *
                       rtb_kxj) + MARS_B.Product4[i + 6] * rtb_ixk) +
      rtb_Sum_hp[i];
  }

  /* Integrator: '<S1>/p,q,r ' */
  if (MARS_DW.pqr_IWORK != 0) {
    MARS_X.pqr_CSTATE[0] = MARS_B.Sum2[0];
    MARS_X.pqr_CSTATE[1] = MARS_B.Sum2[1];
    MARS_X.pqr_CSTATE[2] = MARS_B.Sum2[2];
  }

  /* Sum: '<S1>/Sum4' incorporates:
   *  Integrator: '<S1>/p,q,r '
   *  Product: '<S1>/Product3'
   */
  rtb_Sum2[0] = MARS_X.pqr_CSTATE[0] - rtb_Sum_hp[0];
  rtb_Sum2[1] = MARS_X.pqr_CSTATE[1] - rtb_Sum_hp[1];
  rtb_Sum2[2] = MARS_X.pqr_CSTATE[2] - rtb_Sum_hp[2];

  /* Product: '<S119>/Product' incorporates:
   *  Integrator: '<S4>/q'
   */
  rtb_jxi = MARS_X.q_CSTATE[0] / phi_tmp;

  /* Product: '<S119>/Product1' incorporates:
   *  Integrator: '<S4>/q'
   */
  rtb_kxj = MARS_X.q_CSTATE[1] / phi_tmp;

  /* SignalConversion generated from: '<S4>/q' incorporates:
   *  Fcn: '<S20>/q0dot'
   *  Fcn: '<S20>/q1dot'
   *  Fcn: '<S20>/q2dot'
   *  Fcn: '<S20>/q3dot'
   */
  MARS_B.TmpSignalConversionAtqInport1[0] = ((rtb_kxj * rtb_Sum2[0] + uTmp_idx_0
    * rtb_Sum2[1]) + uTmp_idx_1 * rtb_Sum2[2]) * -0.5;
  MARS_B.TmpSignalConversionAtqInport1[1] = ((rtb_jxi * rtb_Sum2[0] + uTmp_idx_0
    * rtb_Sum2[2]) - uTmp_idx_1 * rtb_Sum2[1]) * 0.5;
  MARS_B.TmpSignalConversionAtqInport1[2] = ((rtb_jxi * rtb_Sum2[1] + uTmp_idx_1
    * rtb_Sum2[0]) - rtb_kxj * rtb_Sum2[2]) * 0.5;
  MARS_B.TmpSignalConversionAtqInport1[3] = ((rtb_jxi * rtb_Sum2[2] + rtb_kxj *
    rtb_Sum2[1]) - uTmp_idx_0 * rtb_Sum2[0]) * 0.5;
  for (i = 0; i < 3; i++) {
    /* Product: '<S5>/Product1' */
    MARS_B.Product1[i] = 0.0;

    /* Product: '<S5>/Product1' incorporates:
     *  Concatenate: '<S139>/Vector Concatenate'
     */
    phi_tmp = rtb_VectorConcatenate_b[i];

    /* Product: '<S5>/Product1' incorporates:
     *  LLA2ECEF: '<S5>/LLA to ECEF Position'
     */
    MARS_B.Product1[i] += phi_tmp * 2.5914855245392346E+6;

    /* Math: '<S5>/Math Function1' incorporates:
     *  Math: '<S4>/Math Function'
     */
    rtb_Product4_g_tmp[3 * i] = phi_tmp;

    /* Product: '<S5>/Product1' incorporates:
     *  Concatenate: '<S139>/Vector Concatenate'
     */
    phi_tmp = rtb_VectorConcatenate_b[i + 3];

    /* Product: '<S5>/Product1' incorporates:
     *  LLA2ECEF: '<S5>/LLA to ECEF Position'
     */
    MARS_B.Product1[i] += phi_tmp * -1.0790178821353326E+6;

    /* Math: '<S5>/Math Function1' incorporates:
     *  Math: '<S4>/Math Function'
     */
    rtb_Product4_g_tmp[3 * i + 1] = phi_tmp;

    /* Product: '<S5>/Product1' incorporates:
     *  Concatenate: '<S139>/Vector Concatenate'
     */
    phi_tmp = rtb_VectorConcatenate_b[i + 6];

    /* Product: '<S5>/Product1' incorporates:
     *  LLA2ECEF: '<S5>/LLA to ECEF Position'
     */
    MARS_B.Product1[i] += phi_tmp * 5.7079969783558194E+6;

    /* Math: '<S5>/Math Function1' incorporates:
     *  Math: '<S4>/Math Function'
     */
    rtb_Product4_g_tmp[3 * i + 2] = phi_tmp;
  }

  for (i = 0; i < 3; i++) {
    /* Product: '<S5>/Product5' incorporates:
     *  Concatenate: '<S114>/Vector Concatenate'
     *  Integrator: '<S6>/ub,vb,wb'
     *  Math: '<S5>/Math Function2'
     */
    MARS_B.Product5[i] = 0.0;
    MARS_B.Product5[i] += rtb_VectorConcatenate_f[3 * i] * MARS_X.ubvbwb_CSTATE
      [0];
    MARS_B.Product5[i] += rtb_VectorConcatenate_f[3 * i + 1] *
      MARS_X.ubvbwb_CSTATE[1];
    MARS_B.Product5[i] += rtb_VectorConcatenate_f[3 * i + 2] *
      MARS_X.ubvbwb_CSTATE[2];

    /* Product: '<S5>/Product4' incorporates:
     *  Math: '<S5>/Math Function1'
     */
    rtb_Sum2[i] = (rtb_Product4_g_tmp[i + 3] * 0.0 + rtb_Product4_g_tmp[i] * 0.0)
      + rtb_Product4_g_tmp[i + 6] * 7.292115E-5;
  }

  /* Integrator: '<S5>/p' */
  if (MARS_DW.p_IWORK != 0) {
    MARS_X.p_CSTATE[0] = MARS_B.Product1[0];
    MARS_X.p_CSTATE[1] = MARS_B.Product1[1];
    MARS_X.p_CSTATE[2] = MARS_B.Product1[2];
  }

  /* Sum: '<S122>/Sum' incorporates:
   *  Integrator: '<S5>/p'
   *  Product: '<S123>/i x j'
   *  Product: '<S123>/j x k'
   *  Product: '<S123>/k x i'
   *  Product: '<S124>/i x k'
   *  Product: '<S124>/j x i'
   *  Product: '<S124>/k x j'
   */
  tmp_0[0] = MARS_X.p_CSTATE[1] * rtb_Sum2[2];
  tmp_0[1] = rtb_Sum2[0] * MARS_X.p_CSTATE[2];
  tmp_0[2] = MARS_X.p_CSTATE[0] * rtb_Sum2[1];
  rtb_Sum2_0[0] = rtb_Sum2[1] * MARS_X.p_CSTATE[2];
  rtb_Sum2_0[1] = MARS_X.p_CSTATE[0] * rtb_Sum2[2];
  rtb_Sum2_0[2] = rtb_Sum2[0] * MARS_X.p_CSTATE[1];
  for (i = 0; i < 3; i++) {
    /* Sum: '<S5>/Sum2' */
    phi_tmp = 0.0;
    for (int32_T i_0{0}; i_0 < 3; i_0++) {
      /* Math: '<S5>/Math Function' incorporates:
       *  Concatenate: '<S114>/Vector Concatenate'
       *  Math: '<S4>/Math Function'
       *  Product: '<S4>/Product1'
       *  Product: '<S5>/Product2'
       */
      Product4_tmp = 3 * i_0 + i;
      rtb_VectorConcatenate_b[Product4_tmp] = 0.0;
      rtb_VectorConcatenate_b[Product4_tmp] += rtb_Product4_g_tmp[3 * i] *
        rtb_VectorConcatenate_f[i_0];
      rtb_VectorConcatenate_b[Product4_tmp] += rtb_Product4_g_tmp[3 * i + 1] *
        rtb_VectorConcatenate_f[i_0 + 3];
      rtb_VectorConcatenate_b[Product4_tmp] += rtb_Product4_g_tmp[3 * i + 2] *
        rtb_VectorConcatenate_f[i_0 + 6];

      /* Sum: '<S5>/Sum2' incorporates:
       *  Integrator: '<S6>/ub,vb,wb'
       *  Product: '<S5>/Product2'
       */
      phi_tmp += rtb_VectorConcatenate_b[Product4_tmp] *
        MARS_X.ubvbwb_CSTATE[i_0];
    }

    /* Sum: '<S5>/Sum2' incorporates:
     *  Product: '<S5>/Product2'
     *  Sum: '<S122>/Sum'
     */
    MARS_B.Sum2_e[i] = phi_tmp - (tmp_0[i] - rtb_Sum2_0[i]);

    /* Sum: '<S6>/Sum2' incorporates:
     *  Integrator: '<S1>/p,q,r '
     *  Product: '<S6>/Product2'
     */
    rtb_Sum_hp[i] += MARS_X.pqr_CSTATE[i];
  }

  /* Sum: '<S127>/Sum' incorporates:
   *  Constant: '<S1>/omega_earth2'
   *  Constant: '<S1>/omega_earth3'
   *  Integrator: '<S5>/p1'
   *  Product: '<S132>/i x j'
   *  Product: '<S132>/j x k'
   *  Product: '<S132>/k x i'
   *  Product: '<S133>/i x k'
   *  Product: '<S133>/j x i'
   *  Product: '<S133>/k x j'
   */
  rtb_Sum2[0] = 0.0 * MARS_X.p1_CSTATE[2] - 7.292115E-5 * MARS_X.p1_CSTATE[1];
  rtb_Sum2[1] = 7.292115E-5 * MARS_X.p1_CSTATE[0] - 0.0 * MARS_X.p1_CSTATE[2];
  rtb_Sum2[2] = 0.0 * MARS_X.p1_CSTATE[1] - 0.0 * MARS_X.p1_CSTATE[0];

  /* Sum: '<S125>/Sum' incorporates:
   *  Integrator: '<S6>/ub,vb,wb'
   *  Product: '<S128>/i x j'
   *  Product: '<S128>/j x k'
   *  Product: '<S128>/k x i'
   *  Product: '<S129>/i x k'
   *  Product: '<S129>/j x i'
   *  Product: '<S129>/k x j'
   */
  tmp_0[0] = MARS_X.ubvbwb_CSTATE[1] * rtb_Sum_hp[2];
  tmp_0[1] = rtb_Sum_hp[0] * MARS_X.ubvbwb_CSTATE[2];
  tmp_0[2] = MARS_X.ubvbwb_CSTATE[0] * rtb_Sum_hp[1];
  rtb_Sum2_0[0] = rtb_Sum_hp[1] * MARS_X.ubvbwb_CSTATE[2];
  rtb_Sum2_0[1] = MARS_X.ubvbwb_CSTATE[0] * rtb_Sum_hp[2];
  rtb_Sum2_0[2] = rtb_Sum_hp[0] * MARS_X.ubvbwb_CSTATE[1];

  /* Sum: '<S126>/Sum' incorporates:
   *  Constant: '<S1>/omega_earth2'
   *  Constant: '<S1>/omega_earth3'
   *  Product: '<S130>/i x j'
   *  Product: '<S130>/j x k'
   *  Product: '<S130>/k x i'
   *  Product: '<S131>/i x k'
   *  Product: '<S131>/j x i'
   *  Product: '<S131>/k x j'
   */
  phi_tmp = 0.0 * rtb_Sum2[2] - 7.292115E-5 * rtb_Sum2[1];
  rtb_q0dot = 7.292115E-5 * rtb_Sum2[0] - 0.0 * rtb_Sum2[2];
  uTmp_idx_0 = 0.0 * rtb_Sum2[1] - 0.0 * rtb_Sum2[0];
  for (i = 0; i < 3; i++) {
    /* Sum: '<S6>/Sum' incorporates:
     *  Concatenate: '<S114>/Vector Concatenate'
     *  Product: '<S6>/Product1'
     *  Sum: '<S125>/Sum'
     */
    uTmp_idx_1 = ((tmp_0[i] - rtb_Sum2_0[i]) - ((rtb_VectorConcatenate_f[i + 3] *
      rtb_q0dot + rtb_VectorConcatenate_f[i] * phi_tmp) +
      rtb_VectorConcatenate_f[i + 6] * uTmp_idx_0)) + MARS_ConstB.Product[i];

    /* DeadZone: '<S6>/Dead Zone' */
    if (uTmp_idx_1 > 2.2204460492503131E-16) {
      /* DeadZone: '<S6>/Dead Zone' */
      MARS_B.DeadZone[i] = uTmp_idx_1 - 2.2204460492503131E-16;
    } else if (uTmp_idx_1 >= -2.2204460492503131E-16) {
      /* DeadZone: '<S6>/Dead Zone' */
      MARS_B.DeadZone[i] = 0.0;
    } else {
      /* DeadZone: '<S6>/Dead Zone' */
      MARS_B.DeadZone[i] = uTmp_idx_1 - -2.2204460492503131E-16;
    }

    /* End of DeadZone: '<S6>/Dead Zone' */

    /* Sum: '<S6>/Sum' incorporates:
     *  Integrator: '<S1>/p,q,r '
     *  Product: '<S135>/Product'
     *  Product: '<S6>/Product1'
     *  Selector: '<S7>/Selector'
     */
    rtb_Sum_hp[i] = (MARS_ConstB.Selector[i + 3] * MARS_X.pqr_CSTATE[1] +
                     MARS_ConstB.Selector[i] * MARS_X.pqr_CSTATE[0]) +
      MARS_ConstB.Selector[i + 6] * MARS_X.pqr_CSTATE[2];
  }

  /* Sum: '<S134>/Sum' incorporates:
   *  Integrator: '<S1>/p,q,r '
   *  Product: '<S137>/i x j'
   *  Product: '<S137>/j x k'
   *  Product: '<S137>/k x i'
   *  Product: '<S138>/i x k'
   *  Product: '<S138>/j x i'
   *  Product: '<S138>/k x j'
   */
  tmp_0[0] = MARS_X.pqr_CSTATE[1] * rtb_Sum_hp[2];
  tmp_0[1] = rtb_Sum_hp[0] * MARS_X.pqr_CSTATE[2];
  tmp_0[2] = MARS_X.pqr_CSTATE[0] * rtb_Sum_hp[1];
  rtb_Sum2_0[0] = rtb_Sum_hp[1] * MARS_X.pqr_CSTATE[2];
  rtb_Sum2_0[1] = MARS_X.pqr_CSTATE[0] * rtb_Sum_hp[2];
  rtb_Sum2_0[2] = rtb_Sum_hp[0] * MARS_X.pqr_CSTATE[1];
  for (i = 0; i <= 0; i += 2) {
    /* Sum: '<S134>/Sum' */
    tmp_3 = _mm_loadu_pd(&tmp_0[i]);
    tmp_4 = _mm_loadu_pd(&rtb_Sum2_0[i]);
    _mm_storeu_pd(&rtb_Sum_hp[i], _mm_sub_pd(tmp_3, tmp_4));

    /* Product: '<S136>/Product' incorporates:
     *  Integrator: '<S1>/p,q,r '
     *  Sum: '<S134>/Sum'
     *  Sum: '<S7>/Sum2'
     */
    _mm_storeu_pd(&rtb_Sum2[i], _mm_set1_pd(0.0 * MARS_X.pqr_CSTATE[2] + (0.0 *
      MARS_X.pqr_CSTATE[1] + 0.0 * MARS_X.pqr_CSTATE[0])));
  }

  for (i = 2; i < 3; i++) {
    /* Sum: '<S134>/Sum' */
    rtb_Sum_hp[i] = tmp_0[i] - rtb_Sum2_0[i];

    /* Product: '<S136>/Product' incorporates:
     *  Integrator: '<S1>/p,q,r '
     *  Sum: '<S7>/Sum2'
     */
    rtb_Sum2[i] = (0.0 * MARS_X.pqr_CSTATE[1] + 0.0 * MARS_X.pqr_CSTATE[0]) +
      0.0 * MARS_X.pqr_CSTATE[2];
  }

  /* Sum: '<S7>/Sum2' incorporates:
   *  Constant: '<Root>/M_x'
   *  Constant: '<Root>/M_y'
   *  Constant: '<Root>/M_z'
   */
  tmp_0[0] = (0.0 - rtb_Sum2[0]) - rtb_Sum_hp[0];
  tmp_0[1] = (0.0 - rtb_Sum2[1]) - rtb_Sum_hp[1];
  tmp_0[2] = (0.01 - rtb_Sum2[2]) - rtb_Sum_hp[2];

  /* Product: '<S7>/Product2' incorporates:
   *  Selector: '<S7>/Selector2'
   */
  rt_mrdivide_U1d1x3_U2d_9vOrDY9Z(tmp_0, MARS_ConstB.Selector2, MARS_B.Product2);
  if (rtmIsMajorTimeStep((&MARS_M))) {
    /* Update for Integrator: '<S5>/p1' */
    MARS_DW.p1_IWORK = 0;

    /* Update for Integrator: '<S4>/q' */
    MARS_DW.q_IWORK = 0;

    /* Update for Integrator: '<S1>/p,q,r ' */
    MARS_DW.pqr_IWORK = 0;

    /* Update for Integrator: '<S5>/p' */
    MARS_DW.p_IWORK = 0;
  }                                    /* end MajorTimeStep */

  if (rtmIsMajorTimeStep((&MARS_M))) {
    rt_ertODEUpdateContinuousStates(&(&MARS_M)->solverInfo);

    /* Update absolute time for base rate */
    /* The "clockTick0" counts the number of times the code of this task has
     * been executed. The absolute time is the multiplication of "clockTick0"
     * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
     * overflow during the application lifespan selected.
     */
    ++(&MARS_M)->Timing.clockTick0;
    (&MARS_M)->Timing.t[0] = rtsiGetSolverStopTime(&(&MARS_M)->solverInfo);

    {
      /* Update absolute timer for sample time: [0.016666666666666666s, 0.0s] */
      /* The "clockTick1" counts the number of times the code of this task has
       * been executed. The resolution of this integer timer is 0.016666666666666666, which is the step size
       * of the task. Size of "clockTick1" ensures timer will not overflow during the
       * application lifespan selected.
       */
      (&MARS_M)->Timing.clockTick1++;
    }
  }                                    /* end MajorTimeStep */
}

/* Derivatives for root system: '<Root>' */
void MARS::MARS_derivatives()
{
  XDot_MARS_T *_rtXdot;
  _rtXdot = ((XDot_MARS_T *) (&MARS_M)->derivs);

  /* Derivatives for Integrator: '<S5>/p1' */
  _rtXdot->p1_CSTATE[0] = MARS_B.Product5[0];
  _rtXdot->p1_CSTATE[1] = MARS_B.Product5[1];
  _rtXdot->p1_CSTATE[2] = MARS_B.Product5[2];

  /* Derivatives for Integrator: '<S4>/q' */
  _rtXdot->q_CSTATE[0] = MARS_B.TmpSignalConversionAtqInport1[0];
  _rtXdot->q_CSTATE[1] = MARS_B.TmpSignalConversionAtqInport1[1];
  _rtXdot->q_CSTATE[2] = MARS_B.TmpSignalConversionAtqInport1[2];
  _rtXdot->q_CSTATE[3] = MARS_B.TmpSignalConversionAtqInport1[3];

  /* Derivatives for Integrator: '<S8>/Integrator' incorporates:
   *  Constant: '<S8>/Constant1'
   */
  _rtXdot->Integrator_CSTATE = 7.292115E-5;

  /* Derivatives for Integrator: '<S6>/ub,vb,wb' */
  _rtXdot->ubvbwb_CSTATE[0] = MARS_B.DeadZone[0];

  /* Derivatives for Integrator: '<S1>/p,q,r ' */
  _rtXdot->pqr_CSTATE[0] = MARS_B.Product2[0];

  /* Derivatives for Integrator: '<S5>/p' */
  _rtXdot->p_CSTATE[0] = MARS_B.Sum2_e[0];

  /* Derivatives for Integrator: '<S6>/ub,vb,wb' */
  _rtXdot->ubvbwb_CSTATE[1] = MARS_B.DeadZone[1];

  /* Derivatives for Integrator: '<S1>/p,q,r ' */
  _rtXdot->pqr_CSTATE[1] = MARS_B.Product2[1];

  /* Derivatives for Integrator: '<S5>/p' */
  _rtXdot->p_CSTATE[1] = MARS_B.Sum2_e[1];

  /* Derivatives for Integrator: '<S6>/ub,vb,wb' */
  _rtXdot->ubvbwb_CSTATE[2] = MARS_B.DeadZone[2];

  /* Derivatives for Integrator: '<S1>/p,q,r ' */
  _rtXdot->pqr_CSTATE[2] = MARS_B.Product2[2];

  /* Derivatives for Integrator: '<S5>/p' */
  _rtXdot->p_CSTATE[2] = MARS_B.Sum2_e[2];
}

/* Model initialize function */
void MARS::initialize()
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&(&MARS_M)->solverInfo, &(&MARS_M)->Timing.simTimeStep);
    rtsiSetTPtr(&(&MARS_M)->solverInfo, &rtmGetTPtr((&MARS_M)));
    rtsiSetStepSizePtr(&(&MARS_M)->solverInfo, &(&MARS_M)->Timing.stepSize0);
    rtsiSetdXPtr(&(&MARS_M)->solverInfo, &(&MARS_M)->derivs);
    rtsiSetContStatesPtr(&(&MARS_M)->solverInfo, (real_T **) &(&MARS_M)
                         ->contStates);
    rtsiSetNumContStatesPtr(&(&MARS_M)->solverInfo, &(&MARS_M)
      ->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&(&MARS_M)->solverInfo, &(&MARS_M)
      ->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&(&MARS_M)->solverInfo, &(&MARS_M)
      ->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&(&MARS_M)->solverInfo, &(&MARS_M)
      ->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&(&MARS_M)->solverInfo, (&rtmGetErrorStatus((&MARS_M))));
    rtsiSetRTModelPtr(&(&MARS_M)->solverInfo, (&MARS_M));
  }

  rtsiSetSimTimeStep(&(&MARS_M)->solverInfo, MAJOR_TIME_STEP);
  (&MARS_M)->intgData.y = (&MARS_M)->odeY;
  (&MARS_M)->intgData.f[0] = (&MARS_M)->odeF[0];
  (&MARS_M)->intgData.f[1] = (&MARS_M)->odeF[1];
  (&MARS_M)->intgData.f[2] = (&MARS_M)->odeF[2];
  (&MARS_M)->contStates = ((X_MARS_T *) &MARS_X);
  rtsiSetSolverData(&(&MARS_M)->solverInfo, static_cast<void *>(&(&MARS_M)
    ->intgData));
  rtsiSetIsMinorTimeStepWithModeChange(&(&MARS_M)->solverInfo, false);
  rtsiSetSolverName(&(&MARS_M)->solverInfo,"ode3");
  rtmSetTPtr((&MARS_M), &(&MARS_M)->Timing.tArray[0]);
  (&MARS_M)->Timing.stepSize0 = 0.016666666666666666;
  rtmSetFirstInitCond((&MARS_M), 1);

  /* Start for If: '<S17>/If' */
  MARS_DW.If_ActiveSubsystem = -1;

  /* Start for If: '<S16>/If' */
  MARS_DW.If_ActiveSubsystem_h = -1;

  /* Start for If: '<S46>/If1' */
  MARS_DW.If1_ActiveSubsystem = -1;

  /* Start for If: '<S69>/If1' */
  MARS_DW.If1_ActiveSubsystem_l = -1;

  /* InitializeConditions for Integrator: '<S5>/p1' incorporates:
   *  Integrator: '<S4>/q'
   */
  if (rtmIsFirstInitCond((&MARS_M))) {
    MARS_X.p1_CSTATE[0] = 0.0;
    MARS_X.p1_CSTATE[1] = 0.0;
    MARS_X.p1_CSTATE[2] = 0.0;
    MARS_X.q_CSTATE[0] = 0.0;
    MARS_X.q_CSTATE[1] = 0.0;
    MARS_X.q_CSTATE[2] = 0.0;
    MARS_X.q_CSTATE[3] = 0.0;
  }

  MARS_DW.p1_IWORK = 1;

  /* End of InitializeConditions for Integrator: '<S5>/p1' */

  /* InitializeConditions for Integrator: '<S4>/q' */
  MARS_DW.q_IWORK = 1;

  /* InitializeConditions for Integrator: '<S8>/Integrator' */
  MARS_X.Integrator_CSTATE = 0.0;

  /* InitializeConditions for Integrator: '<S6>/ub,vb,wb' */
  MARS_X.ubvbwb_CSTATE[0] = 0.0;
  MARS_X.ubvbwb_CSTATE[1] = 0.0;
  MARS_X.ubvbwb_CSTATE[2] = 0.0;

  /* InitializeConditions for Integrator: '<S1>/p,q,r ' incorporates:
   *  Integrator: '<S5>/p'
   */
  if (rtmIsFirstInitCond((&MARS_M))) {
    MARS_X.pqr_CSTATE[0] = 0.0;
    MARS_X.pqr_CSTATE[1] = 0.0;
    MARS_X.pqr_CSTATE[2] = 0.0;
    MARS_X.p_CSTATE[0] = 0.0;
    MARS_X.p_CSTATE[1] = 0.0;
    MARS_X.p_CSTATE[2] = 0.0;
  }

  MARS_DW.pqr_IWORK = 1;

  /* End of InitializeConditions for Integrator: '<S1>/p,q,r ' */

  /* InitializeConditions for Integrator: '<S5>/p' */
  MARS_DW.p_IWORK = 1;

  /* SystemInitialize for Merge: '<S17>/Merge' */
  MARS_B.Merge[0] = 1.0;
  MARS_B.Merge[1] = 0.0;
  MARS_B.Merge[2] = 0.0;
  MARS_B.Merge[3] = 0.0;

  /* SystemInitialize for IfAction SubSystem: '<S16>/AxisRotZeroR3' */
  /* Start for If: '<S51>/If' */
  MARS_DW.If_ActiveSubsystem_c0 = -1;

  /* End of SystemInitialize for SubSystem: '<S16>/AxisRotZeroR3' */

  /* SystemInitialize for IfAction SubSystem: '<S16>/AxisRotDefault' */
  /* Start for If: '<S47>/If' */
  MARS_DW.If_ActiveSubsystem_c = -1;

  /* End of SystemInitialize for SubSystem: '<S16>/AxisRotDefault' */

  /* set "at time zero" to false */
  if (rtmIsFirstInitCond((&MARS_M))) {
    rtmSetFirstInitCond((&MARS_M), 0);
  }
}

/* Model terminate function */
void MARS::terminate()
{
  /* (no terminate code required) */
}

/* Constructor */
MARS::MARS() :
  MARS_B(),
  MARS_DW(),
  MARS_X(),
  MARS_M()
{
  /* Currently there is no constructor body generated.*/
}

/* Destructor */
MARS::~MARS()
{
  /* Currently there is no destructor body generated.*/
}

/* Real-Time Model get method */
RT_MODEL_MARS_T * MARS::getRTM()
{
  return (&MARS_M);
}
