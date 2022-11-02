/*
 * MARS_data.c
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

/* Invariant block signals (default storage) */
const ConstB_MARS_T MARS_ConstB = {
  { -0.82949238588247243, 0.38438233146596967, -0.4052068669455009,
    0.34537608216871168, 0.92317399403188638, 0.1687161480386685,
    0.43892794810628633, 0.0, -0.89852226259075252 },/* '<S61>/Product2' */

  { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 },/* '<S64>/Selector' */

  { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 },/* '<S64>/Selector2' */

  { 1.1102230246251565E-16, 1.3877787807814457E-17, 0.0, 1.3877787807814457E-17,
    1.1102230246251565E-16, 0.0, 0.0, 0.0, 0.0 },/* '<S153>/Abs2' */

  { 0.95540418704665353, 0.044147846096478795, -0.22088421489473919 },/* '<S128>/Product' */

  { 0.19095541539610336, 0.95540418704665309, -0.22088421489473903 },/* '<S130>/Product' */

  { 0.19095541539610345, 0.044147846096478795, -0.22088421489473928 },/* '<S129>/Product' */
  -0.84413481505178722,                /* '<S147>/Add' */
  -0.1687161480386685,                 /* '<S148>/Add' */
  -0.039006249297257989                /* '<S149>/Add' */
};

/* Constant parameters (default storage) */
const ConstP_MARS_T MARS_ConstP = {
  /* Pooled Parameter (Expression: -eye(3))
   * Referenced by:
   *   '<S115>/Bias1'
   *   '<S153>/Bias1'
   */
  { -1.0, -0.0, -0.0, -0.0, -1.0, -0.0, -0.0, -0.0, -1.0 }
};
