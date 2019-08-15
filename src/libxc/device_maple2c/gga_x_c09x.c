/* 
  This file was generated automatically with /nfs/data-012/marques/software/source/libxc/svn/scripts/maple2c.pl.
  Do not edit this file directly as it can be overwritten!!

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.

  Maple version     : Maple 2016 (X86 64 LINUX)
  Maple source      : ../maple/gga_x_c09x.mpl
  Type of functional: work_gga_x
*/

#ifdef CUDA
__device__ void xc_gga_x_c09x_enhance
  (const void *p,  xc_gga_work_x_t *r)
#else
void xc_gga_x_c09x_enhance
  (const xc_func_type *p,  xc_gga_work_x_t *r)
#endif
{
  double t1, t2, t3, t4, t5, t6, t7, t8;
  double t10, t11, t15, t17, t20, t23, t24, t25;
  double t28, t35;


  t1 = M_CBRT6;
  t2 = 0.31415926535897932385e1 * 0.31415926535897932385e1;
  t3 = cbrt(t2);
  t4 = t3 * t3;
  t5 = 0.1e1 / t4;
  t6 = t1 * t5;
  t7 = r->x * r->x;
  t8 = t6 * t7;
  t10 = exp(-0.20125000000000000000e-2 * t8);
  t11 = t7 * t10;
  t15 = exp(-0.10062500000000000000e-2 * t8);
  r->f = 0.2245e1 + 0.25708333333333333333e-2 * t6 * t11 - 0.1245e1 * t15;

  if(r->order < 1) return;

  t17 = r->x * t10;
  t20 = t1 * t1;
  t23 = t20 / t3 / t2;
  t24 = t7 * r->x;
  t25 = t24 * t10;
  t28 = r->x * t15;
  r->dfdx = 0.51416666666666666666e-2 * t6 * t17 - 0.10347604166666666667e-4 * t23 * t25 + 0.25055625000000000000e-2 * t6 * t28;

  if(r->order < 2) return;

  t35 = t7 * t7;
  r->d2fdx2 = 0.51416666666666666666e-2 * t6 * t10 - 0.51738020833333333334e-4 * t23 * t11 + 0.25654139461969691649e-8 * t35 * t10 + 0.25055625000000000000e-2 * t6 * t15 - 0.50424445312500000000e-5 * t23 * t7 * t15;

  if(r->order < 3) return;

  r->d3fdx3 = -0.12417125000000000000e-3 * t23 * t17 + 0.23088725515772722484e-7 * t25 - 0.10325791133442800889e-10 * t35 * r->x * t1 * t5 * t10 - 0.15127333593750000000e-4 * t23 * t28 + 0.62507017639236404079e-9 * t24 * t15;

  if(r->order < 4) return;


}

#ifndef CUDA
#define maple2c_order 3
#define maple2c_func  xc_gga_x_c09x_enhance
#endif