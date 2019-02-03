#ifndef _TECURVECOSTS_INCLUDE_
#define _TECURVECOSTS_INCLUDE_

#include "tecurve.h"

/*
    Copyright 2012 Andrew V. Sutherland

    This file is part of classpoly.

    classpoly is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    classpoly is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with classpoly.  If not, see <http://www.gnu.org/licenses/>.
*/

static double tecurve_fix2_costs[4*TECURVE_MAX_N][2*TECURVE_FIX_K] = {
/* N=1 (2 mod 3) (3 mod 4) */ { 0.84, 0.00, 0.00, 0.00, 0.00, 0.00, 1.41, 0.00, 0.00, 0.00, 0.00, 0.00},
/* N=1 (2 mod 3) (1 mod 4) */ { 1.02, 0.00, 0.00, 0.00, 0.00, 0.00, 1.77, 0.00, 0.00, 0.00, 0.00, 0.00},
/* N=1 (1 mod 3) (3 mod 4) */ { 0.84, 0.00, 0.00, 0.00, 0.00, 0.00, 1.82, 0.00, 0.00, 0.00, 0.00, 0.00},
/* N=1 (1 mod 3) (1 mod 4) */ { 1.02, 0.00, 0.00, 0.00, 0.00, 0.00, 2.09, 0.00, 0.00, 0.00, 0.00, 0.00},
/* N=2 (2 mod 3) (3 mod 4) */ { 0.00, 0.93, 1.57, 2.75, 5.18, 9.94, 0.00, 1.57, 2.75, 5.18, 10.00, 19.64},
/* N=2 (2 mod 3) (1 mod 4) */ { 0.00, 0.93, 1.57, 2.82, 5.27, 10.21, 0.00, 1.59, 2.75, 5.25, 10.18, 19.96},
/* N=2 (1 mod 3) (3 mod 4) */ { 0.00, 0.95, 1.57, 2.84, 5.30, 10.25, 0.00, 1.98, 2.93, 4.93, 8.91, 16.73},
/* N=2 (1 mod 3) (1 mod 4) */ { 0.00, 0.93, 1.57, 2.77, 5.25, 10.14, 0.00, 1.95, 2.91, 4.86, 8.73, 16.48},
/* N=3 (2 mod 3) (3 mod 4) */ { 0.95, 2.11, 3.82, 7.32, 14.21, 28.01, 0.95, 2.11, 3.84, 7.32, 14.25, 28.19},
/* N=3 (2 mod 3) (1 mod 4) */ { 1.18, 2.14, 3.91, 7.46, 14.64, 28.83, 1.18, 2.14, 3.91, 7.46, 14.62, 28.76},
/* N=3 (1 mod 3) (3 mod 4) */ { 0.93, 2.09, 3.77, 7.23, 14.10, 27.85, 0.95, 2.07, 3.80, 7.21, 14.07, 27.69},
/* N=3 (1 mod 3) (1 mod 4) */ { 1.20, 2.05, 3.77, 7.14, 13.87, 27.50, 1.20, 2.07, 3.75, 7.12, 13.94, 27.62},
/* N=4 (2 mod 3) (3 mod 4) */ { 0.00, 0.00, 1.43, 2.50, 4.62, 8.93, 0.00, 0.00, 2.50, 4.64, 8.91, 17.42},
/* N=4 (2 mod 3) (1 mod 4) */ { 0.00, 0.00, 1.43, 2.50, 4.71, 9.05, 0.00, 0.00, 2.48, 4.68, 9.00, 17.76},
/* N=4 (1 mod 3) (3 mod 4) */ { 0.00, 0.00, 1.43, 2.55, 4.75, 9.09, 0.00, 0.00, 2.71, 4.46, 7.93, 14.94},
/* N=4 (1 mod 3) (1 mod 4) */ { 0.00, 0.00, 1.41, 2.48, 4.64, 8.87, 0.00, 0.00, 2.66, 4.39, 7.80, 14.60},
/* N=5 (2 mod 3) (3 mod 4) */ { 0.95, 2.16, 3.96, 7.52, 14.69, 28.97, 1.64, 3.96, 7.48, 14.67, 28.89, 56.86},
/* N=5 (2 mod 3) (1 mod 4) */ { 1.23, 2.18, 4.02, 7.68, 15.05, 29.73, 2.11, 4.05, 7.62, 14.96, 29.74, 58.69},
/* N=5 (1 mod 3) (3 mod 4) */ { 0.98, 2.11, 3.86, 7.43, 14.46, 28.44, 2.05, 3.84, 6.64, 12.32, 23.55, 46.04},
/* N=5 (1 mod 3) (1 mod 4) */ { 1.25, 2.11, 3.86, 7.37, 14.32, 28.32, 2.41, 3.80, 6.59, 12.23, 23.51, 45.44},
/* N=6 (2 mod 3) (3 mod 4) */ { 0.00, 1.41, 2.50, 4.59, 8.80, 17.10, 0.00, 1.43, 2.48, 4.59, 8.82, 17.23},
/* N=6 (2 mod 3) (1 mod 4) */ { 0.00, 1.46, 2.50, 4.68, 9.05, 17.64, 0.00, 1.43, 2.52, 4.68, 9.00, 17.62},
/* N=6 (1 mod 3) (3 mod 4) */ { 0.00, 1.46, 2.50, 4.71, 9.03, 17.69, 0.00, 1.43, 2.52, 4.68, 9.07, 17.78},
/* N=6 (1 mod 3) (1 mod 4) */ { 0.00, 1.41, 2.50, 4.64, 8.96, 17.39, 0.00, 1.43, 2.50, 4.64, 8.91, 17.46},
/* N=7 (2 mod 3) (3 mod 4) */ { 1.00, 2.23, 4.09, 7.78, 15.16, 30.01, 1.70, 4.09, 7.80, 15.25, 29.97, 59.51},
/* N=7 (2 mod 3) (1 mod 4) */ { 1.30, 2.27, 4.14, 7.96, 15.55, 30.76, 2.23, 4.16, 7.91, 15.55, 30.63, 61.26},
/* N=7 (1 mod 3) (3 mod 4) */ { 1.02, 2.18, 4.00, 7.71, 15.09, 29.76, 2.07, 3.93, 6.89, 12.73, 24.44, 48.09},
/* N=7 (1 mod 3) (1 mod 4) */ { 1.32, 2.18, 3.98, 7.64, 14.87, 29.33, 2.50, 3.93, 6.82, 12.66, 24.33, 47.73},
/* N=8 (2 mod 3) (3 mod 4) */ { 0.00, 0.00, 0.00, 1.57, 2.77, 5.18, 0.00, 0.00, 0.00, 2.77, 5.16, 9.96},
/* N=8 (2 mod 3) (1 mod 4) */ { 0.00, 0.00, 0.00, 1.55, 2.80, 5.23, 0.00, 0.00, 0.00, 2.75, 5.18, 10.07},
/* N=8 (1 mod 3) (3 mod 4) */ { 0.00, 0.00, 0.00, 1.57, 2.82, 5.32, 0.00, 0.00, 0.00, 2.93, 4.93, 8.87},
/* N=8 (1 mod 3) (1 mod 4) */ { 0.00, 0.00, 0.00, 1.55, 2.75, 5.18, 0.00, 0.00, 0.00, 2.89, 4.84, 8.66},
/* N=9 (2 mod 3) (3 mod 4) */ { 1.02, 2.23, 4.07, 7.86, 15.32, 30.04, 1.00, 2.25, 4.07, 7.82, 15.30, 30.28},
/* N=9 (2 mod 3) (1 mod 4) */ { 1.32, 2.27, 4.14, 8.00, 15.71, 30.94, 1.32, 2.27, 4.14, 7.98, 15.62, 30.92},
/* N=9 (1 mod 3) (3 mod 4) */ { 1.00, 2.20, 4.05, 7.73, 15.05, 29.78, 1.02, 2.20, 4.05, 7.73, 15.14, 29.87},
/* N=9 (1 mod 3) (1 mod 4) */ { 1.32, 2.20, 4.02, 7.68, 15.05, 29.53, 1.34, 2.18, 4.02, 7.68, 14.96, 29.55},
/* N=10 (2 mod 3) (3 mod 4) */ { 0.00, 1.55, 2.71, 5.07, 9.75, 19.05, 0.00, 2.73, 5.02, 9.73, 19.05, 37.85},
/* N=10 (2 mod 3) (1 mod 4) */ { 0.00, 1.55, 2.75, 5.14, 9.91, 19.50, 0.00, 2.73, 5.12, 9.89, 19.39, 38.65},
/* N=10 (1 mod 3) (3 mod 4) */ { 0.00, 1.57, 2.75, 5.18, 9.98, 19.51, 0.00, 2.93, 4.80, 8.68, 16.30, 31.74},
/* N=10 (1 mod 3) (1 mod 4) */ { 0.00, 1.55, 2.73, 5.09, 9.84, 19.28, 0.00, 2.89, 4.77, 8.57, 16.14, 31.35},
/* N=11 (2 mod 3) (3 mod 4) */ { 1.18, 2.39, 4.21, 7.89, 15.23, 29.90, 1.86, 4.25, 7.84, 15.21, 29.63, 58.73},
/* N=11 (2 mod 3) (1 mod 4) */ { 1.45, 2.43, 4.27, 8.02, 15.51, 30.47, 2.39, 4.27, 7.98, 15.55, 30.42, 60.59},
/* N=11 (1 mod 3) (3 mod 4) */ { 1.18, 2.36, 4.16, 7.78, 14.96, 29.34, 2.25, 4.11, 6.98, 12.73, 24.35, 47.38},
/* N=11 (1 mod 3) (1 mod 4) */ { 1.48, 2.34, 4.14, 7.73, 14.94, 29.14, 2.66, 4.07, 6.93, 12.69, 24.14, 46.95},
/* N=12 (2 mod 3) (3 mod 4) */ { 0.00, 0.00, 1.55, 2.75, 5.09, 9.84, 0.00, 0.00, 1.55, 2.73, 5.12, 9.87},
/* N=12 (2 mod 3) (1 mod 4) */ { 0.00, 0.00, 1.52, 2.75, 5.14, 9.94, 0.00, 0.00, 1.55, 2.73, 5.18, 9.96},
/* N=12 (1 mod 3) (3 mod 4) */ { 0.00, 0.00, 1.55, 2.77, 5.23, 10.07, 0.00, 0.00, 1.55, 2.77, 5.23, 10.03},
/* N=12 (1 mod 3) (1 mod 4) */ { 0.00, 0.00, 1.52, 2.73, 5.14, 9.87, 0.00, 0.00, 1.55, 2.73, 5.12, 9.89},
/* N=13 (2 mod 3) (3 mod 4) */ { 1.52, 3.07, 5.62, 10.59, 20.55, 40.44, 2.55, 5.59, 10.62, 20.49, 40.53, 80.19},
/* N=13 (2 mod 3) (1 mod 4) */ { 2.18, 3.11, 5.71, 10.82, 21.10, 41.67, 3.80, 5.71, 10.82, 20.96, 41.80, 82.92},
/* N=13 (1 mod 3) (3 mod 4) */ { 1.55, 3.07, 5.57, 10.68, 20.69, 40.76, 2.82, 5.23, 9.30, 17.35, 33.37, 65.23},
/* N=13 (1 mod 3) (1 mod 4) */ { 2.18, 3.05, 5.50, 10.53, 20.48, 40.23, 3.77, 5.18, 9.14, 17.09, 33.03, 65.01},
/* N=14 (2 mod 3) (3 mod 4) */ { 0.00, 1.73, 2.93, 5.27, 10.02, 19.37, 0.00, 2.93, 5.25, 10.05, 19.55, 37.96},
/* N=14 (2 mod 3) (1 mod 4) */ { 0.00, 1.77, 2.96, 5.39, 10.18, 19.87, 0.00, 2.96, 5.34, 10.18, 19.75, 39.17},
/* N=14 (1 mod 3) (3 mod 4) */ { 0.00, 1.75, 2.96, 5.37, 10.23, 19.98, 0.00, 3.12, 5.05, 8.91, 16.60, 32.21},
/* N=14 (1 mod 3) (1 mod 4) */ { 0.00, 1.75, 2.93, 5.36, 10.16, 19.66, 0.00, 3.09, 4.98, 8.84, 16.53, 31.70},
/* N=15 (2 mod 3) (3 mod 4) */ { 1.20, 2.46, 4.36, 8.18, 15.84, 31.01, 1.23, 2.46, 4.36, 8.12, 15.80, 30.90},
/* N=15 (2 mod 3) (1 mod 4) */ { 1.52, 2.52, 4.41, 8.34, 16.05, 31.85, 1.52, 2.50, 4.43, 8.37, 16.10, 31.76},
/* N=15 (1 mod 3) (3 mod 4) */ { 1.23, 2.41, 4.30, 8.09, 15.57, 30.76, 1.23, 2.45, 4.30, 8.07, 15.66, 30.62},
/* N=15 (1 mod 3) (1 mod 4) */ { 1.55, 2.43, 4.27, 8.00, 15.48, 30.33, 1.55, 2.43, 4.27, 8.00, 15.42, 30.39},
/* N=16 (2 mod 3) (3 mod 4) */ { 0.00, 0.00, 0.00, 0.00, 2.45, 4.34, 0.00, 0.00, 0.00, 0.00, 4.34, 8.16},
/* N=16 (2 mod 3) (1 mod 4) */ { 0.00, 0.00, 0.00, 0.00, 2.46, 4.43, 0.00, 0.00, 0.00, 0.00, 4.36, 8.25},
/* N=16 (1 mod 3) (3 mod 4) */ { 0.00, 0.00, 0.00, 0.00, 2.50, 4.43, 0.00, 0.00, 0.00, 0.00, 4.30, 7.41},
/* N=16 (1 mod 3) (1 mod 4) */ { 0.00, 0.00, 0.00, 0.00, 2.45, 4.39, 0.00, 0.00, 0.00, 0.00, 4.21, 7.30},
/* N=17 (2 mod 3) (3 mod 4) */ { 3.07, 6.05, 11.53, 22.53, 44.52, 88.34, 5.55, 11.59, 22.44, 44.60, 88.65, 177.25},
/* N=17 (2 mod 3) (1 mod 4) */ { 5.25, 6.23, 11.86, 23.16, 45.72, 90.41, 9.96, 11.89, 23.14, 45.43, 90.96, 181.67},
/* N=17 (1 mod 3) (3 mod 4) */ { 3.02, 5.98, 11.41, 22.40, 44.08, 87.60, 5.18, 10.00, 18.73, 35.96, 71.00, 140.58},
/* N=17 (1 mod 3) (1 mod 4) */ { 5.11, 6.00, 11.30, 22.21, 43.67, 87.37, 8.46, 9.86, 18.44, 35.70, 70.14, 138.54},
/* N=18 (2 mod 3) (3 mod 4) */ { 0.00, 2.43, 4.34, 8.09, 15.53, 30.67, 0.00, 2.46, 4.32, 8.09, 15.64, 30.56},
/* N=18 (2 mod 3) (1 mod 4) */ { 0.00, 2.50, 4.41, 8.27, 16.07, 31.35, 0.00, 2.48, 4.41, 8.30, 16.03, 31.26},
/* N=18 (1 mod 3) (3 mod 4) */ { 0.00, 2.52, 4.41, 8.32, 16.05, 31.49, 0.00, 2.48, 4.43, 8.34, 16.08, 31.66},
/* N=18 (1 mod 3) (1 mod 4) */ { 0.00, 2.48, 4.39, 8.23, 15.89, 31.17, 0.00, 2.48, 4.36, 8.21, 15.80, 31.14},
/* N=19 (2 mod 3) (3 mod 4) */ { 3.91, 7.80, 15.01, 29.53, 58.59, 115.62, 7.27, 15.12, 29.46, 58.02, 116.14, 230.75},
/* N=19 (2 mod 3) (1 mod 4) */ { 6.93, 7.86, 15.14, 29.69, 58.81, 117.51, 13.28, 15.16, 29.69, 59.01, 116.80, 234.99},
/* N=19 (1 mod 3) (3 mod 4) */ { 3.93, 7.80, 15.03, 29.42, 58.38, 116.14, 6.61, 12.78, 24.37, 47.29, 93.18, 185.18},
/* N=19 (1 mod 3) (1 mod 4) */ { 6.93, 7.82, 14.96, 29.44, 58.11, 116.14, 11.30, 12.78, 24.26, 47.51, 94.07, 185.54},
/* N=20 (2 mod 3) (3 mod 4) */ { 0.00, 0.00, 3.50, 6.41, 12.23, 23.89, 0.00, 0.00, 6.36, 12.25, 23.94, 47.38},
/* N=20 (2 mod 3) (1 mod 4) */ { 0.00, 0.00, 3.52, 6.55, 12.57, 24.58, 0.00, 0.00, 6.48, 12.48, 24.42, 48.13},
/* N=20 (1 mod 3) (3 mod 4) */ { 0.00, 0.00, 3.41, 6.32, 12.14, 23.78, 0.00, 0.00, 5.80, 10.48, 19.75, 38.19},
/* N=20 (1 mod 3) (1 mod 4) */ { 0.00, 0.00, 3.43, 6.27, 12.03, 23.46, 0.00, 0.00, 5.75, 10.32, 19.53, 37.97},
/* N=21 (2 mod 3) (3 mod 4) */ { 3.05, 6.09, 11.62, 22.65, 44.72, 88.56, 3.02, 6.09, 11.55, 22.60, 44.44, 88.70},
/* N=21 (2 mod 3) (1 mod 4) */ { 5.23, 6.23, 11.85, 23.05, 45.79, 90.78, 5.25, 6.21, 11.89, 23.10, 45.74, 90.44},
/* N=21 (1 mod 3) (3 mod 4) */ { 3.05, 6.03, 11.50, 22.51, 44.01, 87.47, 3.04, 6.02, 11.53, 22.37, 44.31, 87.52},
/* N=21 (1 mod 3) (1 mod 4) */ { 5.09, 6.00, 11.37, 22.16, 43.83, 87.14, 5.09, 5.98, 11.37, 22.25, 43.56, 87.12},
/* N=22 (2 mod 3) (3 mod 4) */ { 0.00, 6.84, 13.16, 25.46, 50.27, 100.27, 0.00, 13.10, 25.69, 50.61, 100.32, 200.27},
/* N=22 (2 mod 3) (1 mod 4) */ { 0.00, 6.95, 13.44, 26.15, 52.33, 103.21, 0.00, 13.39, 26.01, 51.93, 102.68, 203.91},
/* N=22 (1 mod 3) (3 mod 4) */ { 0.00, 7.11, 13.59, 26.64, 52.84, 105.27, 0.00, 11.71, 22.05, 42.85, 84.97, 168.10},
/* N=22 (1 mod 3) (1 mod 4) */ { 0.00, 7.07, 13.53, 26.44, 52.29, 104.28, 0.00, 11.61, 22.05, 42.63, 84.80, 167.25},
/* N=23 (2 mod 3) (3 mod 4) */ { 5.57, 11.16, 21.60, 42.43, 83.94, 169.17, 10.61, 21.51, 42.36, 83.85, 168.63, 335.98},
/* N=23 (2 mod 3) (1 mod 4) */ { 10.18, 11.18, 21.60, 42.81, 85.07, 169.64, 19.75, 21.78, 42.94, 84.73, 169.03, 338.82},
/* N=23 (1 mod 3) (3 mod 4) */ { 5.59, 11.05, 21.58, 42.61, 84.67, 168.13, 9.23, 18.07, 35.01, 68.46, 135.06, 270.35},
/* N=23 (1 mod 3) (1 mod 4) */ { 10.25, 11.12, 21.57, 42.59, 84.24, 168.45, 16.73, 18.05, 34.78, 68.27, 135.54, 271.37},
/* N=24 (2 mod 3) (3 mod 4) */ { 0.00, 0.00, 0.00, 7.23, 13.87, 27.15, 0.00, 0.00, 0.00, 7.23, 13.87, 27.37},
/* N=24 (2 mod 3) (1 mod 4) */ { 0.00, 0.00, 0.00, 7.27, 13.98, 27.39, 0.00, 0.00, 0.00, 7.28, 13.98, 27.39},
/* N=24 (1 mod 3) (3 mod 4) */ { 0.00, 0.00, 0.00, 7.27, 14.03, 27.37, 0.00, 0.00, 0.00, 7.30, 14.03, 27.51},
/* N=24 (1 mod 3) (1 mod 4) */ { 0.00, 0.00, 0.00, 7.25, 13.98, 27.37, 0.00, 0.00, 0.00, 7.25, 13.98, 27.37},
/* N=25 (2 mod 3) (3 mod 4) */ { 6.48, 13.28, 25.49, 50.40, 101.83, 203.65, 12.77, 26.01, 51.87, 102.80, 202.80, 410.80},
/* N=25 (2 mod 3) (1 mod 4) */ { 12.25, 12.98, 26.01, 51.23, 102.76, 208.68, 24.01, 25.58, 51.53, 103.71, 200.10, 399.81},
/* N=25 (1 mod 3) (3 mod 4) */ { 6.73, 13.24, 25.54, 51.61, 102.60, 200.59, 11.01, 21.54, 41.12, 82.65, 163.59, 322.58},
/* N=25 (1 mod 3) (1 mod 4) */ { 12.51, 13.48, 25.80, 51.42, 101.31, 200.45, 20.04, 21.51, 41.34, 81.01, 162.49, 328.02},
/* N=26 (2 mod 3) (3 mod 4) */ { 0.00, 9.02, 17.26, 34.26, 68.79, 136.82, 0.00, 17.51, 34.01, 67.21, 134.20, 268.83},
/* N=26 (2 mod 3) (1 mod 4) */ { 0.00, 9.02, 17.30, 33.76, 68.30, 134.93, 0.00, 17.30, 34.23, 66.41, 132.28, 269.36},
/* N=26 (1 mod 3) (3 mod 4) */ { 0.00, 8.75, 17.49, 33.32, 67.05, 133.03, 0.00, 14.52, 28.00, 54.51, 107.80, 217.68},
/* N=26 (1 mod 3) (1 mod 4) */ { 0.00, 8.75, 17.27, 33.77, 66.08, 137.17, 0.00, 14.54, 27.25, 53.06, 104.66, 214.32},
};

static double tecurve_free2_costs[4*TECURVE_MAX_N][2*TECURVE_FREE_K] = {
/* N=1 (2 mod 3) (3 mod 4) */ { 0.77 , 0.00, 0.00, 1.25, 0.00, 0.00},
/* N=1 (2 mod 3) (1 mod 4) */ { 0.75 , 0.00, 0.00, 1.25, 0.00, 0.00},
/* N=1 (1 mod 3) (3 mod 4) */ { 0.75 , 0.00, 0.00, 1.68, 0.00, 0.00},
/* N=1 (1 mod 3) (1 mod 4) */ { 0.75 , 0.00, 0.00, 1.68, 0.00, 0.00},
/* N=2 (2 mod 3) (3 mod 4) */ { 0.70 , 0.70, 0.95, 1.16, 1.16, 1.59},
/* N=2 (2 mod 3) (1 mod 4) */ { 0.73 , 0.68, 0.95, 1.16, 1.16, 1.59},
/* N=2 (1 mod 3) (3 mod 4) */ { 0.73 , 0.70, 0.93, 1.64, 1.61, 2.00},
/* N=2 (1 mod 3) (1 mod 4) */ { 0.70 , 0.73, 0.93, 1.61, 1.61, 1.98},
/* N=3 (2 mod 3) (3 mod 4) */ { 0.70 , 0.93, 1.57, 0.73, 0.93, 1.55},
/* N=3 (2 mod 3) (1 mod 4) */ { 0.70 , 0.95, 1.57, 0.73, 0.95, 1.57},
/* N=3 (1 mod 3) (3 mod 4) */ { 0.73 , 0.95, 1.50, 0.68, 0.95, 1.55},
/* N=3 (1 mod 3) (1 mod 4) */ { 0.70 , 0.95, 1.52, 0.70, 0.93, 1.52},
/* N=4 (2 mod 3) (3 mod 4) */ { 0.70 , 0.73, 0.70, 1.16, 1.16, 1.18},
/* N=4 (2 mod 3) (1 mod 4) */ { 0.73 , 0.70, 0.73, 1.16, 1.16, 1.14},
/* N=4 (1 mod 3) (3 mod 4) */ { 0.73 , 0.70, 0.73, 1.61, 1.64, 1.61},
/* N=4 (1 mod 3) (1 mod 4) */ { 0.70 , 0.73, 0.73, 1.61, 1.61, 1.61},
/* N=5 (2 mod 3) (3 mod 4) */ { 0.73 , 0.95, 1.59, 1.18, 1.64, 2.84},
/* N=5 (2 mod 3) (1 mod 4) */ { 0.73 , 0.98, 1.59, 1.18, 1.61, 2.89},
/* N=5 (1 mod 3) (3 mod 4) */ { 0.70 , 0.95, 1.57, 1.64, 2.05, 2.93},
/* N=5 (1 mod 3) (1 mod 4) */ { 0.70 , 0.98, 1.55, 1.64, 2.00, 2.91},
/* N=6 (2 mod 3) (3 mod 4) */ { 0.73 , 0.73, 0.93, 0.73, 0.73, 0.93},
/* N=6 (2 mod 3) (1 mod 4) */ { 0.73 , 0.73, 0.93, 0.73, 0.73, 0.93},
/* N=6 (1 mod 3) (3 mod 4) */ { 0.73 , 0.73, 0.95, 0.73, 0.70, 0.95},
/* N=6 (1 mod 3) (1 mod 4) */ { 0.73 , 0.70, 0.95, 0.73, 0.73, 0.93},
/* N=7 (2 mod 3) (3 mod 4) */ { 0.73 , 1.02, 1.61, 1.20, 1.71, 2.91},
/* N=7 (2 mod 3) (1 mod 4) */ { 0.75 , 1.00, 1.66, 1.20, 1.70, 2.96},
/* N=7 (1 mod 3) (3 mod 4) */ { 0.75 , 1.02, 1.59, 1.66, 2.07, 3.05},
/* N=7 (1 mod 3) (1 mod 4) */ { 0.75 , 1.00, 1.59, 1.64, 2.05, 2.98},
/* N=8 (2 mod 3) (3 mod 4) */ { 0.75 , 0.75, 0.75, 1.23, 1.23, 1.25},
/* N=8 (2 mod 3) (1 mod 4) */ { 0.75 , 0.75, 0.75, 1.23, 1.23, 1.25},
/* N=8 (1 mod 3) (3 mod 4) */ { 0.75 , 0.73, 0.77, 1.68, 1.68, 1.68},
/* N=8 (1 mod 3) (1 mod 4) */ { 0.73 , 0.77, 0.75, 1.66, 1.66, 1.68},
/* N=9 (2 mod 3) (3 mod 4) */ { 0.75 , 1.00, 1.64, 0.75, 1.00, 1.64},
/* N=9 (2 mod 3) (1 mod 4) */ { 0.75 , 1.00, 1.66, 0.75, 1.00, 1.66},
/* N=9 (1 mod 3) (3 mod 4) */ { 0.75 , 1.00, 1.64, 0.73, 1.02, 1.61},
/* N=9 (1 mod 3) (1 mod 4) */ { 0.73 , 1.02, 1.59, 0.73, 1.02, 1.59},
/* N=10 (2 mod 3) (3 mod 4) */ { 0.75 , 0.75, 0.98, 1.25, 1.23, 1.66},
/* N=10 (2 mod 3) (1 mod 4) */ { 0.75 , 0.75, 0.98, 1.23, 1.23, 1.68},
/* N=10 (1 mod 3) (3 mod 4) */ { 0.75 , 0.75, 0.98, 1.68, 1.68, 2.05},
/* N=10 (1 mod 3) (1 mod 4) */ { 0.75 , 0.75, 0.98, 1.68, 1.66, 2.02},
/* N=11 (2 mod 3) (3 mod 4) */ { 0.91 , 1.18, 1.80, 1.41, 1.86, 3.07},
/* N=11 (2 mod 3) (1 mod 4) */ { 0.93 , 1.18, 1.82, 1.41, 1.86, 3.11},
/* N=11 (1 mod 3) (3 mod 4) */ { 0.93 , 1.20, 1.77, 1.86, 2.25, 3.18},
/* N=11 (1 mod 3) (1 mod 4) */ { 0.91 , 1.20, 1.75, 1.84, 2.25, 3.14},
/* N=12 (2 mod 3) (3 mod 4) */ { 0.75 , 0.75, 0.75, 0.75, 0.73, 0.75},
/* N=12 (2 mod 3) (1 mod 4) */ { 0.75 , 0.75, 0.75, 0.73, 0.75, 0.75},
/* N=12 (1 mod 3) (3 mod 4) */ { 0.75 , 0.75, 0.73, 0.75, 0.75, 0.75},
/* N=12 (1 mod 3) (1 mod 4) */ { 0.75 , 0.73, 0.75, 0.75, 0.75, 0.75},
/* N=13 (2 mod 3) (3 mod 4) */ { 1.09 , 1.52, 2.23, 1.73, 2.55, 3.91},
/* N=13 (2 mod 3) (1 mod 4) */ { 1.14 , 1.55, 2.25, 1.73, 2.59, 3.98},
/* N=13 (1 mod 3) (3 mod 4) */ { 1.11 , 1.55, 2.23, 2.14, 2.82, 3.86},
/* N=13 (1 mod 3) (1 mod 4) */ { 1.11 , 1.52, 2.20, 2.11, 2.80, 3.82},
/* N=14 (2 mod 3) (3 mod 4) */ { 0.91 , 0.95, 1.16, 1.41, 1.41, 1.86},
/* N=14 (2 mod 3) (1 mod 4) */ { 0.95 , 0.95, 1.16, 1.43, 1.41, 1.89},
/* N=14 (1 mod 3) (3 mod 4) */ { 0.93 , 0.95, 1.18, 1.86, 1.89, 2.23},
/* N=14 (1 mod 3) (1 mod 4) */ { 0.95 , 0.93, 1.18, 1.86, 1.84, 2.23},
/* N=15 (2 mod 3) (3 mod 4) */ { 0.93 , 1.23, 1.84, 0.95, 1.23, 1.84},
/* N=15 (2 mod 3) (1 mod 4) */ { 0.95 , 1.20, 1.89, 0.93, 1.23, 1.89},
/* N=15 (1 mod 3) (3 mod 4) */ { 0.95 , 1.23, 1.84, 0.93, 1.25, 1.82},
/* N=15 (1 mod 3) (1 mod 4) */ { 0.95 , 1.20, 1.80, 0.95, 1.23, 1.82},
/* N=16 (2 mod 3) (3 mod 4) */ { 1.11 , 1.11, 1.11, 1.73, 1.75, 1.77},
/* N=16 (2 mod 3) (1 mod 4) */ { 1.11 , 1.14, 1.14, 1.77, 1.75, 1.82},
/* N=16 (1 mod 3) (3 mod 4) */ { 1.11 , 1.14, 1.14, 2.14, 2.16, 2.16},
/* N=16 (1 mod 3) (1 mod 4) */ { 1.14 , 1.11, 1.14, 2.14, 2.14, 2.16},
/* N=17 (2 mod 3) (3 mod 4) */ { 1.86 , 3.07, 4.00, 3.25, 5.52, 7.52},
/* N=17 (2 mod 3) (1 mod 4) */ { 1.89 , 3.09, 4.11, 3.27, 5.68, 7.66},
/* N=17 (1 mod 3) (3 mod 4) */ { 1.86 , 3.04, 3.96, 3.32, 5.18, 6.73},
/* N=17 (1 mod 3) (1 mod 4) */ { 1.84 , 3.02, 3.98, 3.27, 5.11, 6.64},
/* N=18 (2 mod 3) (3 mod 4) */ { 1.11 , 1.11, 1.41, 1.11, 1.11, 1.41},
/* N=18 (2 mod 3) (1 mod 4) */ { 1.11 , 1.14, 1.41, 1.14, 1.14, 1.41},
/* N=18 (1 mod 3) (3 mod 4) */ { 1.14 , 1.11, 1.45, 1.14, 1.11, 1.43},
/* N=18 (1 mod 3) (1 mod 4) */ { 1.14 , 1.11, 1.41, 1.14, 1.11, 1.43},
/* N=19 (2 mod 3) (3 mod 4) */ { 2.29 , 3.93, 5.05, 4.14, 7.36, 9.62},
/* N=19 (2 mod 3) (1 mod 4) */ { 2.32 , 3.91, 5.11, 4.18, 7.30, 9.73},
/* N=19 (1 mod 3) (3 mod 4) */ { 2.29 , 3.98, 5.05, 4.09, 6.64, 8.39},
/* N=19 (1 mod 3) (1 mod 4) */ { 2.34 , 3.93, 5.07, 4.02, 6.62, 8.37},
/* N=20 (2 mod 3) (3 mod 4) */ { 1.36 , 1.36, 1.41, 2.25, 2.27, 2.27},
/* N=20 (2 mod 3) (1 mod 4) */ { 1.39 , 1.39, 1.39, 2.30, 2.30, 2.36},
/* N=20 (1 mod 3) (3 mod 4) */ { 1.34 , 1.34, 1.36, 2.55, 2.52, 2.57},
/* N=20 (1 mod 3) (1 mod 4) */ { 1.36 , 1.34, 1.36, 2.50, 2.52, 2.55},
/* N=21 (2 mod 3) (3 mod 4) */ { 1.86 , 3.02, 4.05, 1.86, 3.05, 4.00},
/* N=21 (2 mod 3) (1 mod 4) */ { 1.91 , 3.09, 4.11, 1.89, 3.09, 4.14},
/* N=21 (1 mod 3) (3 mod 4) */ { 1.84 , 3.04, 3.98, 1.86, 3.05, 3.98},
/* N=21 (1 mod 3) (1 mod 4) */ { 1.86 , 2.98, 3.93, 1.86, 3.00, 3.91},
/* N=22 (2 mod 3) (3 mod 4) */ { 2.25 , 2.20, 2.89, 3.98, 3.93, 5.30},
/* N=22 (2 mod 3) (1 mod 4) */ { 2.25 , 2.27, 2.93, 4.05, 4.02, 5.36},
/* N=22 (1 mod 3) (3 mod 4) */ { 2.29 , 2.27, 3.00, 4.00, 4.00, 5.14},
/* N=22 (1 mod 3) (1 mod 4) */ { 2.25 , 2.27, 2.96, 4.02, 3.98, 5.12},
/* N=23 (2 mod 3) (3 mod 4) */ { 3.13 , 5.57, 7.03, 5.82, 10.55, 13.62},
/* N=23 (2 mod 3) (1 mod 4) */ { 3.16 , 5.59, 7.12, 5.82, 10.66, 13.62},
/* N=23 (1 mod 3) (3 mod 4) */ { 3.16 , 5.61, 7.02, 5.41, 9.25, 11.59},
/* N=23 (1 mod 3) (1 mod 4) */ { 3.16 , 5.61, 7.05, 5.36, 9.21, 11.53},
/* N=24 (2 mod 3) (3 mod 4) */ { 2.32 , 2.29, 2.36, 2.29, 2.29, 2.36},
/* N=24 (2 mod 3) (1 mod 4) */ { 2.32 , 2.32, 2.36, 2.32, 2.29, 2.38},
/* N=24 (1 mod 3) (3 mod 4) */ { 2.32 , 2.32, 2.34, 2.32, 2.32, 2.38},
/* N=24 (1 mod 3) (1 mod 4) */ { 2.32 , 2.32, 2.34, 2.34, 2.32, 2.34},
/* N=25 (2 mod 3) (3 mod 4) */ { 3.49 , 6.76, 8.27, 6.99, 12.77, 16.23},
/* N=25 (2 mod 3) (1 mod 4) */ { 3.49 , 6.75, 8.26, 6.75, 12.73, 15.99},
/* N=25 (1 mod 3) (3 mod 4) */ { 3.49 , 6.73, 8.26, 5.98, 11.03, 13.50},
/* N=25 (1 mod 3) (1 mod 4) */ { 3.49 , 6.49, 8.52, 5.97, 11.00, 13.52},
/* N=26 (2 mod 3) (3 mod 4) */ { 2.74 , 2.74, 3.50, 4.74, 5.01, 6.52},
/* N=26 (2 mod 3) (1 mod 4) */ { 2.74 , 2.49, 3.49, 5.00, 4.75, 6.50},
/* N=26 (1 mod 3) (3 mod 4) */ { 2.74 , 2.49, 3.49, 4.51, 4.50, 6.00},
/* N=26 (1 mod 3) (1 mod 4) */ { 2.74 , 2.49, 3.49, 4.74, 4.51, 6.01},
};

#endif