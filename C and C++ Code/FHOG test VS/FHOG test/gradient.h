#pragma once

#ifndef _DSST__SETTINGS_H_
#define _DSST__SETTINGS_H_
//#define SSEv2  // 使用SSE2是否加快工作
#endif

#ifndef GRADIENT_H
#define GRADIENT_H

//#include "DSSTSettings.h"

float* acosTable(void);
void gradMag(float *I, float *M, float *O, int h, int w, int d, bool full);
void gradMagNorm(float *M, float *S, int h, int w, float norm);
void gradHist(float *M, float *O, float *H, int h, int w, int bin, int nOrients, int softBin, bool full);
void hog(float *M, float *O, float *H, int h, int w, int binSize, int nOrients, int softBin, float clip);
void fhog(float *M, float *O, float *H, int h, int w, int binSize, int nOrients, int softBin, float clip);
#endif
