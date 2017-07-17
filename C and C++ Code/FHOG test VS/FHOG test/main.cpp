// Object tracking algorithm using matchTemplate  


#include <opencv2/opencv.hpp>  
#include "stdlib.h"
#include "malloc.h"
#include "string.h"
#include <math.h>
#include <stddef.h>
#include <fstream> //头文件


#define PI 3.14159265f
#define Nmax 100
#define N_h 12
#define N_w 16

using namespace cv;
using namespace std;

void fhog_test(float *Im, float *H1, float(*H2)[N_w/4][N_h/4], int h, int w);
float* acosTable(void);
void gradientMagnitude(float *I, float *M, float *O, int h, int w, int d, bool full);
void gradMagNormalization(float *M, float *S, int h, int w, float norm);
void gradientHist(float *M, float *O, float *H, int h, int w, int bin, int nOrients, int softBin, bool full);
void hog(float *M, float *O, float *H, int h, int w, int binSize, int nOrients, int softBin,  float clip);
void fhog(float *M, float *O, float *H, int h, int w, int binSize, int nOrients, int softBin, float clip);

float H[(N_h/4)*(N_w/4) * 32] = { 0 }, M[N_h * N_w] = { 0 }, O[N_h * N_w] = { 0 }, H_out[32][N_w/4][N_h/4] = { 0 };

int main(int argc, char * argv[])
{

	Mat src_mat = imread("suv_2.jpg");
	Mat gray_mat;
	//int h = src_mat.rows, w = src_mat.cols;
	int h = N_h, w = N_w;
	cvtColor(src_mat, gray_mat, CV_BGR2GRAY);
	float data[N_h][N_w] = { 0 };
	float data_zhuan[N_h*N_w] = { 0 };
	unsigned char* idata = (unsigned char*)gray_mat.data;
	for (int i = 0; i < h; i++)
		for (int j = 0; j < w; j++)
		{
			data[i][j] = (float)idata[i*w+j];
		}
	for (int i = 0; i < w; i++)
		for (int j = 0; j < h; j++)
		{
			data_zhuan[i*h+j] = data[j][i];
		}
	//imshow("image", img);
/****************************************************************************/
	gradientMagnitude(data_zhuan, M, O, h, w, 1, true);
	int binSize = 4; int nOrients = 9; int softBin = -1; float clip = 0.2f;
	int hb = h / binSize; int wb = w / binSize; int nChns = nOrients * 3 + 5;
	//hog(M, O, H, h, w, binSize, nOrients, softBin, clip);
	fhog(M, O, H, h, w, binSize, nOrients, softBin, clip);
	int i, j, k;
	//ofstream fout("H_qw.txt");//默认路径是工程文件下
	for (i = 0; i < nChns; i++)
	{
		for (j = 0; j < wb; j++)
		{
			for (k = 0; k < hb; k++)
			{
				H_out[i][j][k] = H[i*hb*wb + j*hb + k];
				//fout << H2[i][j][k] << endl;
				//memcpy(H_out[k][j][i], H[i*hb*wb + j*hb + k], sizeof(float));
			}
		}
	}
	//fout.close();
	waitKey(0);
	//system("pause");
	return 0;
}

/****************************************************************************/

// compute x and y gradients for just one column 仅计算一列的x和y梯度
void grad1(float *I, float *Gx, float *Gy, int h, int w, int x) {
	int y, y1; float *Ip, *In, r;
	// compute column of Gx
	Ip = I - h; In = I + h; r = .5f;
	if (x == 0) { r = 1; Ip += h; }
	else if (x == w - 1) { r = 1; In -= h; }
	if (h<4 || h % 4>0) {													//此处与源程序不同
		for (y = 0; y<h; y++) *Gx++ = (*In++ - *Ip++)*r;//Gx=G(x+1)-G(x-1)
	}
	else {
		for (y = 0; y<h; y++) Gx[y] = (In[y] - Ip[y])*r;
	}

	// compute column of Gy b
#define GRADY(r) *Gy++=(*In++-*Ip++)*r;
	Ip = I; In = Ip + 1;

	// GRADY(1); Ip--; for(y=1; y<h-1; y++) GRADY(.5f); In--; GRADY(1);
	/* y1=((~((size_t) Gy) + 1) & 15)/4;
	if(y1==0) y1=4;
	if(y1>h-1) y1=h-1;
	GRADY(1); Ip--; for(y=1; y<y1; y++) GRADY(.5f);*/
	y1 = h - 1;
	GRADY(1); Ip--; for (y = 1; y<y1; y++) GRADY(.5f);
	for (; y + 4<h - 1; y++)
		Gy[y] = (In[y] - Ip[y])*0.5f;

	for (; y<h - 1; y++) GRADY(.5f); In--; GRADY(1);
#undef GRADY
}

// compute x and y gradients at each location (uses sse)在每个位置计算x和y梯度
void grad2(float *I, float *Gx, float *Gy, int h, int w, int d) {
	int o, x, c, a = w*h; for (c = 0; c<d; c++) for (x = 0; x<w; x++) {
		o = c*a + x*h; grad1(I + o, Gx + o, Gy + o, h, w, x);
	}
}

// build lookup table a[] s.t. a[x*n]~=acos(x) for x in [-1,1]
float* acosTable(void) {
	const int n = 10000, b = 10; int i;
	static float a[10000 * 2 + 10 * 2];
	static bool init = false;
	float *a1 = a + n + b;
	if (init)
		return a1;
	for (i = -n - b; i<-n; i++)
		a1[i] = PI;
	for (i = -n; i<n; i++)
		a1[i] = (float)acos(i / (float)n);
	for (i = n; i<n + b; i++)
		a1[i] = 0;
	for (i = -n - b; i<n / 10; i++)
		if (a1[i] > PI - 1e-6f)
			a1[i] = PI - 1e-6f;
	init = true; return a1;
}
// compute gradient magnitude and orientation at each location (uses sse)在每个位置计算梯度的大小和方向
void gradientMagnitude(float *I, float *M, float *O, int h, int w, int d, bool full) {
	int x, y, y1, c, h4;//s;
	float *acost = acosTable();
	float acMult = 10000.0f;
	// allocate memory for storing one column of output (padded so h4%4==0)分配内存来存储一列的输出
	h4 = (h % 4 == 0) ? h : h - (h % 4) + 4; //s = d*h4*sizeof(float);//h4如果不是4的整数倍就补足，添加到4的整数倍
	float M2[N_h] = { 0 }, Gx[N_h] = { 0 }, Gy[N_h] = { 0 };
	/*M2 = (float*)malloc(sizeof(float)*d*h4);Gx = (float*)malloc(sizeof(float)*d*h4);Gy = (float*)malloc(sizeof(float)*d*h4);
	memset(M2, 0, s);memset(Gx, 0, s);memset(Gy, 0, s);*/

	// compute gradient magnitude and orientation for each column 计算每个列的梯度大小和方向
	for (x = 0; x<w; x++) {
		// compute gradients (Gx, Gy) with maximum squared magnitude (M2)计算梯度（GX、Gy）的最大幅值的平方（M2），本循环中也计算了梯度方向（O）
		for (c = 0; c<d; c++) {
			grad1(I + x*h + c*w*h, Gx + c*h4, Gy + c*h4, h, w, x);
			for (y = 0; y<h4; y++) {
				y1 = h4*c + y;
				M2[y1] = (Gx[y1] * Gx[y1] + Gy[y1] * Gy[y1]);//求Gx和Gy的平方和 M2
				if (c == 0) continue;
				M2[y] = M2[y1] > M2[y] ? M2[y1] : M2[y];
				Gx[y] = M2[y1] > M2[y] ? Gx[y1] : Gx[y];
				Gy[y] = M2[y1] > M2[y] ? Gy[y1] : Gy[y];//取较大的值
			}

		}
		// compute gradient mangitude (M) and normalize Gx 计算梯度幅值M和范数Gx
		for (y = 0; y<h4; y++) {

			//----因为这段代码的版本所以cce存在小的差异
			float m = 1.0f / sqrt((float)M2[y]);
			m = m < 1e10f ? m : 1e10f;
			M2[y] = 1.0f / m; //此处M2内容变成了sqrt（Gx2+Gy2）
			//----------------
			if (O) Gx[y] = (Gx[y] * m)*acMult; //acMult=10000.0f

			if (O) {
				//bitwise AND on floats按位与
				float zero = -0.f;
				unsigned char *pGy = 0;
				unsigned char *pZero = 0;
				unsigned char *pGand = 0;
				float Gand = 0.f;
				pGy = reinterpret_cast<unsigned char *>(&Gy[y]);
				pZero = reinterpret_cast<unsigned char *>(&zero);
				pGand = reinterpret_cast<unsigned char *>(&Gand);
				for (int i = 0; i<4; i++){
					*pGand = (*pGy & *pZero);
					pGand++;
					pGy++;
					pZero++;
				};

				//bitwise XOR on floats按位异或
				unsigned char *pGx = 0;
				unsigned char *pGxor = 0;
				float Gxor = 0;
				pGx = reinterpret_cast<unsigned char *>(&Gx[y]);
				pGand = reinterpret_cast<unsigned char *>(&Gand);
				pGxor = reinterpret_cast<unsigned char *>(&Gxor);
				for (int i = 0; i<4; i++){
					*pGxor = (*pGx ^ *pGand);
					pGxor++;
					pGx++;
					pGand++;
				};

				Gx[y] = Gxor;
			};
		};

		memcpy(M + x*h, M2, h*sizeof(float));//将M2的值复制到 M 中
		// compute and store gradient orientation (O) via table lookup 计算和存储梯度方向（O）通过查找表
		if (O != 0) for (y = 0; y<h; y++)
			O[x*h + y] = acost[(int)Gx[y]];//O存储的是归一化后的梯度方向
		if (O != 0 && full) {
			int y1 = (((~(int)(O + x*h)) + 1) & 15) / 4;
			//y1 = ((~size_t(O + x*h) + 1) & 15) / 4; 
			y = 0;
			for (; y<y1; y++) O[y + x*h] += (Gy[y]<0)*PI; //如果Gy小于0则将O相应位置方向值+180度，否则不作变换
			for (; y<h - 4; y++)
				O[y + x*h] += (Gy[y]<0)*PI;
			for (; y<h; y++) O[y + x*h] += (Gy[y]<0)*PI;//此处不懂：y为何要从0~y1~h-4~h分成三段来判断
		}
	}
}

// helper for gradHist, quantize O and M into O0, O1 and M0, M1 gradhist 辅助，量化O和M到O0、O1和M0、M1
void gradQuantize(float *O, float *M, int *O0, int *O1, float *M0, float *M1,
	int nb, int n, float norm, int nOrients, bool full)//整个文件里full参数代表方向bin数，如果为1则取0-2PI，0则取0-PI
{//	s=(float)bin, hb=h/bin, wb=w/bin, h0=hb*bin, w0=wb*bin, nb=wb*hb;   
	//norm=1/s/s; n=h0; interpolate=false
	// assumes all *OUTPUT* matrices are 4-byte aligned假设所有*输出*矩阵是4字节对齐
	int i, o0, o1; float o, od, m;
	const float oMult = (float)nOrients / (full ? 2 * PI : PI); const int oMax = nOrients*nb;//oMult=9/2PI或9/PI

	for (i = 0; i <= n - 4; i++) {
		o = O[i] * oMult;
		o0 = (int)(o + 0.5f);
		o0 *= nb;
		o0 = (oMax > o0) ? o0 : 0;
		O0[i] = o0;
		M0[i] = M[i] * norm;
		M1[i] = 0.0f; O1[i] = 0;
		}
	// compute trailing locations without sse 不使用sse计算跟踪位置
	//interpolate为false，因为softBin>=0为假//softBin >= 0使用临近两个bin线性插值，softBin<0则取最近的bin（本程序中为 - 1）
	for (i; i<n; i++) {
		o = O[i] * oMult;				//oMult = 9/2PI
		o0 = (int)(o + .5f);
		o0 *= nb;
		if (o0 >= oMax)
			o0 = 0;
		O0[i] = o0;				//O = wb*hb*((O*9/2PI)+0.5)
		M0[i] = M[i] * norm; //M0为梯度幅度M乘以权值后
		M1[i] = 0; O1[i] = 0; //M1和O1都是0
	}
}

// compute nOrients gradient histograms per bin x bin block of pixels 计算nOrients个梯度直方图每个bin x bin块像素
void gradientHist(float *M, float *O, float *H, int h, int w,
	int bin, int nOrients, int softBin, bool full)
	//H维度[wb*hb*9*3];bin=4;nOrients=2*9;softBin=-1;full=1

	//softBin表示选择梯度幅值加权形式，softBin>=0使用临近两个bin线性插值，softBin<0则取最近的bin（本程序中为-1）。
	//在每个w*h的方向通道的每一binsize * binsize大小的像素，决定了空间的分级。如果“softbin”是奇数像素
	//可以影响多个bin（使用双线性插值），否则每个像素只对一个方向bin有贡献。 经过这一步过后的结果为
	//floor([h / binSize w / binSize nOrients])维大小的特征图，代表了每个图层区域中的梯度直方图。
{
	const int hb = h / bin, wb = w / bin, h0 = hb*bin, w0 = wb*bin, nb = wb*hb;
	const float s = (float)bin, sInv = 1 / s, sInv2 = 1 / s / s;
	float *H0, *H1; int x, y; float  xb, init;
	int   O0[N_h] = { 0 }, O1[N_h] = { 0 };
	float M0[N_h] = { 0 }, M1[N_h] = { 0 };

	// main loop
	for (x = 0; x<w0; x++) {
		// compute target orientation bins for entire column - very fast 在整个列计算目标方向bin-非常快
		gradQuantize(O + x*h, M + x*h, O0, O1, M0, M1, nb, h0, sInv2, nOrients, full);//量化O和M
		if (softBin==-1) { // interpolate using trilinear interpolation 采用三线性插值插值 //softBin为奇数
			float ms[4], xyd, yb, xd, yd;
			bool hasLf, hasRt; int xb0, yb0;
			if (x == 0) { init = (0 + .5f)*sInv - 0.5f; xb = init; }   //初始化xb = init = 0.5*(1/4)-0.5=-0.375
			hasLf = xb >= 0; xb0 = hasLf ? (int)xb : -1; hasRt = xb0 < wb - 1;//hasLf=0，xb0=-1；hasRt=1；
			xd = xb - xb0; xb += sInv; yb = init; y = 0;//xd=xb+1=0.625;xb=-0.125;yb=-0.375;y=0;
			// macros for code conciseness代码简洁宏
			#define GHinit yd=yb-yb0; yb+=sInv; H0=H+xb0*hb+yb0; xyd=xd*yd; \
			   ms[0]=1-xd-yd+xyd; ms[1]=yd-xyd; ms[2]=xd-xyd; ms[3]=xyd;
			// leading rows, no top bin领先的行，没有顶bin
			for (; y<bin / 2; y++) {
				yb0 = -1; GHinit;
				if (hasLf) { H0[O0[y] + 1] += ms[1] * M0[y]; H0[O1[y] + 1] += ms[1] * M1[y]; }
				if (hasRt) { H0[O0[y] + hb + 1] += ms[3] * M0[y]; H0[O1[y] + hb + 1] += ms[3] * M1[y]; }
			}

			// main rows, has top and bottom bins, use SSE for minor speedup主行，有顶部和底部的bin，使用SSE的小加速
				for (;; y++) {
					yb0 = (int)yb; if (yb0 >= hb - 1) break; GHinit;
					if (hasLf) {
						H1 = H0 + O0[y]; 
						*H1 += 0 * M0[y]; *(H1 + 1) += 0 * M0[y];*(H1 + 2) += ms[1] * M0[y]; *(H1 + 3) += ms[0] * M0[y];
						//*H1 += (ms[0] * M0[y] + ms[1] * M0[y]);
					}
					if (hasRt) {
						H1 = H0 + O0[y] + hb;
						*H1 += 0 * M0[y]; *(H1 + 1) += 0 * M0[y];*(H1 + 2) += ms[3] * M0[y];*(H1 + 3) += ms[2] * M0[y];
						//*H1 += (ms[2] * M0[y] + ms[3] * M0[y]);
					}
				}
			// final rows, no bottom bin最后一行，没有底部bin
			for (; y<h0; y++) {
				yb0 = (int)yb; GHinit;
				if (hasLf) { H0[O0[y]] += ms[0] * M0[y]; H0[O1[y]] += ms[0] * M1[y]; }
				if (hasRt) { H0[O0[y] + hb] += ms[2] * M0[y]; H0[O1[y] + hb] += ms[2] * M1[y]; }
			}
			#undef GHinit
		}//三线插值else的末尾
	}

	// normalize boundary bins which only get 7/8 of weight of interior bins归一化边界bin，只有内部bin的7/8权重
	int o;
	if (softBin % 2 != 0) for (o = 0; o<nOrients; o++) {
		x = 0; for (y = 0; y<hb; y++) H[o*nb + x*hb + y] *= 8.f / 7.f;
		y = 0; for (x = 0; x<wb; x++) H[o*nb + x*hb + y] *= 8.f / 7.f;
		x = wb - 1; for (y = 0; y<hb; y++) H[o*nb + x*hb + y] *= 8.f / 7.f;
		y = hb - 1; for (x = 0; x<wb; x++) H[o*nb + x*hb + y] *= 8.f / 7.f;
	}
}

/******************************************************************************/

// HOG helper: compute 2x2 block normalization values (padded by 1 pixel)HOG帮手：计算2x2块归一化值（用1像素）
float* hogNormMatrix(float *N ,float *H, int nOrients, int hb, int wb, int bin) {
	float *N1, *n; int o, x, y, dx, dy, hb1 = hb + 1, wb1 = wb + 1;
	float eps = 1e-4f / 4 / bin / bin / bin / bin; // precise backward equality精确的反向等式
	//N = (float*)malloc(sizeof(float)*hb1*wb1);
	//memset(N, 0, hb1*wb1*sizeof(float));
	N1 = N + hb1 + 1;
	for (o = 0; o<nOrients; o++) for (x = 0; x<wb; x++) for (y = 0; y<hb; y++)
		N1[x*hb1 + y] += H[o*wb*hb + x*hb + y] * H[o*wb*hb + x*hb + y];
	for (x = 0; x<wb - 1; x++) for (y = 0; y<hb - 1; y++) {
		n = N1 + x*hb1 + y; *n = 1 / (float)sqrt(n[0] + n[1] + n[hb1] + n[hb1 + 1] + eps);
	}
	x = 0;     dx = 1; dy = 1; y = 0;						N[x*hb1 + y] = N[(x + dx)*hb1 + y + dy];
	x = 0;     dx = 1; dy = 0; for (y = 0; y<hb1; y++)		N[x*hb1 + y] = N[(x + dx)*hb1 + y + dy];
	x = 0;     dx = 1; dy = -1; y = hb1 - 1;				N[x*hb1 + y] = N[(x + dx)*hb1 + y + dy];
	x = wb1 - 1; dx = -1; dy = 1; y = 0;					N[x*hb1 + y] = N[(x + dx)*hb1 + y + dy];
	x = wb1 - 1; dx = -1; dy = 0; for (y = 0; y<hb1; y++)	N[x*hb1 + y] = N[(x + dx)*hb1 + y + dy];
	x = wb1 - 1; dx = -1; dy = -1; y = hb1 - 1;             N[x*hb1 + y] = N[(x + dx)*hb1 + y + dy];
	y = 0;     dx = 0; dy = 1; for (x = 0; x<wb1; x++)		N[x*hb1 + y] = N[(x + dx)*hb1 + y + dy];
	y = hb1 - 1; dx = 0; dy = -1; for (x = 0; x<wb1; x++)	N[x*hb1 + y] = N[(x + dx)*hb1 + y + dy];
	return N;
}

// HOG helper: compute HOG or FHOG channels HOG辅助：计算HOG或FHOG的多个通道
void hogChannels(float *H, const float *R, const float *N,
	int hb, int wb, int nOrients, float clip, int type)
{
#define GETT(blk) t=R1[y]*N1[y-(blk)]; if(t>clip) t=clip; c++;
	const float r = .2357f; int o, x, y, c; float t;
	const int nb = wb*hb, nbo = nOrients*nb, hb1 = hb + 1;
	for (o = 0; o<nOrients; o++) for (x = 0; x<wb; x++) {
		const float *R1 = R + o*nb + x*hb, *N1 = N + x*hb1 + hb1 + 1;
		float *H1 = (type <= 1) ? (H + o*nb + x*hb) : (H + x*hb);
		if (type == 0) for (y = 0; y<hb; y++) {
			// store each orientation and normalization (nOrients*4 channels) 存储每个方向与归一化（nOrients×4通道）
			c = -1; GETT(0); H1[c*nbo + y] = t; GETT(1); H1[c*nbo + y] = t;
			GETT(hb1); H1[c*nbo + y] = t; GETT(hb1 + 1); H1[c*nbo + y] = t;
		}
		else if (type == 1) for (y = 0; y<hb; y++) {
			// sum across all normalizations (nOrients channels) 求所有的归一化值得和（nOrients通道）
			c = -1; GETT(0); H1[y] += t*.5f; GETT(1); H1[y] += t*.5f;
			GETT(hb1); H1[y] += t*.5f; GETT(hb1 + 1); H1[y] += t*.5f;
		}
		else if (type == 2) for (y = 0; y<hb; y++) {
			// sum across all orientations (4 channels) 求所有方向的总和（4个通道）
			c = -1; GETT(0); H1[c*nb + y] += t*r; GETT(1); H1[c*nb + y] += t*r;
			GETT(hb1); H1[c*nb + y] += t*r; GETT(hb1 + 1); H1[c*nb + y] += t*r;
		}
	}
#undef GETT
}

// compute HOG features
void hog(float *M, float *O, float *H, int h, int w, int binSize,
	int nOrients, int softBin, float clip)
{
	const int hb = h / binSize, wb = w / binSize, nb = hb*wb;
	//compute unnormalized gradient histograms 计算非归一化的梯度直方图
	float R[(N_h / 4)*(N_w / 4) * 9] = { 0 };
	gradientHist(M, O, R, h, w, binSize, nOrients, softBin, true);
	// compute block normalization values 计算块归一化值
	float N[(N_h / 4 + 1)*(N_w / 4 + 1)] = { 0 };
	hogNormMatrix(N,R, nOrients, hb, wb, binSize);
	// perform four normalizations per spatial block 每个空间块执行四次归一化
	hogChannels(H, R, N, hb, wb, nOrients, clip, 0);
}
// compute FHOG features
void fhog(float *M, float *O, float *H, int h, int w, int binSize,
	int nOrients, int softBin, float clip)
{
	const int hb = h / binSize, wb = w / binSize, nb = hb*wb, nbo = nb*nOrients;
	int o, x;
	// compute unnormalized constrast sensitive histograms	计算非归一化对比度敏感的直方图
	//R1 = (float*)malloc(sizeof(float)*wb*hb*nOrients * 2);//为了不超出单元格分配的内存，所以应该是2
	//memset(R1, 0, wb*hb*nOrients * 2 * sizeof(float));
	float R1[(N_h / 4)*(N_w / 4) * 18] = { 0 };	
	float R2[(N_h / 4)*(N_w / 4) * 9 ] = { 0 };
	gradientHist(M, O, R1, h, w, binSize, nOrients * 2, softBin, true);
	// compute unnormalized contrast insensitive histograms	计算非归一化对比度不敏感的直方图
	//R2 = (float*)malloc(sizeof(float)*wb*hb*nOrients);
	//memset(R2, 0, wb*hb*nOrients * sizeof(float));
	for (o = 0; o<nOrients; o++)
		for (x = 0; x<nb; x++)
			R2[o*nb + x] = R1[o*nb + x] + R1[(o + nOrients)*nb + x];
	// compute block normalization values 计算块归一化值
	float N[(N_h / 4 + 1)*(N_w / 4 + 1)] = { 0 };
	hogNormMatrix(N,R2, nOrients, hb, wb, binSize);
	// normalized histograms and texture channels 归一化直方图和纹理通道
	hogChannels(H + nbo * 0, R1, N, hb, wb, nOrients * 2, clip, 1);
	hogChannels(H + nbo * 2, R2, N, hb, wb, nOrients * 1, clip, 1);
	hogChannels(H + nbo * 3, R1, N, hb, wb, nOrients * 2, clip, 2);
}
