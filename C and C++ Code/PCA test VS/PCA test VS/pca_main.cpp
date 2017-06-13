#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "vector"
using namespace std;


#define E 0.0000001
#define INF 99999 
#define dimNum 4     //维数
#define MAXITER 1000   //最大迭代次数


typedef vector<double> doubleVector;
typedef vector<doubleVector> dim2Vector;


vector<doubleVector> getInputSample(char* File);  //获取输入样本
vector<doubleVector> normalizationSPSS(vector<doubleVector> inputTrain);  //采用z-score法标准数据
vector<doubleVector> normalizationMAX_MIN(vector<doubleVector> inputTrain);  //采用最大最小法法标准数据
void PAC(vector<doubleVector> inputTrain);  //主成分分析法PAC
vector<doubleVector> calCovariation(vector<doubleVector> inputTrain);   //计算协方差
vector<dim2Vector> Jacobi(vector<doubleVector> Array);   //使用Jacobi计算协方差的特征值和特征矩阵
bool QueryArray(vector<doubleVector> Array);   //检查是否满足
vector<doubleVector> matTran(vector<doubleVector> Array);   //矩阵转置
vector<doubleVector> matMul(vector<doubleVector> mat1, vector<doubleVector> mat2);   //矩阵相乘


double Input_Meam[dimNum] = { 0 };    //每一维的均值
double Input_Dev[dimNum] = { 0 };     //每一维的标准差


void main()
{
	char *File = "input.txt";


	vector<doubleVector> inputTrain;


	inputTrain = getInputSample(File);


	inputTrain = normalizationSPSS(inputTrain);  //采用z-score法标准数据
	// inputTrain = normalizationMAX_MIN(inputTrain);  //采用最大最小法法标准数据


	PAC(inputTrain);  //主成分分析法PAC

	system("pause");
}




//主成分分析法PAC
void PAC(vector<doubleVector> inputTrain)
{
	int i, j, m, n;
	vector<doubleVector> input_Cov;  //协方差
	vector<dim2Vector> jacobi;   //1为特征值，2为特征矩阵
	double rate;  //贡献率
	double rateSum1 = 0;
	double rateSum2 = 0;
	doubleVector tempVector;
	vector<doubleVector> redTemp;
	vector<doubleVector> reduce_Dim_Mat;  //降维矩阵
	vector<doubleVector> reduce_Dim_Sample;  //降维数据


	input_Cov = calCovariation(inputTrain);   //计算协方差


	jacobi = Jacobi(input_Cov);   //使用Jacobi计算协方差的特征值和特征矩阵


	//计算贡献率
	for (i = 0; i<jacobi[0].size(); i++)
	{
		for (j = 0; j<jacobi[0][i].size(); j++)
			rateSum1 += jacobi[0][i][j];


		for (j = 0; j<jacobi[0][i].size(); j++)
		{
			rateSum2 += jacobi[0][i][j];
			rate = rateSum2 / rateSum1;
			if (rate >= 0.85)
				break;
		}


		//获取将维矩阵
		for (m = 0; m <= j; m++)
		{
			tempVector.clear();
			for (n = 0; n<jacobi[1][m].size(); n++)
				tempVector.push_back(jacobi[1][n][m]);


			reduce_Dim_Mat.push_back(tempVector);
		}
	}


	reduce_Dim_Mat = matTran(reduce_Dim_Mat);


	reduce_Dim_Sample = matMul(inputTrain, reduce_Dim_Mat);  //计算降维结果


	printf("协方差为:\n");
	for (i = 0; i<input_Cov.size(); i++)
	{
		for (j = 0; j<input_Cov[i].size(); j++)
			printf("%lf  ", input_Cov[i][j]);
		printf("\n");
	}


	printf("\n特征值：\n");
	for (i = 0; i<jacobi[0].size(); i++)
	{
		for (j = 0; j<jacobi[0][i].size(); j++)
			printf("%lf  ", jacobi[0][i][j]);
		printf("\n");
	}

	printf("\n特征向量：\n");
	for (i = 0; i<jacobi[1].size(); i++)
	{
		for (j = 0; j<jacobi[1][i].size(); j++)
			printf("%lf  ", jacobi[1][i][j]);
		printf("\n");
	}




	printf("\n降维矩阵：\n");
	for (i = 0; i<reduce_Dim_Mat.size(); i++)
	{
		for (j = 0; j<reduce_Dim_Mat[i].size(); j++)
			printf("%lf  ", reduce_Dim_Mat[i][j]);
		printf("\n");
	}


	printf("\n降维结果：\n");
	for (i = 0; i<reduce_Dim_Sample.size(); i++)
	{
		for (j = 0; j<reduce_Dim_Sample[i].size(); j++)
			printf("%lf  ", reduce_Dim_Sample[i][j]);
		printf("\n");
	}

}




//计算协方差
vector<doubleVector> calCovariation(vector<doubleVector> inputTrain)
{
	int i, j, k;
	doubleVector tempDst(dimNum, 0);
	vector<doubleVector> dst(dimNum, tempDst);

	for (i = 0; i<dimNum; i++)
		Input_Meam[i] = 0;

	//计算均值
	for (i = 0; i<dimNum; i++)
	{
		for (j = 0; j<inputTrain.size(); j++)
			Input_Meam[i] += inputTrain[j][i];

		Input_Meam[i] = Input_Meam[i] / inputTrain.size();
	}

	//计算协方差
	for (i = 0; i<dimNum; i++)
		for (j = 0; j<dimNum; j++)
		{
			for (k = 0; k<inputTrain.size(); k++)
				dst[i][j] += (inputTrain[k][i] - Input_Meam[i])*(inputTrain[k][j] - Input_Meam[j]);

			dst[i][j] = dst[i][j] / (inputTrain.size() - 1);

		}


	return dst;
}






//使用Jacobi计算协方差的特征值和特征矩阵
vector<dim2Vector> Jacobi(vector<doubleVector> Array)
{
	int i, j;
	int count;
	bool flag = false;
	vector<dim2Vector> dst;
	doubleVector tempArray(Array.size(), 0);
	vector<doubleVector> charatMat(Array.size(), tempArray);   //特征向量
	vector<doubleVector> sortArray;  //排序后的特征值
	vector<doubleVector> dim2Jac;
	vector<doubleVector> dim2JacT;
	vector<dim2Vector> dim3Jac;
	double maxArrayNum;
	int laber_j, laber_i;

	double theta;


	//开始迭代
	count = 0;
	tempArray.clear();
	tempArray.resize(Array.size(), 0);
	while (count<MAXITER && !flag)
	{
		count++;
		dim2Jac.clear();
		dim2Jac.resize(Array.size(), tempArray);
		maxArrayNum = 0;
		laber_i = laber_j = 0;


		//寻找非对角元中绝对值最大的A[i][j]
		for (i = 0; i<Array.size(); i++)
			for (j = 0; j<Array.size(); j++)
			{
				if (i == j)
					continue;


				if (maxArrayNum<fabs(Array[i][j]))
				{
					maxArrayNum = fabs(Array[i][j]);
					laber_i = i;
					laber_j = j;
				}
			}


		theta = atanf(Array[laber_i][laber_j] * 2 / (Array[laber_i][laber_i] - Array[laber_j][laber_j] + E));


		//构造雅克比矩阵
		for (i = 0; i<Array.size(); i++)
			dim2Jac[i][i] = 1;


		dim2Jac[laber_i][laber_i] = dim2Jac[laber_j][laber_j] = cosf(theta / 2);
		dim2Jac[laber_i][laber_j] = sinf(theta / 2);
		dim2Jac[laber_j][laber_i] = -sinf(theta / 2);


		dim2JacT = matTran(dim2Jac);  //矩阵转置
		dim3Jac.push_back(dim2JacT);  //保存矩阵


		Array = matMul(matMul(dim2Jac, Array), dim2JacT);


		if (QueryArray(Array))
			flag = true;

	}


	//初始化特征矩阵
	for (i = 0; i<Array.size(); i++)
		charatMat[i][i] = 1;


	//计算特征矩阵
	for (i = 0; i<dim3Jac.size(); i++)
		charatMat = matMul(charatMat, dim3Jac[i]);


	//排序
	doubleVector sortA;
	double tempNum;
	for (i = 0; i<Array.size(); i++)
		sortA.push_back(Array[i][i]);


	for (i = 0; i<sortA.size(); i++)
	{
		maxArrayNum = sortA[i];
		laber_j = i;


		for (j = i; j<sortA.size(); j++)
			if (maxArrayNum<sortA[j])
			{
				maxArrayNum = sortA[j];
				laber_j = j;
			}


		tempNum = sortA[i];
		sortA[i] = sortA[laber_j];
		sortA[laber_j] = tempNum;


		for (j = 0; j<charatMat[laber_j].size(); j++)
			tempArray[j] = charatMat[j][i];


		for (j = 0; j<charatMat[laber_j].size(); j++)
			charatMat[j][i] = charatMat[j][laber_j];


		for (j = 0; j<charatMat[laber_j].size(); j++)
			charatMat[j][laber_j] = tempArray[j];


	}


	sortArray.push_back(sortA);


	dst.push_back(sortArray);
	dst.push_back(charatMat);


	return dst;
}




//检查是否满足
bool QueryArray(vector<doubleVector> Array)
{
	int i, j;


	for (i = 0; i<Array.size(); i++)
		for (j = 0; j<Array.size(); j++)
		{
			if (i == j)
				continue;


			if (fabs(Array[i][j])>E)
				return false;
		}


	return true;
}




//矩阵转置
vector<doubleVector> matTran(vector<doubleVector> Array)
{
	int i, j;
	doubleVector temp(Array.size(), 0);
	vector<doubleVector> dst(Array[0].size(), temp);


	for (i = 0; i<Array.size(); i++)
		for (j = 0; j<Array[0].size(); j++)
			dst[j][i] = Array[i][j];


	return dst;
}


//矩阵相乘
vector<doubleVector> matMul(vector<doubleVector> mat1, vector<doubleVector> mat2)
{
	int i, j, k;
	doubleVector temp(mat2[0].size(), 0);
	vector<doubleVector> dst(mat1.size(), temp);


	for (i = 0; i<mat1.size(); i++)
		for (j = 0; j<mat2[0].size(); j++)
			for (k = 0; k<mat2.size(); k++)
				dst[i][j] += mat1[i][k] * mat2[k][j];


	return dst;
}




//采用最大最小法法标准数据
vector<doubleVector> normalizationMAX_MIN(vector<doubleVector> inputTrain)
{
	int i, j;
	double input_Max[dimNum], input_Min[dimNum];
	vector<doubleVector> dst;
	doubleVector tempDst;


	//初始化
	for (i = 0; i<dimNum; i++)
	{
		input_Max[i] = 0;
		input_Min[i] = INF;
	}


	//寻找最大最小值
	for (i = 0; i<dimNum; i++)
		for (j = 0; j<inputTrain.size(); j++)
		{
			if (input_Max[i]<inputTrain[j][i])
				input_Max[i] = inputTrain[j][i];


			if (input_Min[i]>inputTrain[j][i])
				input_Min[i] = inputTrain[j][i];
		}




	//归一化
	for (i = 0; i<inputTrain.size(); i++)
	{
		tempDst.clear();
		for (j = 0; j<inputTrain[i].size(); j++)
			tempDst.push_back((inputTrain[i][j] - input_Min[j]) / (input_Max[j] - input_Min[j]));

		dst.push_back(tempDst);
	}




	return dst;


}




//采用z-score法标准数据
vector<doubleVector> normalizationSPSS(vector<doubleVector> inputTrain)
{
	int i, j;
	vector<doubleVector> dst;
	doubleVector tempDst;


	//初始化
	for (i = 0; i<dimNum; i++)
	{
		Input_Meam[i] = 0;
		Input_Dev[i] = 0;
	}


	//计算均值
	for (i = 0; i<dimNum; i++)
	{
		for (j = 0; j<inputTrain.size(); j++)
			Input_Meam[i] += inputTrain[j][i];


		Input_Meam[i] = Input_Meam[i] / inputTrain.size();
	}


	//计算标准差
	for (i = 0; i<dimNum; i++)
	{
		for (j = 0; j<inputTrain.size(); j++)
			Input_Dev[i] += (inputTrain[j][i] - Input_Meam[i])*(inputTrain[j][i] - Input_Meam[i]);

		Input_Dev[i] = sqrtf(Input_Dev[i] / (inputTrain.size() - 1));
	}


	//标准化
	for (i = 0; i<inputTrain.size(); i++)
	{
		tempDst.clear();
		for (j = 0; j<inputTrain[i].size(); j++)
			tempDst.push_back((inputTrain[i][j] - Input_Meam[j]) / Input_Dev[j]);


		dst.push_back(tempDst);
	}




	return dst;
}




//获取输入样本
vector<doubleVector> getInputSample(char* File)
{
	vector<doubleVector> dst;
	doubleVector temp;
	int i;
	double num;


	FILE *fp = fopen(File, "r");

	if (fp == NULL)
	{
		printf("OPEN FILE ERROR!!\n");
		exit(0);
	}


	//从文件读取样本
	i = 1;
	temp.clear();
	dst.clear();
	while (fscanf(fp, "%lf", &num) != EOF)
	{
		temp.push_back(num);
		if (i%dimNum == 0)
		{
			dst.push_back(temp);
			temp.clear();
		}
		i++;
	}


	return dst;
}
