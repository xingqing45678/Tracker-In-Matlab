#include <stdio.h>
#include <stdlib.h>
#include <string.h>
void main()
{
	float Gy = -0.f, Gx = 0.f;	char s[100];
	//Gx = Gx*10000.0f;
	float zero = -0.f;
	unsigned char *pGy = 0;
	unsigned char *pZero = 0;
	unsigned char *pGand = 0;
	float Gand = 0.f;
	pGy = reinterpret_cast<unsigned char *>(&Gy);
	pZero = reinterpret_cast<unsigned char *>(&zero);
	pGand = reinterpret_cast<unsigned char *>(&Gand);

	ultoa(*(unsigned long*)(void*)&Gy, s, 2);
	printf("Gy= ");
	for (int j = 0; j<4 * 8 - strlen(s); j++) printf("0");
	printf("%s\n", s);
	ultoa(*(unsigned long*)(void*)&Gx, s, 2);
	printf("Gx= ");
	for (int j = 0; j<4 * 8 - strlen(s); j++) printf("0");
	printf("%s\n", s);
	ultoa(*(unsigned long*)(void*)&zero, s, 2);
	printf("zero= ");
	for (int j = 0; j<4 * 8 - strlen(s); j++) printf("0");
	printf("%s\n", s);
	ultoa(*(unsigned long*)(void*)&Gand, s, 2);
	printf("Gand= ");
	for (int j = 0; j<4 * 8 - strlen(s); j++) printf("0");
	printf("%s\n", s);


	for (int i = 0; i<4; i++){
		*pGand = (*pGy & *pZero);
		pGand++;
		pGy++;
		pZero++;
	};
	ultoa(*(unsigned long*)(void*)&Gand, s, 2);
	printf("Gand= ");
	for (int j = 0; j<4 * 8 - strlen(s); j++) printf("0");
	printf("%s\n", s);
	//bitwise XOR on floats°´Î»Òì»ò
	unsigned char *pGx = 0;
	unsigned char *pGxor = 0;
	float Gxor = 0;
	ultoa(*(unsigned long*)(void*)&Gxor, s, 2);
	printf("Gxor= ");
	for (int j = 0; j<4 * 8 - strlen(s); j++) printf("0");
	printf("%s\n", s);
	pGx = reinterpret_cast<unsigned char *>(&Gx);
	pGand = reinterpret_cast<unsigned char *>(&Gand);
	pGxor = reinterpret_cast<unsigned char *>(&Gxor);
	for (int i = 0; i<4; i++){
		*pGxor = (*pGx ^ *pGand);
		pGxor++;
		pGx++;
		pGand++;
	};
	ultoa(*(unsigned long*)(void*)&Gxor, s, 2);
	printf("Gxor= ");
	for (int j = 0; j<4 * 8 - strlen(s); j++) printf("0");
	printf("%s\n", s);
	printf("Gxor= %f\n", Gxor);
	system("pause");
}