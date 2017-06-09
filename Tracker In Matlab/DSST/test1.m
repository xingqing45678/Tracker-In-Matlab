clc;close all;clear all;
im = imread('suv_1.jpg');
im= rgb2gray(im);%×ª»»Îª»Ò¶ÈÍ¼
im=single(im);
softBin = -1; useHog = 2; binSize=4; clip=.2; nOrients=9;

[M,O]=gradientMex('gradientMag',im,0,1);
H = gradientMex('gradientHist',M,O,binSize,nOrients,softBin,useHog,clip);
