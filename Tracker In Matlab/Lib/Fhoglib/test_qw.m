clc;close all;clear all;
im = imread('suv_2.jpg');
im= rgb2gray(im);%×ª»»Îª»Ò¶ÈÍ¼
im=single(im);
binSize = 4; nOrients = 9; softBin = -1; useHog=2;clip = 0.2;

[M,O]=gradienQW('gradientMag',im,0,1);%
H = gradientQW('gradientHist',M,O,binSize,nOrients,softBin,useHog,clip);%

% [M_qw,O_qw]=gradientMag(im,0,1);%new by qw
% H_qw = fhog_features(im,M,O,binSize,nOrients,softBin,useHog,clip);%new by qw
