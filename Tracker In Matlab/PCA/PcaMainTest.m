clc;close all;clear all;
addpath('..\Lib\Fhoglib');
im = imread('suv_1.jpg');
im= rgb2gray(im);%转换为灰度图
im=single(im);
binSize = 4; nOrients = 9; softBin = -1; useHog=2;clip = 0.2;
    [M,O]=gradientMex('gradientMag',im,0,1);%
    H = gradientMex('gradientHist',M,O,binSize,nOrients,softBin,useHog,clip);%
[M_qw,O_qw]=gradientMag(im,0,1);%new by qw
H_qw = fhog_features(im,M,O,binSize,nOrients,softBin,useHog,clip);%new by qw
[w,h,p] = size(H_qw);
feature = zeros(w*h,p);
for i=1:w
    for j=1:h
        for k=1:p
            feature(i*j,k) = H_qw(i,j,k);
        end
    end
end
[pc,score,latent,tsquare] = pca(feature);%我们这里需要他的pc和latent值做分析
a = cumsum(latent)./sum(latent)
feature_after_PCA=score(:,1:5);
pareto(100*latent);%调用matla画图
% axis([0 20 0 23]);

% %%
% %清屏
% clear
% %初始化数据
% x=[-14.8271317103068,-3.00108550936016,1.52090778549498,3.95534842970601;-16.2288612441648,-2.80187433749996,-0.410815700402130,1.47546694457079;-15.1242838039605,-2.59871263957451,-0.359965674446737,1.34583763509479;-15.7031424565913,-2.53005662064257,0.255003254103276,-0.179334985754377;-17.7892158910100,-3.32842422986555,0.255791146332054,1.65118282449042;-17.8126324036279,-4.09719527953407,-0.879821957489877,-0.196675865428539;-14.9958877514765,-3.90753364293621,-0.418298866141441,-0.278063876667954;-15.5246706309866,-2.08905845264568,-1.16425848541704,-1.16976057326753;];
% %调用princomp函数
% [coef,score,latent,t2] = princomp(x);
% score
% %测试score是否和score_test一样
% score_test=bsxfun(@minus,x,mean(x,1))*coef;
% score_test
% latent=100*latent/sum(latent);%将latent总和统一为100，便于观察贡献率
% pareto(latent);%调用matla画图
