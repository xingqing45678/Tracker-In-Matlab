%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                  跟踪主函数DSST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;clear all;clc;
addpath('../MaltabLibrary_qw','../MaltabLibrary_qw/Fhog_qw','../DSST');%添加上一级目录的文件夹1,2,n
sequence = 'D:\ImageData\Coke';%设置图片路径
[params,im] = Load_image(sequence);%读取源文件picture %%% 
%     [paraams,Image16] = InitCap16(sequence);%读取源文件cap16
%ask the user for the video
params.video_path = [sequence '/img/'];
% ground_truth = [params.init_pos,params.target_sz];

[positions, fps] = dsst(params);%目标跟踪

%% 计算精度
[distance_precision, PASCAL_precision, average_center_location_error] = ...
    compute_performance_measures(positions, ground_truth);

fprintf('Center Location Error: %.3g pixels\nDistance Precision: %.3g %%\nOverlap Precision: %.3g %%\nSpeed: %.3g fps\n', ...
    average_center_location_error, 100*distance_precision, 100*PASCAL_precision, fps);
%%
%%
% 算法流程
% 
% Input: 
% 输入图像patch It 
% 上一帧的位置Pt?1和尺度St?1 
% 位置模型Atranst?1、Btanst?1和尺度模型Ascalet?1、Bscalet?1 
% Output: 
% 估计的目标位置Pt和尺度St 
% 更新位置Atranst、Btranst和尺度模型Ascalet、Bscalet
% 其中, 
% 位置评估： 
% 1.参照模板在前一帧的位置，在当前帧中按照前一帧目标尺度的2倍大小提取一个样本Ztrans 
% 2.利用Ztrans和Atranst?1、Btanst?1，根据公式(6)计算ytrans 
% 3.计算max(ytrans)，得到目标新的位置Pt 
% 尺度评估： 
% 4.以目标当前新位置为中心，提取33种不同尺度的样本Ztrans 
% 5.利用Ztrans和Atranst?1、Btanst?公式(6)计算出yscale 
% 6.计算max(yscale)，得到目标准确的尺度St
% 模型更新： 
% 7.提取样本ftrans和fscale 
% 8.更新位置模型Atranst和Btranst 
% 9.更新尺度模型Ascalet和Bscalet