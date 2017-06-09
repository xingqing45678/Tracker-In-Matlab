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
%%


% run_tracker.m

close all;
% clear all;

%choose the path to the videos (you'll be able to choose one with the GUI)
base_path = 'D:\Image Process\对于FHOG各个版本平台代码对比\DSST_code_MATLAB\code - 副本用于修改\sequences';

%parameters according to the paper
params.padding = 1.0;         			% extra area surrounding the target%超出目标的额外区域
params.output_sigma_factor = 1/16;		% standard deviation for the desired translation filter output%期望位置滤波器输出的标准差
params.scale_sigma_factor = 1/4;        % standard deviation for the desired scale filter output%期望尺度滤波器输出的标准差
params.lambda = 1e-2;					% regularization weight (denoted "lambda" in the paper)%正则化权重（文章中‘lamda’）
params.learning_rate = 0.025;			% tracking model learning rate (denoted "eta" in the paper)%跟踪模型学习率（文章中“eta”）
params.number_of_scales = 33;           % number of scale levels (denoted "S" in the paper)%尺度等级数（记为“S”）
params.scale_step = 1.02;               % Scale increment factor (denoted "a" in the paper)%尺度增量因子（表示“a”）
params.scale_model_max_area = 512;      % the maximum size of scale examples%最大尺度的样例

params.visualization = 1;

%ask the user for the video
video_path = choose_video(base_path);
if isempty(video_path), return, end  %user cancelled
[img_files, pos, target_sz, ground_truth, video_path] = ...
	load_video_info(video_path);

params.init_pos = floor(pos) + floor(target_sz/2);
params.wsize = floor(target_sz);
params.img_files = img_files;
params.video_path = video_path;

[positions, fps] = dsst(params);%目标跟踪

% calculate precisions%计算精度
[distance_precision, PASCAL_precision, average_center_location_error] = ...
    compute_performance_measures(positions, ground_truth);

fprintf('Center Location Error: %.3g pixels\nDistance Precision: %.3g %%\nOverlap Precision: %.3g %%\nSpeed: %.3g fps\n', ...
    average_center_location_error, 100*distance_precision, 100*PASCAL_precision, fps);