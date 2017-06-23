
% This demo script runs the ECO tracker with deep features on the
% included "Crossing" video.
close all;clearvars;clc;
% Add paths
setup_paths();
startframe = 1;
sub_flag = 0;%  '1'--取分割中心图像     '0'--不分割图像

% 读取视频文件
video_path = 'F:\testwot\worldoftanks 2017-06-23 12-45-32-856.avi';%视频
[seq] = load_video_info_qw(video_path,startframe,sub_flag);%视频
% % Run ECO
results = testing_ECO(seq);%读取视频

% 读取图片文件
% video_path = 'F:\testwot\121B_623_1317\';%图片
% [seq] = load_video_info_qw_picture(video_path,startframe);%读图片，用读取groundtruth文件来获取初始定位
% % Run ECO
% results = testing_ECO_picture(seq);%读取图片