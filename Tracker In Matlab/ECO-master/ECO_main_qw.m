
% This demo script runs the ECO tracker
% qw 2017-6-23

close all;clearvars;clc;
% Add paths
setup_paths();

startframe = 1;         %  起始帧
choose     = 'fhog';    %  'fhog'--普通特征       'cnn'--深度学习特征
orign_file = 1;         %  '1 '--源文件为视频     '0 '--源文件为图片  
sub_flag   = 0;         %  '1 '--取分割中心图像   '0 '--不分割图像

%% 读取视频文件
if orign_file==1
    video_path = 'F:\testwot\WorldOfTanks 2017-06-21 13-00-30-988.avi';%视频
    [seq] = load_video_info_qw(video_path,startframe,sub_flag);%视频
else
%% 读取图片文件
    video_path = 'F:\testwot\59huangmo_2';%图片
    [seq] = load_video_info_qw_picture(video_path,startframe);%读图片，用读取groundtruth文件来获取初始定位
end
%% Run ECO

switch choose
    case 'fhog'
        % Run ECO
        if orign_file==1
            results = testing_ECO_HC(seq);%读取视频
        else
            results = testing_ECO_HC_picture(seq);%读取图片
        end
    case 'cnn'
        if orign_file==1
            results = testing_ECO(seq);%读取视频
        else
            results = testing_ECO_picture(seq);%读取图片
        end
end
%%
