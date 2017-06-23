
% This demo script runs the ECO tracker
% qw 2017-6-23

close all;clearvars;clc;
% Add paths
setup_paths();
startframe = 1;

choose = 'fhog';%   'fhog'--手动设计特征hand-crafted      'cnn'--深度学习特征deep features
orign_file = 1;%   '1'--源文件为视频     '0'--源文件为图片  


%% 读取视频文件
if orign_file==1
    video_path = 'F:\testwot\worldoftanks 2017-06-23 12-45-32-856.avi';%视频
    [seq] = load_video_info_qw(video_path,startframe);%视频
else
%% 读取图片文件
    video_path = 'F:\testwot\121B_623_1317';%图片
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
