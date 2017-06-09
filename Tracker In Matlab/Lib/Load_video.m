%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%作者：秦威
%E-mail：285980893@qq.com
%子程序功能是读取视频，读取groundtruth数据
%（videoData存储视频的结构体,ground_truth目标真实位置,videoDir视频文件路径,）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [videoData,ground_truth] = Load_video(videoDir)
    sequence_path = [videoDir,'\'];%文件路径
    videoFile = [sequence_path 'david.mpg']; %文件路径
    videoData = VideoReader(videoFile);%读取视频
% % for i = 1 : 50
% %         videoframe = read(videoData,i);% 读取每一帧
% %         videoframe= rgb2gray(videoframe);%转换为灰度图
% %         imshow(videoframe);%显示每一帧
% %         % imwrite(frame,strcat(num2str(i),'.jpg'),'jpg');% 保存每一帧
% % end
    %% Read files 
    ground_rect = csvread([sequence_path 'groundtruth_rect.txt']);%序列中真实目标位置
    %% get position and boxsize 读取groundtruth数据 
    if(size(ground_rect,2)==1)%一列
        error('please add "," in groundtruth');%x,y,w,h目标框大小
    else if(size(ground_rect,2)==4)%4列
        ground_truth=ground_rect;%x,y,w,h目标框大小
    else
        error('something wrong in groundtruth');
        end
    end
end