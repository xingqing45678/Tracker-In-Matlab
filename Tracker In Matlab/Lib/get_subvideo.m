%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%作者：qw
%E-mail：1406820937@qq.com
%程序功能是将视频一部分子窗口截取，并保存为视频文件
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;clearvars;clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Video_path = 'F:\testwot\121B_623_1317\';%设置路径名称
videoFile = [Video_path 'worldoftanks 2017-06-23 12-49-13-324.avi'];%视频文件名
videoData = VideoReader(videoFile);%读取视频文件
start_frame = 1;%起始帧
end_frame = 100;%videoData.NumberOfFrames;%结束帧CurrentTime
pos = [videoData.Height/2,videoData.Width/2];

out_videoName = 'ourtput_video';
writerObj=VideoWriter(out_videoName);  %// 定义一个视频文件用来存动画
writerObj.Quality=100;
writerObj.FrameRate=60;
open(writerObj);                    %// 打开该视频文件
figure
for frames = start_frame : 1 : end_frame
    videoframe = read(videoData,frames);
    im = get_subimg(videoframe,pos);
    imshow(im);
    writeVideo(writerObj,im);
    disp(frames);
end 
close(writerObj); %// 关闭视频文件句柄
disp('视频截取大小完成');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%