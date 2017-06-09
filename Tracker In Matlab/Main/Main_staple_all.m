%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                  跟踪主函数StapleAll
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    close all;clear all;clc;
% 添加文件路径
%     addpath(genpath(pwd));%添加当前文件夹及子文件夹
    addpath('../Lib','../StapleAll');%添加上一级目录的文件夹1,2,n
%     rmpath('./上级目录中的文件夹1');%去除路径是为了修改文件名等操作，否则Matlab会认为你要改的路径正在使用中，是禁止操作的
    sequence = 'D:\ImageData\Toy';%设置图片路径
    [params,im] = Load_image(sequence);%读取源文件picture %%% 
%     [paraams,Image16] = InitCap16(sequence);%读取源文件cap16
	if params.visualization%为1时画图，0时不画图
		params.videoPlayer = vision.VideoPlayer('Position', [100 100 [size(im,2), size(im,1)]+30]);
	end
%% staple-all
    [params, bg_area, fg_area, area_resize_factor] = initializeAllAreas(im, params);
	% 开始跟踪主程序
    result = trackerMain(params, im, bg_area, fg_area, area_resize_factor);
    %% the end 显示误差
    fclose('all');
%     show_precision(result.pos, params.bb_VOT, sequence);%计算误差值（groundtruth完备情况下）
    disp('程序已正常结束');