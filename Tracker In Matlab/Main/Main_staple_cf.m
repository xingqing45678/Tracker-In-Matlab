%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                  跟踪主函数StapleOnlyCf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    close all;clear all;clc;
% 添加文件路径
%     addpath(genpath(pwd));%添加当前文件夹及子文件夹
    addpath('../Lib','../StapleOnlyCF');%添加上一级目录的文件夹1,2,n
%     rmpath('./上级目录中的文件夹1');%去除路径是为了修改文件名等操作，否则Matlab会认为你要改的路径正在使用中，是禁止操作的
    sequence = 'D:\ImageData\Toy';%设置图片路径
    [params,im] = Load_image(sequence);%读取源文件picture %%% 
%     [paraams,Image16] = InitCap16(sequence);%读取源文件cap16
%% staple-only-cf
    [params, bg_area, fg_area, area_resize_factor] = initializeAllAreas(im, params);
	% 开始跟踪主程序
    result = trackerMain_StapleOnlyCf(params, im, bg_area, area_resize_factor);
    %% the end 显示误差
    fclose('all');
    show_precision(result.pos, params.bb_VOT, sequence);%计算误差值（groundtruth完备情况下）
    disp('程序已正常结束');