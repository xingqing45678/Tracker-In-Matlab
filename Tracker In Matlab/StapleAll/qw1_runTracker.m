function qw1_runTracker()
    close all;clear all;clc;
% RUN_TRACKER  is the external function of the tracker - does initialization and calls trackerMain
    start_frame=1;
    sequence='D:\ImageData\Toy';%设置路径名称
    %% 读取文件 params.txt
    params = readParams('params.txt');%读取文件，初始化参数
	%% 加载视频信息
% 	sequence_path = ['../Sequences/',sequence,'/'];
    sequence_path = [sequence,'/'];
    img_path = [sequence_path 'img/'];%设置图像路径
    %% 读取文件
    text_files = dir([sequence_path '*_frames.txt']);
    f = fopen([sequence_path text_files(1).name]);
    frames = textscan(f, '%f,%f');
    if exist('start_frame')%判断变量或者函数是否存在
        frames{1} = start_frame;
    else
        frames{1} = 1;
    end
    
    fclose(f);
    
    params.bb_VOT = csvread([sequence_path 'groundtruth_rect.txt']);
    region = params.bb_VOT(frames{1},:);%读取groundtruth的第一行8个数据
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % 读取所有图像帧
    dir_content = dir([sequence_path 'img/']);
    % skip '.' and '..' from the count
    n_imgs = length(dir_content) - 2;
    img_files = cell(n_imgs, 1);
    for ii = 1:n_imgs
        img_files{ii} = dir_content(ii+2).name;%imag_files存储所有图像帧文件名
    end
       
    img_files(1:start_frame-1)=[];

    im = imread([img_path img_files{1}]);%将第一幅图像读入
    % 判断是否灰度图像 ?
    if(size(im,3)==1)
        params.grayscale_sequence = true;
    end

    params.img_files = img_files;%将图像名赋给参数表params
    params.img_path = img_path;%将图像路径赋给参数表params

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(numel(region)==8)%返回图像像素或数组元素个数
        % polygon format
        [cx, cy, w, h] = getAxisAlignedBB(region);%从真实数据中提取轴对齐边界框作为旋转图像
    else
        x = region(1);
        y = region(2);
        w = region(3);
        h = region(4);
        cx = x+w/2;
        cy = y+h/2;
    end

    % init_pos 是范围框的中心（is the centre of the initial bounding box）
    params.init_pos = [cy cx];
    params.target_sz = round([h w]);%四舍五入往最近的整数靠拢
    [params, bg_area, fg_area, area_resize_factor] = initializeAllAreas(im, params);
	if params.visualization%为1时画图，0时不画图
		params.videoPlayer = vision.VideoPlayer('Position', [100 100 [size(im,2), size(im,1)]+30]);
	end
    % in runTracker we do not output anything
	params.fout = -1;
	% start the actual tracking
	trackerMain(params, im, bg_area, fg_area, area_resize_factor);
    fclose('all');
end
