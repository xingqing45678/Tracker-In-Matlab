%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%作者：qw
%E-mail：1406820937@qq.com
%子程序功能是读取图片帧，并且读取groundtruth数据
%（ground_truth目标真实位置,img_path图片文件夹,img_files图片文件名存储,imgDir文件路径,）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [params,im]=Load_image(imgDir)
%     %% Read params.txt
    params = readParams('params.txt');
	%% load video info
    sequence_path = [imgDir,'/'];%文件路径
    img_path = [sequence_path 'img/'];
    start_frame = 1;
    %% Read files 
    params.bb_VOT = csvread([sequence_path 'groundtruth_rect.txt']);%序列中真实目标位置
    if size(params.bb_VOT,2)==1
        region = params.bb_VOT(1:4);
    else
        region = params.bb_VOT(start_frame,:);%读取groundtruth的第一行4个数据
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % 读取所有图像帧
    dir_content = dir([sequence_path 'img/']);
    % skip '.' and '..' from the count
    n_imgs = length(dir_content) - 2;
    img_files = cell(n_imgs, 1);
    for ii = 1:n_imgs
        img_files{ii} = dir_content(ii+2).name;
    end
    img_files(1:start_frame-1)=[];
    im = imread([img_path img_files{start_frame}]);%将第一幅图像读入
    % 判断是否灰度图像 ?
    if(size(im,3)==1)
        params.grayscale_sequence = true;
    end
    %% get position and boxsize 读取groundtruth数据 
    params.img_files = img_files;%将图像名赋给参数表params
    params.img_path = img_path;%将图像路径赋给参数表params
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(numel(region)==4)%返回图像像素或数组元素个数
        x = region(1);y = region(2);w = region(3);h = region(4);cx = x+w/2;cy = y+h/2;%
    else
        waring('wrong with groundtruth');
    end
    % init_pos 是范围框的中心
    params.init_pos = [cy cx];
    params.target_sz = round([h w]);%四舍五入往最近的整数靠拢
    params.totalframe = numel(params.img_files);
%     params.imgheight = h;
%     params.imgwidth = w;
%     params.fp = fp;
%     im = imread([img_path img_files{1}]);%读取目标帧
%     im= rgb2gray(im);%转换为灰度图
%     imshow(im);
end