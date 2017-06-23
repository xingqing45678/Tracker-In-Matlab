%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%作者：qw
%E-mail：1406820937@qq.com
%程序功能是将视频文件提取为图片序列，或者相反将图片序列转为视频
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;clearvars;clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 视频转化为图片
Video_path = 'F:\testwot\121B_623_1317\';%设置路径名称
videoFile = [Video_path 'worldoftanks 2017-06-23 12-49-13-324.avi'];%视频文件名
img_path = [Video_path 'img'];%图片保存路径
videoData = VideoReader(videoFile);%读取视频文件
start_frame = 1;%起始帧
end_frame = videoData.NumberOfFrames;%结束帧CurrentTime
pos = [videoData.Height/2,videoData.Width/2];
figure
for i = start_frame : 1 : end_frame
    videoframe = read(videoData,i);
    im = get_subimg(videoframe,pos);
    imshow(im);
%     text(5, 18, strcat('#',num2str(i)), 'Color','y', 'FontWeight','bold', 'FontSize',20);
%     imwrite(im,fullfile(img_path,[num2str(i,'%06d') '.jpg']));
    disp(i);
end
disp('视频转换为图片完成');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% 图片转化为视频
% %framesPath :图像序列所在路径，同时要保证图像大小相同
% framesPath = 'D:\ImageData\David_qw\img\';
% %videoName:  表示将要创建的视频文件的名字
% videoName = 'David_video';
% %quality:    生成视频的质量 0-100
% quality = 100;
% %Compressed: 压缩类型， 'Indeo3'（默认）, 'Indeo5', 'Cinepak', 'MSVC', 'RLE' or 'None'
% % Compressed = 'Motion JPEG AVI';
% %FrameRate: 帧率
% FrameRate = 50;
% if(exist('videoName','file'))
%     delete videoName.avi
% end
% %生成视频的参数设定
% aviobj = VideoWriter(videoName);  %创建一个avi视频文件对象，开始时其为空
% aviobj.Quality=quality;
% aviobj.FrameRate=FrameRate;
% % aviobj.VideoCompressionMethod=Compressed;
% %读入图片
% % 读取所有图像帧
% dir_content = dir(framesPath);
% % skip '.' and '..' from the count
% n_imgs = length(dir_content) - 2;
% img_files = cell(n_imgs, 1);
% for ii = 1:n_imgs
%     img_files{ii} = dir_content(ii+2).name;%imag_files存储所有图像帧文件名
% end
% %startFrame ,endFrame ;表示从哪一帧开始，哪一帧结束
% startFrame = 1;
% endFrame = n_imgs;
% open(aviobj);
% for i=startFrame:n_imgs   %
%     frames=imread([framesPath,img_files{i}]);
%     writeVideo(aviobj,frames);
%     disp(i);
% end
% close(aviobj); % 关闭创建视频
% disp('图片转换为视频完成');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%