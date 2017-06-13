%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                  跟踪主函数StapleOnlyCf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;clearvars;clc;
% 添加文件路径
%path to the videos (you'll be able to choose one with the GUI).
addpath('../Lib','../KCF DCF');%添加上一级目录的文件夹1,2,n
sequence = 'D:\ImageData\box';%'./data/Benchmark/';
%
[params,im] = Load_image(sequence);%读取源文件picture %%% 
%     [paraams,Image16] = InitCap16(sequence);%读取源文件cap16
params.video_path = [sequence '/img/'];
kernel.type = params.kernel_type;
features.gray = false;
features.hog = false;

%% 
switch params.feature_type
case 'gray'
    params.interp_factor = 0.075;  %linear interpolation factor for adaptation
    kernel.sigma = 0.2;  %gaussian kernel bandwidth	
    kernel.poly_a = 1;  %polynomial kernel additive term
    kernel.poly_b = 7;  %polynomial kernel exponent	
    features.gray = true;
    params.cell_size = 1;		
case 'hog'
    params.interp_factor = 0.02;		
    kernel.sigma = 0.5;		
    kernel.poly_a = 1;
    kernel.poly_b = 9;		
    features.hog = true;
    features.hog_orientations = 9;
    params.cell_size = 4;		
otherwise
    error('Unknown feature.')
end
assert(any(strcmp(params.kernel_type, {'linear', 'polynomial', 'gaussian'})), 'Unknown kernel.')
%% 
    %call tracker function with all the relevant parameters
    [positions, time] = tracker(params.video_path, params.img_files, params.init_pos, params.target_sz, ...
                                params.padding, kernel, params.lambda, params.output_sigma_factor, params.interp_factor, ...
                                params.cell_size, features, params.show_visualization);
%% 
%     %calculate and show precision plot, as well as frames-per-second
%     precisions = precision_plot(positions, ground_truth, video, show_plots);
%     fps = numel(img_files) / time;
% 
%     fprintf('%12s - Precision (20px):% 1.3f, FPS:% 4.2f\n', video, precisions(20), fps)
% 
%     if nargout > 0
%         %return precisions at a 20 pixels threshold
%         precision = precisions(20);
%     end
