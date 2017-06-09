function [positions, fps] = dsst(params)

% [positions, fps] = dsst(params)

% parameters
padding = params.padding;                         	%extra area surrounding the target%超出目标的额外区域
output_sigma_factor = params.output_sigma_factor;	%spatial bandwidth (proportional to target)%空间带宽（与目标成比例）%期望位置滤波器输出的标准差
lambda = params.lambda;                             % regularization weight (denoted "lambda" in the paper)%正则化权重（文章中‘lamda’）
learning_rate = params.learning_rate;               % tracking model learning rate (denoted "eta" in the paper)%跟踪模型学习率（文章中“eta”）
nScales = params.number_of_scales;                  % number of scale levels (denoted "S" in the paper)%尺度等级数（记为“S”）
scale_step = params.scale_step;                     % Scale increment factor (denoted "a" in the paper)%尺度增量因子（表示“a”）
scale_sigma_factor = params.scale_sigma_factor;     % standard deviation for the desired scale filter output%期望尺度滤波器输出的标准差
scale_model_max_area = params.scale_model_max_area; % the maximum size of scale examples%最大尺度的样例

video_path = params.video_path; %视频路径
img_files = params.img_files;   %文件路径
pos = floor(params.init_pos);   %初始化目标位置
target_sz = floor(params.target_sz);%目标大小

visualization = params.visualization;%可视化参数（等于1）

num_frames = numel(img_files);%p.totalframe;%图片帧的数量

init_target_sz = target_sz;%目标初始大小

% target size att scale = 1
base_target_sz = target_sz;

% window size, taking padding into account
sz = floor(base_target_sz * (1 + padding));

% desired translation filter output (gaussian shaped), bandwidth
% proportional to target size
output_sigma = sqrt(prod(base_target_sz)) * output_sigma_factor;
[rs, cs] = ndgrid((1:sz(1)) - floor(sz(1)/2), (1:sz(2)) - floor(sz(2)/2));
y = exp(-0.5 * (((rs.^2 + cs.^2) / output_sigma^2)));
yf = single(fft2(y));


% desired scale filter output (gaussian shaped), bandwidth proportional to
% number of scales
scale_sigma = nScales/sqrt(33) * scale_sigma_factor;
ss = (1:nScales) - ceil(nScales/2);
ys = exp(-0.5 * (ss.^2) / scale_sigma^2);
ysf = single(fft(ys));

% store pre-computed translation filter cosine window
cos_window = single(hann(sz(1)) * hann(sz(2))');

% store pre-computed scale filter cosine window
if mod(nScales,2) == 0
    scale_window = single(hann(nScales+1));
    scale_window = scale_window(2:end);
else
    scale_window = single(hann(nScales));
end;

% scale factors
ss = 1:nScales;
scaleFactors = scale_step.^(ceil(nScales/2) - ss);

% compute the resize dimensions used for feature extraction in the scale
% estimation
scale_model_factor = 1;
if prod(init_target_sz) > scale_model_max_area
    scale_model_factor = sqrt(scale_model_max_area/prod(init_target_sz));
end
scale_model_sz = floor(init_target_sz * scale_model_factor);

currentScaleFactor = 1;

% to calculate precision
positions = zeros(numel(img_files), 4);

% to calculate FPS
time = 0;

% find maximum and minimum scales
im = imread([video_path img_files{1}]); %%% im = loding_cap(p.fp,p.ImgHeight,p.ImgWidth);
min_scale_factor = scale_step ^ ceil(log(max(5 ./ sz)) / log(scale_step));
max_scale_factor = scale_step ^ floor(log(min([size(im,1) size(im,2)] ./ base_target_sz)) / log(scale_step));

for frame = 1:num_frames,
    %load image
    im = imread([video_path img_files{frame}]); %%% im = loding_cap(p.fp,p.ImgHeight,p.ImgWidth);

    tic;
    
    if frame > 1
        
        %提取特征测试输入F
         %样本中每个像素点计算28维融合特征(1维原始灰度+27维fhog)
         %乘以二维hann后作为输入F
         %用于位置相关滤波器
        % extract the test sample feature map for the translation filter
        xt = get_translation_sample(im, pos, sz, currentScaleFactor, cos_window);
        
        %计算响应值y=F-1{(A*Z)/(B+lambda)}
        % calculate the correlation response of the translation filter
        xtf = fft2(xt);
        response = real(ifft2(sum(hf_num .* xtf, 3) ./ (hf_den + lambda)));
        
        %找到max(y)得到目标新位置
        % find the maximum translation response
        [row, col] = find(response == max(response(:)), 1);
        
        % 更新目标位置
        % update the position
        pos = pos + round((-sz/2 + [row, col]) * currentScaleFactor);
        
        %把每个样本resize成固定大小，分别提取31维fhog特征，每个样本的所有fhog再
        %串联成一个特征向量构成33层金字塔特征，乘以一维hann窗后作为输入F
        % 用于尺度相关滤波器
        % extract the test sample feature map for the scale filter
        xs = get_scale_sample(im, pos, base_target_sz, currentScaleFactor * scaleFactors, scale_window, scale_model_sz);
        
        %得到尺度变换的响应最大值y=F-1{(A*Z)/(B+lambda)}
        % calculate the correlation response of the scale filter
        xsf = fft(xs,[],2);
        scale_response = real(ifft(sum(sf_num .* xsf, 1) ./ (sf_den + lambda)));
        
        %找到max(y)得到当前的尺度
        % find the maximum scale response
        recovered_scale = find(scale_response == max(scale_response(:)), 1);
        
         % 更新当前尺度
        % update the scale
        currentScaleFactor = currentScaleFactor * scaleFactors(recovered_scale);
        if currentScaleFactor < min_scale_factor
            currentScaleFactor = min_scale_factor;
        elseif currentScaleFactor > max_scale_factor
            currentScaleFactor = max_scale_factor;
        end
    end
    %提取特征训练样本输入X
    %样本中每个像素点计算28维融合特征(1维原始灰度+27维fhog)
    %乘以二维hann后作为输入X
    %提取特征用于位置相关滤波器
    % extract the training sample feature map for the translation filter
    xl = get_translation_sample(im, pos, sz, currentScaleFactor, cos_window);
    
    %获取分子A=GF;分母B=F*F;此时没有lambda
    % calculate the translation filter update
    xlf = fft2(xl);
    new_hf_num = bsxfun(@times, yf, conj(xlf));
    new_hf_den = sum(xlf .* conj(xlf), 3);
    
    %把每个样本resize成固定大小，分别提取31维fhog特征，每个样本的所有fhog再
    %串联成一个特征向量构成33层金字塔特征，乘以一维hann窗后作为输入X
    % 提取特征用于尺度相关滤波器
    % extract the training sample feature map for the scale filter
    xs = get_scale_sample(im, pos, base_target_sz, currentScaleFactor * scaleFactors, scale_window, scale_model_sz);
    
    %同样的获取分子A=GF;分母B=F*F;此时没有lambda
    % calculate the scale filter update
    xsf = fft(xs,[],2);
    new_sf_num = bsxfun(@times, ysf, conj(xsf));
    new_sf_den = sum(xsf .* conj(xsf), 1);
    
    
    if frame == 1
        % first frame, train with a single image
        hf_den = new_hf_den;
        hf_num = new_hf_num;
        
        sf_den = new_sf_den;
        sf_num = new_sf_num;
    else
        % subsequent frames, update the model
        hf_den = (1 - learning_rate) * hf_den + learning_rate * new_hf_den;
        hf_num = (1 - learning_rate) * hf_num + learning_rate * new_hf_num;
        sf_den = (1 - learning_rate) * sf_den + learning_rate * new_sf_den;
        sf_num = (1 - learning_rate) * sf_num + learning_rate * new_sf_num;
    end
    
    % calculate the new target size
    target_sz = floor(base_target_sz * currentScaleFactor);
    
    %save position
    positions(frame,:) = [pos target_sz];
    
    time = time + toc;
    
    
    %visualization
    if visualization == 1
        rect_position = [pos([2,1]) - target_sz([2,1])/2, target_sz([2,1])];
        if frame == 1,  %first frame, create GUI
            figure('Number','off', 'Name',['Tracker - ' video_path]);
            im_handle = imshow(im);%im_handle = imshow(uint8(im), 'Border','tight', 'InitialMag', 100 + 100 * (length(im) < 500));
            title('DSST');
            rect_handle = rectangle('Position',rect_position, 'EdgeColor','w');
            text_handle = text(5, 18, strcat('#',num2str(frame)), 'Color','w', 'FontWeight','bold', 'FontSize',15);
%             set(text_handle, 'color', [0 1 1]);
        else
            try  %subsequent frames, update GUI
                set(im_handle, 'CData', im)
                set(rect_handle, 'Position', rect_position)
                set(text_handle, 'string', strcat('#',num2str(frame)))
            catch
                return
            end
        end
        
        drawnow
%         pause
    end
end

fps = num_frames/time;