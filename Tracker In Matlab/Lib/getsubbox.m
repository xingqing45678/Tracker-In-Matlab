%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%作者：秦威
%E-mail：285980893@qq.com
%子程序功能是从原图像中获取感兴趣的目标区域（子框），并加窗处理减小边缘效应
%（pos目标中心位置，target_sz感兴趣子区域的大小，im原图像，windows加窗处理）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% 加窗方式1
% [rs, cs] = ndgrid((1:target_sz(1)) - floor(target_sz(1)/2), (1:target_sz(2)) - floor(target_sz(2)/2));
% dist = rs.^2 + cs.^2;
% hamming_window = hamming(sz(1)) * hann(sz(2))';
% sigma = mean(target_sz);% initial \sigma_1 for the weight function w_{\sigma} in Eq.(11)
% window = hamming_window.*exp(-0.5 / (sigma^2) *(dist));% use Hamming window to reduce frequency effect of image boundary
% window = window/sum(sum(window));%normalization
%% 加窗方式2
% window = hann(target_sz(1)) * hann(target_sz(2))';
%%
function target_box = getsubbox(pos,target_sz,im)
	%get and process the context region
	xs = floor(pos(2) + (1:target_sz(2)) - (target_sz(2)/2));
	ys = floor(pos(1) + (1:target_sz(1)) - (target_sz(1)/2));
	
	%check for out-of-bounds coordinates, and set them to the values at
	%the borders
	xs(xs < 1) = 1;
	ys(ys < 1) = 1;
	xs(xs > size(im,2)) = size(im,2);
	ys(ys > size(im,1)) = size(im,1);	
	%extract image in context region
	target_box = im(ys, xs, :);	
	%pre-process window
    target_box = double(target_box);
    target_box = (target_box-mean(target_box(:)));%normalization
%     target_box = window .* target_box;  %加窗处理
end