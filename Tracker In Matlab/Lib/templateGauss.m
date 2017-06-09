%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%作者：秦威
%E-mail：285980893@qq.com
%子程序功能是产生高斯理想响应模板
%（sz模板的像素大小,im原图像）
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STC算法的模板生成
% alapha = 2.25;                    %parmeter \alpha in Eq.(6)
% [rs, cs] = ndgrid((1:sz(1)) - floor(sz(1)/2), (1:sz(2)) - floor(sz(2)/2));
% dist = rs.^2 + cs.^2;
% conf = exp(-0.5 / (alapha) * sqrt(dist));%confidence map function Eq.(6)
% conf = conf/sum(sum(conf));% normalization
% conff = fft2(conf); %transform conf to frequencey domain
%% KFC算法的模板生成
% %window size, taking padding into account
% sz = floor(target_sz * (1 + padding));
% %desired output (gaussian shaped), bandwidth proportional to target size
% output_sigma = sqrt(prod(target_sz)) * output_sigma_factor;
% [rs, cs] = ndgrid((1:sz(1)) - floor(sz(1)/2), (1:sz(2)) - floor(sz(2)/2));
% y = exp(-0.5 / output_sigma^2 * (rs.^2 + cs.^2));
% yf = fft2(y);
%% CF算法基本高斯模板
function F_response=templateGauss(sz,im)
    [rs, cs] = ndgrid((1:sz(1)) - floor(sz(1)/2), (1:sz(2)) - floor(sz(2)/2));
    dist = rs.^2 + cs.^2;
    conf = exp(-0.5 / (2.25) * sqrt(dist));%生成二维高斯分布
    conf = conf/sum(sum(conf));% normalization
    if(size(im,3)==1)%灰度图像
        response=conf;
    else            %彩色图像
        response(:,:,1)=conf;response(:,:,2)=conf;response(:,:,3)=conf;    
    end       
%         figure
%         imshow(256.*response);
%         mesh(response);
        F_response=fft2(response);%傅里叶变换
end