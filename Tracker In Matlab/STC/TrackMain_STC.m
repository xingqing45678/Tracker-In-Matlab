% Demo for paper "Fast Tracking via Spatio-Temporal Context Learning,"Kaihua Zhang, Lei Zhang,Ming-Hsuan Yang and David Zhang
% Paper can be available from http://arxiv.org/pdf/1311.1939v1.pdf 
% Implemented by Kaihua Zhang, Dept.of Computing, HK PolyU.
% Email: zhkhua@gmail.com
% Date: 11/24/2013
%%
clc;close all;clear all;
%%
fftw('planner','patient')
%% set path
imgDir='D:\ImageData\Coke\';%图片文件夹路径名
addpath([imgDir '\img']);
img_dir = dir([imgDir 'img\*.jpg']);%图片
%% load video info
%% Read files 
ground_rect = csvread([imgDir 'groundtruth_rect.txt']);%序列中真实目标位置
%% initialization get groundtruth and x y w h
if(size(ground_rect,2)==1)%一列
    error('please add "," in groundtruth');%x,y,w,h目标框大小
else if(size(ground_rect,2)==4)%4列
    ground_truth=ground_rect;%x,y,w,h目标框大小
else
    error('something wrong in groundtruth');
    end
end
%% set initial position and size
target_sz = [ground_truth(1,4), ground_truth(1,3)];
pos = [ground_truth(1,2), ground_truth(1,1)] + floor(target_sz/2);
%% parameters according to the paper
padding = 1;					%extra area surrounding the target
rho = 0.075;			        %the learning parameter \rho in Eq.(12)
sz = floor(target_sz * (1 + padding));% size of context region
time = 0;  %to calculate FPS
positions = zeros(numel(img_dir), 2);  %to calculate precision
%% parameters of scale update. See Eq.(15)
scale = 1;%initial scale ratio
lambda = 0.25;% \lambda in Eq.(15)
num = 5; % number of average frames
%% store pre-computed confidence map
alapha = 2.25;                    %parmeter \alpha in Eq.(6)
[rs, cs] = ndgrid((1:sz(1)) - floor(sz(1)/2), (1:sz(2)) - floor(sz(2)/2));
dist = rs.^2 + cs.^2;
conf = exp(-0.5 / (alapha) * sqrt(dist));%confidence map function Eq.(6)
conf = conf/sum(sum(conf));% normalization
conff = fft2(conf); %transform conf to frequencey domain
%% store pre-computed weight window
hamming_window = hamming(sz(1)) * hann(sz(2))';
sigma = mean(target_sz);% initial \sigma_1 for the weight function w_{\sigma} in Eq.(11)
window = hamming_window.*exp(-0.5 / (sigma^2) *(dist));% use Hamming window to reduce frequency effect of image boundary
window = window/sum(sum(window));%normalization
%%
maxconf=zeros(numel(img_dir)-1);
%%
for frame = 1:numel(img_dir),
    sigma = sigma*scale;% update scale in Eq.(15)
    window = hamming_window.*exp(-0.5 / (sigma^2) *(dist));%update weight function w_{sigma} in Eq.(11)
    window = window/sum(sum(window));%normalization
 	%load image
    img = imread(img_dir(frame).name);	    
	if size(img,3) > 1,
		im = rgb2gray(img);
    else
        im=img;
	end
   	contextprior = get_context(im, pos, sz, window);% the context prior model Eq.(4)
    tic()
    %%
    if frame > 1,
		%calculate response of the confidence map at all locations
	    confmap = real(ifft2(Hstcf.*fft2(contextprior))); %Eq.(11) 
       	%target location is at the maximum response
		[row, col] = find(confmap == max(confmap(:)), 1);
        pos = pos - sz/2 + [row, col]; 
        %%
        contextprior = get_context(im, pos, sz, window);
        conftmp = real(ifft2(Hstcf.*fft2(contextprior))); 
        maxconf(frame-1)=max(conftmp(:));
        %% update scale by Eq.(15)
        if (mod(frame,num+2)==0)
            scale_curr = 0;
            for kk=1:num
               scale_curr = scale_curr + sqrt(maxconf(frame-kk)/maxconf(frame-kk-1));
            end            
            scale = (1-lambda)*scale+lambda*(scale_curr/num);%update scale
        end  
        %%
    end	
    %% 
    %save position and calculate FPS
    positions(frame,:) = pos;
    time = time + toc();
	%% update the spatial context model h^{sc} in Eq.(9)
   	contextprior = get_context(im, pos, sz, window); 
    hscf = conff./(fft2(contextprior)+eps);% Note the hscf is the FFT of hsc in Eq.(9)
    %% update the spatio-temporal context model by Eq.(12)
    if frame == 1,  %first frame, initialize the spatio-temporal context model as the spatial context model
		Hstcf = hscf;
	else
		%update the spatio-temporal context model H^{stc} by Eq.(12)
		Hstcf = (1 - rho) * Hstcf + rho * hscf;% Hstcf is the FFT of Hstc in Eq.(12)
    end
    %% visualization
    target_sz([2,1]) = target_sz([2,1])*scale;% update object size
	rect_position = [pos([2,1]) - (target_sz([2,1])/2), (target_sz([2,1]))]; 
        if frame == 1,  %first frame, create GUI
            figure
            im_handle = imagesc(uint8(img));
            rect_handle = rectangle('Position',rect_position,'LineWidth',2,'EdgeColor','r');
            tex_handle = text(5, 18, strcat('#',num2str(frame)), 'Color','y', 'FontWeight','bold', 'FontSize',20);
            drawnow;
        else
            try  %subsequent frames, update GUI
                set(im_handle, 'CData', img)
                set(rect_handle, 'Position', rect_position)
                set(tex_handle, 'string', strcat('#',num2str(frame)))
%                 pause(0.001);
                drawnow;
            catch  % #ok, user has closed the window
                return
            end
        end
%     imagesc(uint8(img))
%     colormap(gray)
%     rectangle('Position',rect_position,'LineWidth',4,'EdgeColor','r');
%     hold on;
%     text(5, 18, strcat('#',num2str(frame)), 'Color','y', 'FontWeight','bold', 'FontSize',20);
%     set(gca,'position',[0 0 1 1]); 
%     pause(0.001); 
%     hold off;
%     drawnow;    
end
disp(['Frames-per-second: ' num2str(numel(img_dir) / time)])

%show the precisions plot
show_precision(positions, ground_truth, imgDir)