function [fig_handle] = run_imshow_qw(im_to_show,rect_position_vis,frame,fig_handle,...
                                        img_support_sz,currentScaleFactor,scaleFactors,scale_ind,det_sample_pos,scores_fs,output_sz)
                                        
    if size(im_to_show,3)>1
        im_to_show = rgb2gray(im_to_show);
    end
    
    if nargin <4  %first frame, create GUI
            fig_handle = figure('Name', 'Tracking');
%             set(fig_handle, 'Position', [100, 100, size(im,2), size(im,1)]);
%             imagesc(im_to_show);
            imshow(im_to_show);
            hold on;
            rectangle('Position',rect_position_vis, 'EdgeColor','w', 'LineWidth',1);
            text(10, 15, int2str(frame), 'Color','w', 'FontSize',15);
            hold off;
            axis off;axis image;set(gca, 'Units', 'normalized', 'Position', [0 0 1 1])
            
        else
            % Do visualization of the sampled confidence scores overlayed
            resp_sz = round(img_support_sz*currentScaleFactor*scaleFactors(scale_ind));
            xs = floor(det_sample_pos(2)) + (1:resp_sz(2)) - floor(resp_sz(2)/2);
            ys = floor(det_sample_pos(1)) + (1:resp_sz(1)) - floor(resp_sz(1)/2);           
            % To visualize the continuous scores, sample them 10 times more
            % dense than output_sz. 
            sampled_scores_display = fftshift(sample_fs(scores_fs(:,:,scale_ind), 10*output_sz));
            
            figure(fig_handle);
%                 set(fig_handle, 'Position', [100, 100, 100+size(im,2), 100+size(im,1)]);
%             imagesc(im_to_show);
            imshow(im_to_show);
            hold on;
%             resp_handle = imagesc(xs, ys, sampled_scores_display); colormap hsv;
%             alpha(resp_handle, 0.5);
            rectangle('Position',rect_position_vis, 'EdgeColor','w', 'LineWidth',1);
            text(10, 15, int2str(frame), 'Color','w', 'FontSize',15);
            hold off;
            
%                 axis off;axis image;set(gca, 'Units', 'normalized', 'Position', [0 0 1 1])
    end
end