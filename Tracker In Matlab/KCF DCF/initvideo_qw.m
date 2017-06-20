function [handle_box] = initvideo_qw(startframe,im,frame,rect_pos,pos)
    if frame==startframe
        [height,width] = size(im);
        figure
        im_handle = imshow(im);
        title('DCF');
%         set(gcf,'Position',[500 300 width height]);
        axis normal;
        
        rectt_handle = rectangle('Position',rect_pos, 'EdgeColor','w');
        tex_handle = text(5, 18, strcat('Ö¡Êý:',num2str(frame)), 'Color','w', 'FontWeight','bold', 'FontSize',10);
        hold on;
        pos_handle = plot(pos(2),pos(1),'w+','MarkerSize',8);

        drawnow;
        handle_box = [im_handle,rectt_handle,tex_handle,pos_handle];
    end
end