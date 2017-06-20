function [handle_box] = showvideo_qw(im,frame,rect_pos,pos,handle_box)
        
        im_handle = handle_box(1);
        rectt_handle = handle_box(2);
        tex_handle = handle_box(3);
        pos_handle = handle_box(4);

        set(im_handle, 'CData', im)
        set(rectt_handle, 'Position', rect_pos)
        set(tex_handle, 'string', strcat('帧数：',num2str(frame)))
        hold on
        set(pos_handle,'XData',pos(2),'YData',pos(1),'MarkerSize',8)
%         pause(0.001);
        drawnow;
         %保存分析数据
%     save(fullfile(output_path),'target_sz','target_sz_saved');
end