function [seq] = load_video_info_qw(video_path,startframe,sub_flag)
    
    videoData = VideoReader(video_path);%读取视频文件
    seq.len = videoData.NumberOfFrames;%结束帧CurrentTime
    seq.startframe = startframe;
    seq.videoData = videoData;
    seq.sub_flag = sub_flag;
    pos = [videoData.Height/2,videoData.Width/2];
    
    videoframe = read(videoData,startframe);
    if sub_flag==1
        im = get_subimg(videoframe,pos);%截取720*576中间图像
    else
        im = videoframe;
    end
    
    [~,seq.init_rect] = imcrop(im);%分割图像
    
end

