function [seq, ground_truth] = load_video_info_qw(video_path,startframe)
    
    img_lib = 1; %0 = 'OTB' , 1 ='VOT';
    if img_lib==0
        if video_path(end) ~= '/', video_path(end+1) = '/'; end

        ground_truth = dlmread([video_path '/groundtruth_rect.txt']);%
%         ground_truth = dlmread([video_path '/groundtruth.txt']);%

        seq.init_rect = ground_truth(startframe,:);
        sequence_path = [video_path,'/'];
        %%%%%%%%%%%%%%%%%%%%%%%%%
            % 读取所有图像帧
            dir_content = dir([sequence_path 'img/']);
%             dir_content = dir(sequence_path);%
            % skip '.' and '..' from the count
            n_imgs = length(dir_content) - 2;
            img_files = cell(n_imgs, 1);
            for ii = 1:n_imgs
                img_files{ii} = dir_content(ii+2).name;
            end
        %%%%%%%%%%%%%%%%%%%%%%%%%
        seq.len = numel(img_files);
        seq.s_frames = cellstr(img_files);
        seq.path = [video_path,'/img/'];
%          seq.path = video_path ;
        seq.startframe = startframe;
    else
        if video_path(end) ~= '/', video_path(end+1) = '/'; end

        ground_truth = dlmread([video_path '/groundtruth.txt']);%

        rect = ground_truth(startframe,:);
        [cx, cy, w, h] = getAxisAlignedBB(rect);
        init_rect(1) = cx-h/2;init_rect(2) = cy-w/2;init_rect(3) = w;init_rect(4) = h;
        seq.init_rect = init_rect;
        
        sequence_path = [video_path,'/'];
        %%%%%%%%%%%%%%%%%%%%%%%%%
            % 读取所有图像帧
        %     dir_content = dir([sequence_path 'img/']);
            dir_content = dir(sequence_path);%
            % skip '.' and '..' from the count
            n_imgs = length(dir_content) - 2;
            img_files = cell(n_imgs, 1);
            for ii = 1:n_imgs
                img_files{ii} = dir_content(ii+2).name;
            end
        %%%%%%%%%%%%%%%%%%%%%%%%%
        seq.len = numel(img_files);
        seq.s_frames = cellstr(img_files);
        %seq.path = [video_path,'/img/'];
        seq.path = video_path ;
        seq.startframe = startframe;
    end
end

