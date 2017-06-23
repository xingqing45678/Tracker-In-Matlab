function [seq, ground_truth] = load_video_info_qw(video_path,startframe)

if video_path(end) ~= '/', video_path(end+1) = '/'; end

ground_truth = dlmread([video_path '/groundtruth_rect.txt']);%

seq.init_rect = ground_truth(1,:);
sequence_path = [video_path,'/'];
%%%%%%%%%%%%%%%%%%%%%%%%%
    % ∂¡»°À˘”–ÕºœÒ÷°
    dir_content = dir([sequence_path 'img/']);
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
seq.startframe = startframe;
end

