
% This demo script runs the ECO tracker with deep features on the
% included "Crossing" video.
close all;clearvars;clc;
% Add paths
setup_paths();

% Load video information
startframe = 500;
video_path = 'F:\testwot\121B_623_1317\worldoftanks 2017-06-23 12-49-13-324.avi';
[seq] = load_video_info_qw(video_path,startframe);
% Run ECO
results = testing_ECO(seq);