
% This demo script runs the ECO tracker with hand-crafted features on the
% included "Crossing" video.
close all;clearvars;clc;
% Add paths
setup_paths();

% Load video information
startframe = 1;
video_path = 'F:\testwot\worldoftanks 2017-06-23 12-45-32-856.avi';%121B_623_1317\worldoftanks 2017-06-23 12-49-13-324.avi';
[seq] = load_video_info_qw(video_path,startframe);

% Run ECO
results = testing_ECO_HC(seq);