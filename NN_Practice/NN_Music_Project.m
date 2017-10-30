% TODO:
% Scaling of frequency bin sizes (e.g. #Hz at f = floor(log(f))+1)
% Automate data loading
% Validation
% Momentum: if fit improves, add c1 to momentum, else subtract c2; if
% momentum is greater than threshold k, keep going that direction
% 
% DONE:
% Simplify to single channel

%%
% % Get Data (run only if needed)%
% SongList = dir('C:\Users\OliWhail\Desktop\Road Trip Music\Road Trip Music (wav)\NN_Training\*.wav');
% SongData = cell(1,length(SongList));
% for index = 1:length(SongData)
%     vals = importdata(strcat('C:\Users\OliWhail\Desktop\Road Trip Music\Road Trip Music (wav)\NN_Training\',SongList(index,1).name),1);
%     % Mono compression
%     val = mean(vals.data,2);
%     SongData{index} = val;
% end

%%
% Pull in time domain song data
% load('SongData_Full.mat')
% Truncate to a manageable size
% SongData = SongData(1:5);

%%
clear all
close all
clc
rng(100*sum(clock))

% Load appropriate data subset
load('SongData_1.mat')

% Assume 44.1k sampling rate
fs = 44100;

UnsongData = cell(1,length(SongData));
for index = 1:length(SongData)
    song = SongData{index};
    unsong = abs(fft(song));
    UnsongData{index} = unsong;
    
    figure
    plot(real(unsong))
end

% Create a vector of weights to change stochastically
w = randn(1,fs);