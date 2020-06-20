    clear; clc;
addpath(genpath('mmread'));
s = pwd;
V = mmread(strcat(s,'/trevor_qcif.y4m'));
all_frames = V.frames;
vidframes = double(cat(4,all_frames.cdata));
[H,W,C,F] = size(vidframes);
vidframes_filtered= zeros(H,W,C,F);
indices=  zeros(H,W,C,F);
doFrames=15;
searchArea=7;
Fsel=5;
%Make sure the frame dimensions are multiples of patchSize and refInt 
patchSize=8;
refInt=4;
neighbourhood=5;
noise_key = {'gaussian','impulsive','poisson'};
noise_value = {30,0.1,0.1};
M = containers.Map(noise_key,noise_value);
vidframes_noisy = ((1-M('poisson'))*vidframes+poissrnd(M('poisson').*vidframes) +randn(size(vidframes)).*M('gaussian'))/255;
vidframes_noisy = imnoise(vidframes_noisy,'salt & pepper',M('impulsive'));
for i = 1:doFrames
    i
    [vidframes_filtered(:,:,:,i), indices(:,:,:,i)]= Med_Filter(vidframes_noisy(:,:,:,i),5);
end
vidframes_o= vidframes(:,:,:, 1:doFrames)/255;
vidframes_n= vidframes_noisy(:,:,:, 1:doFrames);
vidframes_f= vidframes_filtered(:,:,:, 1: doFrames);
tic;
vidframes_a = PatchFinding(vidframes_f, indices, Fsel, patchSize, refInt, searchArea, neighbourhood, 'exhaustive');
toc;

MSE_n= sum((vidframes_o-vidframes_n).*(vidframes_o-vidframes_n), 'all');
MSE_f= sum((vidframes_o-vidframes_f).*(vidframes_o-vidframes_f), 'all');
MSE_a= sum((vidframes_o-vidframes_a).*(vidframes_o-vidframes_a), 'all');
