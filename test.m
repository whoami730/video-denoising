clear; clc;
addpath(genpath('mmread'));
s = pwd;
V = mmread(strcat(s,'/bus.y4m'));
vidframes = double(cat(4,V.frames.cdata));
[H,W,C,F] = size(vidframes);
vidframes_filtered= zeros(H,W,C,F);
indices=  zeros(H,W,C,F);
doFrames=10;
searchArea=7;
Fsel=3;
patchSize=8;
refInt=4;
neighbourhood=7;
noise_key = {'gaussian','impulsive','poisson'};
noise_value = {10,0.3,0.05};
M = containers.Map(noise_key,noise_value);
vidframes_noisy = (vidframes+poissrnd(M('poisson').*vidframes) +randn(size(vidframes)).*M('gaussian'))/255;
vidframes_noisy = imnoise(vidframes_noisy,'salt & pepper',M('impulsive'));
for i = 1:doFrames
i
    [vidframes_filtered(:,:,:,i), indices(:,:,:,i)]= Med_Filter(vidframes_noisy(:,:,:,i), 5);
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
