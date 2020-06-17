clear; clc;
addpath(genpath('mmread'));
% I = imread('lena512.bmp');
% 
% I_noisy = imnoise(I,'salt & pepper',0.7);
% 
% [filtI,indices] = Med_Filter(I_noisy,11);
% imshow(filtI);
% figure;
% imshow(I_noisy);
% 
% psnr(I_noisy,I)
% psnr(filtI,I)
s = pwd;
V = mmread(strcat(s,'/bus.y4m'));
vidframes = cat(4,V.frames.cdata);
[H,W,C,F] = size(vidframes);
vidframes_filtered= zeros(H,W,C,F);
indices=  zeros(H,W,C,F);
doFrames=3;
searchArea=7;
Fsel=3;
patchSize=32;
refInt=16;
neighbourhood=7;
gnoise = randn(H,W,C,F)*20;
vidframes_noisy = (double(vidframes)+gnoise)/255;
vidframes_impnoisy = imnoise(vidframes_noisy,'salt & pepper',0.4);
for i = 1:doFrames
i
    [vidframes_filtered(:,:,:,i), indices(:,:,:,i)]= Med_Filter(vidframes_impnoisy(:,:,:,i), 11);
end
vidframes_o= double(vidframes(:,:,:, 1:doFrames));
vidframes_o= vidframes_o/255;
vidframes_n= vidframes_impnoisy(:,:,:, 1:doFrames);
vidframes_f= vidframes_filtered(:,:,:, 1: doFrames);
tic;
vidframes_a = PatchFinding(vidframes_f, indices, Fsel, patchSize, refInt, searchArea, neighbourhood, 'exhaustive');
toc;

MSE_n= sum((vidframes_o-vidframes_n).*(vidframes_o-vidframes_n), 'all');
MSE_f= sum((vidframes_o-vidframes_f).*(vidframes_o-vidframes_f), 'all');
MSE_a= sum((vidframes_o-vidframes_a).*(vidframes_o-vidframes_a), 'all');
