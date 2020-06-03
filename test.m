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
V = mmread(strcat(s,'/test.mp4'));
vidframes = cat(4,V.frames.cdata);
[H,W,C,F] = size(vidframes);
vidframes_gray = zeros(H,W,1,F);
for i = 1:F
    vidframes_gray(:,:,1,i) = rgb2gray(vidframes(:,:,:,i));
end
vidframes_gray = vidframes_gray/256;
tic;
PatchFinding(vidframes_gray(:,:,:,1:2),5,8,4,'fast');
toc;

tic;
PatchFinding(vidframes_gray(:,:,:,1:2),5,8,4,'exhaustive');
toc;