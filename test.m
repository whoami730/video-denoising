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
vidframes_impnoisy= zeros(H,W,C,F);
vidframes_noisy=zeros(H,W,C,F);
vidframes_filtered= zeros(H,W,C,F);
doFrames=4;
searchArea=10;
Fsel=5;
patchSize=8;
refInt=4;
for i = 1:doFrames
    i
    for j= 1:C
        frame= vidframes(:,:,j,i);
        
        gnoise= randn(H,W)*15;
        vidframes_noisy(:,:,j,i)= double(frame)+gnoise;
        vidframes_noisy(:,:,j,i)=vidframes_noisy(:,:,j,i)/255;
        vidframes_impnoisy(:,:,j,i)= imnoise(vidframes_noisy(:,:,j,i), 'salt & pepper', 0.2);
        %vidframes_noisy(:,:,j,i)= imnoise(vidframes_impnoisy(:,:,j,i),'gaussian', 0);
        vidframes_filtered(:,:,j,i)= Med_Filter(vidframes_impnoisy(:,:,j,i), 11);
    end
end

%tic;
%Ansf = PatchFinding(vidframes_filtered(:,:,:,1:2),5,8,4,'fast');
%toc;
vidframes_o= double(vidframes(:,:,:, 1:doFrames));
vidframes_o= vidframes_o/255;
vidframes_n= vidframes_impnoisy(:,:,:, 1:doFrames);
vidframes_f= vidframes_filtered(:,:,:, 1: doFrames);
tic;
vidframes_a = PatchFinding(vidframes_f, Fsel, patchSize, refInt, searchArea, 'exhaustive');
toc;

MSE_n= sum((vidframes_o-vidframes_n).*(vidframes_o-vidframes_n), 'all');
MSE_f= sum((vidframes_o-vidframes_f).*(vidframes_o-vidframes_f), 'all');
MSE_a= sum((vidframes_o-vidframes_a).*(vidframes_o-vidframes_a), 'all');
