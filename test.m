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
for i = 1:10
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

tic;
Ansf = PatchFinding(vidframes_filtered(:,:,:,1:2),5,8,4,'fast');
toc;

tic;
AnsE = PatchFinding(vidframes_filtered(:,:,:,1:2),5,8,4,'exhaustive');
toc;
