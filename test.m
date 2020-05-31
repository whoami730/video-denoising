clear; clc;
% I = imread('lena512.bmp');
% 
% I_noisy = imnoise(I,'salt & pepper',0.7);
% 
% filtI = Med_Filter(I_noisy,11);
% imshow(filtI);
% figure;
% imshow(I_noisy);
% 
% indices = (filtI == I_noisy)
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
% PatchFinding(vidframes_gray);
