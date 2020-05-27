clear; clc;
I = imread('lena512.bmp');

I_noisy = imnoise(I,'salt & pepper',0.7);

filtI = Med_Filter(I_noisy,11);
imshow(filtI);
figure;
imshow(I_noisy);

indices = (filtI == I_noisy)

psnr(I_noisy,I)
psnr(filtI,I)