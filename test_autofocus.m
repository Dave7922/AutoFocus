%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference:
% Sun, Y., Duthaler, S., & Nelson, B. J. (2004). Autofocusing in 
% computer microscopy: Selecting the optimal focus algorithm. 
% Microsc Res Tech, 65(3), 139-149. doi:10.1002/jemt.20118
%
% All Rights Reserved.
% Created by
% Mingzhai Sun
% Lewis-Sigler Institute for Integrative Genomics
% Princeton University
% 5-16-2011


close all; 
clear all;
clc;

dirName = uigetdir; 
fileNames = dir(fullfile(dirName, '*.tif'));
frameNum = size(fileNames, 1); 

%select ROI if wanted
% figure,imagesc(imread(fullfile(dirName, fileNames(32).name)));
% rect = getrect;

%gaussian kernel with sigma as the std.
[X, Y] = meshgrid(-2:1:2);
sigma = 3.;
Pi = 3.1415926;
Z = 1/(2*Pi*sigma^2)*exp(-(X.^2 + Y.^2)/(2*sigma^2));

%Sobel filter
hy = fspecial('sobel');
hx = hy';

%Laplacian
L = [-1, -4, -1; -4, 20, -4; -1, -4, -1]; %equavelent with fspecial('laplacian')*6

%KH filter size
filter_size = 3;

normVar = [];%normalized variance
non_normVar = [];%non-normalized variance.
Fscore = [];%Fscore using the sqrt(fx^2+fy^2)
Fscore_sobel = [];
Fscore_energy = [];
Fscore_lap = [];
Fscore_KH = [];
Fscore_cor = [];

step = 1;
for i=1:step:frameNum
    im = double(imread(fullfile(dirName, fileNames(i).name)));
    im_sm = imfilter(im, Z, 'replicate', 'conv');
    
    if(exist('rect')~=0)
        im_sm = imcrop(im_sm, rect);
    end
    im = im_sm;
    
    %correlation method F-12 in the paper.
    temp1 = im_sm.*circshift(im_sm, 1) - im_sm.*circshift(im_sm,2);
    temp2 = im_sm.*circshift(im_sm, [0,1]) - im_sm.*circshift(im_sm, [0, 2]);
    Fscore_cor_temp = sum(sum(temp1+temp2));
    Fscore_cor = cat(1, Fscore_cor, Fscore_cor_temp);
    
    %Fscore using the sqrt(fx^2+fy^2)
    [f_x, f_y] = gradient(im);
    temp = f_x.^2 + f_y.^2;
    Fscore = cat(1, Fscore, sum(temp(:)));
    
    %normalized and non-normalized standard deviation
    im_var = var(im(:));
    im_norm_var = im_var/(mean(im(:)));
    normVar = cat(1, normVar, im_norm_var);
    non_normVar = cat(1, non_normVar, im_var);
    
    %sobel operator
    f_x_sobel = imfilter(im_sm, hy, 'replicate');
    f_y_sobel = imfilter(im_sm, hx, 'replicate');
    temp_sobel = sqrt(f_x_sobel.^2+f_y_sobel.^2);
    Fscore_sobel = cat(1, Fscore_sobel, sum(temp_sobel(:)));
    
    %Laplacian
    F_energy = imfilter(im_sm, L, 'replicate');
    F_energy_temp = sum(sum(F_energy.*F_energy));
    Fscore_energy = cat(1, Fscore_energy, F_energy_temp);
    
    %ridge measure
    [Laplacian,H,KH,KL,VH,VL]=RidgeMeasures(im_sm,filter_size); 
    Fscore_KH_temp = sum(sum(KH.*KH));
    Fscore_KH = cat(1, Fscore_KH, Fscore_KH_temp);
end

Fscore = Fscore/max(Fscore(:));
Fscore_sobel = Fscore_sobel/max(Fscore_sobel(:));
normVar = normVar/max(normVar(:));
non_normVar = non_normVar/max(non_normVar(:));
Fscore_energy = Fscore_energy/max(Fscore_energy(:));
Fscore_KH = Fscore_KH/max(Fscore_KH(:));
Fscore_cor = Fscore_cor/max((Fscore_cor(:)));

colors = distinguishable_colors(7);
figure, plot(normVar, 'o--', 'Color', colors(1, :)), 
hold on
plot(Fscore, '*--', 'Color', colors(2,:))
plot(Fscore_sobel, 'o--', 'Color', colors(3,:))
plot(non_normVar, '*--', 'Color', colors(4,:))
plot(Fscore_energy, 'd-','Color', colors(5, :));
plot(Fscore_KH, '*--', 'Color', colors(6,:));
plot(Fscore_cor, '+--', 'Color', colors(7,:));
legend('normalized variance', 'sqrt(fx^2+fy^2)', ...
    'Sobel operator', 'non-normalized variance', 'Laplacian', 'KH', 'Cor');

[~, ind] = min(Fscore_energy);
if(ind==60)
    ind = 18;
end
Fscore_energy_sub = Fscore_energy(ind:end);
Fscore_energy_z = ind:(numel(Fscore_energy_sub)+ind-1);

cfun = fit(Fscore_energy_z', Fscore_energy_sub, 'poly9');
xx = min(Fscore_energy_z):0.1:max(Fscore_energy_z);
yy = feval(cfun, xx);
[~, yyind] = max(yy);

figure, plot(Fscore_energy_z, Fscore_energy_sub, '*--');
hold on
plot(cfun, 'r-');
plot(Fscore_energy, 'cd--');
vline(xx(yyind))
hold off
fprintf('%s: %0.2f\n', 'maximum is at: ', xx(yyind))
