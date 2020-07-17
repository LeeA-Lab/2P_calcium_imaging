%%
% ROI_select_Main.m
% Jae Sung Lee
%
% This code allows manual selection of ROIs.  The user clicks near the
% center of a cell and an ROI is created from the surrounding pixels that
% consists of a subset of pixels.  See example image ROI_example.jpg.  Note
% that masks do not try to contain all pixels of a cell.  The user needs to
% set the initial candidate mask size (with par_R) so that it doesn't 
% invade the neighboring cells for most cells, then can choose to shrink
% the candidate mask in individual cases if the nearby cells are too close.
%
% This code was developed and used for 2-photon in vivo calcium imaging 
% data (GCaMP6f) from mouse dorsal hippocampal area CA1 acquired at 30Hz
% using Scanimage 2015-2018.  The raw image resolution was 512 x 512 pixels
% covering 500 x 500 um.  For this data, the parameter par_R was generally
% set to 4 in trs (i.e. 9 x 9 pixel ROI candidate area per cell) and to 3
% (7 x 7 pixels) in trs_small. 
%
% Before running this ROI selection process, motion correction was done
% using the MOCO (Yuste lab) ImageJ plugin.  We then cut out 6 pixels 
% around the 4 edges after motion correction, giving a 500 x 500 pixel
% image.  We then made a standard deviation image (std_img.tif) using the
% ImageJ Z projection of a time series of the imaging data.  This standard
% deviation image was created using at least ~3500 frames (at 30Hz).  For
% computing correlations between pixels, which was one feature used to 
% determine the ROI mask, we again used a time series of at least ~3500 
% frames (at 30Hz): mov_tif.tif.  Using fewer frames made the correlation
% values noisy so that's not advised.  In the case where a given ROI 
% candidate area did not have a minimum number of pairs of correlated 
% pixels, the pixels with highest standard deviation were taken as the ROI
% mask.

%% Load STD image of Tiff stack (STD z-projection of time series of imaging) 
std_img=imread('std_img.tif'); % load standard deviation image
  % image is assumed to be square

%% Load Tiff mov & Image Normalization
% load time series of calcium imaging to get temporal correlation of pixels
% for ROI selection
numfr = 3500; 
% numfr = number of frames of time series of calcium imaging
% We used a minimum of 3500 frames for this (more frame numbers are better 
% for computing temporal correlations)
for idx = 1 : numfr; Tifmov(:,:,idx) = double(imread('mov_tif.tif',idx)); end; clear idx; clear numfr;
  % each Tifmov frame should be same dimensions as std_img
Tifmov = Tifmov-min(min(min(Tifmov))); 
  % note imagesc(std(Tifmov,ones(size(Tifmov,3),1),3)); colormap(gray); axis equal
  %   should look similar to (but not necessarily exactly like) std_img 
  %   when plotted below by Show_IMG.m (e.g. orientations should match or
  %   else need to change orienatations so that they do)

%% Making empty ROI mask
pxn=200; % for image data of pxn x pxn pixels
Mask_ROI = zeros(pxn,pxn); % make empty Mask_ROI 
Mask_Cent = zeros(pxn,pxn); % make empty Mask_Cent

%%
% ROI selection by temporal correlation
% low_threshold = 0.45; % temporal correlation lower limit
% high_threshold = 0.85; % temporal correlation upper limit (must be < 1)
% group_option : include only the largest group of correlated pixels in ROI
%   (1: use option / 0: don't use option)
%   (but also see the other condition in ADD_ROI.m regarding the minimum
%   number of correlated pairs)
% par_R = 3; % 3 -> default size of ROI candidate field is 7 x 7 pixels
% par_R = 4; % 4 -> default size of ROI candidate field is 9 x 9 pixels
% par_R = 5; % 5 -> default size of ROI candidate field is 11 x 11 pixels

% trs= [low_threshold, high_threshold, group_option, par_R]
trs = [0.45, 0.85, 0, 4]; 
  % note: 2nd value must be < 1
% same as trs, except smaller size of ROI candidate field
trs_small = [0.45, 0.85, 0, 3]; 
  % note: 2nd value must be < 1
  % note: for consistency, the first 3 values of trs and trs_small should
  %       match
  
figure
LBWH = get(gcf,'Position'); close
figure('Position',[LBWH(1)+LBWH(3)/2-LBWH(4) LBWH(2) LBWH(4)*2 LBWH(4)])
  % make figure initially 2x as wide as high, with same height, same center
ax = subplot('Position',[0.025 0.05 0.45 0.9]);
ax2 = subplot('Position',[0.525 0.05 0.45 0.9]);
axes_lim = [];
cont_f = 2; % increase contrast of image (x2), or 1 for original std image
Show_IMG(std_img, pxn, cont_f, Mask_ROI, Mask_Cent, ax, ax2, axes_lim); 

% ROI color: HSV colormap was used (25 steps of colormap, ROI#1-25 have  
%   different colors, ROI#26 will have same color as ROI#1, etc)
% User can zoom in/out using figure control tools during selection process
%   and add (ADD_ROI) or delete (DEL_ROI) ROI in zoomed figure
% User can add ROI in either panel but recommend using right panel, and use
%   left panel to assess ROI relative to original underlying image
% All ROIs will be shown in the right panel if user re-runs this section 
% Solid colors : ROI mask pixels
% Open Squares : location that was clicked for each ROI
% When running this section:
%   If user wants to add ROI, press 1 (then enter) and locate cursor near
%     the center (or the region you want) of the cell that you want to 
%     select then click
%       ROI of the cell will then show up in the panel
%   If user wants delete the last ROI, press 2 (then enter)
%   If user wants a smaller ROI candidate field, press -1 (then enter) and
%     select the cell you want (this is for a small cell or region with a
%     higher cell density)
%   If user wants to delete a specific ROI, press 3 (then enter) and select
%     the target ROI by clicking on one of the pixels in the mask (that are
%     marked by color); if any ROI other thatn the last is deleted, the
%     cell IDs numbers of each ROI will be updated to fill in the gap left
%     by the deleted ROI, and the ROI colors displayed will be changed
%     correspondingly
%   If user wants to end this process, press any other number (then enter)

inp = input('1:add ROI / 2:delete last ROI / 3:del ROI(select) / -1:add ROI(smaller par_R) / else:out_');
while inp == -1 || inp == 1 || inp == 2 || inp == 3
    if     inp ==  1 ;[Mask_ROI, Mask_Cent] = ADD_ROI(std_img, Mask_ROI, Mask_Cent, Tifmov, trs);
    elseif inp == -1 ;[Mask_ROI, Mask_Cent] = ADD_ROI(std_img, Mask_ROI, Mask_Cent, Tifmov, trs_small);
    elseif inp ==  2 ;[Mask_ROI, Mask_Cent] = DEL_ROI(Mask_ROI, Mask_Cent,0); cla(ax); cla(ax2); clear axes_lim; axes_lim([1 2]) = get(ax2,'XLim'); axes_lim([3 4]) = get(ax2,'YLim'); Show_IMG(std_img, pxn, cont_f, Mask_ROI, Mask_Cent, ax, ax2, axes_lim); 
    elseif inp ==  3; [Mask_ROI, Mask_Cent] = DEL_ROI(Mask_ROI, Mask_Cent,1); cla(ax); cla(ax2); clear axes_lim; axes_lim([1 2]) = get(ax2,'XLim'); axes_lim([3 4]) = get(ax2,'YLim'); Show_IMG(std_img, pxn, cont_f, Mask_ROI, Mask_Cent, ax, ax2, axes_lim); 
    end;   inp = input('1:add ROI / 2:delete last ROI / 3:del ROI(select) / -1:add ROI(smaller par_R) / else:out_');
end; clear inp; 

%% Just Show image and ROI mask 
figure; LBWH = get(gcf,'Position'); close
figure('Position',[LBWH(1)+LBWH(3)/2-LBWH(4) LBWH(2) LBWH(4)*2 LBWH(4)])
ax = subplot('Position',[0.025 0.05 0.45 0.9]);
ax2 = subplot('Position',[0.525 0.05 0.45 0.9]); axes_lim = [];
Show_IMG(std_img, pxn, cont_f, Mask_ROI, Mask_Cent, ax, ax2, axes_lim); 

%% Option: Eliminate ROI pixels that are close to other ROIs then show image and modified ROI mask
% Running this results in at least a 1 pixel gap with the neighboring ROIs
% It does nothing if ROIs are already separated from each other
% Do not need to run if user think ROIs are not overlapping
% Note that this process reduces the number of pixels in some ROIs
[ Mask_ROI2 ] = ROI_edge_cut( Mask_ROI );
figure; LBWH = get(gcf,'Position'); close
figure('Position',[LBWH(1)+LBWH(3)/2-LBWH(4) LBWH(2) LBWH(4)*2 LBWH(4)])
ax = subplot('Position',[0.025 0.05 0.45 0.9]);
ax2 = subplot('Position',[0.525 0.05 0.45 0.9]); axes_lim = [];
Show_IMG(std_img, pxn, cont_f, Mask_ROI2, Mask_Cent, ax, ax2, axes_lim); 
