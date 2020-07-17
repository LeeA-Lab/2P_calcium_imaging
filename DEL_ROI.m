function [ROIM, ROCM] = DEL_ROI(ROIM, ROCM, op)

% This function allows the user to eliminate a single ROI using one of 2
%   options.
%
% ROIM : Mask matrix (should be same dimensions as IMG size, one entry per
%        pixel) that marks each mask with the corresponding ROI numbers 
%        (cell ID), with zeros at unassigned pixels.
% ROCM : Same as ROIM except only have nonzero entries at the originally
%        clicked centers of each ROI (number is the cell ID).
% op : Option:
%      if set to 1, delete ROI that user clicks on (must click on a pixel
%        in the ROI mask), then reassign cell ID numbers to fill in any gap
%        in ID numbers left by deletion
%      if set to 0, delete ROI with the highest cell ID number (i.e. the
%        last ROI added)

if op==1
    [y,x] = ginput(1);
    xp = round(x);
    yp = round(y);
    delnum = ROIM(xp, yp);
    if delnum ~= 0
        ROIM(find(ROIM==delnum)) = 0;
        ROCM(find(ROCM==delnum)) = 0;
        if max(max(ROIM)) > delnum
            ROIM(find(ROIM > delnum)) = ROIM(find(ROIM > delnum)) - 1;
        end
        if max(max(ROCM)) > delnum
            ROCM(find(ROCM > delnum)) = ROCM(find(ROCM > delnum)) - 1;
        end
    end    
elseif op==0
    delnum = max(max(ROIM));
    ROIM(find(ROIM==delnum)) = 0;
    ROCM(find(ROCM==delnum)) = 0;
end
