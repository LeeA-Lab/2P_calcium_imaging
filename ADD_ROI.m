function [ROIM, ROCM] = ADD_ROI(IMG, ROIM, ROCM, Tifmov, trs)

% This function determines the pixels for an ROI mask around a center point
%   that was clicked by user.
%
% IMG : Should be averaged or standard deviation image of tiff image stack.
%       The code assumes a square image.
%       Note imagesc(IMG) displays IMG in the orientation of the matrix
%       itself, so IMG(1,1) is shown at the top left and
%       IMG(length(IMG(:,1)),1) is shown at the bottom left.
% ROIM : Mask matrix (should be same dimensions as IMG size, one entry per
%        pixel) that marks each mask with the corresponding ROI numbers 
%        (cell ID), with zeros at unassigned pixels.
% ROCM : Same as ROIM except only have nonzero entries at the originally
%        clicked centers of each ROI (number is the cell ID).
% Tifmov : To select pixels for ROI, this function needs a sample Tiffmovie
%          image.  The time series of the tiff movie is used for computing
%          the temporal correlation of pixels.
% trs(1) : Low limit of correlation.
% trs(2) : High limit of correlation (must be < 1).
% trs(3) : Option:
%          if set to 1, selects as ROI mask only the largest group of
%            pixels where each pixel is linked to all other pixels via
%            a chain of correlations between pixels (where correlated 
%            pixels means the correlation value between the low and high
%            limits); note that if there's a tie for the largest group, the
%            pixels of all of those tied groups will be included in the
%            mask
%          if set to 0 (default), selects as ROI mask all pixels that are
%            correlated with at least one other pixel in the candidate area
%          note that if the number of pairs of correlated pixels is < 5, 
%            the mask will be set to be the 20 pixels with highest 
%            intensity (in terms of the IMG, e.g. if it's a standard 
%            deviation image then it's the 20 pixels with largest standard 
%            deviation)
% trs(4) : Determines size of the ROI candidates area.  This was generally
%          set to 4 for our dorsal hippocampal CA1 data where 1 pixel was
%          ~1 um.  Adjust appropriately for the data (and if much
%          different, this may also require changing the number of
%          correlated pairs conditions described for trs(3)).
%            e.g. if trs(4) = 4 then the candidate are is 9 x 9 pixels 
%                 centered on the user chosen pixel, i.e. candidate area
%                 contains all pixels 4 units away from the center pixel in
%                 either dimension
%
% Note: This code is not vectorized since the speed didn't matter relative
%       to the manual part of the selection process; however it could be.

bs=trs(4);  % bs is same as par_R in ROI_select_Main.m
ref = 1:((2*bs+1)^2)^2;  % 2*bs+1 = length of side of ROI candidate area
refmtx = reshape(ref,(2*bs+1)^2,(2*bs+1)^2);
  % note: if np = (2*bs+1)^2 = total # of pixels in candidate area, then
  %
  %                [1   np+1  . . .    .
  %                 2   np+2
  %       refmtx =  .    .             .
  %                 .    .             .
  %                 np        . . .   np^2]
  
[y,x] = ginput(1);
  % note: flipped x & y so that IMG(row,col) represented as IMG(x,y)
  %       This means for example the Matlab data cursor will, when pointing
  %       on the imagesc-displayed IMG itself, indicate a pixel in the 8th 
  %       row from the top and 2nd row from the left in the image as [X,Y] 
  %       = [2 8], but this will be recorded as xp = 8 and yp = 2, so that
  %       IMG(row,col) can be treated as IMG(x,y) below, i.e. x indicates 
  %       the row# in IMG and increases when going down in the IMG, and y
  %       indicates the col# in IMG and increases when going to the right
  %       in the IMG.
xp = round(x);
yp = round(y);
    
	% if area based on the chosen pixel is too close to the edge of the 
    % image, shift area towards inside of image
    if xp+bs > length(IMG(1,:))
        xp = length(IMG(1,:))-bs;
    elseif xp <= bs
        xp = bs+1;
    end
    if yp <= bs
        yp = bs+1;
    elseif yp+bs > length(IMG(1,:))
        yp = length(IMG(1,:))-bs;
    end

    % mark center of ROI candidate area with next unassigned cell ID# by
    % adding to existing ROCM that stores previously added ROI centers
    rocnum = max(max(ROCM))+1;
    ROCM(xp, yp) = rocnum;
    
    % get info of pixels in candidates area   
    fdx = 1;
    for yidx = yp-bs : yp+bs
        for xidx = xp-bs : xp+bs 
            FS(fdx,1) = xidx;
            FS(fdx,2) = yidx;
            FS(fdx,3) = IMG(xidx, yidx);
            fdx = fdx + 1;
        end
    end
      % note: for instance if the center is (5,5) and bs=4, then have 81
      %       pixels in a 9 x 9 array, and they and their values in IMG
      %       will be represented in FS as:
      %
      %             x   y  val
      %            [1   1   i1    % for pixel 1
      %             2   1   i2    % for pixel 2
      %             3   1   i3    % for pixel 3
      %             .   .   .
      %       FS =  1   2   i10   % for pixel 10
      %             2   2   i11   % for pixel 11
      %             3   2   i12   % for pixel 12
      %             .   .   .
      %             9   9   i81]  % for pixel 81
      %
      %       and in terms of the IMG displayed by imagesc, the pixel
      %       numbers j (x,y) are:
      %
      %             1 (1,1)   10 (1,2)    .      .
      %             2 (2,1)   11 (2,2)    .      .
      %             .
      %             9 (9,1)      .        .   81 (9,9)
    
    % get fluorescence traces from Tifmov to calculate temporal correlation	
    FStraces = zeros(length(Tifmov(1,1,:)),length(FS(:,1)));
      % e.g. for bs=4 (9 x 9 pixel ROI candidate area)
      %
      %      FStraces = [pixel 1   pixel 2  ...  pixel 81
      %
      %                     |         |             |
      %                     v         v             v
      %                                                  ]
      %
      %      where the time series trace for each pixel is given in each
      %      col
      %
    for idx = 1 : length(FS(:,1))
        FStraces(:,idx) = Tifmov(FS(idx,1),FS(idx,2),:); 
    end
    FScorr = corr(FStraces); % get correlation values 
      % e.g. for bs=4 (9 x 9 pixel ROI candidate area)
      %
      %                            pixel 1   pixel 2  ...  pixel 81
      %                pixel 1   [
      %                pixel 2
      %      FScorr =     .          pairwise correlation matrix
      %                   .
      %                pixel 81                                    ]
      %
    FSidx = find(FScorr > trs(1) & FScorr < trs(2)); % filtering
      % note: FSidx gives a single col of indices that correspond to the
      %       correlation "hits" in FScorr, columnwise, e.g. for the
      %       example FScorr above, if there is a hit in 
      %       FScorr(pixel 1, pixel 2), this will be given as a 2 in FSidx,
      %       and for FScorr(pixel 2, pixel 1) which would have the same
      %       value, this will be given as an 82 in FSidx.
      %       Importantly, these indices corresponds to refmtx above.

% beginning of options %

% find qualifying pixel positions

if trs(3)==0
    fdx = 1; FSC=[];
    for idx = 1 : length(FSidx)
        [a, b] = find(refmtx == FSidx(idx));
          % because FSidx and refmtx are "in register", this line gives
          %   values a and b such that a corresponds to the pixel # of one
          %   of the pixels in the pair of correlated pixels represented by
          %   FSidx(idx), and b corresponds to the pixel # of the other
          %   pixel in the pair
        if idx == 1 
            FSC(fdx) = a; fdx = fdx + 1;
            FSC(fdx) = b; fdx = fdx + 1;
        elseif idx > 1 && isempty(find(FSC == a)) == 1
            FSC(fdx) = a; fdx = fdx + 1;
        elseif idx > 1 && isempty(find(FSC == b)) == 1
            FSC(fdx) = b; fdx = fdx + 1;
        else
        end
    end

elseif trs(3)==1 && isempty(FSidx)==0
    
    FSCC = zeros((2*bs+1)^2,(2*bs+1)^2);
    for i = 1 : length(FScorr(1,:))
        tmp=find(FScorr(:,i)> trs(1) & FScorr(:,i) < trs(2));
        FSCC(1:length(tmp),i) = tmp; clear tmp;
    end
      % note: FSCC has size #pixels x #pixels, e.g.
      %
      %              [3  0  1  0  0  1  .
      %               6
      %               0
      %       FSCC =  .
      %               .
      %               .
      %               0                 0]
      %
      %       for the case where pixel 1 is correlated with pixels 3 & 6

    FSgr=[]; gn=1; nmg=0;  % nmg = no more groups
    while nmg==0
        % Find first pixel in new group
        cg=[]; i=1;
        while isempty(cg) & i<=length(FSCC(1,:))
            if ~isempty(find(FSCC(:,i)>0)) & isempty(intersect(i,FSgr)) 
                cg = i;  % seed the new group
            else
                i=i+1;
            end
        end
        if isempty(cg)
            nmg=1;  % if no new groups
        else
            % Complete new group
            gnd=0;  % gnd = new group not done
            while gnd==0
                % Check if all pixels correlated to pixels currently in
                %   group are in the group, in which case the group is
                %   done, otherwise update group
                tmp = unique(FSCC(:,cg));  % pixels correlated to pixels in group
                if size(tmp,2)>1
                    tmp = tmp';  % make into column vector if not already
                end
                tmp = [cg ; tmp];
                tmp = unique(tmp(tmp>0));
                if ~isequal(tmp,cg)  % if group not "closed"
                    cg = tmp;  % update group
                else  % else group is complete
                    gnd=1;
                    FSgr(1:length(cg),gn)=cg; gn=gn+1;  % update group list
                end
            end
        end
    end

    for i = 1: length(FSgr(1,:)); FSgrn(i)=length(find(FSgr(:,i)>0)); end
    m=find(FSgrn==max(FSgrn)); FSC=FSgr(:,m);

end

% end of options %
        
	% if number of pairs of pixels (that have qualifying high temporal 
    % correlation) is less than 5 (e.g. a cell that has very low calcium
    % activity), then select the 20 pixels with highest intensity (as in 
    % highest values in IMG, which would mean the highest standard 
    % deviation values if that's what's input as IMG)
    if  length(FSidx) >= 10
        % since every correlated pair of pixels is counted twice in FSidx, 
        % 10 means 5 pairs
        for index = 1 : length(FSC)
            FS2(index,:) = FS(FSC(index),:);
        end
    elseif isempty(FSidx) == 1 || length(FSidx) < 10 
        FSsort = sortrows(FS,-3);
        FS2 = FSsort(1:20,:);
    end
    
    roinum = max(max(ROIM))+1;
    for index = 1 : length(FS2(:,1))
        ROIM(FS2(index,1),FS2(index,2)) = roinum;
    end
    
% update figure
c=hsv(25);
[y,x]=find(ROIM==roinum); m=rem(roinum,25)+1; 
    for k=1:length(y)
        hold on; plot(x(k),y(k),'MarkerSize',8,'Marker','.','LineStyle','none','Color',c(m,:));
    end; clear x; clear y;
[y,x]=find(ROCM==roinum); 
plot(x,y,'MarkerSize',8,'Marker','s','LineStyle','none','color','b');  
