function Show_IMG(std_img, pxn, cont_f, Mask_ROI, Mask_Cent, ax, ax2, axes_lim)

pf=500/pxn; npx=round(pf*pxn);

%% Main Window, showing image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Note imagesc(std_img) displays std_img in the orientation of the matrix
%       itself, so std_img(1,1) is shown at the top left and
%       std_img(length(std_img(:,1)),1) is shown at the bottom left.

subplot(ax)
clims = [min(min(std_img))/cont_f ,max(max(std_img))/cont_f]; 
imagesc(std_img, clims); colormap(gray); axis equal; hold all;
set(ax,'xtick',[],'ytick',[]);  % Remove ticks.
if ~isempty(axes_lim); axis(axes_lim); end

subplot(ax2)
imagesc(std_img, clims); colormap(gray); axis equal; hold all;
set(ax2,'xtick',[],'ytick',[]);  % Remove ticks.
if ~isempty(axes_lim); axis(axes_lim); end

linkaxes([ax ax2],'xy');
c=hsv(25);
if max(max(Mask_Cent))>0
    for i=1:max(max(Mask_Cent))
        [y,x]=find(Mask_ROI==i); m=rem(i,25)+1; 
        for k=1:length(y)
            plot(x(k),y(k),'MarkerSize',8,'Marker','.','LineStyle','none','Color',c(m,:));
        end; clear x; clear y;
    end
    for i=1:max(max(Mask_Cent))
        [y,x]=find(Mask_Cent==i); 
        plot(x,y,'MarkerSize',8,'Marker','s','LineStyle','none','color','b');
    end
end
