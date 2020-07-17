function [ RM2 ] = ROI_edge_cut( RM )
RM2 = RM;
for i = 2 : length(RM(:,1))-1
    for j = 2 : length(RM(1,:))-1
        
        if RM(i,j)>0 && RM(i+1,j-1)>0 && RM(i+1,j-1) ~= RM(i,j); RM2(i+1,j-1)=0; RM2(i,j)=0; end
        if RM(i,j)>0 && RM(i  ,j-1)>0 && RM(i  ,j-1) ~= RM(i,j); RM2(i  ,j-1)=0; RM2(i,j)=0; end
        if RM(i,j)>0 && RM(i-1,j-1)>0 && RM(i-1,j-1) ~= RM(i,j); RM2(i-1,j-1)=0; RM2(i,j)=0; end
        
        if RM(i,j)>0 && RM(i+1,j  )>0 && RM(i+1,j  ) ~= RM(i,j); RM2(i+1,j  )=0; RM2(i,j)=0; end
       %if RM(i,j)>0 && RM(i  ,j  )>0 && RM(i  ,j  ) ~= RM(i,j); RM2(i  ,j  )=0; RM2(i,j)=0; end
        if RM(i,j)>0 && RM(i-1,j  )>0 && RM(i-1,j  ) ~= RM(i,j); RM2(i-1,j  )=0; RM2(i,j)=0; end

        if RM(i,j)>0 && RM(i+1,j+1)>0 && RM(i+1,j+1) ~= RM(i,j); RM2(i+1,j+1)=0; RM2(i,j)=0; end
        if RM(i,j)>0 && RM(i  ,j+1)>0 && RM(i  ,j+1) ~= RM(i,j); RM2(i  ,j+1)=0; RM2(i,j)=0; end
        if RM(i,j)>0 && RM(i-1,j+1)>0 && RM(i-1,j+1) ~= RM(i,j); RM2(i-1,j+1)=0; RM2(i,j)=0; end
        
    end
end
