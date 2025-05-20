for i=1:36
    if rem(i, 2) == 0
        tNewEPI(:,:,i,1,1,:) = fliplr(trajCorrectedEPI(:,:,i,1,1,:));
        kSPaceNewRARE(1,:,i,:,1,:)= fliplr(kSpaceRARE(1,:,i,:,1,:));
    else
        tNewEPI(:,:,i,1,1,:) = (trajCorrectedEPI(:,:,i,1,1,:));
        kSPaceNewRARE(1,:,i,:,1,:)= (kSpaceRARE(1,:,i,:,1,:));
    end
end