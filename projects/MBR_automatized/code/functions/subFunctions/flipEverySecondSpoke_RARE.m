for i=1:36
    if rem(i, 2) == 0
        anglesRAREnew(i,:) = anglesRARE(i,:) + pi;
        kSpaceRAREnew(1,:,i,:,1,:)= fliplr(kSpaceRARE(1,:,i,:,1,:));
    else
        anglesRAREnew(i,:) = anglesRARE(i,:);
        kSpaceRAREnew(1,:,i,:,1,:)= (kSpaceRARE(1,:,i,:,1,:));
    end
end