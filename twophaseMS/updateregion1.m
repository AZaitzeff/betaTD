
smallest=10;
for i=1:M
    part=current(i);
    save(['results/here' num2str(part) '.mat'],'map');
    region=(map==part);
    g1=dict(part);
    [u,~,~]=oneregion(region,EBSD,CI,betas,0,0,fid,region,g1,[1,0,0,0],0);
    newregion1=region.*(u<.5);
    newregion2=region.*(u>=.5);
    diff1=sum(sum(newregion1));
    diff2=sum(sum(newregion2));
    if diff1>smallest && diff2>smallest
        newpart=part;
        CC1 = bwconncomp(newregion1>.5);
        CC2 = bwconncomp(newregion2>.5);
        PixelIdxList=[CC1.PixelIdxList CC2.PixelIdxList];
        numPixels = cellfun(@numel,PixelIdxList);
        [~,indexd] = sort(numPixels,'descend');
        [~,indexa] = sort(numPixels);
        for ind=indexa
            if numPixels(ind)<smallest
                region1=zeros(size(map));
                region1(PixelIdxList{ind})=1;
                se = strel('disk',2);
                for ind2=indexd
                    if ind2==ind
                        warning('this code should not be run')
                    end
                    region2=zeros(size(map));
                    region2(PixelIdxList{ind2})=1;
                    dreg1=imdilate(region1,se);
                    dreg2=imdilate(region2,se);
                    overlap=sum(sum(dreg1.*dreg2));
                    if overlap
                        PixelIdxList{ind2}=[PixelIdxList{ind2}; PixelIdxList{ind}];
                        break
                    end
                    
                end
                    
            else
                while sum(sum(newpart==map))>0
                    newpart=newpart+1;
                end
                map(PixelIdxList{ind})=newpart;
                new(z)=newpart;
                z=z+1;
            end
        end
            

    end
    
end