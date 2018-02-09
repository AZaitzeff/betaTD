
for i=1:M
    for j=1:N
        part1=current(i);
        part2=clusterlist(j);
        region1=(map==part1);
        region2=(map==part2);
        se = strel('disk',1);
        dreg1=imdilate(region1,se);
        dreg2=imdilate(region2,se);
        overlap=sum(sum(dreg1.*dreg2));

        %if ~sum((find(part2==current)<=i)) && overlap>0
        if part1~=part2 && overlap>0

            region=region1+region2;
            [u]=oneregion(region,fEBSDq,fBSE,0,fidEBSD,fidBSE,region1,0);
            newregion1=region.*(u>.5);
            newregion2=region.*(u<=.5);

            CC1 = bwconncomp(newregion1>.5);
            CC2 = bwconncomp(newregion2>.5);
            PixelIdxList=[CC1.PixelIdxList CC2.PixelIdxList];
            numPixels = cellfun(@numel,PixelIdxList);
            if numel(PixelIdxList)==1
                map(PixelIdxList{1})=part1;
            else
            
            
                [~,indexd] = sort(numPixels,'descend');

                for ind=indexd(3:end)
                    newregion1=zeros(size(map));
                    newregion1(PixelIdxList{ind})=1;
                    se = strel('disk',2);
                    for ind2=indexd(1:2)
                        newregion2=zeros(size(map));
                        newregion2(PixelIdxList{ind2})=1;
                        dreg1=imdilate(newregion1,se);
                        dreg2=imdilate(newregion2,se);
                        overlap=sum(sum(dreg1.*dreg2));
                        if overlap
                            PixelIdxList{ind2}=[PixelIdxList{ind2}; PixelIdxList{ind}];
                            break
                        end

                    end
                end
                region=zeros(size(map));
                region(PixelIdxList{1})=1;
                origover1=sum(sum(abs(region1-region)))/sum(sum(region1));
                origover2=sum(sum(abs(region2-region)))/sum(sum(region2));
                if origover1>0 && origover2>0
                    new(z)=part1;
                    z=z+1;
                    new(z)=part2;
                    z=z+1;
                    if origover1<origover2
                        map(PixelIdxList{1})=part1;
                        map(PixelIdxList{2})=part2;
                    else
                        map(PixelIdxList{1})=part2;
                        map(PixelIdxList{2})=part1;
                     end
                end
            end
        end
    end
end