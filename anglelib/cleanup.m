[clusterlist,~,labels] = unique(map);
N=numel(clusterlist);
smallest=16;
for i=1:N
    part1=clusterlist(i);
    region1=(map==part1);
    if sum(region1(:))<smallest && part1~=0
        biggestoverlap=0;
        part=part1;
        for j =1:N
            part2=clusterlist(j);
            region2=(map==part2);
            se = strel('disk',1);
            dreg1=imdilate(region1,se);
            dreg2=imdilate(region2,se);
            overlap=sum(sum(sum(dreg1.*dreg2)));
            if part1~=part2 && overlap>biggestoverlap
                biggestoverlap=overlap;
                part=part2;
            end
        end
        map(region1)=part;
    end
    
end