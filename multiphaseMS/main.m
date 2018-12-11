names={'sim','AFone','RX','mapcenter','mapedge','hardbot','hardmid','hardtop'};
for i=1:8
    filename=names{i};
    name=names{i};
    runMStd(filename,name,12,24,50);
end


 names={'rand','sim2f','sim3f','sim4f','AFbig'};
 for i=1:5
     filename=names{i};
     name=names{i};
     runMStd(filename,name,12,24,50);
 end