names={'sim','AFone','RX','mapcenter','mapedge','hardbot','hardmid','hardtop'};
for i=1:8
    filename=names{i};
    name=names{i};
    runMStd(filename,name,12,24,50);
end


% names={'sim','AFone','RX','mapcenter','mapedge','hardbot','hardmid','hardtop'};
% for i=1:9
%     filename=names{i};
%     name=names{i};
%     runMStd(filename,name,12,24,50);
% end