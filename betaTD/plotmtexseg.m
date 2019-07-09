
m=611;
n=1083;
map=reshape(ID10, [1083,611]);
[colorsbnd,bnds]=  bndmap(map,[],[]);
imagesc(alphacolors);hold on
imagesc(permute(colorsbnd, [2,1,3]),'AlphaData',bnds')
M=round(n/50);
N=round(m/50);

f=gca;
f.PlotBoxAspectRatio=[M,N,1];


% axis equal
print(['figures/AFbetaMTEX10'],'-dpng')
close


m=611;
n=1083;
map=reshape(ID5, [1083,611]);
[colorsbnd,bnds]=  bndmap(map,[],[]);
imagesc(alphacolors);hold on
imagesc(permute(colorsbnd, [2,1,3]),'AlphaData',bnds')
M=round(n/50);
N=round(m/50);

f=gca;
f.PlotBoxAspectRatio=[M,N,1];


% axis equal
print(['figures/AFbetaMTEX5'],'-dpng')
close


m=611;
n=1083;
map=reshape(ID15, [1083,611]);
[colorsbnd,bnds]=  bndmap(map,[],[]);
imagesc(alphacolors);hold on
imagesc(permute(colorsbnd, [2,1,3]),'AlphaData',bnds')
M=round(n/50);
N=round(m/50);

f=gca;
f.PlotBoxAspectRatio=[M,N,1];


% axis equal
print(['figures/AFbetaMTEX15'],'-dpng')
close

m=611;
n=1083;
map=reshape(ID2, [1083,611]);
[colorsbnd,bnds]=  bndmap(map,[],[]);
imagesc(alphacolors);hold on
imagesc(permute(colorsbnd, [2,1,3]),'AlphaData',bnds')
M=round(n/50);
N=round(m/50);

f=gca;
f.PlotBoxAspectRatio=[M,N,1];


% axis equal
print(['figures/AFbetaMTEX2'],'-dpng')
close