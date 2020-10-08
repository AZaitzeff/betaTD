function [EBSD,CI,IQ,betas,scale]=ang2mat(ebsd,ind)
addpath('..')
%pname = '/Users/azaitzeff/Documents/Research/Ti64/';
% which files to be imported
%fname = [pname fileload];
if numel(ebsd.indexedPhasesId)==2
    alpha=ebsd('Titanium (Alpha)');
    if ind

        beta=ebsd('Titanium - Beta');
    else
        beta=ebsd('Titanium (Beta)');
    end

    LOGbeta=isempty(beta);
    LOGalpha=isempty(alpha);
else
    LOGalpha=0;
    LOGbeta=1;
end

if ~LOGbeta && ~LOGalpha


if size(ebsd.unitCell,1)~= 6
    ext=zeros(1,4);
    bext=beta.extend;
    aext=alpha.extend;
    for i=1:2
        ext(2*i-1)=min(aext(2*i-1),bext(2*i-1));
        ext(2*i)=max(aext(2*i),bext(2*i));
    end
    [newalpha,~]=gridifyz(alpha,ext);
    [newbeta,~]=gridifyz(beta,ext);
    scale=1;
else
    ext=zeros(1,4);
    bext=beta.extend;
    aext=alpha.extend;
    for i=1:2
        ext(2*i-1)=min(aext(2*i-1),bext(2*i-1));
        ext(2*i)=max(aext(2*i),bext(2*i));
    end
    [newalpha,~,~,~,~]=hexify(alpha,ext);
    [newbeta,~,~,axis,~]=hexify(beta,ext);
    scale=axis*sqrt(3)/2+(~axis)*2/sqrt(3);
end


CIa=newalpha.prop.ci;
CIa(isnan(CIa))=0;
IQa=newalpha.prop.iq;
IQa(isnan(IQa))=0;


CIb=newbeta.prop.ci;
CIb(isnan(CIb))=0;
IQb=newbeta.prop.iq;
IQb(isnan(IQb))=0;
fac=5;
betas=(CIb>CIa)*1;
[m,n]=size(newalpha);
EBSD=zeros(m,n,3);
b1=fillmissing(newbeta.orientations.phi1,'movmedian',fac);
b1(isnan(b1))=0;
a1=fillmissing(newalpha.orientations.phi1,'movmedian',fac);
a1(isnan(a1))=0;
EBSD(:,:,1)=a1.*(1-betas)+b1.*betas;
b2=fillmissing(newbeta.orientations.Phi,'movmedian',fac);
b2(isnan(b2))=0;
a2=fillmissing(newalpha.orientations.Phi,'movmedian',fac);
a2(isnan(a2))=0;
EBSD(:,:,2)=a2.*(1-betas)+b2.*betas;
b3=fillmissing(newbeta.orientations.phi2,'movmedian',fac);
b3(isnan(b3))=0;
a3=fillmissing(newalpha.orientations.phi2,'movmedian',fac);
a3(isnan(a3))=0;
EBSD(:,:,3)=a3.*(1-betas)+b3.*betas;
CI=CIa.*(1-betas)+CIb.*betas;
IQ=IQa.*(1-betas)+IQb.*betas;
else
    if size(ebsd.unitCell,1)~= 6
        [newalpha,~]=gridify(ebsd,ebsd.extend);
        scale=1;
    else
        [newalpha,~,~,axis,~]=hexify(ebsd,ebsd.extend);
        scale=axis*sqrt(3)/2+(~axis)*2/sqrt(3);
    end
    

%CIa=newalpha.prop.confidenceindex;%for osc files
CIa=newalpha.prop.ci;
CIa(isnan(CIa))=0;

%IQa=newalpha.prop.imagequality;%for osc files
IQa=newalpha.prop.iq;
IQa(isnan(IQa))=0;

if LOGbeta
    betas=zeros(size(CIa));
else
    betas=ones(size(CIa));
end
[n,m]=size(newalpha);
EBSD=zeros(n,m,3);
EBSD(:,:,1)=fillmissing(newalpha.orientations.phi1,'nearest',1);
EBSD(:,:,2)=fillmissing(newalpha.orientations.Phi,'nearest',1);
EBSD(:,:,3)=fillmissing(newalpha.orientations.phi2,'nearest',1);
CI=CIa;
IQ=IQa;
end
%EBSD = fillmissing(EBSD,'nearest',1);

end