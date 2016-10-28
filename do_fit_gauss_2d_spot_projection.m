function [A,sx,sy,Intensity,bk,res] = do_fit_gauss_2d_spot_projection(it,sx,sy,sz,zmax,fname_ms2,Lx,Ly,window)
%[A,sx,sy,bk,res,zmean,nz]
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%fprintf(1,'hello!! %d\n',nuc);
%%
set(0,'DefaultTextFontname', 'Comic')
set(0,'DefaultTextFontSize', 12)

% Change default axes fonts.
set(0,'DefaultAxesFontName', 'Comic')
set(0,'DefaultAxesFontSize', 12)


stack=zeros(Lx,Ly,zmax);
zmean=0;
nz=0;

sz=round(sz);



Im=zeros(Lx,Ly);
Im(round(sy),round(sx))=1;
Imd=imdilate(Im,strel('disk',window));
%figure,imshow(Imd);
if sz>1 && sz<zmax
    for iz=sz-1:sz+1
        iis=sz+(it-1)*zmax;
        II=imread(fname_ms2,iis);
        stack(:,:,iz)=II.*uint8(Imd);
    end
end

if sz==1
    for iz=1:sz+1
        iis=sz+(it-1)*zmax;
        II=imread(fname_ms2,iis);
        stack(:,:,iz)=II.*uint8(Imd);
    end
end
if sz==zmax
    for iz=sz-1:zmax
        iis=sz+(it-1)*zmax;
        II=imread(fname_ms2,iis);
        stack(:,:,iz)=II.*uint8(Imd);
    end
end

Imtofit=uint8(max(stack,[],3));


%%% compute the bacground


xf=round(sx);
yf=round(sy);

d(1)=(xf-window/2);
d(2)=(yf-window/2);
d(3)=Ly-(xf+window/2);
d(4)=Lx-(yf+window/2);



if d(1)>=0 && d(2)>=0 && d(3)>=0 && d(4)>=0
    Ifit_raw = imcrop(Imtofit,[xf-window/2 yf-window/2 window window]);
    %figure,subplot(2,2,1),imshow(Ifit_raw);
    
else
    dmin=min(d(:));
    iii=find(d==dmin);
    window=window-d(iii);
    Ifit_raw = imcrop(Imtofit,[xf-window/2 yf-window/2 window window]);
    %figure,subplot(2,2,1),imshow(Ifit_raw);
end
    
%%

if window>=9
    Intensity=0;
    bkgr=0;
    normr=0;
    for i=1:size(Ifit_raw,1)
        for j=1:size(Ifit_raw,2)
            if sqrt( (i-window/2)^2+(j-window/2)^2 )>=8
                bkgr=bkgr+double(Ifit_raw(i,j));
                normr=normr+1;
            end
            if sqrt( (i-window/2)^2+(j-window/2)^2 )<8
                Intensity=Intensity+double(Ifit_raw(i,j));
                
            end
        end
    end
    
    bkgr=bkgr/normr;
    
    
    
    
    
    
    
    Ifit=Ifit_raw-bkgr;
    
    %subplot(2,2,2),imshow(Ifit);
    %%
    %do the fit with a gaussian
    
    [ny,nx] = size(Ifit);
    [mx,my] = meshgrid(1:nx,1:ny);
    
    
    %%% 5 parameters 
    p0 = [100,   nx/2,   nx/5,     ny/2,     ny/5];
   
    %upper and lower limits
    lb = [0,    nx/2*0.05,  0,    ny/2*0.05,  0  ];
    ub = [255, nx/2*1.95, nx/3,  ny/2*1.95,  ny/3];
    
    
    options = optimset('lsqnonlin');
    %options = optimset;
    %options = optimset('Algorithm',{'levenberg-marquardt',0.005});
    % %options.Display = 'iter';
    %options = optimset('Jacobian','on')
    
    
    myfun = @(p) p(1)*exp(-(((mx-p(2)).^2/(2*p(3)^2))+((my-p(4)).^2/(2*p(5)^2))))- double(Ifit);
    [p,resnorm] = (lsqnonlin(myfun,p0,lb,ub,options));
    %[p,resnorm] = (lsqnonlin(myfun,p0,lb,ub));
    A=p(1);
    sx=p(3);
    sy=p(5);
    bk=bkgr;
    res=resnorm;
    
    %     %%% 6 parameters
    %
    %     % starting values:
    %     %p0=[amplitude,meanx,sigx,meany,sigy,angle]
    %     p0 = [100,nx/2,10,ny/2,10,0.7];
    %     %upper and lower limits
    %     lb = [0,  nx/2*0.1,  1,  ny/2*0.1,  1,   0];
    %     ub = [250, nx/2*1.9, nx,  ny/2*1.9,  ny, 2*pi];
    %     myfun =@(p) p(1) *exp( -0.5/(1-p(6)^2) * ( (mx-p(2)).^2/p(3)^2 + (my-p(4)).^2/p(5)^2 - 2*p(6).*(mx-p(2)).*(my-p(4))./(p(3)*p(5)) ) ) - double(Ifit);
    
    %options = optimset('Algorithm',{'levenberg-marquardt',0.005});
    
    
    % % %% just plots
    %
    TwoDGaussian = p(1)*exp(-(((mx-p(2)).^2/(2*p(3)^2))+((my-p(4)).^2/(2*p(5)^2))));
    
%      subplot(2,2,4);
%     % scatter3(mx(:),my(:),TwoDGaussian(:),20,[0 0 0],'filled');
%     hold on;
%     m=mesh(mx,my,TwoDGaussian,'LineWidth',0.7);
%     
%     title( {sprintf('A= %g  sx=%g  sy=%g',p(1),p(3),p(5) ) ; sprintf('I2d=%g  Act=%g',Intensity,p(1)*p(3)*p(5) )} ) ;
%     
%     
%     set(m,'facecolor','none'); %makes transparent the mesh
%     
%     hold on
%     %II = double(Ifit);
%     II = double(Ifit_raw);
%     scatter3(mx(:),my(:),II(:),'filled');
    
    fprintf(1,' A = %f sx = %f sy = %f I2d= %d  bk=%g \n',p(1),p(3),p(5),Intensity,bkgr);
    
else
    A=-999;
    sx=-999;
    sy=-999;
    bk=-999;
    res=-999;
end


end

