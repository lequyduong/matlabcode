%% This file is derived from Bi2Te3
clear all
close all


NLayer = 2;
Mz=0.32;

kr=linspace(0,.1,31)';
kth=linspace(0,2*pi,25);

% kr=linspace(0,1.3,132)';
% kth=linspace(0,2*pi,61);
kx=kr*cos(kth);
ky=kr*sin(kth);
%[kx ky]=meshgrid(kkone,kkone);
kz=zeros(size(kx));

kk=[kx(:) ky(:) kz(:)];

if 0==0
%%%%lattice
abc=[8.560462  8.560462 57.882334]*0.529177; %% angstrom
BR=[sqrt(3)/2 0 0
    -1/2      1 0
    0         0 1];
BRrec=diag(abc)*BR;
BZrec=inv(BRrec)'*2*pi; %% reciprocal lattice (colunm vectors)
BZ=diag(abc)*BZrec/2/pi;
%%%%%%%%%%%%%

krec=kk/diag(abc)*2*pi; 
kxrec=reshape(krec(:,1),size(kx));
kyrec=reshape(krec(:,2),size(kx));
kzrec=reshape(krec(:,3),size(kx));
krh=krec/(BZrec');
krh=kk/BZ'; %% same as above

krh1=mod(krh,1);

Eplane1=zeros(size(kxrec));
Eplane2=zeros(size(kxrec));
wso1=zeros(size(kxrec));
wso2=zeros(size(kxrec));
wso2vecx=zeros(size(kxrec));
wso2vecy=zeros(size(kxrec));
wso2vecz=zeros(size(kxrec));
wso1vecx=zeros(size(kxrec));
wso1vecy=zeros(size(kxrec));
wso1vecz=zeros(size(kxrec));
end
% kF = 0.1;
% Nk=31;
% kxrec=linspace(-kF,kF,Nk);
% kyrec=linspace(-kF,kF,Nk);

sigx=[0 0 1 0;0 0 0 1;1 0 0 0;0 1 0 0 ];
sigy=[0 0 -1i 0;0 0 0 -1i;1i 0 0 0;0 1i 0 0 ];
sigz=[1 0 0 0; 0 1 0 0; 0 0 -1 0; 0 0 0 -1];

ssx = zeros(4*NLayer,4*NLayer);
ssy = zeros(4*NLayer,4*NLayer);
ssz = zeros(4*NLayer,4*NLayer);
firstQL=1:4;
ssx(firstQL,firstQL)=sigx;
ssy(firstQL,firstQL)=sigy;
ssz(firstQL,firstQL)=sigz;

H0 = genH(0,0,NLayer,Mz );
[E0, Wave0] = eig(H0);
EgammaD=E0(2*NLayer,2*NLayer);
EgammaU=E0(2*NLayer+1,2*NLayer+1);


for ii=1:prod(size(kxrec))
    H=genH(kxrec(ii),kyrec(ii),NLayer,Mz ); %% kx,ky interchanged!!
[eigvec, eigenvalues]=eig(H);
Eplane1(ii)=eigenvalues(2*NLayer,2*NLayer);
%Eplane2(ii)=eigenvalues(2,2)+tpara(5);
Eplane2(ii)=eigenvalues(2*NLayer+1,2*NLayer+1);
 if (Mz == 0.2424)&&(kxrec(ii) == 0)&&(kyrec(ii) == 0)
    Eplane2(ii)=Eplane1(ii); 
 end
% wso1(ii)=abs(eigvec(1,1)).^2;
% wso2(ii)=abs(eigvec(1,2)).^2;
%Eplane2(ii)=funEk(tpara,kxrec(ii),kyrec(ii));
if (Mz ~=0)
    wso2vecx(ii)=-eigvec(:,2*NLayer+1)'*ssy*eigvec(:,2*NLayer+1);
    wso2vecy(ii)=-eigvec(:,2*NLayer+1)'*ssx*eigvec(:,2*NLayer+1);
    wso2vecz(ii)=eigvec(:,2*NLayer+1)'*ssz*eigvec(:,2*NLayer+1);
    wso1vecx(ii)=-eigvec(:,2*NLayer)'*ssy*eigvec(:,2*NLayer);
    wso1vecy(ii)=-eigvec(:,2*NLayer)'*ssx*eigvec(:,2*NLayer);
    wso1vecz(ii)=eigvec(:,2*NLayer)'*ssz*eigvec(:,2*NLayer);
end
if (Mz == 0)
    wso2vecx(ii)=-eigvec(:,2*NLayer+1)'*ssy*eigvec(:,2*NLayer+1)-eigvec(:,2*NLayer+2)'*ssy*eigvec(:,2*NLayer+2);
    wso2vecy(ii)=-eigvec(:,2*NLayer+1)'*ssx*eigvec(:,2*NLayer+1)-eigvec(:,2*NLayer+2)'*ssx*eigvec(:,2*NLayer+2);
    wso2vecz(ii)=eigvec(:,2*NLayer+1)'*ssz*eigvec(:,2*NLayer+1)+eigvec(:,2*NLayer+2)'*ssz*eigvec(:,2*NLayer+2);
    wso1vecx(ii)=real(-eigvec(:,2*NLayer)'*ssy*eigvec(:,2*NLayer)-eigvec(:,2*NLayer-1)'*ssy*eigvec(:,2*NLayer-1));
    wso1vecy(ii)=-eigvec(:,2*NLayer)'*ssx*eigvec(:,2*NLayer)-eigvec(:,2*NLayer-1)'*ssx*eigvec(:,2*NLayer-1);
    wso1vecz(ii)=eigvec(:,2*NLayer)'*ssz*eigvec(:,2*NLayer)+eigvec(:,2*NLayer-1)'*ssz*eigvec(:,2*NLayer-1);
end
end
wso2vec(:,:,1)=wso2vecx;
wso2vec(:,:,2)=wso2vecy;
wso2vec(:,:,3)=wso2vecz;
wso1vec(:,:,1)=wso1vecx;
wso1vec(:,:,2)=wso1vecy;
wso1vec(:,:,3)=wso1vecz;

%Esurfhigh=Eplane2-min(Eplane2(:));
%Esurfhigh=(Esurfhigh+fliplr(Esurfhigh))/2*1000;
%Esurfhigh=(Esurfhigh+(Esurfhigh))/2*1000;
Esurfhigh=Eplane2; %% meV
Esurflow=Eplane1; %% meV

wso2(1,:)=0.5;
wsohigh=(wso2-0.5)*2*100;
wso1(1,:)=0.5;
wsolow=(wso1-0.5)*2*100;

if 1==0
figure;
%contourf(kx,ky,Eplane)
pcolor(kx,ky,Eplane2);shading flat
grid on
axis equal
figure;
contourf(kxrec,kyrec,Eplane2)
%pcolor(kx,ky,Eplane);shading flat
grid on
axis equal
end
%----------------------------------
if 1==0
figure;
quiver(kxrec,kyrec,wso2vec(:,:,1),wso2vec(:,:,2),0.1)
return
end
%----------------------------------
if 1==0
figure;hold on;axis equal

krange=-0.2:0.01:0.2;
[kxx,kyy]=meshgrid(krange,krange);
spinxx=griddata(kxrec,kyrec,wso2vec(:,:,1),kxx,kyy);
spinyy=griddata(kxrec,kyrec,wso2vec(:,:,2),kxx,kyy);
quiver(kxx,kyy,spinxx,spinyy,0.5)
return
end

%--------------------------------------------------------------------------
% Ecut=0.25;
% wsohigh(Esurfhigh>Ecut)=NaN;
% Esurfhigh(Esurfhigh>Ecut)=NaN;
kabs=sqrt(kxrec.^2+kyrec.^2);
kcut=0.9;
wsohigh(kabs>kcut)=NaN;
Esurfhigh(kabs>kcut)=NaN;
% wcut=30;
% Esurfhigh(wsohigh>wcut)=NaN;
% wsohigh(wsohigh>wcut)=NaN;
% Esurfhighmesh=Esurfhigh;
% kcutmesh=0.05;
% Esurfhighmesh(kabs>kcutmesh)=NaN;

figure;
%h1=surf(kxrec,kyrec,Esurfhigh,wsohigh); shading interp
h1=surf(kxrec,kyrec,Esurfhigh,wso2vecz); shading interp
alpha(0.5);
%set(h1,'FaceColor',[1 0 0],'FaceAlpha',1.0);

%colormap hot;
% h1=surf(kxrec,kyrec,Esurfhigh); shading interp
% alpha(0.25);
hold on
%h2=surf(kxrec,kyrec,Esurflow,wsolow); shading interp
h2=surf(kxrec,kyrec,Esurflow,wso1vecz); shading interp
alpha(0.5);
%mesh(kxrec,kyrec,Esurfhighmesh,wsohigh); shading interp
Esurfhigh=(Esurfhigh+fliplr(Esurfhigh))/2;

%contour(kxrec,kyrec,Esurfhigh,[0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.38]*1000,'k')
%contour(kxrec,kyrec,Esurfhigh,'k')
%mesh(kxrec,kyrec,Esurfhigh,wsohigh); shading interp
       daspect([1 1 1.2])
        xlim([-1 1]*0.3);
        ylim([-1 1]*0.3);
       zlim([-0.4 0.25]);
       view([-2 9])
%       view([15 24])
     % camlight; 
lighting phong
%pcolor(kxrec,kyrec,Esurfhigh); shading interp
%pcolor(kxrec,kyrec,wsohigh); shading interp
%xlabel('k_x [1/angstrom]');
%ylabel('k_y [1/angstrom]');
%zlabel('Energy [eV]');
% cc=flipud(pink);
% cc(:,3)=0.6;
% colormap(flipud(cc))
caxis([-1 1]*.95);
%colormap
%hc=COLORBAR('horiz');
%set(hc,'position',[0.25    0.15    0.4    0.01]);
%set(hc,'xtick',[0 10 20 30]);
% load datacolormapa mapa
% colormap(mapa);
%COLORBAR('horiz')
%grid on
set(gcf,'Color','w')

if 1==0
savefname=['fig2conev4'];
%print(gcf,'-depsc',[savefname '.eps']);
print(gcf,'-djpeg90','-r600',[savefname '.jpg']);
%saveas(gcf,[savefname '.fig'],'fig');
end

%--------------------------------
%EBselect=[0.05 0.15 0.2]*1000;
EBselect=([0.04 0.14 0.26]); %% relative to Dirac point after TPT
%EBselect=([0.05 0.10 0.20]);
expEBeV=[0 60 120 230];

%figure('Position',[50 50 1050 300],'paperposition',[0.2500    2.5000    6    2.5]);
%figure('Position',[50 50 600 600]);
scolor='m';
for ieb=[1:3]
%openfig(['..\FromSuYang20100327FS/FS' num2str(expEBeV(ieb)) '.fig']);
%    subplot(1,4,ieb);
    eb=EBselect(ieb);
hold on;
[atemp hcon]=contour(kxrec,kyrec,Esurfhigh,[eb eb],scolor);
delete(hcon);
%grid on
if ieb==1
    fs0rec=atemp(:,2:end);
end
fsrec{ieb}=atemp(:,2:end);
    fsrectemp=fsrec{ieb};
%    fsrectemp(3,:)=0;
    fsrectemp(3,:)=eb;
%    fsrectemp(:,fsrectemp(2,:)<0)=[];
plot3(fsrectemp(1,:),fsrectemp(2,:),fsrectemp(3,:),'linewidth',2,'color',[0.5 0.5 1]);
clear spin
spin(1,:)=griddata(kxrec,kyrec,wso2vec(:,:,1),fsrectemp(1,:),fsrectemp(2,:));
spin(2,:)=griddata(kxrec,kyrec,wso2vec(:,:,2),fsrectemp(1,:),fsrectemp(2,:));
spin(3,:)=griddata(kxrec,kyrec,wso2vec(:,:,3),fsrectemp(1,:),fsrectemp(2,:));
cmap=colormap;
czz=linspace(-1,1,64);
%spincolor=interp1(czz,cmap,spin(3,:)/max(abs(spin(3,:))));
spincolor=interp1(czz,cmap,spin(3,:));
if ieb==1
    for iplot=1:3:length(fsrectemp(1,:))
    
arrow3D(fsrectemp(:,iplot),spin(:,iplot)*0.08,'r');
    end
end %ieb==1
if ieb==2
for iplot=1:3:length(fsrectemp(1,:))
    
arrow3D(fsrectemp(:,iplot),spin(:,iplot)*0.12,'k');
end
end%%ieb>1
%
if ieb==3
    for iplot=1:3:length(fsrectemp(1,:))
    
arrow3D(fsrectemp(:,iplot),spin(:,iplot)*0.15,'k');
    end
end %ieb==1
end %ieb
%%%%%%%%%%%%%%%%%%%%%%%%lower cone spin arrow
hold on;

if 0==0
    spinU=[0 0 -1];
    spinD=[0 0 1];
end

if 0==0
    EBselect=([-0.1 -0.13 -0.17]); %% relative to Dirac point after TPT
    %EBselect=([-0.12 -0.17 -0.14]); % zero
    
expEBeV=[0 60 120 230];

%figure('Position',[50 50 1050 300],'paperposition',[0.2500    2.5000    6    2.5]);
%figure('Position',[50 50 600 600]);
scolor='m';
for ieb=[1:3]
%openfig(['..\FromSuYang20100327FS/FS' num2str(expEBeV(ieb)) '.fig']);
%    subplot(1,4,ieb);
    eb=EBselect(ieb);
hold on;
[atemp hcon]=contour(kxrec,kyrec,Esurflow,[eb eb],scolor);
delete(hcon);
%grid on
if ieb==1
    fs0rec=atemp(:,2:end);
end
fsrec{ieb}=atemp(:,2:end);
    fsrectemp=fsrec{ieb};
%    fsrectemp(3,:)=0;
    fsrectemp(3,:)=eb;
%    fsrectemp(:,fsrectemp(2,:)<0)=[];
plot3(fsrectemp(1,:),fsrectemp(2,:),fsrectemp(3,:),'linewidth',2,'color',[0.5 0.5 1]);
clear spin
spin(1,:)=griddata(kxrec,kyrec,wso1vec(:,:,1),fsrectemp(1,:),fsrectemp(2,:));
spin(2,:)=griddata(kxrec,kyrec,wso1vec(:,:,2),fsrectemp(1,:),fsrectemp(2,:));
spin(3,:)=griddata(kxrec,kyrec,wso1vec(:,:,3),fsrectemp(1,:),fsrectemp(2,:));
cmap=colormap;
czz=linspace(-1,1,64);
%spincolor=interp1(czz,cmap,spin(3,:)/max(abs(spin(3,:))));
spincolor=interp1(czz,cmap,spin(3,:));
if ieb==1
    for iplot=1:3:length(fsrectemp(1,:))
    
arrow3D(fsrectemp(:,iplot),spin(:,iplot)*0.09,'k');
    end
end %ieb==1
if ieb==2
for iplot=1:3:length(fsrectemp(1,:))
    
arrow3D(fsrectemp(:,iplot),spin(:,iplot)*0.15,'r');
end

end%%ieb>1
if ieb==3
for iplot=1:3:length(fsrectemp(1,:))
    
arrow3D(fsrectemp(:,iplot),spin(:,iplot)*.3,'r');
end

end%%ieb>1
end %ieb


end

if ieb==1
    fs0rec=atemp(:,2:end);
end
fsrec{ieb}=atemp(:,2:end);
    fsrectemp=fsrec{ieb};
%    fsrectemp(3,:)=0;
    fsrectemp(3,:)=eb;
%    fsrectemp(:,fsrectemp(2,:)<0)=[];
plot3(fsrectemp(1,:),fsrectemp(2,:),fsrectemp(3,:),'linewidth',2,'color',[0.5 0.5 1]);
clear spin
spin(1,:)=griddata(kxrec,kyrec,wso1vec(:,:,1),fsrectemp(1,:),fsrectemp(2,:));
spin(2,:)=griddata(kxrec,kyrec,wso1vec(:,:,2),fsrectemp(1,:),fsrectemp(2,:));
spin(3,:)=griddata(kxrec,kyrec,wso1vec(:,:,3),fsrectemp(1,:),fsrectemp(2,:));
cmap=colormap;
czz=linspace(-1,1,64);
%spincolor=interp1(czz,cmap,spin(3,:)/max(abs(spin(3,:))));
spincolor=interp1(czz,cmap,spin(3,:));
% for iplot=1:3:length(fsrectemp(1,:))
%     
% arrow3D(fsrectemp(:,iplot),spin(:,iplot)*0.1,spincolor(iplot,:));
% end


%set(hcon,'linewidth',2,'color',[0.5 0.5 1]);
%axis equal
% xlim([-1 1]*0.25);
% ylim([-1 1]*0.25);
%xlim([-1 1]*0.3);
%ylim([-1 1]*0.3);
% set(gca,'YAxisLocation','right');
%set(gca,'xtick',[-0.2 -0.1 0 0.1 0.2]);
%set(gca,'ytick',[-0.2 -0.1 0 0.1 0.2]);
%set(gca,'ztick',[]);
% set(gca,'xticklabel','');
% set(gca,'yticklabel','');
% axespos=get(gca,'position');
% axespos(3:4)=axespos(3:4)*1.05;
% set(gca,'position',axespos);
%set(gca,'ticklength',[0.01 0.025]*5)
%box on
%set(gca,'linewidth',2)
%view(-18,76);
%ha=arrow([0 0],[0.1 0.1])
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gca,'Visible','off')
if 1==0
savefname=['figSuYangspinarrow'];
%print(gcf,'-depsc',[savefname '.eps']);
print(gcf,'-djpeg90','-r600',[savefname '.jpg']);
%saveas(gcf,[savefname '.fig'],'fig');
end
return


figure;
scolor='b';
EBselect=[-0.05 0.05 0.15 0.2];
for ieb=1:4
    subplot(3,1,4-ieb);
    eb=EBselect(ieb);
hold on;
[atemp hcon]=contour(kxrec,kyrec,Esurfhigh,[eb eb],scolor);
set(hcon,'linewidth',4);
axis equal
% xlim([-1 1]*0.15);
% ylim([-1 1]*0.15);
box on
set(gca,'YAxisLocation','right');
set(gca,'xtick',[-0.1 0 0.1]);
set(gca,'ytick',[-0.1 0 0.1]);
set(gca,'xticklabel','');
set(gca,'yticklabel','');
axespos=get(gca,'position');
axespos(3:4)=axespos(3:4)*1.25;
set(gca,'position',axespos);
%ha=arrow([0 0],[0.1 0.1])
end

% figure('Position',[50 50 200 200]);
% kytemp=[0:0.001:0.3];
% kxtemp=zeros(size(kytemp));
% spinz=griddata(kxrec,kyrec,wsohigh,kxtemp,kytemp);
% ebz=griddata(kxrec,kyrec,Esurfhigh,kxtemp,kytemp);
% plot(ebz,spinz,'xx','linewidth',3)

if 1==0
savefname=['fig2FSv4'];
%print(gcf,'-depsc',[savefname '.eps']);
print(gcf,'-djpeg90','-r600',[savefname '.jpg']);
%saveas(gcf,[savefname '.fig'],'fig');
end
