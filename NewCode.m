% mseed to m file
% This code describes the STA/LTA method for removing the blasting
% activities and Erathquake noises from the original signal
%Written by Mostafa Ebrahimi
clc
close all
clear all
sps=100;
VH=[];
vohs=[];
cd('C:\Users\mebrahi\Desktop\New folder (2)');
Zdata=dir('C:\Users\mebrahi\Desktop\New folder (2)\*EHZ.mseed');
Ndata=dir('C:\Users\mebrahi\Desktop\New folder (2)\*EHN.mseed');
Edata=dir('C:\Users\mebrahi\Desktop\New folder (2)\*EHE.mseed');

for k=1:length(Zdata)
    clear  spZ spN spE l1 l 
    kz=Zdata(k).name;
    kn=Ndata(k).name;
    ke=Edata(k).name;

z=rdmseed(kz);
n=rdmseed(kn);
e=rdmseed(ke);
[z1 z2]=size(z);
[n1 n2]=size(n);
[e1 e2]=size(e);
sz=[];
sn=[];
se=[];
for i=1:z2
    sz=[sz;z(1,i).d];
end

for j=1:n2
    sn=[sn;n(1,j).d];
 end

for q=1:e2
    se=[se;e(1,q).d];
end

      l1=10000;  
         l=size(sz)/l1-1;
     for j3=1:floor(l(:,1));
%          %---------------------
        Fs=sps; Wn=[0.5 12]/(100);n=4;
        [b,a]=butter(n,Wn);
        
        dcn=mean(sn((l1*(j3-1)+1):(j3*l1)));
        dce=mean(se((l1*(j3-1)+1):(j3*l1)));
        dcz=mean(sz((l1*(j3-1)+1):(j3*l1)));
        sn2=sn-dcn;
        se2=se-dce;
        sz2=sz-dcz;
% %         velz=filter(b,a,sz2((l1*(j3-1)+1):(j3*l1)));
% %         vele=filter(b,a,se2((l1*(j3-1)+1):(j3*l1)));
% %         veln=filter(b,a,sn2((l1*(j3-1)+1):(j3*l1)));
% % 
% % 
% %         nfft=2^nextpow2(length(fft(velz)));
% %         psdz = fft(velz,nfft)/(length(fft(velz)));
% %         psdn = fft(veln,nfft)/(length(fft(veln)));
% %         psde = fft(vele,nfft)/(length(fft(vele)));
% %         f = sps/2*linspace(0,1,nfft/2);
% %         psdz=2*abs(psdz(1:nfft/2));
% %         psdn=2*abs(psdn(1:nfft/2));
% %         psde=2*abs(psde(1:nfft/2));
% %         %---------------------------------------------------------
% %         spZ=psdz;
% %         spN=psdn;
% %         spE=psde;
% %         vohs(:,(k+j3-1))=(spZ)./sqrt((spN.^2+spE.^2)./2);
% %         vohs(:,(k+j3-1))=smooth(vohs(:,(k+j3-1)),100);
% %         fL= f;
% %      end
% % end
% % 
% % figure(40)
% % title('SKHZ')
% % hold on
% % VH=[VH,vohs];
% % vhmean=mean(VH,2);
% % z=ones(size(fL));
% % plot(fL,z);
% % k2=(std((vohs)'))';
% % plot(fL,vhmean-k2,'r--');
% % plot(fL,vhmean+k2,'r--');
% % plot(fL,vhmean,'k','linewidth',2);
% % xlim([.5,7]);
% % ylim([0,3]);
% % xlabel('Frq (Hz)');
% % ylabel('V/H')
% % figure(41)
% % plot(fL,vhmean,'k','linewidth',2);
% % hold on
% % plot(fL,z);
% % xlim([.75,7]);
% % ylim([0,3]);
% % figure(42)
% % plot(fL,vhmean)
     end
end


       sz2=sz2(70000:100000);
         se2=se2(70000:100000);
          sn2=sn2(70000:100000);
% % % % % % % % % % pol
t=0:0.01:(size(sz2)/100)-0.01;
t=t';
z=sz2;
n=sn2;
e=se2;

flp=0.01;
fhi=0.07;
if flp > 0;
   w = [flp fhi];
   [b,a]=butter(4,w);
%   
%    e1=filter(b,a,e); 
%    n1=filter(b,a,n); 
%    z1=filter(b,a,z); 
%    e=e1; 
%    n=n1; 
%    z=z1; 
   % 
 f(1)=figure('name','Filterd seismograms');  % plot the filtered data 
 subplot(3,1,1); 
 plot(t,e); 
 xlabel('time (s)'); 
 title('E')
%  ylabel(strcat('EW Comp. of ',IDe)); 
 subplot(3,1,2); 
 plot(t,n); 
 xlabel('time (s)'); 
  title('N')
%  ylabel(strcat('NS Comp. of ',IDn)); 
 subplot(3,1,3); 
 plot(t,z); 
 xlabel('time (s)'); 
  title('Z')
%  ylabel(strcat('Z comp. of ',IDz)); 
else; 
end; 

% 
%   Moving window loop 
%
% delt=1/sps;
delt=1/sps;
ttot=600;
twin=1;
npts=ttot/delt + 1; 

npts1=fix(((size(z)-1)*1/sps)/delt) + 1;   %  total number of samples to analyze 
nwin=fix(twin/delt) + 1 ;   %  number of samples in a time window 
npshift=fix(twin/(2*delt))+1 ; % number of samples to shift over 
kfin=fix((npts1-nwin)/(npshift+1))+1; % number of time windows considered 
mxde1=0; 
mxde2=0; 
mxde3=0; 
for k=1:kfin; 
   nwinst=(k-1)*(npshift-1)+1;  % start of time window 
   nwinfn=nwinst+nwin-1;    % end of time window 
   a=[];
   a(:,1)=e(nwinst:nwinfn);
   a(:,2)=n(nwinst:nwinfn);
   a(:,3)=z(nwinst:nwinfn);   
   %a 
   c=(1/length(e))*(a'*a);         % correlation matrix 
   %c 
   [v,d]=eig(c);      % eigen vectors and eigen values
   d ;
   v;
   
          % sort the eigenvalues and eigenvectors 
          % it seems matlab sort the eign values itself
    
   ang1(k)=atan2(v(1,1),v(2,1)) * 180/pi; % Azimuth for first  eigenvalue 
   ang2(k)=atan2(v(1,2),v(2,2)) * 180/pi; % Azimuth for second eigenvalue 
   ang3(k)=atan2(v(1,3),v(2,3)) * 180/pi; % Azimuth for third  eigenvalue
   dip1(k)=atan2(v(3,1),((v(1,1)^2+v(2,1)^2)^0.5)) * 180/pi; % Dip for first  eigenvalue 
   dip2(k)=atan2(v(3,2),((v(1,2)^2+v(2,2)^2)^0.5)) * 180/pi; % Dip for second eigenvalue 
   dip3(k)=atan2(v(3,3),((v(1,3)^2+v(2,3)^2)^0.5)) * 180/pi; % Dip for third  eigenvalue
   
   
   de1(k)=d(1) ;
   de2(k)=d(2); 
   de3(k)=d(3); 
   % find the maximum values 
   mxde1=max(mxde1,de1(k));   
   mxde2=max(mxde2,de2(k)); 
   mxde3=max(mxde3,de3(k)); 
   %angle from the vertical
 vang1(k)=acos(abs(v(3,1)))* 180/pi;   
 vang2(k)=acos(abs(v(3,2)))* 180/pi;  
 vang3(k)=acos(abs(v(3,3)))* 180/pi;
 
 t2(k)=delt*(nwinst-1);    % assign time for this window to the window start 
end; 
% 
f(2)=figure('name','Eigenvalues and Inferred Azimuth'); 
subplot(3,1,1); 
% plot(t2,de1,'-or',t2,de2,'-dg',t2,de3,'-+b');
plot(t2,de1,'-+b');
xlabel('time sec'); 
ylabel('eigenvalues'); 
subplot(3,1,2); 
% plot(t2,ang1,'-or',t2,ang2,'-dg',t2,ang3,'-+b'); 
plot(ang3,'-+b'); 
xlabel('time sec'); 
ylabel('Azimuth '); 
subplot(3,1,3); 
% plot(t2,vang1,'-or',t2,vang2,'-dg',t2,vang3,'-+b'); 
plot(t2,vang3,'-+b');
xlabel('time sec'); 
ylabel('incidence angle '); 
% 


%  Rose plot 
f(3)=figure('name','Azimuth Distribution'); 
subplot(2,3,1); 
title('Az- Largest Eig'); 
rose(ang1*pi/180,100); 
subplot(2,3,2); 
title('Az- Intermediate Eig'); 
rose(ang2*pi/180,100); 
subplot(2,3,3); 
title('Az- Smallest Eig'); 
rose(ang3*pi/180,100); 
subplot(2,3,4); 
% title('Azimuth - Largest Eigenvalue,50% Threshold'); 
rose(ang1*pi/180,100); 
subplot(2,3,5); 
% title('Azimuth - Intermediate Eigenvalue,50% Threshold'); 
rose(ang2*pi/180,100); 
subplot(2,3,6); 
% title('Azimuth - Smallest Eigenvalue,50% Threshold'); 
rose(ang3*pi/180,100); 

r_e = stalta(e,50,500);
r_n = stalta(n,50,500);
r_z = stalta(z,50,500);

save pqfile.mat r_e r_n r_z t

 f(4)=figure('name','Filterd seismograms');  % plot the filtered data 
 subplot(3,1,1); 
 plot(t,r_e); 
 xlabel('time (s)'); 
 title('E')
%  ylabel(strcat('EW Comp. of ',IDe)); 
 subplot(3,1,2); 
 plot(t,r_n); 
 xlabel('time (s)'); 
  title('N')
%  ylabel(strcat('NS Comp. of ',IDn)); 
 subplot(3,1,3); 
 plot(t,r_z); 
 xlabel('time (s)'); 
  title('Z')
  
  for k=1:5
  saveas(f(k),sprintf('figure_%d.jpg',k))
end
  
  
  
  
  
  
  
  
  
 