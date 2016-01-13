clear all; clc
% LEO Rain att profile using Rec. ITU-R P.1853-1 

% Authors: Paulo Ferreira & Alexander Wyglinski, Worcester Polytechnic Institute (WPI), Worcester - MA
% January 2016
% Please cite this code using its DOI number
% Any quetions: email <prferreira@wpi.edu> <alexw@wpi.edu>

%Based on the ITU-R P.1853-1 implementation for GEO by Laurent CASTANET, Guillaume CARRIE & Nicolas Jeannin, ONERA, France



%ITU-R P.1853-1 input parameters :
% ****************
% 1) ccdf :            input CCDF, a suggested set of time percentages is  
%          [0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10] (%)
% 2) P0 :              Probability of the rain on the path. P0 can be well 
%                   approximated as P0(Lat,Lon) deriverd in ITU-R P.837 (%)
% 3) Ar :              Attenuation levels exceeded for percentages of time stored in ccdf (dB)
% ****************
%%
%Use ITU-R P.837 to obtain 2) P0 and rain rate CCDF. Using the rain rates, then use
%ITU P.618-12 to obtain 3) Ar. 

%ITU-R P.837 input parameters :
%1) and (Lat,Lon)

%Atwater Kent Laboratories coordinates (Worcester, MA, USA): 
%42ï¿½16'30.8"N 71ï¿½48'25.2"W
lat=42+((16+(30.8/60))/60);
lon=360-(71+((48+(25.2/60))/60));
%If negative or W lon=360-71.807030= 289.1250 (0-360)

%18ï¿½55'23.6"S 48ï¿½13'41.3"W

% y = [0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10];
y = [0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1, 2, 3, 5, 10];
[rr, P0]=itur_p837_5(y,lat,lon);   

%%
% %Convert rain rates CCDF into attenuation levels (dB) CCDF using ITU P.618-12
% 
% %Need all geometrical info for each point:
% %Given by AER report from STK (.csv Excel file imported with correction)
% load('ISS_angle.mat') %Elevation angle (degrees)
% 
% %rain height by ITU P839
% hr_prime=5;
% 
% %Slant path (distance from ground antenna to rain layer)
% d_rain=hr_prime./tand(angle);
% 
% %From ITU-R P.838-3
% %Specific rain attenuation coefficients for CIRCULAR POLARIZATION by ITU P838 
% %26GHz
% freq=26*10^9;
% k_h=0.1724;
% k_v=0.1669;
% a_h=0.9884;
% a_v=0.9421;
% 
% %circular polarization
% tau=45;
% 
% kay= (k_h+k_v+((k_h-k_v)*((cosd(angle)).^2)*cosd(2*tau)))/2;
% alpha_= ((k_h*a_h)+(k_v*a_v)+(((k_h*a_h)-(k_v*a_v))*((cosd(angle).^2))*cosd(2*tau)))/(2.*kay);
% 
% %Specif att
% gamma=kay.*((rr).^alpha_);
% 
% att=gamma.*d_rain;

%%
%Convert rain rates into attenuation levels (dB)using ITU P.618-12 

%For LEO use sectorization of elevation angles; and consider the first of
%each interval:
%[5-7;7-10;10-15;15-22;22-37;37-90]
% theta=[5 7 10 15 22 37];
% theta=[5 19 33 47 61 75];
% theta=[19];
theta=[5 8.869 14 21.37 33.48 56.48];
theta2=[5 8.869 14 21.37 33.48 56.48 83.41];

%Sec. 2.2.1.1
%Rrain height by ITU P839
hr_=5;%hr-hs

Ls=hr_./sind(theta);%Slant path

Lg=Ls.*(cosd(theta));

%From ITU-R P.838-3
%Specific rain attenuation coefficients for CIRCULAR POLARIZATION by ITU P838 
%26GHz
freq=26; %GHz
k_h=0.1724;%26GHz
k_v=0.1669;%26GHz
a_h=0.9884;%26GHz
a_v=0.9421;%26GHz

%circular polarization
tau=45;
kay= (k_h+k_v+((k_h-k_v).*((cosd(theta)).^2).*cosd(2*tau)))./2;
alpha_= ((k_h*a_h)+(k_v*a_v)+(((k_h*a_h)-(k_v*a_v)).*((cosd(theta).^2)).*cosd(2*tau)))./(2.*kay);
%Specif att for rr(0.01%) [for all elevation angles]
%FOR CIRCULAR POLARIZATION, GAMMA_R IS NOT FUNCTION OF ELEVATION ANGLE!!!
if tau==45
    kay=kay(1);
    alpha_=alpha_(1);
    gamma_R=(kay.*((rr(1)).^alpha_));%should be a matrix (% vs elev.angle)
end

%Chi (X)
if abs(lat)<36
    chi=36-abs(lat);
else
    chi=0;
end

% %Compute beta_ for each y. Repeat for all elev. angles:
beta_=zeros(length(y),length(theta));
att1=beta_;
% att(1,:)=A_001;
for jj=1:length(theta)%cols
    % Compute r0.01% for all % in ccdf; and for all elev angles  [f in GHz]
    r_001(jj)=1./(1+(0.78.*sqrt((Lg(jj).*gamma_R)/(freq)))-(0.38.*(1-exp(-2.*Lg(jj)))));      
    %csi
    csi(jj)=atand(hr_./(Lg(jj).*r_001(jj)));
    %Lr
    if csi(jj)>theta(jj)
        Lr(jj)=(Lg(jj).*r_001(jj))./cosd(theta(jj));
    else
        Lr(jj)=hr_/sind(theta(jj));
    end
    % v0.01%
    v_001(jj)=1./(1+((sqrt(sind(theta(jj)))).*((31.*(1-exp(-theta(jj)./(1+chi))).*((sqrt(Lr(jj).*gamma_R))./((freq)^2)))-0.45)));
    %Le
    Le(jj)=Lr(jj).*v_001(jj);
    % A 0.01% [dB]
    A_001(1,jj)=gamma_R.*Le(jj);
    att1(1,jj)=A_001(1,jj);
    for ii=2:1:length(y)%rows %y=p%(% exceeded)        
        %beta_
        if y(ii)>=1 || abs(lat)>=36
            beta_(ii,jj)=0;
        elseif y(ii)<1 && abs(lat)<36 && theta(jj)>=25
            beta_(ii,jj)=(-0.005*(abs(lat)-36));
        else
            beta_(ii,jj)=(-0.005*(abs(lat)-36))+1.8-(4.25*sind(theta(jj)));
        end
        %In this case all beta_ are zeros because of abs(lat)>=36
        att1(ii,jj)=A_001(1,jj)*( (y(ii)/0.01)^(-0.655+(0.033*log(y(ii)))-(0.045*log(A_001(1,jj)))-(beta_(ii,jj)*(1-y(ii))*sind(theta(jj)))) );                
    end
end


% % Plots
% y_ccdf=y./100;
% 
% semilogy(att(:,1),y_ccdf)
% hold on
% semilogy(att(:,2),y_ccdf)
% hold on
% semilogy(att(:,3),y_ccdf)
% hold on
% semilogy(att(:,4),y_ccdf)
% hold on
% semilogy(att(:,5),y_ccdf)
% hold on
% semilogy(att(:,6),y_ccdf)
% grid on

% 
% %% Compute Rain Att time series for each elevation angle
% 
% sampling_rate = 1; %Hz
% duration = 86400-1; %sec
% ccdf=y;
% Arain2=zeros(duration+1,length(theta));
% for jj=1:1:length(theta)
%     Ar=att(:,jj)';
%     [Arain, Time]=itur_p1853_0_annex1_s2(ccdf,Ar,P0,sampling_rate,duration);
%      Arain2(:,jj)=Arain'; 
% end

%%
%ITU P1853
sampling_rate = 1; %Hz
duration = 84600; %sec
% sampling_rate = (0.6503e-6)/2; %Hz
% duration = 512-1; %sec
D=duration;
Ts=1/sampling_rate;

itr=2;
att3=zeros((duration*sampling_rate)+1,6,itr);
for ij=1:itr
    rng('shuffle')
%Noise filter
% --------------------------------------------------------
if exist('seed','var') && ~isempty(seed)
   randn('state',seed);
end


beta = 2e-4; % (s-1) value recommended by ITU-R P.1853-0 annex 1 section 2.2 B.

Nadd = 200000;  

% Number of samples to synthesise
Ns = floor(D/Ts)+1;

n1 = randn(1,Ns+Nadd);		% Random synthesis of the WGN

% Definition of the low-pass filter
ro = exp(-beta*Ts);
b = sqrt(1-ro^2);
a = [1 -ro]; 

% Low-pass filtering of the WGN
Y = filter(b,a,n1); 
Y(1:Nadd) = [];

att2=att1;
%Iterate over all elev.angles:
for jj=1:length(theta)

% Estimation of  the Log-normal and Aoffset parameters of the ONERA-CNES channel model :
% ------------------------------------------------------------------------------------
ccdf=y;
Aexc=att2(:,jj)';
ccdf = ccdf(:);
Aexc = Aexc(:);
ccdf = ccdf(ccdf<=100);
Aexc = Aexc(ccdf<=100);

% Log-normal fit of the input CCDF between pinf and psup
ind_perc = (ccdf>=0 & ccdf<=P0);
coeff = polyfit(log(Aexc(ind_perc)),-sqrt(2).*erfinv(1-2*ccdf(ind_perc)/100),1);
lnm = -coeff(2)/coeff(1);
lns = -1/coeff(1);
% Assessment of the attenuation offset Aoffset correponding to P0 in the input CCDF
Aoffset = exp(lnm + lns*sqrt(2)*erfinv(1-2*P0/100));


% %% ITU-R  P.1057-4  Appendix 2 (finding sigma and m)
% % Least Squares fit to linear function ln(xi)=sigma*Zi + m
% % Transform (Gi,xi) to (Zi,ln(xi)); Zi=invQ(Gi)
% %G=ccdf/100; x=Att;
% idx=sum(ind_perc); %find the index
% Att=Aexc(1:idx);
% ccdf2=ccdf(1:idx);
% Z=sqrt(2).*erfinv(1-2*(ccdf(ind_perc)./100));
% %Ax=B solve for x; A=[X 1] B=Y
% A=[Z ones(idx,1)];
% B=log(Att(ind_perc));
% x = mldivide(A,B)
% new=100*0.5*(1-erf((log(Aexc)-x(2))./(sqrt(2)*x(1))));
% semilogx(ccdf,Ar,new,Ar);

%%
% Execution of the ONERA-CNES rain attenuation synthesiser :

% Log-normal ponderation of the WGN
att = exp(lnm+lns.*Y);

% Subtraction of the attenuation offset to the time series
att = max(att - Aoffset,0);

att3(:,jj,itr)=att;
end
end
% 
% time = (0:Ns-1)*Ts;
% %% Plots
% % Synthesized Time-Series
% figure ;
% plot(time,att3(:,1)) ; hold on 
% plot(time,att3(:,2)) ; hold on 
% plot(time,att3(:,3)) ; hold on 
% plot(time,att3(:,4)) ; hold on 
% plot(time,att3(:,5)) ; hold on 
% plot(time,att3(:,6)) ; grid on 
% title('\fontsize{14}\bfitur-p1853-0-annex1-s2 model') ;
% xlabel('\fontsize{14}Time (s)');
% ylabel('\fontsize{14}Rain attenuation (dB)');
% 
% 
% % CCDF Fit
% figure
% ccdfth = 100*0.5*(1-erf((log(att1(:,3))-lnm)/(lns*sqrt(2))));
% semilogx(ccdf,att1(:,3),ccdfth,att1(:,3),'LineWidth',2)
% title(['Original and approximated CCDF']);
% ylabel('Attenuation Exceeded (dB)')
% xlabel('Time (%)')
% legend('Original','LogN')

%%
%Composed time series

[maxval, maxind]=max(att3(:));
[maxxidx, maxyidx, maxzidx]=ind2sub(size(att3), maxind);

%Choose worst case matrix
Att=att3(:,:,maxzidx);

%Get 512 samples around the max in the worst case
if maxxidx>=512
    x1=maxxidx-256;
    x2=maxxidx+255;
else
    x1=1;
    x2=512;
end

Att=Att(x1:x2,:);
%%
%Conveting all fixed elev.angle time series for GEO to one for LEO

% load('Att.mat')
load('ISS_angle.mat')
y1=1;
sector_=[1:1:6 6:-1:1];
idx_=1;
att_final=[];
xx=1:1:512;

%Elev.angle sector
% Splitting 512 samples into 12 ectors
theta2=[angle(1) angle(43.*(1:11)) angle(end)];
for jj=1:size(sector_,2)
    if jj<6
        y2=(43*jj);
        
%         %Transition technique 1) 
%         %Use 50% of current slice, then 50% of this are samples from the
%         %current elev.angle interval and the other 50% belongs to the next
%         %interval.
%         %Find the line connecting these edge points and then subtract this
%         %value from the 50% samples of the current interval and add to this
%         %value the other 50% from the next interval.
%         %
%         %Use the last 21 samples of current interval
%         y3=y2-21+1;
%         z1=Att(y3,sector_(jj));
%         z2=Att(y2+1,sector_(jj+1));        
%         m=(z2-z1)/((y2+1)-y3);        
%         line_=m.*xx(1:20)+z1;        
%         Att(y3+1:y3+10,sector_(jj))=Att((y3+1):(y3+10),sector_(jj))-(z1-line_(1:10))';
%         Att(y3+11:y3+20,sector_(jj+1))=Att(y3+11:y3+20,sector_(jj+1))+(line_(11:20)-z2)';
        
        a=sector_(jj)
        %Transition technique 2)        
        z1=Att(y1+1:y2,sector_(jj));
        z2=Att(y1+1:y2,sector_(jj+1));
        
        angle_ref=angle(y1+1:y2);%ISS angle
        angle_d=abs(theta2(jj)-theta2(jj+1));% Angle difference for current sector
                
        w1=abs(theta2(jj)-angle_ref)./angle_d;
        w2=abs(angle_ref-theta2(jj+1))./angle_d;
        
        interp_=(z1.*w2')+(z2.*w1');%should be inverted so that at the first ISS angles (beggining) z1 has bigger weights (w2).
        slice_=[Att(y1,sector_(jj)) interp_']; 
        
       
%         slice_=Att(y1:y2,sector_(jj));%concatenates att within the elev.angle interval (each column has data from a different elev.angle att)
        att_final=[att_final slice_];
        
    elseif jj==6 || jj==7
        a=sector_(jj)
        y2=(43*jj);
        slice_=Att(y1:y2,sector_(jj));
        att_final=[att_final slice_'];
        
    elseif jj>=8 && jj<=11
        y2=(43*jj);
        a=sector_(jj)
        %Transition technique 2)        
        z1=Att(y1:y2,sector_(jj));
        z2=Att(y1:y2,sector_(jj-1));
        
        angle_ref=angle(y1:y2);%ISS angle
        angle_d=abs(theta2(jj)-theta2(jj+1));% Angle difference for current sector
                
        w1=abs(theta2(jj)-angle_ref)./angle_d;
        w2=abs(angle_ref-theta2(jj+1))./angle_d;
        
        interp_=(z1.*w1')+(z2.*w2');        
        slice_=[interp_' Att(y2+1,sector_(jj))];
        att_final=[att_final slice_];
        
        y2=y2+1;
        
    else %jj==12
        a=sector_(jj)
        %Transition technique 2) 
        
        z1=Att(y1:end-1,sector_(jj));
        z2=Att(y1:end-1,sector_(jj-1));
        
        angle_ref=angle(y1:end-1);%ISS angle
        angle_d=abs(theta2(jj)-theta2(jj+1));% Angle difference for current sector
                
        w1=abs(theta2(jj)-angle_ref)./angle_d;
        w2=abs(angle_ref-theta2(jj+1))./angle_d;
        
        interp_=(z1.*w1')+(z2.*w2');
        slice_=[interp_' Att(end,sector_(jj))];
        att_final=[att_final slice_];
        
    end
    y1=y2+1;
end

plotyy(xx,Att,xx,angle)
hold on
plot(xx,att_final)
grid on

%%
%Theoretical Coherence Time
freq=26e9;
v=7500; %27,000mph
fd_=(7500/((3e8)/(26e9)).*cosd(angle));
fd_max=max(fd_);

Td_=0.423/fd_max;

rate_threshold=1/Td_  %data_rate>1/Td ->Slow fading

%%
%0-600 sec
% time_series=

%Test
% Ar=att(:,6)';    
% [Arain, Time]=itur_p1853_0_annex1_s2(ccdf,Ar,P0,sampling_rate,duration);
% 
% 
% ccdf = [kron(10.^(-2:0),[1 2 3 5]),10];
% Ar = [16.0243, 12.1351, 10.1633, 8.0040, 5.6312, 3.8381, 3.0224, 2.2024, ...
%       1.3946, 0.8555, 0.6334, 0.4271, 0.2434];
% P0 = 4.9095;
% sampling_rate = 1;
% duration = 86400;
% [Arain, Time]=itur_p1853_0_annex1_s2(ccdf,Ar,P0,sampling_rate,duration);
