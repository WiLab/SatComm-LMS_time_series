% Land Mobile Satellite time series synthesizer based on ITU-R P.681-12 [1] for  
%             Ka-band and GEO satellite with fixed ground receiver
% 
% Authors: Paulo Ferreira, Alex M. Wyglinski
% Wireless Innovation Laboratory
% Worcester Polytechnic Institute 
% Worcester, MA
%
% December 2015
%
% Please cite this code as: Ferreira P. and Wyglinski A.M. Satellite-Communications: LMS time series gen GEO v1.1. GitHub, 2015. doi: 10.5281/zenodo.44695.
% or
% Ferreria, P.V.R., Paffenroth, R. and Wyglinski, A.M. "IMM performance analysis for variable measurement delay for satellite communications at Ka-band for LMS channel under rain fading." [Under review].

% For Rain time series simulation for LEO/GEO please refer to:
% P. Ferreira and A. Wyglinski, Satellite-Communications: 3D rain field and LEO time series v1.1. GitHub, 2016, doi: 10.5281/zenodo.44881. [Online]. Available: https://zenodo.org/badge/latestdoi/19596/WiLab/SatComm-3D rain field

clear all
clc
%Time series synthethizer ITU P.681-12 [1]

%LMS Ka-band GEO

%Parameters
u_GOOD=1.0125;
sigma_GOOD=1.6944;
u_BAD=-0.8026;
sigma_BAD=1.288;
dur_min_GOOD=1.5;
dur_min_BAD=1.1;

uMA_GOOD=-0.02;
sigmaMA_GOOD=0;
uMA_BAD=-5.4;
sigmaMA_BAD=7.3;
pB_min=0.1;
pB_max=0.6;

h1_GOOD=0;
h2_GOOD=-38.17;
h1_BAD=0.69;
h2_BAD=-15.97;
g1_GOOD=0;
g2_GOOD=0.39;
g1_BAD=-0.21;
g2_BAD=0;

f1=0.036;
f2=0.8;

L_corr_G=0.5;
L_corr_B=0.5;

Loo_range_G_min=(uMA_GOOD-(1.645*sigmaMA_GOOD));
Loo_range_G_max=(uMA_GOOD+(1.645*sigmaMA_GOOD));
Loo_range_B_min=uMA_BAD+(sigmaMA_BAD*(sqrt(2))*(1/(erf((2*pB_min)-1))));
Loo_range_B_max=uMA_BAD+(sigmaMA_BAD*(sqrt(2))*(1/(erf((2*pB_max)-1))));

%%
%State series generator
freq=26e9;%Carrier frequency
vm=round(13.88);%[m/s] vehicle speed 50 km/h
% vsat=0;%sat speed
LEO_dur=600; %[secs] pass duration of a LEO
dis_=vm*LEO_dur; %distance traveled [m] 
fs=10000; %DLR sampling rate based on 26 GHz and vm speed [2]
theta=30;%Sat elev.angle
phi=0;%angle between mobile direction and sat. azimuth
fd=freq*(vm/(3e8))*cosd(phi)*cosd(theta);
fc=uint64(vm*freq/3e8);




time_series=zeros(1,(fs*LEO_dur));%oversampled
state2=time_series;
time_series_dB=time_series;
time_series_dB_filter=time_series;
time_series_cplx=time_series;
time_series_multi=time_series;
time_series_multi_filter=time_series;
time_series_LMS_cplx=time_series;
time_series_LMS=time_series;
time_series_LMS_dB=time_series;

pointer_state2=1;

dur_G=0;
dur_B=0;
res_=1;%distance resolution
series_=dis_*res_;%series length in [m] with resolution of 0.1
state_series=zeros(1,series_);%logical vector
pointer_=1;
pointer_2=0;
flagGB=1;
flag_tran=0;%flag for transition (used once)
rng default
while pointer_<=length(state_series)    
    if flagGB==1%GOOD
        dur_G=0;
        while dur_G<dur_min_GOOD 
            rng shuffle
            dur_G = round(lognrnd(u_GOOD,sigma_GOOD,1,1),1);%duration of Good state in [m]
        end
        flagGB=0;
        state_series(1,pointer_:pointer_+(dur_G*res_))=1; %Good state series in [m]       
        pointer_=pointer_+(dur_G*res_)+1;%undersampled
        pointer_tran=pointer_+(dur_G*res_);%End of Good/start of Bad undersampled
        state2(1,pointer_state2:pointer_state2+(round((dur_G/vm),2)*fs))=1;%Series conversion from [m] to [sec]
        idx_G_start=pointer_state2;%start of GOOD oversampled
        pointer_state2=pointer_state2+(round((dur_G/vm),2)*fs)+1;
        idx_G_stop=pointer_state2-1;%stop of GOOD oversampled
        idx_B_start=pointer_state2;%start of BAD oversampled
        pointer_2=pointer_2+1;
        
        if flag_tran==1
            MA_1=MA_B;%saves the previous MA_B for BAD-GOOD transition length 
        end
        

        
    else%BAD 
        dur_B=0;
        while dur_B<dur_min_BAD 
            rng shuffle
            dur_B = round(lognrnd(u_BAD,sigma_BAD,1,1),1);
        end
        flagGB=1;
        state_series(1,pointer_:pointer_+(dur_B*res_))=0;        
        pointer_=pointer_+(dur_B*res_)+1;
        state2(1,pointer_state2:pointer_state2+(round((dur_B/vm),2)*fs))=0;
        pointer_state2=pointer_state2+(round((dur_B/vm),2)*fs)+1;
        idx_B_stop=pointer_state2-1;%stop of BAD oversampled
        pointer_2=pointer_2+1;
        
        
    end
    
    if pointer_2==2 %Already generated two consecutive different states and Loo triplet parameters
        pointer_2=0;
        if Loo_range_G_min~=Loo_range_G_max% when sigmaMA_GOOD=0, Loo_range_G_min==Loo_range_G_max
            while MA_G<Loo_range_G_min || MA_G>Loo_range_G_max
                rng shuffle
                MA_G=normrnd(uMA_GOOD,sigmaMA_GOOD,1,1);
            end
        else
            MA_G=uMA_GOOD;%when sigmaMA_GOOD=0
        end
        sum_MA_G=(g1_GOOD*MA_G)+g2_GOOD;%g1_GOOD=0, sum_MA_G==g2_GOOD
        sum_MP_G=(h1_GOOD*MA_G)+h2_GOOD;        
        
        MA_B=Loo_range_B_min-1;%This lures the algorithm for it to start
        rng default
        while MA_B<Loo_range_B_min || MA_B>0
            rng shuffle
            MA_B=normrnd(uMA_BAD,sigmaMA_BAD,1,1);
        end
        sum_MA_B=(g1_BAD*MA_B)+g2_BAD;
        sum_MP_B=(h1_BAD*MA_B)+h2_BAD;
        
        %Transition
        idx_trans_start_B=idx_B_start;
        L_trans_GB=round((f1*abs(MA_G-MA_B))+f2,1);
        idx_trans_stop_B=idx_B_start+(round((L_trans_GB/vm),2)*fs);
        weights_B=zeros(1,uint64(round((L_trans_GB/vm),2)*fs));
        for ii=1:length(weights_B)
            weights_B(1,ii)=ii*(1/(length(weights_B)));%Vector containing weights for the transition
        end
        weights_BB=[0, weights_B];
        if flag_tran==1
            idx_trans_start_G=idx_G_start;
            L_trans_BG=round((f1*abs(MA_1-MA_G))+f2,1);
            idx_trans_stop_G=idx_G_start+(round((L_trans_BG/vm),2)*fs);
            weights_G=zeros(1,uint64(round((L_trans_BG/vm),2)*fs));
            for ii=1:length(weights_G)
                weights_G(1,ii)=ii*(1/(length(weights_G)));%Vector containing weights for the transition
            end
            weights_GG=[0, weights_G];
        end
            
        
        %Gaussian series direct signal (slow variation)        
        %GOOD direct signal
        rng shuffle
        dir_sig_G=[];
        dir_sig_LP=[];
        n_G=idx_G_stop-idx_G_start+1;
        dir_sig_G=normrnd(MA_G,(sum_MA_G),1,uint64(n_G));
     
        time_series_dB(idx_G_start:idx_G_stop)=dir_sig_G;%[dB]
        %Transition BAD-GOOD

        if flag_tran==1 
            L_trans_BG_=idx_trans_stop_G-idx_trans_start_G;
                    w1=weights_GG;
                w2=fliplr(w1);
%                         time_series_dB(idx_trans_start_G-round(L_trans_BG/2):idx_trans_start_G)=(time_series_dB(idx_trans_start_G-round(L_trans_BG/2):idx_trans_start_G).*w2(1:round(L_trans_BG/2)))+(w1(1:round(L_trans_BG/2))*MA_1);
%              time_series_dB(idx_trans_start_G+1:idx_trans_start_G+(L_trans_BG-round(L_trans_BG/2)))=(time_series_dB(idx_trans_start_G+1:idx_trans_start_G+(L_trans_BG-round(L_trans_BG/2))).*w1(1:(L_trans_BG-round(L_trans_BG/2))-1))+(fliplr(weights_GG)*MA_1);
             
            time_series_dB(idx_trans_start_G-L_trans_BG_:idx_trans_start_G)=(time_series_dB(idx_trans_start_G-L_trans_BG_:idx_trans_start_G).*w2)+(MA_1.*w1);
             
        end
        

        %GOOD multipath
        rng shuffle
        multi_G_real=normrnd(0,1,1,uint64(n_G));
        rng shuffle
        multi_G_im=1i.*(normrnd(0,1,1,uint64(n_G)));
        multi_G=multi_G_real+multi_G_im;
%         %Multipath Butterworth filter                
%         [zb,pb,kb] = butter(3,fc/(n_G/2));
%         [bb,ab] = zp2tf(zb,pb,kb);
%         multi_G_b=filter(bb,ab,multi_G);%Doppler spread through Butterworth filter
        multi_G_total=multi_G.*(sqrt(0.5*(10^(sum_MP_G/10))));
        time_series_multi(idx_G_start:idx_G_stop)=multi_G_total;
       
        
        %BAD direct signal
        rng shuffle
        dir_sig_B=[];
        dir_sig_LP=[];
        n_B=idx_B_stop-idx_B_start+1;
        dir_sig_B=normrnd(MA_B,(sum_MA_B),1,uint64(n_B));
     
        time_series_dB(idx_B_start:idx_B_stop)=dir_sig_B;%[dB] 
        %Transition GOOD-BAD
        time_series_dB(idx_trans_start_B:idx_trans_stop_B)=(time_series_dB(idx_trans_start_B:idx_trans_stop_B).*weights_BB)+(fliplr(weights_BB)*MA_G);
        
        %BAD multipath
        rng shuffle
        multi_B_real=normrnd(0,1,1,uint64(n_B));
        rng shuffle
        multi_B_im=1i.*(normrnd(0,1,1,uint64(n_B)));
        multi_B=multi_B_real+multi_B_im;
%         %Multipath Butterworth filter 
%         multi_B_b=filter(bb,ab,multi_B);%Doppler spread through Butterworth filter
        multi_B_total=multi_B.*(sqrt(0.5*(10^(sum_MP_B/10))));
        time_series_multi(idx_B_start:idx_B_stop)=multi_B_total;       
        
        flag_tran=1;        
    end
end

%%Butterworth filter
[zb,pb,kb] = butter(10,(double(fc)/(double(fs)/2)));
[bb,ab] = zp2tf(zb,pb,kb);
% freqz(bb,ab)
for ii=1:length(time_series_multi)/double(fs)
    time_series_multi_filter(((ii-1)*double(fs))+1:((ii-1)*double(fs))+double(fs))=filter(bb,ab,time_series_multi(((ii-1)*double(fs))+1:((ii-1)*double(fs))+double(fs)));
end

%Low-pass filter direct signal
ro_s_B=exp((-vm*(1/fs))/(L_corr_B));
a=[1 -ro_s_B];
b=[sqrt(1-(ro_s_B^2))];
for ii=1:length(time_series_dB)/double(fs)
    time_series_dB_filter(((ii-1)*double(fs))+1:((ii-1)*double(fs))+double(fs))=filter(b,a,time_series_dB(((ii-1)*double(fs))+1:((ii-1)*double(fs))+double(fs)));
end

time_series=10.^(time_series_dB_filter./20);%[Linear]
time_series_cplx=time_series.*(exp(1i*2*pi*fd*(1/fs)));%[Linear complex]

time_series_LMS_cplx=time_series_cplx+time_series_multi_filter;
time_series_LMS=abs(time_series_LMS_cplx);
time_series_LMS_dB=20.*log10(time_series_LMS);

    if Loo_range_G_min~=Loo_range_G_max % when sigmaMA_GOOD=0, Loo_range_G_min==Loo_range_G_max
        rng default
        while MA_G_2<Loo_range_G_min || MA_G_2>Loo_range_G_max
            rng shuffle
            MA_G_2=normrnd(uMA_GOOD,sigmaMA_GOOD,1,1);
        end
    else
        MA_G_2=uMA_GOOD;%when sigmaMA_GOOD=0
    end
    while dur_G<dur_min_GOOD
        rng shuffle
        dur_G = round(lognrnd(u_GOOD,sigma_GOOD,1,1));%duration [m] with resolution at 0.1
    end

%% 
% REFERENCES
% [1] “Rec. ITU-R P.618-12 - Propagation data and prediction methods required for the design of earth-space telecommunication systems.” International Telecommunication Union, Tech. Rep., 2015. 
% [2] “Technical note on the implementation of the land mobile satellite channel model, software usage.” German Aerospace Center DLR, Tech. Rep., 2007. 
