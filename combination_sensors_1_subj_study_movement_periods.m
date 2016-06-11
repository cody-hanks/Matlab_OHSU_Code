
clear;
close all;
clc;

addpath('/home/cody/Documents/for_cody_11_30/Helpersfunc/');


wtbar = waitbar(0,'Running Calulation...');
DataSheet=zeros(1000,100);
positions =10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cords of the sensors for COP X and COP Y 
co_ords = [4,165;
    4,86;
    4,27;
    149,168;
    149,83;
    149,29;
    5,18;
    148,18;
    76,89;
    5,161;
    148,161];

positionstring = ['L calb';
'R calb';
'L back';
'L left';
'L rigt';
'L stom';
'R back';
'R left';
'R rigt';
'R stom'];

    labels = ['BCK_C_L';'BCK_C_R';'BCK_L  ';'LFT_L  ';'RHT_L  ';'STM_L  '; ...
        'BCK_R  ';'LFT_R  ';'RHT_R  ';'STM_R  '];
    
positionval = [1;1;1;2;3;4;1;2;3;4];
positionval_all = [1;2;3;4;5;6;7;8;9;10];
positionval_FS = [1;1;1;2;3;1;1;2;3;1];
positionval_LR = [1;2;1;3;4;5;2;6;7;8];
positionval_FBLR = [1;2;1;3;4;1;2;5;6;2];
positionval_FB_LR = [1;1;1;2;2;1;1;2;2;1];
positionval_F_B_L_R = [1;1;1;-1;1;-1;1;-1;1;-1];

%%%%%%% clear old data 
clear estimate_corr lc_aligned pressure_aligned clear edf1_plot_cell clear edf2_plot_cell ...
t_edf1cell t_edf2cell cepstral_coeff lpc_coeff clear data_cell rr_cell ...
N_distr x_val_hist ar_coeff wdth prm rr_psg_cell ALL_SUM;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%flags for code
prdplt_flg =0;
plt_filt = 0;
plt_xy =1;
pltfft =0;


winoverlap = 0;
pltplt =0;
% subindarr = [4:9];
filename ='DataSheet2-2.mat';
%subindarr = [4:9, 13:14, 16, 18:20]
%subindarr = [subindarr 2:3,10:12,15,21:22,24]
subindarr = 4:9;
%subindarr = 10:17;
dec_factor = 10;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% loop for subjects 
for kkkk = 1:length(subindarr) 
    
    clear data
    waitbar(kkkk/length(subindarr),wtbar,strcat('Running',num2str(subindarr(kkkk))));
    %close all;
%         subind = two_subjects_time_indices(kkkk,1);
    % get the subject number 
    subind = subindarr(kkkk);
    
    ipos =1;
    % get data sheet list of rows for subject 
    strings  = char(positionstring);
    dataset =[ ((kkkk-1)*positions+1):((kkkk)*positions) ];
    for j = dataset(1:end)
        DataSheet(j,1) = subind;% subject # added to sheed 
        DataSheet(j,2) = ipos; % add the string of position to data sheet 
        DataSheet(j,3) = positionval(ipos); % add the position value to the sheet 
        DataSheet(j,25) = positionval_all(ipos); % add the 
        DataSheet(j,26) = positionval_FS(ipos);
        DataSheet(j,27) = positionval_LR(ipos);
        DataSheet(j,28) = positionval_FBLR(ipos);
        DataSheet(j,29) = positionval_FB_LR(ipos);
        DataSheet(j,30) = positionval_F_B_L_R(ipos);
        ipos = ipos +1;
        if (ipos >10) 
            ipos = 1;
        end
         
    end
    % read data from subject 
    read_1subject_on_bed;

    %%%%%%%%%%%%%%%%%%
    % read XLS file 
    nums = xlsread(['/home/cody/Documents/for_cody_11_30/sp/sp/SP_ID_' num2str(subind) '/SP_ID_' num2str(subind) '_Info.xlsx']);    
    
    % sensors selected 
    sensors_selected = [7,8,10,11,9];
    cop_x = zeros(ceil(size(data,1)/5/5/5/2),1);
    cop_y = zeros(ceil(size(data,1)/5/5/5/2),1);
    X_tot = zeros(ceil(size(data,1)/5/5/5/2),1);
    Y_tot = zeros(ceil(size(data,1)/5/5/5/2),1);
    
    for j = sensors_selected(1:end)
        decsig = decimate(data(:,j),5);
        decsig = decimate(decsig,5);
        decsig = decimate(decsig,5);
        decsig = decimate(decsig,2);
       
        cop_x =cop_x + decsig*co_ords(j,1);
        cop_y =cop_y + decsig*co_ords(j,2);
        X_tot = X_tot+ decsig;
        Y_tot = Y_tot+ decsig;
        
    end
    
    
    
    
    cop_x = cop_x ./ X_tot;
    cop_y = cop_y ./ Y_tot;
    
    clear X_tot Y_tot decsig innersig outersig;
    
    
    

    
    fs = 25;
    %%%%% filter the data to window of pass <10Hz and  pass > .5 Hz 
     [N, Ws] = cheb2ord(0.1/(fs/2),0.06/(fs/2), 3, 40, 's');
     [B,A] = cheby2(N,40,Ws,'high');
     [B2,A2] = cheby2(7,40,(10/(fs/2)),'low');
     
    devid = 10;
    %decimate to a lower level and combine for filtering 
    %decimate to 20hz IE remove high freq details 
    cop_x25 = zeros(ceil(size(data,1)/devid),1);
    cop_y25 = zeros(ceil(size(data,1)/devid),1);
    X_tot = zeros(ceil(size(data,1)/devid),1);
    Y_tot = zeros(ceil(size(data,1)/devid),1);
    outersig = zeros(ceil(size(data,1)/devid),length(sensors_selected));
    column = 1;
    for j = sensors_selected(1:end)
        decsig = decimate(data(:,j),5);
        decsig = decimate(decsig,2);
        %decsig = filtfilt(B,A,decsig);
        decsig = filtfilt(B2,A2,decsig);
        %filter the signals and then decimate 
        %decsig = filtfilt(B,A,data(:,j));
        %decsig = filtfilt(B2,A2,data(:,j));
        %decsig = decimate(decsig,5);
        %decsig = decimate(decsig,2);
        %decsig = data(:,j);
        
        cop_x25 =cop_x25 + decsig*co_ords(j,1);
        cop_y25 =cop_y25 + decsig*co_ords(j,2);
        X_tot = X_tot+ decsig;
        Y_tot = Y_tot+ decsig;
        if j ~=9
            outersig(:,column)  = decsig;
        elseif j==9 
            innersig = decsig;
        end
        column = column+1;
    end
    
    %ratio = innersig'./outersig';
    
    
    cop_x25 = cop_x25 ./ X_tot;
    cop_y25 = cop_y25 ./ Y_tot;
    
    clear X_tot Y_tot decsig;
    
%    [N,Ws] = cheb2ord(0.75/250,1.3/250,3,40);
%    [B,A] = cheby2(N,40,Ws);
   
    % get start and stop time for each segment get a combo of all segs and
    % each seg ind.
    stoptime = 200;%2*60*2-1-(2*15);
    back_cal_L = floor(nums(2,1)/250); %1
    back_cal_R =  floor(nums(3,1)/250); %2
    bck_L =  floor(nums(4,1)/250)+6;     %3
    lft_L = floor(nums(5,1)/250)+6;     %4
    rht_L = floor(nums(6,1)/250)+6;     %5
    stm_L = floor(nums(7,1)/250)+6;     %6
    bck_R = floor(nums(8,1)/250)+6;     %7
    lft_R = floor(nums(9,1)/250)+6;     %8
    rht_R = floor(nums(10,1)/250)+6;     %9
    stm_R = floor(nums(11,1)/250)+6;     %10
    
    % get start and stop time for each segment get a combo of all segs and
    % each seg ind.
    stoptime = 100*(250/devid);%2*60*2-1-(2*15);
    seg_back_cal_L_25 = floor(nums(2,1)/devid):floor(nums(2,1)/devid)+stoptime; %1
    seg_back_cal_R_25 = floor(nums(3,1)/devid):floor(nums(3,1)/devid)+stoptime; %2
    seg_bck_L_25 = floor(nums(4,1)/devid)+75:floor(nums(4,1)/devid)+75+stoptime;     %3
    seg_lft_L_25 = floor(nums(5,1)/devid)+75:floor(nums(5,1)/devid)+75+stoptime;     %4
    seg_rht_L_25 = floor(nums(6,1)/devid)+75:floor(nums(6,1)/devid)+75+stoptime;     %5
    seg_stm_L_25 = floor(nums(7,1)/devid)+75:floor(nums(7,1)/devid)+75+stoptime;     %6
    seg_bck_R_25 = floor(nums(8,1)/devid)+75:floor(nums(8,1)/devid)+75+stoptime;     %7
    seg_lft_R_25 = floor(nums(9,1)/devid)+75:floor(nums(9,1)/devid)+75+stoptime;     %8
    seg_rht_R_25 = floor(nums(10,1)/devid)+75:floor(nums(10,1)/devid)+75+stoptime;     %9
    seg_stm_R_25 = floor(nums(11,1)/devid)+75:floor(nums(11,1)/devid)+75+stoptime;     %10
        
    %reset the stop time 
    stoptime = 200;

    
    seg_bck_cal_L = [back_cal_L:back_cal_L+stoptime];
    seg_bck_cal_R = [back_cal_R:back_cal_R+stoptime];
    seg_bck_L = [bck_L:bck_L+stoptime];
    seg_bck_R = [bck_R:bck_R+stoptime];
    seg_lft_L = [lft_L:lft_L+stoptime];
    seg_lft_R = [lft_R:lft_R+stoptime];
    seg_rht_L = [rht_L:rht_L+stoptime];
    seg_rht_R = [rht_R:rht_R+stoptime];
    seg_stm_L = [stm_L:stm_L+stoptime];
    seg_stm_R = [stm_R:stm_R+stoptime];
      
    all_periods = [ back_cal_L:back_cal_L+stoptime, ...
        back_cal_R:back_cal_R+stoptime, ...
        bck_L:bck_L+stoptime,lft_L:lft_L+stoptime,rht_L:rht_L+stoptime,stm_L:stm_L+stoptime, ...
        bck_R:bck_R+stoptime,lft_R:lft_R+stoptime,rht_R:rht_R+stoptime,stm_R:stm_R+stoptime];
    all_periods_it = [ back_cal_L:back_cal_L+stoptime; ...
        back_cal_R:back_cal_R+stoptime; ...
        bck_L:bck_L+stoptime;lft_L:lft_L+stoptime;rht_L:rht_L+stoptime;stm_L:stm_L+stoptime; ...
        bck_R:bck_R+stoptime;lft_R:lft_R+stoptime;rht_R:rht_R+stoptime;stm_R:stm_R+stoptime];
   
    all_periods_it_25 = [ seg_back_cal_L_25; ...
        seg_back_cal_R_25; ...
        seg_bck_L_25;seg_lft_L_25;seg_rht_L_25;seg_stm_L_25; ...
        seg_bck_R_25;seg_lft_R_25;seg_rht_R_25;seg_stm_R_25];
    
    if prdplt_flg ==1
    first = figure;
    figure(first);
    plot(cop_x,'b');
    hold on;
    plot(all_periods,cop_x(all_periods),'.r');
    title(strcat('subject ',num2str(subind),' cop x periods of intereset'));
    second =figure;
    figure(second);
    plot(cop_y,'b');
    hold on;
    plot(all_periods,cop_y(all_periods),'.r');
    title(strcat('subject ',num2str(subind),' cop y periods of intereset'));
    
    figure;
    plot(cop_x25,'b');
    hold on;
    plot(cop_y25,'g');
        for j =1:10
            hold on;
            plot(all_periods_it_25(j,:),cop_x25(all_periods_it_25(j,:)),'.r');
            hold on;
            plot(all_periods_it_25(j,:),cop_y25(all_periods_it_25(j,:)),'.r');
        end
    title(strcat('subject ',num2str(subind),' 25 Hz period of interest'));
    end
    % end of plot for all periods of intrest 
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% add data to data sheet pull the axis and pca data from the 2hz signal 
    dec_cop_x = zeros(201,10);
    dec_cop_y = zeros(201,10);
    for j =1:10
        % get the median filtered cop x and y %8 = 4 seconds per mahesh 
        temp_dec_cop_x = medfilt1(cop_x(all_periods_it(j,:)),8);
        temp_dec_cop_y = medfilt1(cop_y(all_periods_it(j,:)),8);
        % subtract the median filter from the signal ... 
        dec_cop_x(2:end,j) = cop_x(all_periods_it(j,2:end))-temp_dec_cop_x(2:end);
        dec_cop_y(2:end,j) = cop_y(all_periods_it(j,2:end))-temp_dec_cop_y(2:end);
        % get the pca 
        coeff = pca([dec_cop_x(:,j),dec_cop_y(:,j)]);
        temp = zeros(2,4);% create container 
        temp(:,2:2:end) = coeff;% get the coeffs
        xtemp(j,:) = temp(1:2:end);
        ytemp(j,:) = temp(2:2:end);
        
        DataSheet(dataset(j),4) = atan2(xtemp(j,4)-xtemp(j,3),xtemp(j,2)-xtemp(j,1))*180/pi;
        DataSheet(dataset(j),5) = atan2(ytemp(j,4)-ytemp(j,3),ytemp(j,2)-ytemp(j,1))*180/pi;
        %DataSheet(dataset(j),24) = 
        outersigmix = mean([...
            rms(outersig(all_periods_it(j,:),1))...
            ,rms(outersig(all_periods_it(j,:),2))...
            ,rms(outersig(all_periods_it(j,:),3))...
            ,rms(outersig(all_periods_it(j,:),4))]); 
        
        ratio = rms(innersig(all_periods_it(j,:)))'/outersigmix';
        DataSheet(dataset(j),24) = ratio;
        xtemp(j,:) = xtemp(j,:)+ mean(dec_cop_x(:,j));
        ytemp(j,:) = ytemp(j,:)+ mean(dec_cop_y(:,j));
        
        if prdplt_flg ==1 
            figure(first);
            hold on;
            plot(all_periods_it(j,:),dec_cop_x(:,j),'.g');
            figure(second);
            hold on;
            plot(all_periods_it(j,:),dec_cop_y(:,j),'.g');
            
        end
        
    end
    
    

    strings = char(labels);
    if plt_xy ==1
        figure;
         for j = 1:10
            xmax = mean(dec_cop_x(:,j))+.06;
            xmin = mean(dec_cop_x(:,j))-.06;
            ymax = mean(dec_cop_y(:,j))+.06;
            ymin = mean(dec_cop_y(:,j))-.06;  
           subplot(4,3,j);
           plot(dec_cop_x(:,j),dec_cop_y(:,j));
           hold on; 
           plot(xtemp(j,:),ytemp(j,:),'r','LineWidth',2);
           axis([xmin xmax ymin ymax]);
           title(strings(j,:));
         end
     end
    
    
    
    
    clear all_periods dec_sig dec_cop_x dec_cop_y xtemp ytemp temp xmax ymax xmin ymin
    
    
    
    %%%%% filter the data to window of pass <10Hz and  pass > .5 Hz 
%      [N, Ws] = cheb2ord(0.1/(fs/2),0.06/(fs/2), 3, 40, 's');
%      [B,A] = cheby2(N,40,Ws,'high');
%      [B2,A2] = cheby2(7,40,(10/(fs/2)),'low');
%     
%     if(plt_filt ==1)
%         figure;
%         freqz(B,A,1000,500);
%        	figure;
%         freqz(B2,A2,1000,500);
%     end
    
    %%%% filter cop_y 
    % filters generated by mahesh
 %   cop_y25 = filtfilt(B,A,cop_y25);
    
%     cop_y = filtfilt(B,A,cop_y);
%     cop_y = filtfilt(B2,A2,cop_y);
%     cop_y = cop_y - mean(cop_y);

%     cop_y = cop_y/(max(abs(cop_y)));

%     cop_x = filtfilt(B,A,cop_x);
%     cop_x = filtfilt(B2,A2,cop_x);
%     cop_x = cop_x - mean(cop_x);
%     
    
    
    %%%%% 
    % get a med filter for bck copy 
    
      Fs =25;
      T = 1/Fs;
     for p = 1:10 
        tempx = medfilt1(cop_x25(all_periods_it_25(p,:)),((250/devid)*4));
        tempy = medfilt1(cop_y25(all_periods_it_25(p,:)),((250/devid)*4));
        
        sigfiltx = cop_x25(all_periods_it_25(p,2:end))-tempx(2:end);
        sigfilty = cop_y25(all_periods_it_25(p,2:end))-tempy(2:end);
        

        if pltplt == 1 
            figure;
            plot(sigfiltx(2:end));
            figure;
            plot(sigfilty(2:end));
        end
        
        DataSheet(dataset(p),6) = rms(sigfiltx);
        DataSheet(dataset(p),7) = rms(sigfilty);
        
        [envx fsx] = hilbert2(sigfiltx);
        [envy fsy] = hilbert2(sigfilty); 
        
        DataSheet(dataset(p),8) = rms(envx);
        DataSheet(dataset(p),9) = rms(envy);


         L = size(sigfiltx,1);
         f = Fs*(0:(L/2))/L;
         p1 = abs(fft(sigfiltx)/L);
         p2 = p1(1:L/2+1);
         fft_sig = 2*p2;
         if pltfft ==1
             figure;
             plot(f,fft_sig);
         end
         
         p1y = abs(fft(sigfilty)/L);
         p2y = p1y(1:L/2+1);
         fft_sigy = 2*p2y;
         if pltfft ==1
             figure;
             plot(f,fft_sigy);
         end
         
         [pks locs] = findpeaks(fft_sig);
         [pksy locsy] = findpeaks(fft_sigy);
         
         DataSheet(dataset(p),10) = max(pks);
         DataSheet(dataset(p),11) = max(pksy);
         DataSheet(dataset(p),12) = f(locs(find(pks==max(pks))));
         DataSheet(dataset(p),13) = f(locsy(find(pksy==max(pksy))));
         
         DataSheet(dataset(p),18) = sum(fft_sig);
         DataSheet(dataset(p),19) = sum(fft_sigy);
         DataSheet(dataset(p),20) = std(fft_sig);
         DataSheet(dataset(p),21) = std(fft_sigy);
         DataSheet(dataset(p),22) = entropy(fft_sig);
         DataSheet(dataset(p),23) = entropy(fft_sigy);
        
%         DataSheet(dataset(p),9) = f(locs(find(pks==max(pks))));
%         [pxx,f] =periodogram(cop_y(all_periods_it(p,:)),[],[],500);
%         [pks locs] = findpeaks(pxx); 
%         DataSheet(dataset(p),10) = max(pks);
%         DataSheet(dataset(p),11) = f(locs(find(pks==max(pks))));
%         DataSheet(dataset(p),12) = sum(pxx);    
%         
         sigmaxx = max(sigfiltx);
         sigmaxy = max(sigfilty);
         sigminx = min(sigfiltx);
         sigminy = min(sigfilty);
%         DataSheet(dataset(p),13) = sigmax;
%         DataSheet(dataset(p),14) = sigmin;
         stepsizex = (sigmaxx-sigminx)/50;
         stepsizey = (sigmaxy-sigminy)/50;
         edgesx = sigminx:stepsizex:sigmaxx;
         edgesy = sigminy:stepsizey:sigmaxy;
         
         histogrmx = histc(sigfiltx,edgesx);
         histogrmy = histc(sigfilty,edgesy);
         DataSheet(dataset(p),14) = kurtosis(histogrmx);
         DataSheet(dataset(p),15) = skewness(histogrmy);
         DataSheet(dataset(p),16) = kurtosis(histogrmx);
         DataSheet(dataset(p),17) = skewness(histogrmy);
%         %bar(edges,histogrm,'histc');
%         DataSheet(dataset(p),15) = max(histogrm);
%         DataSheet(dataset(p),16) = min(histogrm);
%         DataSheet(dataset(p),17) = kurtosis(histogrm);
%         DataSheet(dataset(p),18) = skewness(histogrm);
        
        

        clear sigfiltx sigfilty
    end
    clear env fs pks locs p1 p2 all_periods_it L f sigmax sigmin edges histogrm
    
    
    
    

    
    
    
    
    
    
%     last = 1;
%     L = size(seg_bck_L,2)/10;
%     for j = 1:10
%         subseg = seg_bck_L(last:last+L-1);
%         cop_y(subseg) = medfilt1(cop_y(subseg),500);
%         last = L*j;
%     end
%     
    %plot(seg_bck_L,cop_y(seg_bck_L));
    
    %plot(bck_seg_L);
    
    
    
    %[m_ind,m_seg] = detect_movement(ALL_SUM,7,2.5,500);
%     for mov_indx = 1:size(m_seg,1)
% %      cop_y(m_seg(mov_indx,1):m_seg(mov_indx,2)) = mean(cop_y([(m_seg(mov_indx,1) - 1000):(m_seg(mov_indx,1) - 500) (m_seg(mov_indx,2) + 500):(m_seg(mov_indx,2) + 1000)]));
%         cop_y((m_seg(mov_indx,1):(m_seg(mov_indx,2)))) = cop_y(((m_seg(mov_indx,1))):(m_seg(mov_indx,2))).*(1-hamming(length((m_seg(mov_indx,1)):(m_seg(mov_indx,2))))');
%         
%     end;
%     cop_y = filter(B,A,cop_y);
%     cop_y = cop_y_nm(t_lc < t_end & t_lc > t_start);
%     cop_x = cop_x_nm(t_lc < t_end & t_lc > t_start);
    
    
%     fig_ind = 1;
%     close all;
%     for t1t = 1:size(window_loc,2)
%         
% %         t_start = t_lc(1) + ((window_loc(t1t)*60)/86400);
% %         t_end = t_start + ((winwidth*60)/86400);
% 
%         t_start = t_lc(window_loc(subindarr(kkkk),t1t))+10/86400;
%         
%         if(t1t == size(window_loc,2))
%             t_end = t_lc(window_loc(subindarr(kkkk),t1t))+100/86400;
%         else 
%             t_end = t_lc(window_loc(subindarr(kkkk),t1t+1))-20/86400;
%         end
%         
%         t_lc_window = t_lc(t_lc < t_end & t_lc > t_start);
%         t_edf1_window = t_edf1{edf_sensor_ind}(t_edf1{edf_sensor_ind} < t_end & t_edf1{edf_sensor_ind} > t_start);
%         %     t_edf2_window = t_edf2{edf_sensor_ind}(t_edf2{edf_sensor_ind} < t_end & t_edf2{edf_sensor_ind} > t_start);
%         
%         ind_lc_start = find(t_lc == t_lc_window(1));
%         ind_lc_end = find(t_lc == t_lc_window(end));
%         
%         ind_edf1_start = find(t_edf1{edf_sensor_ind} == t_edf1_window(1));
%         ind_edf1_end = find(t_edf1{edf_sensor_ind} == t_edf1_window(end));
%         
%         %     ind_edf2_start = find(t_edf2{edf_sensor_ind} == t_edf2_window(1));
%         %     ind_edf2_end = find(t_edf2{edf_sensor_ind} == t_edf2_window(end));
%         edf1_windowed = sig1{edf_sensor_ind}(ind_edf1_start:ind_edf1_end);
%         %     edf2_windowed = sig2{edf_sensor_ind}(ind_edf2_start:ind_edf2_end);
%         edf1_windowed_interp{kkkk} = interp1(t_edf1{edf_sensor_ind}(ind_edf1_start:ind_edf1_end),edf1_windowed,t_lc(ind_lc_start:ind_lc_end));
%         %     edf2_windowed_interp{kkkk} = interp1(t_edf2{4}(ind_edf2_start:ind_edf2_end),edf2_windowed,t_lc(ind_lc_start:ind_lc_end));
%         
%         for tt = 1:length(edf1_windowed_interp{kkkk})
%             if(isnan(edf1_windowed_interp{kkkk}(tt)) == 1)
%                 edf1_windowed_interp{kkkk}(tt) = 0;
%             end
%        end
        
        
%         pressure_aligned{kkkk,t1t} = edf1_windowed_interp{kkkk};
%         
%         t_edf1cell{kkkk,t1t} = t_edf1{4}(ind_edf1_start:ind_edf1_end);
%         
%         %     t_edf2cell{kkkk,t1t} = t_edf2{4}(ind_edf2_start:ind_edf2_end);
%         edf1_plot_cell{kkkk,t1t} = sig1{edf_sensor_ind}(ind_edf1_start:ind_edf1_end);
%         %     edf2_plot_cell{kkkk,t1t} = sig2{edf_sensor_ind}(ind_edf2_start:ind_edf2_end);
%         
%         
%         data_trunc = data(ind_lc_start:ind_lc_end,sensor_indx);
%         
%         sensor_indx
%         
%         cop_data_trunc_y = cop_y(ind_lc_start:ind_lc_end);
        
        
            
        
        %
        %
        %         [empty_ind,empty_seg] = LC_empty_bed(ALL_SUM,10,m_ind);
        %
        %
        %         if ~isempty(m_seg) && isempty(empty_seg)
        %             empty_seg = m_seg;
        %             empty_ind = find(m_ind);
        %         end
        
        %         cop_y_nm = LC_remove_data(COP_Y',empty_seg,1);
        %         cop_x_nm = LC_remove_data(COP_X',empty_seg,1);
        
        
        
%        lc_aligned{kkkk,t1t} =  cop_data_trunc_y;
%         lc_aligned{kkkk,t1t} = smooth(cop_data_trunc_y,500);
%        LC_sig = {lc_aligned{kkkk,t1t}-mean(lc_aligned{kkkk,t1t})'};
%        LC_sig{1} = LC_sig{1}.*(1-m_ind(ind_lc_start:ind_lc_end));
%        LC_sig{1} = LC_sig{1} - smooth(LC_sig{1},2000)';
%        is_there_movement{kkkk,t1t} = prod((1-m_ind(ind_lc_start:ind_lc_end)));
%         if(is_there_movement{kkkk,t1t} == 0)
%             LC_sig{1}(abs(LC_sig{1}) > 2.5*std(LC_sig{1})) = 0;
%             pressure_aligned{kkkk,t1t}(abs(LC_sig{1}) > 2*std(LC_sig{1})) = 0;
%         end
%        lc_aligned{kkkk,t1t}  = LC_sig{1};
        
        % %         LC_sig{1}(abs(LC_sig{1}) > 2*std(LC_sig{1})) = 0;
        %                 NM_SEG = zeros(1,length(t_lc));
        %
        % %         Set movement indices to '1'.
        %                 NM_SEG(empty_ind) = 1;
        %
        % %         Shorten 'NM_SEG' to only include the time period between 'start_time'
        % %         and 'stop_time'.
        %                 nm_seg = NM_SEG(t_lc < t_end & t_lc > t_start);
        %
        % %         Non-movement segments.
        %                 seg = Binary_2_SEG(nm_seg,0);
%        seg = [1 59999];
        %         seg1 = [1 floor((ind_edf1_end-ind_edf1_start)/2); floor((ind_edf1_end-ind_edf1_start)/2) ind_edf1_end-ind_edf1_start];
%        seg1 = [1 ind_edf1_end-ind_edf1_start];
%        pressure_aligned_rr{1} = sig1{edf_sensor_ind}(ind_edf1_start:ind_edf1_end)/max(abs(sig1{edf_sensor_ind}(ind_edf1_start:ind_edf1_end)));
        %         [Xsum, rr, T, Hz ] = Smart_LC_2_RR_Combo(LC_sig,seg,fs,t_lc(t_lc < t_end & t_lc > t_start),1,1);
        %         [Xsum_psg, rr_psg, T_psg, Hz_psg ] = Smart_LC_2_RR_Combo(pressure_aligned_rr,seg1,info1.FS(edf_sensor_ind),t_edf1{edf_sensor_ind}(ind_edf1_start:ind_edf1_end),1,1);
        %         [pks_lc,locs_lc,wid_lc,prom_lc] = findpeaks(LC_sig{1});
        %         [pks_psg,locs_psg,wid_psg,prom_psg] = findpeaks(pressure_aligned_rr{1});
        %         prom_lc_norm = prom_lc/max(prom_lc);clc
        %         prom_psg_norm = prom_psg/max(prom_psg);
        %
        %         rr = sum((prom_lc_norm > 0.15))/2
        %         rr_psg = sum((prom_psg_norm > 0.15))/2
        %if(plot_flag == 1)
%        plot(t_lc_window,LC_sig{1}/max(abs(LC_sig{1}))); hold on; plot(t_edf1{edf_sensor_ind}(ind_edf1_start:ind_edf1_end),pressure_aligned_rr{1});
%        datetick;
        %pause; 
        %clf;
        %end;
%         [Pxx,F] = periodogram(LC_sig{1},[],[],fs);
%         [maxPxx,indmaxPxx] = max(Pxx(F > 0.1));
%         rr = F(indmaxPxx+max(find(F < 0.1)))*60;
%         
%         
%         [Pxx1,F1] = periodogram(pressure_aligned_rr{1},[],[],info1.FS(edf_sensor_ind));
%         [maxPxx1,indmaxPxx1] = max(Pxx1(F1>0.1));
%         rr_psg = F1(indmaxPxx1+max(find(F1 < 0.1)))*60;
        
        %     lc_aligned{kkkk} = cop_y;
%         lc_aligned{kkkk,t1t} = lc_aligned{kkkk,t1t}-mean(lc_aligned{kkkk,t1t});
        %     [lc_aligned, pressure_aligned] = align_data_segment(cop_data_trunc_y,pressure_signal);
%         kkkk
%         estimate_corr(kkkk,t1t) = max(abs(xcorr(lc_aligned{kkkk,t1t},pressure_aligned{kkkk,t1t},'coeff')));
%         pressure_aligned_norm{kkkk,t1t} = pressure_aligned{kkkk,t1t}/max(abs((pressure_aligned{kkkk,t1t})));
%         
%         
%         if(plot_flag == 2)
%             [Pxx1,F1] = periodogram(pressure_aligned_rr{1},[],[],info1.FS(edf_sensor_ind));
%             [Pxx,F] = periodogram(LC_sig{1},[],[],fs);
%             plot(F1,Pxx1./norm(Pxx1));
%             hold on;
%             plot(F,Pxx./norm(Pxx));
%             axis([0 10 0 0.5])
%             pause;
%             clf;
%         end
        
%         if(plot_flag == 3 )%&& is_there_movement{kkkk,t1t} == 1)
%             
%             cepstral_coeff{kkkk,t1t} = ifft(log(abs(fft(decimate(lc_aligned{kkkk,t1t},dec_factor)))));
%             cepstral_coeff{kkkk,t1t} = (cepstral_coeff{kkkk,t1t}(10:end/2)-mean(cepstral_coeff{kkkk,t1t}(10:end/2)))./max(abs(cepstral_coeff{kkkk,t1t}(10:end/2)-mean(cepstral_coeff{kkkk,t1t}(10:end/2))));
% %             cepstral_coeff{kkkk,t1t} = hist(Pxx);
%             lpc_coeff{kkkk,t1t} = lpc(lc_aligned{kkkk,t1t},15);
%             if(sum(isnan(cepstral_coeff{kkkk,t1t})) > 0)
%                 display('error');
%                 continue;
%             end
%             h = figure('units','normalized','outerposition',[0 0 1 1]);
%             h1 = zeros(1,4);
%             h1(1) = subplot(511)
%             
%              plot(t_lc_window,LC_sig{1}/max(abs(LC_sig{1}))); hold on; plot(t_edf1{edf_sensor_ind}(ind_edf1_start:ind_edf1_end),pressure_aligned_rr{1});
%             datetick;
%             title(['Subject ' num2str(subind) ' Position ' positionstring(t1t,:)])
%             h2(2) = subplot(512)
% %             bins_vec = linspace(0,1,20);
% %             hist(cepstral_coeff{kkkk,t1t},bins_vec);
%             stem(F,Pxx,'.');
% %             hist()
%             axis([0,2,0,max(Pxx)])
%             title('Power Spectrum')
%                 bins_vec = linspace(-1,1,500/dec_factor);
%             h1(3) = subplot(513)
% %             mdl = ar(LC_sig{1}/max(abs(LC_sig{1})),5,'ls','Ts',1/500);
% %             ar_coeff{kkkk,t1t} = mdl.a;
%             [PKS,LOCS,wdth{kkkk,t1t},prm{kkkk,t1t}] = findpeaks(LC_sig{1});
%             bins1 = linspace(-1,1,10);
%             hist((wdth{kkkk,t1t}-mean(wdth{kkkk,t1t}))/max(abs(wdth{kkkk,t1t})),bins1)
%             title('Histogram of widths of peaks')
%             hold on
%             
% 
% 
%             subplot(515)
%             hist(LC_sig{1}/max(abs(LC_sig{1})),bins_vec)
%             
%             title('Histogram of raw time series')
%             h1(4) = subplot(514);
% %             stem(cepstral_coeff{kkkk,t1t}(10:end/2),'.')
% %             yl = ylim;r
%             whiten_LC = LC_sig{1}./abs(fft(LC_sig{1}));
%             whiten_LC = whiten_LC/max(abs(whiten_LC));
% %             [N_distr{kkkk,t1t},x_val_hist{kkkk,t1t}] = hist(whiten_LC,bins_vec);
% %             bar(x_val_hist{kkkk,t1t},N_distr{kkkk,t1t})
%             hist((prm{kkkk,t1t}-mean(prm{kkkk,t1t}))./max(abs((prm{kkkk,t1t}))),bins1)
%             title('Histogram of prominences')
%             
% %             axis([-0.3 0.3 yl])
% %             pause
% %             clf
% %             fig_ind = fig_ind+1;
% %             linkaxes(h1,'x');
%             savefig(h,['.\histograms_12subj_studies\sub_' num2str(subind) 'position_indx' num2str(t1t)])
%             saveas(h,['.\histograms_12subj_studies\sub_' num2str(subind) 'position_indx' num2str(t1t) '.png'])
%         end
%         
%         %     [Xsum,rr,T,Hz] = Smart_LC_2_RR_Combo(data_cell,,fs,t_lc,0,0);
%         %         rr((rr~=0));
%         %     Hz
%         rr1 = rr(rr~=0 & ~isnan(rr));
%         rr_psg1 = rr_psg(rr_psg~=0 & ~isnan(rr_psg));
%         rr_cell{kkkk,t1t} = mean(rr1)
%         %         Hz_cell{kkkk,t1t} = Hz;
%         rr_psg_cell{kkkk,t1t} = mean(rr_psg1)
%         %     plot(t_edf1cell{kkkk},edf1_plot_cell{kkkk}/norm(edf1_plot_cell{kkkk}));
%         %     subplot(211)
%         %     plot(data(start_indx:(start_indx+60000),1));
%         %         A1 = figure(1);
%         %         set(A1, 'Units', 'normalized', 'Position', [0,0,1,1]);
%         %         % %     subplot(212)
%         %         plot(t_edf1{edf_sensor_ind}(ind_edf1_start:ind_edf1_end),sig1{edf_sensor_ind}(ind_edf1_start:ind_edf1_end)/max(abs(sig1{edf_sensor_ind}(ind_edf1_start:ind_edf1_end))),'linewidth',1.5);
%         %         hold on
%         %         plot(t_lc(ind_lc_start:ind_lc_end),(lc_aligned{kkkk,t1t})/max((abs(lc_aligned{kkkk,t1t}))),'linewidth',1.5);
%         %         datetick
%         %         title(['Blue: nasal cannula, Red: load cell CoP_y study '  num2str(subind) '\newline Subj 1 Sensors ' num2str(sensor_indx)],'fontsize',16);
%         %
%         %         saveas(A1,['.\subj_1_figs_report\figure_comparing_single_person_on_bed_study_1st_study_'  num2str(subind) '_' num2str(sensor_indx)],'png')
%         %         saveas(A1,['.\subj_1_figs_report\figure_comparing_single_person_on_bed_study_1st_study_' num2str(subind) '_' num2str(sensor_indx)],'pdf')
%         %         saveas(A1,['.\subj_1_figs_report\figure_comparing_single_person_on_bed_study_1st_study_' num2str(subind) '_' num2str(sensor_indx), '.fig'])
%         %
%         %         close all
%         %     rr_lc(kkkk,:) = estimate_rr(lc_aligned);
%         %     rr_pressure(kkkk,:) = estimate_rr(pressure_aligned);
%         
%         
%         
    for p = 4:24
        
        DataSheet(dataset,p) = (DataSheet(dataset,p)- min(DataSheet(dataset,p)))/ ...
                (max(DataSheet(dataset,p)) -min(DataSheet(dataset,p)));
        
    end
end
% strings= char(labels);
% figure;
% high = 150;
% low = -100;
% stepsize = (high-low)/100;
for j = 1:10 
% plot rms value  
% 
% edges = low:stepsize:high;
% 
% n= histc(histo(j,:),edges);
% subplot(5,2,j);
% bar(edges,n,'histc');
% title(strings(j,:));
end

    
    

%%%%%% Normalize table 
%     for j = 4:17
%         max_col = max(DataSheet(1:(kkkk)*positions,j));
%         min_col = min(DataSheet(1:(kkkk)*positions,j));
%         DataSheet(1:(kkkk)*positions,j) = (DataSheet(1:(kkkk)*positions,j)-min_col)/(max_col-min_col);
%         
% 
%     end
%end

save sensor_1_SP_NEW_SENSORS
save (filename,'DataSheet')