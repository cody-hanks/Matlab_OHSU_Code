%%%
% Cody Hanks 
% 6/7/2016  
% Step 2 uses data from step 1 to generate a box and kernel for each item
% Step 2 aprox 2.1 hours on 16 cores 
%%%
clc
clear
close all

tic

surf_plt =1;

% subjects 
subjects = 21;
% results
resultscolumn = [29,3,3];
lines = 8;


%load data from step 1 
load('minute1.mat');
load('pcadatamin1.mat');
load('endpoint.mat');

%box and kern setup 
boxsize = 200;
kernelsize = 150;

fb = [1:4:(subjects-1)*8];
fb = [fb 4:4:(subjects-1)*8];
lr = [2:4:(subjects-1)*8];
lr = [lr 3:4:(subjects-1)*8];

test_fb = [1:4:8,4:4:8];
test_lr = [2:4:8,3:4:8];

TempDataset = zeros(168,32);
box=0;kern=0;


% warning has to do with tsquared we dont use and its only due to
% parforloop 
id = ('stats:pca:ColRankDefX');
warning('off',id);
per_cor = zeros(3,subjects+1,boxsize,kernelsize);
% ok now add the plots for leave one out 
for level =1:3
 for subject= 1:subjects
     % get the Temp data set with all rows except current subject 
    TempDataset = Datasheet2(~ismember(1:subjects*lines,[((subject-1)*lines)+1:(subject*lines)]),:);    
    % get PCA parameters 
    [tco,tsc,tlat,ttsq,texp,tmu]=pca(TempDataset(:,4:24));
    % get the subject Database 
    Singleperson = Datasheet2(~~ismember(1:subjects*lines,[((subject-1)*lines)+1:(subject*lines)]),:);
    % apply the coeff transform to the data cols 4:24
    Singleperson(:,4:24) = (tco'*(Singleperson(:,4:24)-repmat(tmu,lines,1))')';
    %-apply scored matrix to labels 
    TempDataset(1:end,4:24) = tsc;
    %iterate acros box and kernel for each subject and level 
    for box = 1:boxsize 
     parfor kern = 1:kernelsize  % parallel to make this take reasonable amounts of time
        
        switch level
            case 1
               % evaluate the score for each 
               mdl = fitcsvm(TempDataset(:,4:endpoint),TempDataset(:,resultscolumn(1)),'KernelFunction','rbf','Standardize',true,'BoxConstraint',box,'KernelScale',kern);
               prediction = predict(mdl,Singleperson(:,4:endpoint));
               % save percent correct
               per_cor(level,subject,box,kern) = sum(Singleperson(:,resultscolumn(1))==prediction)/lines;
            case 2
               % evaluate the score for each 
               mdl_fb = fitcsvm(TempDataset(fb,4:endpoint),TempDataset(fb,resultscolumn(2)),'KernelFunction','rbf','Standardize',true,'BoxConstraint',box,'KernelScale',kern);
               prediction = predict(mdl_fb,Singleperson(test_fb,4:endpoint));
               % save the percent correct 
               per_cor(level,subject,box,kern) = sum(Singleperson(test_fb,resultscolumn(2))==prediction)/length(test_fb);
            case 3
               %evaluate for left/right
               mdl_lr = fitcsvm(TempDataset(lr,4:endpoint),TempDataset(lr,resultscolumn(2)),'KernelFunction','rbf','Standardize',true,'BoxConstraint',box,'KernelScale',kern);
               prediction = predict(mdl_lr,Singleperson(test_lr,4:endpoint));
               % save % correct 
               per_cor(level,subject,box,kern) = sum(Singleperson(test_lr,resultscolumn(2))==prediction)/length(test_lr);
        end 
        
     end %kernparam
    end %boxparam
 end% subject  
 %%%%%%%%%%%%%%%%%%%%
 for box = 1:boxsize
    for kern = 1:kernelsize
        % combine indevidual % correct into one average 
        per_cor(level,22,box,kern)  =  mean(per_cor(level,1:21,box,kern) );
    end
 end
 %%%%%%%%%%%%%%%%%%%
end% level 
save 'Box_Kern_Search.mat' 'per_cor';

toc

% get from the different levels the plot of %correct and find max for each.
% 

tmp(1:boxsize,1:kernelsize) = per_cor(1,22,:,:);
highlv = max(tmp(:));
[box1 kern1] = find(tmp==highlv,1);
save 'Parameterslv1.mat' 'box1' 'kern1'
if surf_plt
    figure;
    surf(tmp);
    savefig('fig1.fig');
end
clear tmp
tmp(1:boxsize,1:kernelsize) = per_cor(2,22,:,:);
highlv = max(tmp(:));
[box2 kern2] = find(tmp==highlv,1);
save 'Parameterslv2.mat' 'box2' 'kern2'
if surf_plt
    figure;
    surf(tmp);
    savefig('fig2.fig');
end
clear tmp
tmp(1:boxsize,1:kernelsize) = per_cor(3,22,:,:);
highlv = max(tmp(:));
[box3 kern3] = find(tmp==highlv,1);
save 'Parameterslv3.mat' 'box3' 'kern3'
if surf_plt
    figure;
    surf(tmp);
    savefig('fig3.fig');
end
clear tmp
