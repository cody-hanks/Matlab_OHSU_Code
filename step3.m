%%%
% Cody Hanks 
% 6/7/2016  
% Step 3 use folds to determine % correct acros the number of features for
% each subject should have a % correct from each # of features 
% use error bar to plot std and mean 
%%%
clc 
clear 
close all

tic
% load the Data 
load('minute1.mat');
load('Parameterslv1.mat');
load('Parameterslv2.mat');
load('Parameterslv3.mat');

% subjects 
subjects = 21;
% results
resultscolumn = [29,3,3];
lines = 8;



%row identifier setup 
fb = [1:4:(subjects-1)*8];
fb = [fb 4:4:(subjects-1)*8];
lr = [2:4:(subjects-1)*8];
lr = [lr 3:4:(subjects-1)*8];

test_fb = [1:4:8,4:4:8];
test_lr = [2:4:8,3:4:8];


%Data holder 
TempDataset = zeros(168,32);
per_cor = zeros(3,subjects+2,21);
% warning has to do with tsquared we dont use and its only due to
% parforloop 
id = ('stats:pca:ColRankDefX');
warning('off',id);


% need to run by leave 1 out across the subjects
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
    %iterate acros # of features to get % correct for each # of features
    %with given box and kernel 
    parfor features = 1:21
        switch level
            case 1
               % evaluate the score for each 
               mdl = fitcsvm(TempDataset(:,4:features+3),TempDataset(:,resultscolumn(1)),'KernelFunction','rbf','Standardize',true,'BoxConstraint',box1,'KernelScale',kern1);
               prediction = predict(mdl,Singleperson(:,4:features+3));
               % save percent correct
               per_cor(level,subject,features) = sum(Singleperson(:,resultscolumn(1))==prediction)/lines;
            case 2
               % evaluate the score for each 
               mdl_fb = fitcsvm(TempDataset(fb,4:features+3),TempDataset(fb,resultscolumn(2)),'KernelFunction','rbf','Standardize',true,'BoxConstraint',box2,'KernelScale',kern2);
               prediction = predict(mdl_fb,Singleperson(test_fb,4:features+3));
               % save the percent correct 
               per_cor(level,subject,features) = sum(Singleperson(test_fb,resultscolumn(2))==prediction)/length(test_fb);
            case 3
               %evaluate for left/right
               mdl_lr = fitcsvm(TempDataset(lr,4:features+3),TempDataset(lr,resultscolumn(2)),'KernelFunction','rbf','Standardize',true,'BoxConstraint',box3,'KernelScale',kern3);
               prediction = predict(mdl_lr,Singleperson(test_lr,4:features+3));
               % save % correct 
               per_cor(level,subject,features) = sum(Singleperson(test_lr,resultscolumn(2))==prediction)/length(test_lr);
        end 
    end% features  
 end% subject
 %%%%%%%%%%%%%
 for features = 1:21 
    per_cor(level,22,features) = std(per_cor(level,1:21,features));
    per_cor(level,23,features) = mean(per_cor(level,1:21,features));
 end
 %%%%%%%%%%%%%
end% level 
toc

h=figure;
tmp(1:21) = per_cor(1,23,:);
tmp2(1:21) = per_cor(1,22,:);
mx1 = find(tmp==max(tmp),1);
save 'feat1.mat' 'mx1';
errorbar(1:21,tmp,tmp2);
hold on ;
errorbar(mx1,tmp(mx1),tmp2(mx1),'r');
savefig('errorbar1.fig');
saveas(h,'errorbar1.jpg');

h=figure;
tmp(1:21) = per_cor(2,23,:);
tmp2(1:21) = per_cor(2,22,:);
mx2 = find(tmp==max(tmp),1);
save 'feat2.mat' 'mx2';
errorbar(1:21,tmp,tmp2);
hold on ;
errorbar(mx2,tmp(mx2),tmp2(mx2),'r');
savefig('errorbar2.fig')
saveas(h,'errorbar2.jpg');

h=figure;
tmp(1:21) = per_cor(3,23,:);
tmp2(1:21) = per_cor(3,22,:);
mx3= find(tmp==max(tmp),1);
save 'feat3.mat' 'mx3';
errorbar(1:21,tmp,tmp2);
hold on ;
errorbar(mx3,tmp(mx3),tmp2(mx3),'r');
savefig('errorbar3.fig')
saveas(h,'errorbar3.jpg');





