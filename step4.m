%%%
% Cody Hanks 
% 6/7/2016  
% Step 4 use determined features and determined kernel/box along with coeff
% from step 1 to apply to data minute 2 then test  for accuracy 
%
%
%%%
clc 
clear
close all;


% load varialbes
subjects =21;
datalength = 10;
lines = 8;


%load minute 2 data 
load('DataSheet4.mat');

%creat the row starts for each subject 
rows = 1:datalength:(subjects*datalength);


%Get the Dataset for Testing
Datasheet_tst = zeros(subjects*lines,32);
for subject = 1:subjects
    Datasheet_tst(((subject-1)*lines)+1:(subject*lines),1:32) ...
    = DataSheet(rows(subject)+2:rows(subject)+datalength-1,1:32);
end
save('minute2.mat','Datasheet_tst');

clear

% reaload the data variables 
subjects =21;
datalength = 10;
lines = 8;
% results
resultscolumn = [29,3,3];
fb = [1:4:(subjects-1)*8];
fb = [fb 4:4:(subjects-1)*8];
lr = [2:4:(subjects-1)*8];
lr = [lr 3:4:(subjects-1)*8];

tic
% load the Data 
load('minute1.mat');% Datasheet2
load('minute2.mat');% Datasheet_tst
load('Parameterslv1.mat');% box1 kern1
load('Parameterslv2.mat');% box2 kern2
load('Parameterslv3.mat');% box3 kern3
load('feat1.mat'); % mx1 Max % number of features
load('feat2.mat'); % mx2
load('feat3.mat'); % mx3
load('pcadatamin1.mat'); % 'mu','coeff','explained','score'

% apply the transform from step 1 to Data for test 
Datasheet_tst(:,4:24) = (coeff'*(Datasheet_tst(:,4:24)-repmat(mu,(subjects*lines),1))')';
% re apply transform from step 1 to Train Data 
%[coeff,score,latent,tsquared,explained,mu] = pca(Datasheet2(1:subjects*lines,4:24));
Datasheet2(:,4:24)    = (coeff'*(Datasheet2(:,4:24)   -repmat(mu,(subjects*lines),1))')';



%Train models 
mdl = fitcsvm(Datasheet2(:,4:mx1+3),Datasheet2(:,resultscolumn(1)),'KernelFunction','rbf','Standardize',true,'BoxConstraint',box1,'KernelScale',kern1);
mdl_fb = fitcsvm(Datasheet2(fb,4:mx2+3),Datasheet2(fb,resultscolumn(2)),'KernelFunction','rbf','Standardize',true,'BoxConstraint',box2,'KernelScale',kern2);
mdl_lr = fitcsvm(Datasheet2(lr,4:mx3+3),Datasheet2(lr,resultscolumn(3)),'KernelFunction','rbf','Standardize',true,'BoxConstraint',box3,'KernelScale',kern3);

% ok test models 
prediction = predict(mdl,Datasheet_tst(:,4:mx1+3));
percent_correct1 = sum(Datasheet_tst(:,resultscolumn(1))==prediction)/length(Datasheet_tst);
lbl = Datasheet_tst(:,resultscolumn(1));
save 'predict1.mat' 'prediction' 'lbl'

prediction = predict(mdl_fb,Datasheet_tst(fb,4:mx2+3));
percent_correct2 = sum(Datasheet_tst(fb,resultscolumn(2))==prediction)/length(fb);
lbl = Datasheet_tst(fb,resultscolumn(2));
save 'predict2.mat' 'prediction' 'lbl'

prediction = predict(mdl_lr,Datasheet_tst(lr,4:mx3+3));
percent_correct3 = sum(Datasheet_tst(lr,resultscolumn(3))==prediction)/length(lr);
lbl = Datasheet_tst(lr,resultscolumn(2));
save 'predict3.mat' 'prediction' 'lbl'



