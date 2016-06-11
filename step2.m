%%%
% Cody Hanks 
% 6/7/2016  
% Step 1 Use Data minute 1 to get the Coeff and get the Explained Values 
%
%%%
clc
clear
close all

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
boxsize = linspace(1,201,201);%200
kernelsize = linspace(1,151,151);%150

fb = [1:4:(subjects-1)*8];
fb = [fb 4:4:(subjects-1)*8];
lr = [2:4:(subjects-1)*8];
lr = [lr 3:4:(subjects-1)*8];

test_fb = [1:4:8,4:4:8];
test_lr = [2:4:8,3:4:8];



% ok now add the plots for leave one out 
for level =1:3
 for boxparam = 1:length(boxsize)
     box=boxsize(boxparam);
  for kernparam = 1:length(kernelsize)
      kern = kernelsize(kernparam);
    for subject = 1:subjects
        TempDataset = Datasheet2(~ismember(1:subjects*lines,[((subject-1)*lines)+1:(subject*lines)]),:);
        [tcoeff,tscore,tlatent,ttsquared,texplained,tmu]=pca(TempDataset(:,4:24));
        Singleperson = Datasheet2(~~ismember(1:subjects*lines,[((subject-1)*lines)+1:(subject*lines)]),:);
        Singleperson(:,4:24) = (tcoeff'*(Singleperson(:,4:24)-repmat(tmu,lines,1))')';
        %-get the data set ready for this. 
        TempDataset(1:end,4:24) = tscore;
        switch level
            case 1
               mdl = fitcsvm(TempDataset(:,4:endpoint),TempDataset(:,resultscolumn(1)),'KernelFunction','rbf','Standardize',true,'BoxConstraint',box,'KernelScale',kern);
               prediction = predict(mdl,Singleperson(:,4:endpoint));
               per_cor_lv1(subject)  = sum(Singleperson(:,resultscolumn(1))==prediction)/lines;
            case 2
               mdl_fb = fitcsvm(TempDataset(fb,4:endpoint),TempDataset(fb,resultscolumn(2)),'KernelFunction','rbf','Standardize',true,'BoxConstraint',box,'KernelScale',kern);
               prediction = predict(mdl,Singleperson(test_fb,4:endpoint));
               per_cor_lv2(subject)  = sum(Singleperson(test_fb,resultscolumn(2))==prediction)/length(test_fb);
            case 3
               mdl_fb = fitcsvm(TempDataset(lr,4:endpoint),TempDataset(lr,resultscolumn(2)),'KernelFunction','rbf','Standardize',true,'BoxConstraint',box,'KernelScale',kern);
               prediction = predict(mdl,Singleperson(test_lr,4:endpoint));
               per_cor_lv3(subject)  = sum(Singleperson(test_lr,resultscolumn(2))==prediction)/length(test_fb);
        end       
        clear TempDataset tcoeff tscore tlatent ttsquared texplained tmu tendpoint
    end  % end subjects 
        switch level
            case 1
                svm_params_lv1(box,kern) = mean(per_cor_lv1(:));
            case 2
                svm_params_lv2(box,kern) = mean(per_cor_lv2(:));
            case 3
                svm_params_lv3(box,kern) = mean(per_cor_lv3(:));
        end
  end% kernel 
 end% box size 
end % level 


if surf_plt
    figure;
    surf(svm_params_lv1);
    figure;
    surf(svm_params_lv2);
    figure
    surf(svm_params_lv3);
end


