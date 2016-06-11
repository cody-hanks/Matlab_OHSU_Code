%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cody hanks - 4/7/2016 - Run clasifier on determined Data for result 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear;
close all;
clc;
for endpoint=[4:24]

load('DataSheet3.mat');
tic

 [coeff,score,latent,tsquared,explained,mu] = pca(DataSheet(1:210,4:24));
 DataSheet(1:210,4:24)=score;
 save('pcadat.mat','mu','coeff');
 %disp(latent)



svmfeat  = [4:endpoint];
svmfeat1 = [4:endpoint];
svmfeat2 = [4:endpoint];

%svmfeat  = [7,6,23];
%svmfeat1 = [20,18,8];
%svmfeat2 = [5,4,13,10 ];


bayfeat  = [18,4,5,6,7,8];
bayfeat1 = [7,4];
bayfeat2 = [11,4,5,6,7,8,9,10,12,13];

knnfeat  = [23,4,5,6,7,8,9,10,11,12];
knnfeat1 = [17,4,5,6,7,8,9,10,11,12,13];
knnfeat2 = [7,4,5];


featurel1 = svmfeat;
featurel2 = svmfeat1;
featurel3 = svmfeat2;

totalrows = 210;
folds = 21;
datalength = totalrows/folds;
rows = 1:datalength:totalrows;
resultscolumn = [29,3,3];
box1 = 118;
kern1 = 29;
box2 = 172;
kern2 = 6;
box3 = 92;
kern3 = 12;
numbs1 = 20;
numbs2 = 1;
numbs3 = 3;

boxsize = linspace(1,200,200);%130
kernel = linspace(1,150,150);%57
numneighbors = linspace(1,150,150);
Resultssheet = zeros(length(boxsize),length(kernel));
%Resultssheet = zeros(length(numneighbors));
 %for box = 1:length(boxsize)
     %parfor kern = 1:length(kernel)
%    parfor numbs = 1:length(numneighbors)
correctlv1 = 0;
correctlv2 =0;
correctlv3 = 0;
back =0;
front =0;
left=0;
right=0;
fb = 0;
lr = 0;
    FBIND =[ ];
    for rowid = 1:8:(folds-1)*8
      FBIND = [FBIND rowid];
      FBIND = [FBIND rowid+3];
      FBIND = [FBIND rowid+4];
      FBIND = [FBIND rowid+7];
    end 
    LRIND =[];
    for rowid = 1:8:(folds-1)*8
      LRIND = [LRIND rowid+ 1];
      LRIND = [LRIND rowid+ 2];
      LRIND = [LRIND rowid+ 5];
      LRIND = [LRIND rowid+ 6];
    end
    
    
    for fold = 1:folds
        dataset = [];
        % build the Dataset rows (these are all rows except the current
        % fold that is the test data ... note skipping the 2 calib rows
        for datarows = 1:folds
            if(datarows < fold || datarows > fold) 
                dataset = [dataset rows(datarows)+2:rows(datarows)+datalength-1];
            end
        end
        % generate other tables of Data for training 
        datasetFB = dataset(FBIND);
        datasetLR = dataset(LRIND);
        testset = [rows(fold)+2:rows(fold)+(datalength-1)];
    
        mdl = fitcsvm(DataSheet(dataset,featurel1),DataSheet(dataset,resultscolumn(1)),'KernelFunction','rbf','Standardize',true,'BoxConstraint',box1,'KernelScale',kern1);
        mdlfb = fitcsvm(DataSheet(datasetFB,featurel2),DataSheet(datasetFB,resultscolumn(2)),'KernelFunction','rbf','Standardize',true,'BoxConstraint',box2,'KernelScale',kern2);
        mdllr = fitcsvm(DataSheet(datasetLR,featurel3),DataSheet(datasetLR,resultscolumn(3)),'KernelFunction','rbf','Standardize',true,'BoxConstraint',box3,'KernelScale',kern3);
        %mdl = fitcknn(DataSheet(dataset,featurel1),DataSheet(dataset,resultscolumn(1)),'NumNeighbors',numbs1);
        %mdlfb = fitcknn(DataSheet(datasetFB,featurel2),DataSheet(datasetFB,resultscolumn(2)),'NumNeighbors',numbs2);
        %mdllr = fitcknn(DataSheet(datasetLR,featurel3),DataSheet(datasetLR,resultscolumn(3)),'NumNeighbors',numbs3);
        predictionl1 = predict(mdl,DataSheet(testset,featurel1));
        correctlv1 =correctlv1 + sum(DataSheet(testset,resultscolumn(1))==predictionl1);
        
     for pred = 1:length(predictionl1)
            if (predictionl1(pred) ==1 )
                if predictionl1(pred) == DataSheet(testset(pred),resultscolumn(1))
                    fb =  fb +1;
                end
                pred1 = predict(mdlfb,DataSheet(testset(pred),featurel2));
                if(pred1 == DataSheet(testset(pred),resultscolumn(2)))
                    correctlv2 = correctlv2 +1;
                    if pred1 ==1
                        back =back+1;
                    else
                        front = front +1;
                    end

                end
            elseif (predictionl1(pred) ~= 1)
                if predictionl1(pred) == DataSheet(testset(pred),resultscolumn(1))
                    lr =  lr +1;
                end
                pred2 = predict(mdllr,DataSheet(testset(pred),featurel3));
                if(pred2 == DataSheet(testset(pred),resultscolumn(3)))
                    correctlv3 = correctlv3 +1;
                    if pred2 ==2
                        left = left +1;
                    else
                        right = right +1;
                    end
                end
            end
      end% end for each pred
        
      
    end%end folds 
    % Resultssheet(numbs) = correctlv3;
    %Resultssheet(box,kern) = correctlv3;
    %end% kernel 
%end%boxconst
 toc
 %save('output3_pca.mat','kernel','boxsize','Resultssheet');
 %figure('visible','off');
 %surf(kernel,boxsize,Resultssheet);
 %saveas(gcf,'file3_pca.fig','fig');
 
 %clear
 %openfig('file.fig','new','visible')
 disp('[4:')
 disp(endpoint)
 disp(']')
 disp('[correctlv1 correctlv2 correctlv3 fb lr front back left right]')
[correctlv1 correctlv2 correctlv3 fb lr front back left right]    
 %figure();
 %surf(Resultssheet);
 clear;
end  
    
    
    
    
    
    
    
    
    
    