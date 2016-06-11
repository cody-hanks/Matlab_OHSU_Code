%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cody hanks - 4/7/2016 - Run clasifier on determined Data for result 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear;
close all;
clc;
for endpoint=[5:24]
load('DataSheet3.mat');

 [coeff,score,latent,tsquared,explained,mu] = pca(DataSheet(1:210,4:24));
 DataSheet(1:210,4:24)=score;
 %disp(explained);
 save('pcadat.mat','mu','coeff','explained');
 load('pcadat.mat');
 %cntremar = linspace(5,1,5);
 %for cntleft = 1:100
 %   for cntremv =1:5
 %       cntrem=cntremar(cntremv);
 %leave_out = randi(21,1,cntrem);
 %disp([ 'leave out:' ])
 %disp([leave_out])
%svmfeat  = [7,24,12,18,13,10];
%svmfeat1 = [11,13,15,24,21,10];
%svmfeat2 = [5,4,23];


svmfeat  = [4:endpoint];
svmfeat1 = [4:endpoint];
svmfeat2 = [4:endpoint];

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

% boxsize = linspace(1,200,200);%130
% kernel = linspace(1,150,150);%57
numneighbors = linspace(1,150,150);
% Resultssheet = zeros(length(boxsize),length(kernel));
Resultssheet = zeros(length(numneighbors));
% for box = 1:length(boxsize)
%     parfor kern = 1:length(kernel)
%    parfor numbs = 1:length(numneighbors)
correctlv1 = 0;
correctlv2 =0;
correctlv3 = 0;
back =0;
front =0;
left=0;
right=0;
fb=0;
lr=0;
    FBIND =[ ];
    for rowid = 1:8:96
      FBIND = [FBIND rowid];
      FBIND = [FBIND rowid+3];
      FBIND = [FBIND rowid+4];
      FBIND = [FBIND rowid+7];
    end 
    LRIND =[];
    for rowid = 1:8:96
      LRIND = [LRIND rowid+ 1];
      LRIND = [LRIND rowid+ 2];
      LRIND = [LRIND rowid+ 5];
      LRIND = [LRIND rowid+ 6];
    end
    
    
    dataset = [];
    
        % build the Dataset rows (these are all rows except the current
        % fold that is the test data ... note skipping the 2 calib rows
        for datarows = 1:folds
            %if( ~ ismember([datarows],leave_out))
                dataset = [dataset rows(datarows)+2:rows(datarows)+datalength-1];
            %end
        end
        % generate other tables of Data for training 
        datasetFB = dataset(FBIND);
        datasetLR = dataset(LRIND);
        
    
        mdl = fitcsvm(DataSheet(dataset,featurel1),DataSheet(dataset,resultscolumn(1)),'KernelFunction','rbf','Standardize',true,'BoxConstraint',box1,'KernelScale',kern1);
        mdlfb = fitcsvm(DataSheet(datasetFB,featurel2),DataSheet(datasetFB,resultscolumn(2)),'KernelFunction','rbf','Standardize',true,'BoxConstraint',box2,'KernelScale',kern2);
        mdllr = fitcsvm(DataSheet(datasetLR,featurel3),DataSheet(datasetLR,resultscolumn(3)),'KernelFunction','rbf','Standardize',true,'BoxConstraint',box3,'KernelScale',kern3);
        
        %trained on data now run test 
        load('DataSheet4.mat');
        %DataSheet(1:210,4:24) = (coeff'*(DataSheet(1:210,4:24)-repmat(mu,210,1))')';
        
        coeff = pca(DataSheet(1:210,4:24));
        newdatasheet = DataSheet(1:210,4:24)*coeff';
        DataSheet(1:210,4:24) = newdatasheet;
        
        
        
        testset =[];
        for datarows = 1:folds
                testset = [testset rows(datarows)+2:rows(datarows)+datalength-1];
        end
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
        
      
    %end%end folds 
    % Resultssheet(numbs) = correctlv3;
    %end% kernel 
% end%boxconst[correctlv1 correctlv2 correctlv3 fb lr front back left right]
disp(['[4:' num2str(endpoint) ']' ]);
disp('[correctlv1 correctlv2 correctlv3 fb lr front back left right]')
disp([correctlv1 correctlv2 correctlv3 fb lr front back left right])

%fiveline(cntrem,cntleft,:) = [correctlv1 correctlv2 correctlv3 fb lr front back left right];

 %   end%
 %clear
   end
    
 %minus=[mean(fiveline(5,:,1)) mean(fiveline(4,:,1)) mean(fiveline(3,:,1)) mean(fiveline(2,:,1)) mean(fiveline(1,:,1))];
 %save('fiveline.mat','fiveline');
    
    
    
    
    
    
    
    
    
    
    
    