%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cody Hanks 
% 3/29/2016
% classifier function for inclusion into a search algorithm 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Greedy_pca(DataSheet,features,featuresl1,featuresl2,featuresl3,Classifier,Optimize)
    [coeff,score,latent,tsquared,explained,mu] = pca(DataSheet(1:210,4:24));
    DataSheet(1:210,4:24)=score;
    save('pcadat.mat','mu','coeff');
    disp(latent)
    %newdatasheet = DataSheet(1:210,4:24)*coeff';
    %DataSheet(1:210,4:24) = newdatasheet;
    if Optimize(2)==1
        featuersleft = setdiff(features,featuresl1);
    elseif Optimize(2)==2
        featuersleft = setdiff(features,featuresl2);
    elseif Optimize(2)==3
        featuersleft = setdiff(features,featuresl3);
    end
    ResultsFeatures = zeros(13,21);
    
    for ind = 1:length(featuersleft)
        Feature = featuersleft(ind);
    if Optimize(2)==1
        Newfeaturesl1 = [featuresl1 Feature];
    elseif Optimize(2)==2
        Newfeaturesl1 = featuresl1;
        Newfeaturesl2 = [featuresl2 Feature];
    elseif Optimize(2)==3
        Newfeaturesl1 = featuresl1;
        Newfeaturesl3 = [featuresl3 Feature];
    end    
        
            
    box1 = 118;
    kern1 = 29;
    box2 = 172;
    kern2 = 6;
    box3 = 55;
    kern3 = 41;


    totalrows = 210;
    folds = 21;
    datalength = totalrows/folds;
    rows = 1:datalength:totalrows;
    resultscolumn = [29,3,3];
    correctlv1 = 0;
    correctlv2 =0;
    correctlv3 = 0;
    % FB data set === 13 folds and 10 rows per subject ignoring row 0 1 due
    % in order to even out the number of Data in each group IE 2 rows of
    % each position in each test/train
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
        
        switch(Classifier)
            case 'knn'
                mdl = fitcknn(DataSheet(dataset,featuresl1),DataSheet(dataset,resultscolumn(1)));
                if(Optimize(2)==2)
                mdlfb = fitcknn(DataSheet(datasetFB,featuresl2),DataSheet(datasetFB,resultscolumn(2)));
                elseif (Optimize(2)==3)
                mdllr = fitcknn(DataSheet(datasetLR,featuresl3),DataSheet(datasetLR,resultscolumn(3)));
                end
            case 'svm'
                mdl = fitcsvm(DataSheet(dataset,Newfeaturesl1),DataSheet(dataset,resultscolumn(1)),'KernelFunction','rbf','Standardize',true,'BoxConstraint',box1,'KernelScale',kern1);
                if(Optimize(2)==2)
                mdlfb = fitcsvm(DataSheet(datasetFB,Newfeaturesl2),DataSheet(datasetFB,resultscolumn(2)),'KernelFunction','rbf','Standardize',true,'BoxConstraint',box2,'KernelScale',kern2);
                elseif (Optimize(2)==3)
                mdllr = fitcsvm(DataSheet(datasetLR,Newfeaturesl3),DataSheet(datasetLR,resultscolumn(3)),'KernelFunction','rbf','Standardize',true,'BoxConstraint',box3,'KernelScale',kern3);
                end
            case 'bay'
                mdl = NaiveBayes.fit(DataSheet(dataset,Newfeaturesl1),DataSheet(dataset,resultscolumn(1)));
                if(Optimize(2)==2)
                mdlfb = NaiveBayes.fit(DataSheet(datasetFB,Newfeaturesl2),DataSheet(datasetFB,resultscolumn(2)));
                elseif (Optimize(2)==3)
                mdllr = fitcnb(DataSheet(datasetLR,featuresl3),DataSheet(datasetLR,resultscolumn(3)));
                end
        end % end switch 
        
        predictionl1 = predict(mdl,DataSheet(testset,Newfeaturesl1));
        correctlv1 = sum(DataSheet(testset,resultscolumn(1))==predictionl1);
        
        for pred = 1:length(predictionl1)
            if (predictionl1(pred) ==1 && Optimize(2)==2)
                pred1 = predict(mdlfb,DataSheet(testset(pred),Newfeaturesl2));
                if(pred1 == DataSheet(testset(pred),resultscolumn(2)))
                    correctlv2 = correctlv2 +1;
                end
            elseif (predictionl1(pred) ~= 1 && Optimize(2)==3)
                pred2 = predict(mdllr,DataSheet(testset(pred),Newfeaturesl3));
                if(pred2 == DataSheet(testset(pred),resultscolumn(3)))
                    correctlv3 = correctlv3 +1;
                    
                end
            end
        end% end for each pred
        
        switch(Optimize(2))
        case 1
            ResultsFeatures(fold,ind) = correctlv1;
        case 2
           ResultsFeatures(fold,ind) = correctlv2;
           correctlv2 =0;
        case 3
            ResultsFeatures(fold,ind) = correctlv3;
            correctlv3 =0;
         end% end switch case optimize 
    end %% end for folds
     ResultsFeatures(24,ind)=Feature;
    end  %% end of the features 
  
   for ind = 1:length(featuersleft)
        ResultsFeatures(22,ind)=sum(ResultsFeatures(1:21,ind));
        ResultsFeatures(23,ind)=mean(ResultsFeatures(1:21,ind));
       
   end
  save('ResultsFeatures.mat','ResultsFeatures');
end%% end function 
