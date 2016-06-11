%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cody Hanks 
% 3/29/2016
% classifier function for inclusion into a search algorithm 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Class_Feat(DataSheet,features,featuresl1,featuresl2,featuresl3,Classifier,Optimize,filename,best,feats,classif)
    if Optimize(2)==1
        featuersleft = setdiff(features,featuresl1);
    elseif Optimize(2)==2
        featuersleft = setdiff(features,featuresl2);
    elseif Optimize(2)==3
        featuersleft = setdiff(features,featuresl3);
    end

    box1 = 130;
    kern1 = 57;
    box2 = 139;
    kern2 = 5;
    box3 = 200;
    kern3 = 45;
    
    totalrows = 130;
    folds = 13;
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
                mdl = fitcsvm(DataSheet(dataset,featuresl1),DataSheet(dataset,resultscolumn(1)),'KernelFunction','rbf','Standardize',true,'BoxConstraint',box1,'KernelScale',kern1);
                if(Optimize(2)==2)
                mdlfb = fitcsvm(DataSheet(datasetFB,featuresl2),DataSheet(datasetFB,resultscolumn(2)),'KernelFunction','rbf','Standardize',true,'BoxConstraint',box2,'KernelScale',kern2);
                elseif (Optimize(2)==3)
                mdllr = fitcsvm(DataSheet(datasetLR,featuresl3),DataSheet(datasetLR,resultscolumn(3)),'KernelFunction','rbf','Standardize',true,'BoxConstraint',box3,'KernelScale',kern3);
                end
            case 'bay'
                mdl = fitcnb(DataSheet(dataset,featuresl1),DataSheet(dataset,resultscolumn(1)));
                if(Optimize(2)==2)
                mdlfb = fitcnb(DataSheet(datasetFB,featuresl2),DataSheet(datasetFB,resultscolumn(2)));
                elseif (Optimize(2)==3)
                mdllr = fitcnb(DataSheet(datasetLR,featuresl3),DataSheet(datasetLR,resultscolumn(3)));
                end
        end % end switch 
        
        predictionl1 = predict(mdl,DataSheet(testset,featuresl1));
        correctlv1 =correctlv1 + sum(DataSheet(testset,resultscolumn(1))==predictionl1);
        
        for pred = 1:length(predictionl1)
            if (predictionl1(pred) ==1 && length(featuresl2)>0 && Optimize(2)==2)
                pred1 = predict(mdlfb,DataSheet(testset(pred),featuresl2));
                if(pred1 == DataSheet(testset(pred),resultscolumn(2)))
                    correctlv2 = correctlv2 +1;
                end
            elseif (predictionl1(pred) ~= 1 && length(featuresl3)>0 && Optimize(2)==3)
                pred2 = predict(mdllr,DataSheet(testset(pred),featuresl3));
                if(pred2 == DataSheet(testset(pred),resultscolumn(3)))
                    correctlv3 = correctlv3 +1;
                    
                end
            end
        end% end for each pred
        
        
    end %% end for folds
    
    switch(Optimize(2))
        case 1
            if correctlv1 > Optimize(1)&&~isempty(featuersleft) && correctlv1 > best
                Optimize(1) = correctlv1;
                for feature = featuersleft
                    NewFeaturesl1 = [featuresl1 feature];
                    Class_Feat(DataSheet,features,NewFeaturesl1,featuresl2,featuresl3,Classifier,Optimize,filename,best,feats,classif);
                end
            else
                if correctlv1 > best
                    %best,feats,classif
                    classif = Classifier;
                    feats = featuresl1;
                    best = correctlv1;
                    fid = fopen(filename,'a+');
                    fprintf(fid,Classifier);
                    fprintf(fid,'[');
                    fprintf(fid,'%d,',featuresl1);
                    fprintf(fid,'];%d\n',Optimize(1));
                    fclose(fid);
                end
                
            end
        case 2
            if (correctlv2 > Optimize(1)&&~isempty(featuersleft)) && correctlv2 > best
                Optimize(1) = correctlv2;
                for feature = featuersleft
                    NewFeaturesl2 = [featuresl2 feature];
                    Class_Feat(DataSheet,features,featuresl1,NewFeaturesl2,featuresl3,Classifier,Optimize,filename,best,feats,classif);
                end
            else
                if correctlv2 >best
                classif = Classifier;
                feats = featuresl2;
                best = correctlv2;
                fid = fopen(filename,'a+');
                fprintf(fid,Classifier);
                fprintf(fid,'[');
                fprintf(fid,'%d,',featuresl2);
                fprintf(fid,'];%d\n',Optimize(1));
                fclose(fid);
                end
            end
        case 3
            if correctlv3 > Optimize(1)&&~isempty(featuersleft)&& correctlv3 > best
                Optimize(1) = correctlv3;
                for feature = featuersleft
                    Newfeaturesl3 = [featuresl3 feature];
                    Class_Feat(DataSheet,features,featuresl1,featuresl2,Newfeaturesl3,Classifier,Optimize,filename,best,feats,classif);
                end
            else
                if correctlv3 >best
                classif = Classifier;
                feats = featuresl3;
                best = correctlv3;
                fid = fopen(filename,'a+');
                fprintf(fid,Classifier);
                fprintf(fid,'[');
                fprintf(fid,'%d,',featuresl3);
                fprintf(fid,'];%d\n',Optimize(1));
                fclose(fid);
                end
            end
    end% end switch case optimize 
    
    

end%% end function 