%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cody Hanks 
% 3/7/2016 
% classifier
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clear all
close all
clc;
load('DataSheet.mat');





totalrows = 130;
folds = 13;
datalength = totalrows/folds;
rows = 1:datalength:totalrows;
Resultscolumn = 28;

% FB data set === 13 folds and 10 rows per subject
FBIND =[ ];
for rowid = 1:10:120
  FBIND = [FBIND rowid];
  FBIND = [FBIND rowid+1];
  FBIND = [FBIND rowid+2];
  FBIND = [FBIND rowid+5];
  FBIND = [FBIND rowid+6];
  FBIND = [FBIND rowid+9];
end 
LRIND =[];
for rowid = 1:10:120
  LRIND = [LRIND rowid+ 3];
  LRIND = [LRIND rowid+ 4];
  LRIND = [LRIND rowid+ 7];
  LRIND = [LRIND rowid+ 8];
end

features = 4:23;
% 
% %permfeatbay = [5 9 7 12];
% % for all positions features 5 and 10 gave 43 correct 
% % for front back Vs sides 5,9,7,12 gave 110 correct
 permfeatbay = [ 5 23 15 18 4 8 ];
% front =0;
% back =0;
% left=0;
% right=0;
% bestfeat = 0;
% bestfeatcnt = 0;
% %for feat = features
% %%%%% using folds and greedy approach find the least # of metrics for front
% %%%%% back
% %      if sum(find(permfeatbay==feat))>0
% %          continue
% %      end
%     fb=0;
%     lr=0;
%     front =0;
%     back =0;
%     left=0;
%     right=0;
%     correct =0;
%     correct2 = 0;
%     for fold = 1:folds 
%         subjcor = 0;
%         dataset = []; 
%         if (rows(fold) >1) 
%             dataset = [dataset 1:rows(fold)-1];
%         end
%         if ((rows(fold)+(datalength-1)) <totalrows)
%             dataset = [dataset rows(fold)+(datalength):totalrows];
%         end
%         datasetFB = dataset(FBIND);
%         datasetLR = dataset(LRIND);
%         testset = rows(fold):rows(fold)+(datalength-1);
%         
%         nb = fitcnb(DataSheet(dataset,[permfeatbay]),DataSheet(dataset,Resultscolumn));
%         pred = predict(nb,DataSheet(testset,[permfeatbay]));
%         
%         tempdsfb = DataSheet(datasetFB,1:end);
%         tempdslr = DataSheet(datasetLR,1:end);
%         
%         nbfb = fitcnb(DataSheet(datasetFB,[permfeatbay]),DataSheet(datasetFB,3));
%         nblr = fitcnb(DataSheet(datasetLR,[permfeatbay]),DataSheet(datasetLR,3));
%         
% 
%         
%         disp(strcat('Fold ',num2str(fold)))
%         for predind = 1:datalength
%             if pred(predind) ==1
%                 if pred(predind) == DataSheet(testset(predind),Resultscolumn)
%                     subjcor = subjcor +1;
%                     correct = correct +1;
%                     fb = fb +1;
%                 end 
%                 predictnext = predict(nbfb,DataSheet(testset(predind),[permfeatbay]));
%                 if predictnext == DataSheet(testset(predind),3)
%                     if(DataSheet(testset(predind),3)==1)
%                            back = back +1;
%                     else 
%                         front = front +1;
%                     end
%                     disp(strcat('+Pred 1 FB, ACT-',num2str(DataSheet(testset(predind),Resultscolumn)),' Pred 2 FB:', ...
%                     num2str(predictnext),' ACT-',num2str(DataSheet(testset(predind),3))))
%                     correct2 = correct2+1;
%                 else
%                     disp(strcat('Pred 1 FB, ACT-',num2str(DataSheet(testset(predind),Resultscolumn)),' Pred 2 FB:', ...
%                     num2str(predictnext),' ACT-',num2str(DataSheet(testset(predind),3))))
%                 end
%             else
%                 if pred(predind) == DataSheet(testset(predind),Resultscolumn)
%                     subjcor = subjcor +1;
%                     correct = correct +1;
%                     lr = lr +1;
%                 end 
%                 predictnext = predict(nblr,DataSheet(testset(predind),[permfeatbay]));
%                 if predictnext == DataSheet(testset(predind),3)
%                     if(DataSheet(testset(predind),3)==2)
%                            left = left +1;
%                     else 
%                         right = right +1;
%                     end
%                     disp(strcat('+Pred 1 LR, ACT-',num2str(DataSheet(testset(predind),Resultscolumn)),' Pred 2 LR:', ...
%                     num2str(predictnext),' ACT-',num2str(DataSheet(testset(predind),3))))
%                     correct2 = correct2+1;
%                 else
%                     disp(strcat('Pred 1 LR, ACT-',num2str(DataSheet(testset(predind),Resultscolumn)),' Pred 2 LR:', ...
%                     num2str(predictnext),' ACT-',num2str(DataSheet(testset(predind),3))))
%                 end
%             end
%             
%             
%         end
%         
%         if subjcor ==0
%             disp(strcat('Bad Subj Data',num2str(DataSheet(testset(1,1),1))))
%         end
%     end
% %     if(correct2 > bestfeatcnt)
% %         bestfeat = feat;
% %         bestfeatcnt = correct2;
% %     end
%     disp(strcat('correct ',num2str(correct),' correct 2-',num2str(correct2)))
%     disp(strcat(' front ',num2str(front),' left ',num2str(left) ...
%                 ,' right ',num2str(right),' back ',num2str(back)))
%     disp(strcat('FB ',num2str(fb),' LR ',num2str(lr)))
%             
% %end 
% %    disp(strcat('Bestfeat ',num2str(bestfeat),' feat correct = ',num2str(bestfeatcnt)))
%     
%     


% disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
% front =0;
% back =0;
% left=0;
% right=0;
%  bestfeat = 0;
%  bestfeatcnt = 0;
% %for all positions featurs 5,8,11,12,4 gave 43 correct 
% % for front back Vs sides 4 5 13 gave 104 correct 
% %permfeatknn = [4 5 13];
 permfeatknn = [ 5 4 9 15];
% %for feat = features
%     fb=0;
%     lr=0;
%     front =0;
%     back =0;
%     left=0;
%     right=0;
%     correct =0;
%     correct2 = 0;
% %%%%% using folds and greedy approach find the least # of metrics for front
% %%%%% back
% %      if sum(find(permfeatknn==feat))>0
% %          continue
% %      end
%     correct = 0;
%     correct2 = 0;
%     for fold = 1:folds 
%         subjcor = 0;
%         dataset = []; 
%         if (rows(fold) >1) 
%             dataset = [dataset 1:rows(fold)-1];
%         end
%         if ((rows(fold)+(datalength-1)) <totalrows)
%             dataset = [dataset rows(fold)+(datalength):totalrows];
%         end
%         datasetFB = dataset(FBIND);
%         datasetLR = dataset(LRIND);
%         testset = rows(fold):rows(fold)+(datalength-1);
% 
%         knn = fitcknn(DataSheet(dataset,[permfeatknn]),DataSheet(dataset,Resultscolumn));
%         pred = predict(knn,DataSheet(testset,[permfeatknn]));
%         
%         
%         knnfb = fitcknn(DataSheet(datasetFB,[permfeatknn]),DataSheet(datasetFB,3));
%         knnlr = fitcknn(DataSheet(datasetLR,[permfeatknn]),DataSheet(datasetLR,3));
%        
%         disp(strcat('Fold ',num2str(fold)))
%         for predind = 1:datalength
%             if pred(predind) ==1
%                 if pred(predind) == DataSheet(testset(predind),Resultscolumn)
%                     subjcor = subjcor +1;
%                     correct = correct +1;
%                     fb = fb +1;
%                 end 
%                 predictnext = predict(knnfb,DataSheet(testset(predind),[permfeatknn]));
%                 if predictnext == DataSheet(testset(predind),3)
%                     if(DataSheet(testset(predind),3)==1)
%                            back = back +1;
%                     else 
%                         front = front +1;
%                     end
%                     disp(strcat('+Pred 1 FB, ACT-',num2str(DataSheet(testset(predind),Resultscolumn)),' Pred 2 FB:', ...
%                     num2str(predictnext),' ACT-',num2str(DataSheet(testset(predind),3))))
%                     correct2 = correct2+1;
%                 else
%                     disp(strcat('Pred 1 FB, ACT-',num2str(DataSheet(testset(predind),Resultscolumn)),' Pred 2 FB:', ...
%                     num2str(predictnext),' ACT-',num2str(DataSheet(testset(predind),3))))
%                 end
%             else
%                 if pred(predind) == DataSheet(testset(predind),Resultscolumn)
%                     subjcor = subjcor +1;
%                     correct = correct +1;
%                     lr = lr +1;
%                 end 
%                 predictnext = predict(knnlr,DataSheet(testset(predind),[permfeatknn]));
%                 if predictnext == DataSheet(testset(predind),3)
%                    if(DataSheet(testset(predind),3)==2)
%                         left = left +1;
%                     else 
%                         right = right +1;
%                     end
%                     disp(strcat('+Pred 1 LR, ACT-',num2str(DataSheet(testset(predind),Resultscolumn)),' Pred 2 LR:', ...
%                     num2str(predictnext),' ACT-',num2str(DataSheet(testset(predind),3))))
%                     correct2 = correct2+1;
%                 else
%                     disp(strcat('Pred 1 LR, ACT-',num2str(DataSheet(testset(predind),Resultscolumn)),' Pred 2 LR:', ...
%                     num2str(predictnext),' ACT-',num2str(DataSheet(testset(predind),3))))
%                 end
%             end
% 
%         end
%          if subjcor ==0
%             disp(strcat('Bad Subj Data',num2str(DataSheet(testset(1,1),1))))
%         end
%     end
% %     if(correct2 > bestfeatcnt)
% %          bestfeat = feat;
% %          bestfeatcnt = correct2;
% %     end
%     disp(strcat('correct ',num2str(correct),' correct 2-',num2str(correct2)))
%     disp(strcat(' front ',num2str(front),' left ',num2str(left) ...
%                 ,' right ',num2str(right),' back ',num2str(back)))
%     disp(strcat('FB ',num2str(fb),' LR ',num2str(lr)))
% % end 
% % disp(strcat('Bestfeat ',num2str(bestfeat),' feat correct = ',num2str(bestfeatcnt)))



% disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
% front =0;
% back =0;
% left=0;
% right=0;
%  bestfeat = 0;
%  bestfeatcnt = 0;
% %for all positions featurs 5,8,11,12,4 gave 43 correct 
% % for front back Vs sides 4 5 13 gave 104 correct 
% %permfeatknn = [4 5 13];
 permfeatsvm = [ 23 5 4 6 12 7 10 8 20 13 ];
% % for feat = features
%     fb=0;
%     lr=0;
%     front =0;
%     back =0;
%     left=0;
%     right=0;
%     correct =0;
%     correct2 = 0;
% %%%%% using folds and greedy approach find the least # of metrics for front
% %%%%% back
% %      if sum(find(permfeatsvm==feat))>0
% %          continue
% %      end
%     correct = 0;
%     correct2 = 0;
%     for fold = 1:folds 
%         subjcor = 0;
%         dataset = []; 
%         if (rows(fold) >1) 
%             dataset = [dataset 1:rows(fold)-1];
%         end
%         if ((rows(fold)+(datalength-1)) <totalrows)
%             dataset = [dataset rows(fold)+(datalength):totalrows];
%         end
%         datasetFB = dataset(FBIND);
%         datasetLR = dataset(LRIND);
%         testset = rows(fold):rows(fold)+(datalength-1);
% 
%         svm = fitcsvm(DataSheet(dataset,[ permfeatsvm]),DataSheet(dataset,Resultscolumn));
%         pred = predict(svm,DataSheet(testset,[ permfeatsvm]));
%         
%         
%         svmfb = fitcsvm(DataSheet(datasetFB,[ permfeatsvm]),DataSheet(datasetFB,3));
%         svmlr = fitcsvm(DataSheet(datasetLR,[ permfeatsvm]),DataSheet(datasetLR,3));
%        
%         disp(strcat('Fold ',num2str(fold)))
%         for predind = 1:datalength
%             if pred(predind) ==1
%                 if pred(predind) == DataSheet(testset(predind),Resultscolumn)
%                     subjcor = subjcor +1;
%                     correct = correct +1;
%                     fb = fb +1;
%                 end 
%                 predictnext = predict(svmfb,DataSheet(testset(predind),[ permfeatsvm]));
%                 if predictnext == DataSheet(testset(predind),3)
%                     if(DataSheet(testset(predind),3)==1)
%                            back = back +1;
%                     else 
%                         front = front +1;
%                     end
%                     disp(strcat('+Pred 1 FB, ACT-',num2str(DataSheet(testset(predind),Resultscolumn)),' Pred 2 FB:', ...
%                     num2str(predictnext),' ACT-',num2str(DataSheet(testset(predind),3))))
%                     correct2 = correct2+1;
%                 else
%                     disp(strcat('Pred 1 FB, ACT-',num2str(DataSheet(testset(predind),Resultscolumn)),' Pred 2 FB:', ...
%                     num2str(predictnext),' ACT-',num2str(DataSheet(testset(predind),3))))
%                 end
%             else
%                 if pred(predind) == DataSheet(testset(predind),Resultscolumn)
%                     subjcor = subjcor +1;
%                     correct = correct +1;
%                     lr = lr +1;
%                 end 
%                 predictnext = predict(svmlr,DataSheet(testset(predind),[ permfeatsvm]));
%                 if predictnext == DataSheet(testset(predind),3)
%                    if(DataSheet(testset(predind),3)==2)
%                         left = left +1;
%                     else 
%                         right = right +1;
%                     end
%                     disp(strcat('+Pred 1 LR, ACT-',num2str(DataSheet(testset(predind),Resultscolumn)),' Pred 2 LR:', ...
%                     num2str(predictnext),' ACT-',num2str(DataSheet(testset(predind),3))))
%                     correct2 = correct2+1;
%                 else
%                     disp(strcat('Pred 1 LR, ACT-',num2str(DataSheet(testset(predind),Resultscolumn)),' Pred 2 LR:', ...
%                     num2str(predictnext),' ACT-',num2str(DataSheet(testset(predind),3))))
%                 end
%             end
% 
%         end
%          if subjcor ==0
%             disp(strcat('Bad Subj Data',num2str(DataSheet(testset(1,1),1))))
%         end
%     end
% %     if(correct2 > bestfeatcnt)
% %          bestfeat = feat;
% %          bestfeatcnt = correct2;
% %     end
%     disp(strcat('correct ',num2str(correct),' correct 2-',num2str(correct2)))
%     disp(strcat(' front ',num2str(front),' left ',num2str(left) ...
%                 ,' right ',num2str(right),' back ',num2str(back)))
%     disp(strcat('FB ',num2str(fb),' LR ',num2str(lr)))
% % end 
% % disp(strcat('Bestfeat ',num2str(bestfeat),' feat correct = ',num2str(bestfeatcnt)))


disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
front =0;
back =0;
left=0;
right=0;
 bestfeat = 0;
 bestfeatcnt = 0;
%for all positions featurs 5,8,11,12,4 gave 43 correct 
% for front back Vs sides 4 5 13 gave 104 correct 
%for feat = features
    fb=0;
    lr=0;
    front =0;
    back =0;
    left=0;
    right=0;
    correct =0;
    correct2 = 0;
%%%%% using folds and greedy approach find the least # of metrics for front
%%%%% back
%      if sum(find(permfeatsvm==feat))>0
%          continue
%      end
    correct = 0;
    correct2 = 0;
    parfor fold = 1:folds 
        subjcor = 0;
        dataset = []; 
        if (rows(fold) >1) 
            dataset = [dataset 1:rows(fold)-1];
        end
        if ((rows(fold)+(datalength-1)) <totalrows)
            dataset = [dataset rows(fold)+(datalength):totalrows];
        end
        datasetFB = dataset(FBIND);
        datasetLR = dataset(LRIND);
        testset = rows(fold):rows(fold)+(datalength-1);

        nb = fitcnb(DataSheet(dataset,[ permfeatbay]),DataSheet(dataset,Resultscolumn));
        pred = predict(nb,DataSheet(testset,[ permfeatbay]));
        
        
        knnfb = fitcknn(DataSheet(datasetFB,[ permfeatknn]),DataSheet(datasetFB,3));
        svmlr = fitcsvm(DataSheet(datasetLR,[ permfeatsvm]),DataSheet(datasetLR,3));
       
        disp(strcat('Fold ',num2str(fold)))
        for predind = 1:datalength
            if pred(predind) ==1
                if pred(predind) == DataSheet(testset(predind),Resultscolumn)
                    subjcor = subjcor +1;
                    correct = correct +1;
                    fb = fb +1;
                end 
                predictnext = predict(knnfb,DataSheet(testset(predind),[ permfeatknn]));
                if predictnext == DataSheet(testset(predind),3)
                    if(DataSheet(testset(predind),3)==1)
                           back = back +1;
                    else 
                        front = front +1;
                    end
                    disp(strcat('+Pred 1 FB, ACT-',num2str(DataSheet(testset(predind),Resultscolumn)),' Pred 2 FB:', ...
                    num2str(predictnext),' ACT-',num2str(DataSheet(testset(predind),3))))
                    correct2 = correct2+1;
                else
                    disp(strcat('Pred 1 FB, ACT-',num2str(DataSheet(testset(predind),Resultscolumn)),' Pred 2 FB:', ...
                    num2str(predictnext),' ACT-',num2str(DataSheet(testset(predind),3))))
                end
            else
                if pred(predind) == DataSheet(testset(predind),Resultscolumn)
                    subjcor = subjcor +1;
                    correct = correct +1;
                    lr = lr +1;
                end 
                predictnext = predict(svmlr,DataSheet(testset(predind),[ permfeatsvm]));
                if predictnext == DataSheet(testset(predind),3)
                   if(DataSheet(testset(predind),3)==2)
                        left = left +1;
                    else 
                        right = right +1;
                    end
                    disp(strcat('+Pred 1 LR, ACT-',num2str(DataSheet(testset(predind),Resultscolumn)),' Pred 2 LR:', ...
                    num2str(predictnext),' ACT-',num2str(DataSheet(testset(predind),3))))
                    correct2 = correct2+1;
                else
                    disp(strcat('Pred 1 LR, ACT-',num2str(DataSheet(testset(predind),Resultscolumn)),' Pred 2 LR:', ...
                    num2str(predictnext),' ACT-',num2str(DataSheet(testset(predind),3))))
                end
            end

        end
         if subjcor ==0
            disp(strcat('Bad Subj Data',num2str(DataSheet(testset(1,1),1))))
        end
    end
%     if(correct2 > bestfeatcnt)
%          bestfeat = feat;
%          bestfeatcnt = correct2;
%     end
    disp(strcat('correct ',num2str(correct),' correct 2-',num2str(correct2)))
    disp(strcat(' front ',num2str(front),' left ',num2str(left) ...
                ,' right ',num2str(right),' back ',num2str(back)))
    disp(strcat('FB ',num2str(fb),' LR ',num2str(lr)))
% end 
% disp(strcat('Bestfeat ',num2str(bestfeat),' feat correct = ',num2str(bestfeatcnt)))

