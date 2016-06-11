clear
close all
clc;
load('DataSheet3.mat');
%delete(gcp('nocreate'))
features = 4:24;



filename = 'Lv3_WCenter.txt';
mail = 'cody';
password = 'pk#5QBD0';

setpref('Internet','E_mail','cody@byterule.com');
setpref('Internet','SMTP_Server','smtp.byterule.com');

%sendmail('cody.hanks2@gmail.com','Script running started','started');

tic
% parfor Feat=features
% Class_Feat(DataSheet,features, [ Feat ], [ ], [ ],'knn',[ 0 1 ],filename);
% end
% sendmail('cody.hanks2@gmail.com','Script running-KNN Done','KNN Done');
% parfor Feat=features
% Class_Feat(DataSheet,features, [ Feat ], [ ], [ ],'bay',[ 0 1 ],filename);
% end
%sendmail('cody.hanks2@gmail.com','Script running-BAY Done','Bay Done');
%65-22,73-9,74-20,75-16

    Greedy_pca(DataSheet,features, [4:loop],[ ], [ ],'svm',[ 0 1]);

%sendmail('cody.hanks2@gmail.com','Script running-SVM Done','SVM Done');
toc
load('ResultsFeatures.mat');
ResultsFeatures