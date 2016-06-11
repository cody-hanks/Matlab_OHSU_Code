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

Greedy(DataSheet,features, [ 23,18,12,11,10 ], [ 20,4,10,24,5 ], [ 5,4,21  ],'svm',[ 0 3 ]);

%sendmail('cody.hanks2@gmail.com','Script running-SVM Done','SVM Done');
toc
load('ResultsFeatures.mat');
ResultsFeatures