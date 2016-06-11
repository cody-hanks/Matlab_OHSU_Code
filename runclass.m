clear all
close all
clc;
load('DataSheet2.mat');
delete(gcp('nocreate'))
features = 4:24;



filename = 'Lv3_WCenter.txt';
mail = 'cody';
password = 'pk#5QBD0';

setpref('Internet','E_mail','cody@byterule.com');
setpref('Internet','SMTP_Server','smtp.byterule.com');

%sendmail('cody.hanks2@gmail.com','Script running started','started');
global best
best = 0;
tic
% parfor Feat=features
% Class_Feat(DataSheet,features, [ Feat ], [ ], [ ],'knn',[ 0 1 ],filename);
% end
% sendmail('cody.hanks2@gmail.com','Script running-KNN Done','KNN Done');
% parfor Feat=features
% Class_Feat(DataSheet,features, [ Feat ], [ ], [ ],'bay',[ 0 1 ],filename);
% end
%sendmail('cody.hanks2@gmail.com','Script running-BAY Done','Bay Done');
parfor Feat=features
Class_Feat(DataSheet,features, [ Feat ], [  ], [  ],'svm',[ 0 1 ],filename,best,[],'');
end
%sendmail('cody.hanks2@gmail.com','Script running-SVM Done','SVM Done');
toc
