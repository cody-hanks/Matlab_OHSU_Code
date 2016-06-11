%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cody hanks 4/25/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;
close all;
load('DataSheet2.mat');

totalrows = 100;
folds = 10;
datalength = totalrows/folds;
rows = 1:datalength:totalrows;

dataset = [];
% build the Dataset rows (these are all rows except the current
% fold that is the test data ... note skipping the 2 calib rows
for datarows = 1:folds
    %if(datarows < fold || datarows > fold) 
        dataset = [dataset rows(datarows)+2:rows(datarows)+datalength-1];
    %end
end

FeatureTable = DataSheet(dataset,4:24);
[Coeff] = pca(FeatureTable');

DataTable=[DataSheet(dataset,29) DataSheet(dataset,3) Coeff];
