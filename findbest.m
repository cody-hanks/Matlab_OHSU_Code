
function [svmfeat,svmline,corr3] = findbest(filename)
fid = fopen(filename);
bestresult =[0 0 0];
while ~feof(fid)
    tline = fgetl(fid);
    switch(tline(1:3))
        case 'knn'
            corr = str2num(tline(strfind(tline,';')+1:length(tline)));
            if corr >= bestresult(1)
                knnfeat =[];
                try 
                    featurestr = tline(strfind(tline,'[')+1:strfind(tline,']')-2);
                    splits = strsplit(featurestr,',');
                    for split = splits
                        knnfeat = [knnfeat  str2num(char(split))];
                    end
                    
                    bestresult(1) = corr;
                    corr1 = corr;
                    knnline = tline;
                catch 
                    continue
                end 
            end
        case 'bay'
            corr = str2num(tline(strfind(tline,';')+1:length(tline)));
            if corr >= bestresult(2)
              bayfeat =[];
                try 
                    featurestr = tline(strfind(tline,'[')+1:strfind(tline,']')-2);
                    splits = strsplit(featurestr,',');
                    for split = splits
                        bayfeat = [bayfeat  str2num(char(split))];
                    end
                    bestresult(2) = corr;
                    corr2 = corr;
                    bayline = tline;
                catch 
                    continue
                end 
            end
        case 'svm'
            corr = str2num(tline(strfind(tline,';')+1:length(tline)));
             if corr >= bestresult(3)
                svmfeat =[];
                try 
                    featurestr = tline(strfind(tline,'[')+1:strfind(tline,']')-2);
                    splits = strsplit(featurestr,',');
                    for split = splits
                        svmfeat = [svmfeat  str2num(char(split))];
                    end
                    %if (length(svmfeat) < currlen)
                        %svmfeat = featurs
                        bestresult(3) = corr;
                        corr3 = corr;
                        svmline = tline;
                   % end 
                catch 
                    continue
                end 
            end
    end
end

