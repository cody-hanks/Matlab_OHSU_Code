% info1 = zinfo('D:\mahesh\Sleep_Project\Legacy_Stuff_Zach\zachs_archive\DATA\BSS Data\BSS_Pete and Zach 5-26-2015\MMP_ID_1_1.edf');
% info2 = zinfo('D:\mahesh\Sleep_Project\Legacy_Stuff_Zach\zachs_archive\DATA\BSS Data\BSS_Pete and Zach 5-26-2015\MMP_ID_1_2.edf');

% subind = 7;

% load(['D:\Projects\zachs_archive\DATA\POCL_6308_Data\POCL_ID_6308_' num2str(subind) '\POCL_ID_6308_' num2str(subind) '_data' '.mat']);
sub1_file = ['/home/cody/Documents/for_cody_11_30/sp/sp/SP_ID_' num2str(subind) '/SP_ID_' num2str(subind) '.edf'];
% sub2_file = ['D:\mahesh\Sleep_Project\Data\Two_Subjects\MMP_ID_' num2str(subind) '\MMP_ID_' num2str(subind) '_2.edf'];
info1 = zinfo(sub1_file);
% info2 = zinfo(sub2_file);


N_labels1 = length(info1.Labels);
for label_index = 1:N_labels1
    label_cell{label_index} = info1.Labels(label_index);
    %         clear label
    %         label = [label_temp repmat([' '],1,16-length(label_temp))];
    [sig1{label_index}, t01{label_index}, edf_FS1{label_index}, Pdim1{label_index}] = zload(sub1_file,find(strcmp(info1.Labels,label_cell{label_index})));
%     [sig1{label_index}, t01{label_index}, edf_FS1{label_index}, Pdim1{label_index}] = zload(sub1_file,find(strcmp(info1.Labels,label_cell{label_index})));
    sig1{label_index} = sig1{label_index} - mean(sig1{label_index});
    
end
load(['/home/cody/Documents/for_cody_11_30/sp/sp/SP_ID_' num2str(subind) '/LC_Data_SP_ID_' num2str(subind) '.mat']);


% N_labels2 = length(info2.Labels);
% for label_index = 1:N_labels2
%     label_cell{label_index} = info2.Labels(label_index);
%     %         clear label
%     %         label = [label_temp repmat([' '],1,16-length(label_temp))];
% %     [sig2{label_index}, t02{label_index}, edf_FS2{label_index}, Pdim2{label_index}] = zload('D:\mahesh\Sleep_Project\Legacy_Stuff_Zach\zachs_archive\DATA\BSS Data\BSS_Pete and Zach 5-26-2015\MMP_ID_1_2.edf',find(strcmp(info2.Labels,label_cell{label_index})));
%     [sig2{label_index}, t02{label_index}, edf_FS2{label_index}, Pdim2{label_index}] = zload(sub2_file,find(strcmp(info2.Labels,label_cell{label_index})));
% end

