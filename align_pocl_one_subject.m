
fs = 1/(time(2) - time(1));
% t_lc = datenum(abstime)/86400 + [0:(size(data,1)-1)]/fs_lc;
T_LC = 0:(1/fs):((size(data,1)-1)/fs);  
t_lc = datenum(abstime) + T_LC/86400;

N_labels1 = length(info1.Labels);
t_edf1 = cell(1,N_labels1);
for i = 1:N_labels1
    t_edf1{i} = t01{i} + (0:(1/edf_FS1{i}):((length(sig1{i})-1)/edf_FS1{i}))/86400;  
%     t_edf1_segment_start{i} = find(floor(t_edf1{i}*40) == floor(t_lc(end/2)*40));
%     t_edf1_segment_end{i} = find(floor(40*t_edf1{i}) == floor(40*t_lc(floor(5*length(t_lc)/6))));
    sig_align_scale1{i} = interp1(t_edf1{i},sig1{i},t_lc);
end


% t_edf2 = cell(1,N_labels2);
% N_labels2 = length(info2.Labels);
% for i = 1:N_labels2
%     t_edf2{i} = t02{i} + (0:(1/edf_FS2{i}):((length(sig2{i})-1)/edf_FS2{i}))/86400;
%     
% %     t_edf2_segment_start{i} = find(floor(t_edf2{i}*40) == floor(t_lc(end/2)*40));
% %     t_edf2_segment_end{i} = find(floor(40*t_edf2{i}) == floor(40*t_lc(floor(5*end/6))));
%     sig_align_scale2{i} = interp1(t_edf2{i},sig2{i},t_lc);
% end