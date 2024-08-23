clc
clear all
basepath = 'H:/depression_repair/sanzu/zhu/';
outpath = 'H:/depression_repair/sanzu/zhu/outputs/';
Nsub=192;
Nroi = 1001;
Nstep = 7;




%% 1) Binarize connectivity matrix with threshold 10
th = 10;

for sidx = 1 : Nsub
    disp(['subject = ', num2str(sidx)])
    
    
    csv_file = [basepath, '1.connqugsi/sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.csv'];
    
    
    conn = csvread(csv_file);
    
   
    binconn = binarize_conn(conn,th);
    
    
    save([outpath, '1.bincoon/sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat'], 'binconn');
end


%% 2) Construct SFC matrix
Nstep = 7;
for sidx = 1 : Nsub
    disp(['subject = ', num2str(sidx)])
   
    load([outpath, '1.bincoon/sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat'])
    sfc = compute_sfc(binconn, Nstep)
    save([outpath, '2.sfc/sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat'], 'sfc');
end


%% 3) Save DC values per ROI
seed_idx=1001
load([outpath, 'a_group.mat'])
Nstep=7
Nsub=192
roi_dc = zeros(Nsub, Nroi, Nstep);

for sidx = 1 : Nsub
    load([outpath, '2.sfc/sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat'])
    for step = 1 : Nstep
        dc = sum(sfc(seed_idx,:,step), 1);
        dc(isinf(dc)|isnan(dc)) = 0;
        roi_dc(sidx, :, step) = dc;
    end
end
save([outpath, 'wholesub_ROI'], 'roi_dc');


%% 4) Compute group average SFC
seed_idx=1001;
grpmean_SFC = cell(2,1);
grpmean_SFC{1} = zeros(length(seed_idx), Nroi, Nstep);
grpmean_SFC{2} = zeros(length(seed_idx), Nroi, Nstep);



for sidx = 1 : Nsub

    disp(['subject = ', num2str(sidx)])
    if group(sidx) == 0
        load([outpath, '2.sfc/sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat']);
        grpmean_SFC{1} = grpmean_SFC{1} + sfc(1001,:,:);
    elseif group(sidx) == 1
        load([outpath, '2.sfc/sub', pad(num2str(sidx, '%d'), 3, 'left', '0'), '.mat']);
        grpmean_SFC{2} = grpmean_SFC{2} + sfc(1001,:,:);
    end
end
grpmean_SFC{1} = grpmean_SFC{1} / sum(group == 0);
grpmean_SFC{2} = grpmean_SFC{2} / sum(group == 1);

grpmean_DC = zeros(2, Nroi, Nstep);
grpmean_DC(1,:,:) = squeeze(sum(grpmean_SFC{1}(:,:,1:Nstep), 1));
grpmean_DC(2,:,:) = squeeze(sum(grpmean_SFC{2}(:,:,1:Nstep), 1));
save([outpath, 'groupmean_SFC.mat'], 'grpmean_SFC','grpmean_DC')

%% 5) Group difference test: roi-level
load([outpath, 'wholesub_ROI.mat']);
Nroi = 1000;
H = zeros(Nroi, Nstep);
P = zeros(Nroi, Nstep);
T = zeros(Nroi, Nstep);

for step = 1 : Nstep
    %roi_dc(:,:,step)=(roi_dc(:,:,step) - mean(roi_dc(:,:,step))/std(roi_dc(:,:,step)));
    jixian = [roi_dc(group==1,1:1000,step)];
    suifang = [roi_dc(group==0,1:1000,step)];
    for roi = 1 : Nroi
        [~,p,~,stats] = ttest(suifang(:,roi), jixian(:,roi));
        P(roi, step) = p;
        T(roi, step) = stats.tstat;
    end
   
    
    corrected_P = mafdr(P(:, step), 'BHFDR', true);
    
    
    H(:, step) = corrected_P < 0.05;
    
    
    HT(:, step) = H(:, step) .* T(:, step);
end
H=logical(P<0.05)
HT=H.*T




%% 6£©Extracting the pairwise difference

load([outpath, 'wholesub_ROI.mat']);
Nroi = 1000;

diff_values_1=zeros(23,Nroi)
diff_values_2=zeros(23,Nroi)
diff_values_3=zeros(23,Nroi)
diff_values_4=zeros(23,Nroi)



step=1
jixian = [roi_dc(group==1,1:1000,step)];
suifang = [roi_dc(group==0,1:1000,step)];
for roi = 1 : Nroi
        
    diff_values_1(:,roi) = jixian(:,roi) - suifang(:,roi);
end
save([outpath, 'diff_values_1.mat'], 'diff_values_1'); 


step=2
jixian = [roi_dc(group==1,1:1000,step)];
suifang = [roi_dc(group==0,1:1000,step)];
for roi = 1 : Nroi
        
    diff_values_2(:,roi) = jixian(:, roi) - suifang(:, roi);
end
save([outpath, 'diff_values_2.mat'], 'diff_values_2'); 


step=3
jixian = [roi_dc(group==1,1:1000,step)];
suifang = [roi_dc(group==0,1:1000,step)];
for roi = 1 : Nroi
        
    diff_values_3(:,roi) = jixian(:, roi) - suifang(:, roi);
end
save([outpath, 'diff_values_3.mat'], 'diff_values_3'); 



step=4
jixian = [roi_dc(group==1,1:1000,step)];
suifang = [roi_dc(group==0,1:1000,step)];
for roi = 1 : Nroi
     
    diff_values_4(:,roi) = jixian(:, roi) - suifang(:, roi);
end
save([outpath, 'diff_values_4.mat'], 'diff_values_4'); 




%% Select Features



selected_features_1 = [];
selected_features_2 = [];
selected_features_3 = [];
selected_features_4 = [];
Nstep=4
for step = 1:Nstep
    
    significant_indices = find(H(:, step));  
    
    
    if step == 1
        selected_features_1 = diff_values_1(:, significant_indices);
    elseif step == 2
        selected_features_2 = diff_values_2(:, significant_indices);
    elseif step == 3
        selected_features_3 = diff_values_3(:, significant_indices);
    elseif step == 4
        selected_features_4 = diff_values_4(:, significant_indices);
    end
end


