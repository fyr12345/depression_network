clc
clear
resmaple=table2array(t3);%%Drag and drop loading

load(['D:\pls\','Schaefer201810003mmexpression.mat'],'X');

Y=T(:,3);%%Drag and drop loading





resmaple=resmaple'

nanRows = any(isnan(X),2); 

X(nanRows,:) = []; Y(nanRows,:) = [];
resmaple(nanRows,:) = [];
X=zscore(X);Y=zscore(Y);
resmaple=zscore(resmaple)


rho = corr(X,Y,'type','spearman');

%%%%%%%%%%%%%%Spatial autocorrelation resampling

rep=10000

tem=zeros(10000, 15630);
dim=15630



delete(gcp('nocreate'));

parpool(24); 

parfor j=1:rep
    j
    Yp=resmaple(:,j);

    rho = corr(X,Yp,'type','spearman');


    tem(j,:)=rho;
   
end



rho = corr(X,Y,'type','spearman');

for k=1:dim
    p_singlezheng(k)=length(find(tem(:,k)>=rho(k)))/10000;
    p_singlefu(k)=length(find(tem(:,k)<=rho(k)))/10000;
end




%%%%%%%%%%%
%%%%%%%%%%%%
%%%%%%%%%

filename = 'C:\Users\Administrator\Desktop\expression.csv';
K = fopen(filename);
data = textscan(K , '%s', 15630, 'Delimiter', ',');

fclose(K);


geneNames = data{1};
genes  = geneNames;
geneindex=1:15630;


result = cell(dim,3); 


for i = 1:dim
    
    geneName = geneNames{i};
    
    rho = corr(X(:,i),Y,'type','spearman');
    
    if rho > 0 
        p = p_singlezheng(i); 
    else 
        p = p_singlefu(i); 
    end

    result{i,1} = geneName; 
    result{i,2} = rho; 
    result{i,3} = p; 
end

xlswrite('H:\depression_repair\sanzu\zhu\jiyin\t3_jyxg_nmae_p_r.xlsx',result)





