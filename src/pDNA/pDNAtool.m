clear all
close all

prompt = 'Please put the first directory:';
str1 = input(prompt,'s');
prompt = 'Please put the second directory:';
str2 = input(prompt,'s');

prompt = 'Pleaee choose your lambda:';
lambda = input(prompt);
% load data

% read data
x1 = readtable(str1, 'Delimiter', '\t', 'FileType','text');
x2 = readtable(str2,'Delimiter', '\t', 'FileType','text');

% get class names
class1= extractBetween(str1, 'out_', '.tsv');
class2= extractBetween(str2, 'out_', '.tsv');

% extract genes list
geneNames = x1{:,1};

% exclude gene column from data
x1 = x1(:, 2:end);
x2= x2(:, 2:end);

% create algorithm input in struct format 
x1=table2array(x1);
x2=table2array(x2);
OV=struct('X', {{x1',x2'}}, 'Gene', {geneNames},'class',{{class1,class2}});


% compute sigma and sigma_svd
[OV.Sigma, OV.Sigma_svd] = Sigma_compute(OV.X);

% estimate the differential networks using pDNA (with lambda = 0.45)
[Delta_hat, V_hat] = pDNA(OV.Sigma, lambda, 'Sigma_svd', OV.Sigma_svd);

% non-zero values for genes interaction
geneNames=string(OV.Gene);
[rows, cols] = find(Delta_hat{1});
non_zero_values = Delta_hat{1} (Delta_hat{1} ~= 0);


condition = repmat(OV.class{1}, size(non_zero_values));
condition(non_zero_values > 0) = OV.class{2};

result = [ geneNames(cols), geneNames(rows), string(condition), non_zero_values];

%%% Save the result file
T = array2table(result, 'VariableNames', {'target','regulator','condition','weight'});
writetable(T, 'network.tsv', 'FileType', 'text', 'Delimiter', '\t');



%%
p = size(Delta_hat{1},1);
K = length(Delta_hat);

% compue the weighted differential network
W = zeros(p,p);
for k = 1:K
    W = W + double(Delta_hat{k}~=0);
end


% compute the degrees of nodes
Degree = sum(W,2);

[~,ID] = sort(Degree,'descend');

% show the weighted differential network accoring to the degrees of nodes
W = W(ID, ID);
imagesc(W)
Gene_sorted = OV.Gene(ID);
