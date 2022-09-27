%% load data
function [distance] = Calculate_Similarity(HyperCube, param2)

% compute pairwise distances
[Height,Width,Bands] = size(HyperCube);
ALL_data = reshape(HyperCube, Width*Height, Bands); 
[N,num_feature] = size(ALL_data); 

tic
distance = zeros(N, N); % norm matrix
if strcmp(param2,'DIST')
    X2=repmat(sum(ALL_data.^2,2),[1 size(ALL_data,1)]);
    R2=ALL_data*ALL_data';
%     distance= abs(sqrt(X2+X2'-2*R2));
%     WeightMat(nodes_to_retain) = exp( -distance(nodes_to_retain).^2 / (2*sigma^2) );
    distance= X2+X2'-2*R2;
    
elseif strcmp(param2,'DIST2')               %先去均值处理，再计算距离
    ALL_data_n=ALL_data-repmat(mean(ALL_data,2),[1 size(ALL_data,2)]);
    X2=repmat(sum(ALL_data_n.^2,2),[1 size(ALL_data_n,1)]);
    R2=ALL_data_n*ALL_data_n';
%     distance= abs(sqrt(X2+X2'-2*R2));
%     WeightMat(nodes_to_retain) = exp( -distance(nodes_to_retain).^2 / (2*sigma^2) );
    distance= X2+X2'-2*R2;
    
elseif strcmp(param2,'DIST3')
%     ALL_data=ALL_data-repmat(mean(ALL_data,1),[size(ALL_data,1) 1]);
    % compute pairwise distances
    distance = zeros(N);
    % Assign lower triangular part only in loop, saves time
    for i = 1:N
        for j = 1:i-1
            distance(i,j) = sum(abs(ALL_data(i,:) - ALL_data(j,:)));
        end
    end
    % Complete upper triangular part
    distance = distance + distance';
    distance = distance.^2;
    
elseif strcmp(param2,'CORR')  %Feature Mean-Residual Normalized Correlation
    ALL_data_n=ALL_data-repmat(mean(ALL_data,2),[1 size(ALL_data,2)]);
    X1=abs(sqrt(sum(ALL_data_n.^2,2)));
    distance= abs(ALL_data_n*ALL_data_n'./(X1*X1'));
    
% elseif strcmp(param2,'CORR2')
%     distance= abs(corr(ALL_data'));
elseif strcmp(param2,'CORR2')
    distance= abs(ALL_data*ALL_data');
    
elseif strcmp(param2,'CORR3')
    X1=abs(sqrt(sum(ALL_data.^2,2)));
    distance= abs(ALL_data*ALL_data'./(X1*X1'));
    
elseif strcmp(param2,'CORR4')  %Sample Mean-Residual Normalized Correlation
    ALL_data_n=ALL_data-repmat(mean(ALL_data,1),[size(ALL_data,1) 1]);
    X1=abs(sqrt(sum(ALL_data_n.^2,2)));
    distance= abs(ALL_data_n*ALL_data_n'./(X1*X1'));
    
elseif strcmp(param2,'COVA')
%     distance= abs(cov(ALL_data'));
    ALL_data_n=ALL_data-repmat(mean(ALL_data,1),[size(ALL_data,1) 1]);
%     X1=abs(sqrt(sum(ALL_data_n.^2,2)));
    distance= abs(ALL_data_n*ALL_data_n'/(num_feature-1));
    
elseif strcmp(param2,'COVA2')
    distance= abs(cov(ALL_data'));
    
end
