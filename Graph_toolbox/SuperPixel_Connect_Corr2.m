%% load data
function [nodes_to_retain,distance,connect_representative_pixel,distance2] = SuperPixel_Connect_Corr2(HyperCube, superpix_img , param2)

sp_num=max(superpix_img);
% compute pairwise distances
distance = zeros(sp_num);
connect_representative_pixel = zeros(sp_num);
% connect_representative_pixel2 = zeros(sp_num);

%% method1:median/mean
[Height,Width,Bands] = size(HyperCube);
ALL_data = reshape(HyperCube, Width*Height, Bands); %3D±ä 2D

Superpixel_Mean_Features=[];
represent_coord_vec=[];
for i=1:sp_num
    sp_idx1= find(superpix_img==i);
    sp_data_i=ALL_data(sp_idx1,:);
    %     X1=mean(sp_data_i);
    X1 = median(sp_data_i);
    Superpixel_Mean_Features=[Superpixel_Mean_Features;X1];
    
    sp_data_i = sp_data_i-repmat(X1,size(sp_data_i,1),1);
    diffs = sqrt(sum(sp_data_i.^2,2));
    usedidx_tmp = find(diffs==min(diffs));
    representative_pixel = sp_idx1(usedidx_tmp(1));
    
    represent_coord_vec=[represent_coord_vec;representative_pixel];
end

for i=1:sp_num
    for j=1:i-1
        
        connect_representative_pixel(i,j) = represent_coord_vec(i);
        connect_representative_pixel(j,i) = represent_coord_vec(j);
        
    end
end

% ALL_data_n=Superpixel_Mean_Features-repmat(mean(Superpixel_Mean_Features,2),[1 size(Superpixel_Mean_Features,2)]);
% X1=abs(sqrt(sum(ALL_data_n.^2,2)));
% if ~isempty(find(X1==0))
%     X1(X1==0)=1;
% end
% distance= abs(ALL_data_n*ALL_data_n'./(X1*X1'));

[distance] = Calculate_Similarity_basic(Superpixel_Mean_Features, param2);


N=sp_num;

% Number of nearest neighbors
knn_param = 10;

% Calculating distances of k-nearest neighbors
knn_distance = zeros(N,1);
nn_distance = zeros(N,1);
for i = 1:N
    % sort all possible neighbors according to distance
    temp = sort(distance(i,:), 'ascend');
    % select k-th neighbor: knn_param+1, as the node itself is considered
    knn_distance(i) = temp(knn_param + 1);
    nn_distance(i) = temp(1 + 1);
end

% sparsification matrix
nodes_to_retain = true(N);
nodes_to_nn = true(N);
for i = 1:N
    nodes_to_retain(i, distance(i,:) > knn_distance(i) ) = false;
    nodes_to_retain(i,i) = false; % diagonal should be zero
    nodes_to_nn(i, distance(i,:) > nn_distance(i) ) = false;
    nodes_to_nn(i,i) = false; % diagonal should be zero
end
nodes_to_retain( nodes_to_retain ~= nodes_to_retain' ) = false;
nodes_to_nn( nodes_to_nn ~= nodes_to_nn' ) = true;

distance2 = zeros(sp_num);
distance2(nodes_to_retain)=distance(nodes_to_retain);
end
