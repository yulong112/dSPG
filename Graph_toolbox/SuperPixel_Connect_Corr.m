
function [nodes_to_retain,distance,connect_representative_pixel] = SuperPixel_Connect_Corr(HyperCube, Corr_Mat, superpix_img)

%
%  Calculating most representative pixels and weights
%
%  see paper:
%
%  Yu, L.; Li, J.; et al. "dSPG: A New Discriminant Superpixel Graph Regularizer for Semi-Supervised Image Classification"
%
% -- Input Parameters ----------------
%
% HyperCube     ->   data cube
% Corr_Mat      ->   precomputed contract weight matrix between superpixels
% superpix_img  ->   superpixel segmentation results
%
% Author: Long Yu, September. 2022


sp_num=max(superpix_img);
% compute pairwise distances
distance = zeros(sp_num);
connect_representative_pixel = zeros(sp_num);
% connect_representative_pixel2 = zeros(sp_num);

%% method1:median/mean
% [Height,Width,Bands] = size(HyperCube);
% ALL_data = reshape(HyperCube, Width*Height, Bands); 
% 
% Superpixel_Mean_Features=[];
% represent_coord_vec=[];
% for i=1:sp_num
%     sp_idx1= find(superpix_img==i);
%     sp_data_i=ALL_data(sp_idx1,:);
% %     X1=mean(sp_data_i);
%     X1 = median(sp_data_i);
%     Superpixel_Mean_Features=[Superpixel_Mean_Features;X1];
% end
% 
% ALL_data_n=Superpixel_Mean_Features-repmat(mean(Superpixel_Mean_Features,2),[1 size(Superpixel_Mean_Features,2)]);
% X1=abs(sqrt(sum(ALL_data_n.^2,2)));
% if ~isempty(find(X1==0))
%     X1(X1==0)=1;
% end
% distance= abs(ALL_data_n*ALL_data_n'./(X1*X1'));
% 
% 
% for i=1:sp_num
%     for j=1:i-1
%         sp_idx1=find(superpix_img==i);
%         sp_idx2=find(superpix_img==j);
%         distance_temp = Corr_Mat(sp_idx1,sp_idx2);
%         median_dist = distance(i,j);
%         res_dist = abs(distance_temp-median_dist);
% 
%         [maxdistant1,maxid1] =min(res_dist);
%         [maxdistant2,maxid2] =min(maxdistant1);
%         representative_pixel1_idx=maxid1(maxid2);
%         representative_pixel2_idx=maxid2;
%         
%         representative_pixel1 = sp_idx1(representative_pixel1_idx);
%         representative_pixel2 = sp_idx2(representative_pixel2_idx);
%         
%         connect_representative_pixel(i,j) = representative_pixel1;
%         connect_representative_pixel(j,i) = representative_pixel2;
%         
%         %             distance(i,j_idx2) = sum(sum(distance_temp,1),2);
%     end
% end


%% method2:max
for i=1:sp_num
    for j=1:i-1
        sp_idx1=find(superpix_img==i);
        sp_idx2=find(superpix_img==j);
        distance_temp = []; % norm matrxi
        distance_temp = Corr_Mat(sp_idx1,sp_idx2);
%         distance(i,j) = mean(mean(distance_temp(distance_temp<0.9),1),2);
%         distance(i,j) = mean(mean(distance_temp,1),2);

        [maxdistant1,maxid1] =max(distance_temp);
        [maxdistant2,maxid2] =max(maxdistant1);
        representative_pixel1_idx=maxid1(maxid2);
        representative_pixel2_idx=maxid2;
        distance(i,j) = maxdistant2;
        
        representative_pixel1 = sp_idx1(representative_pixel1_idx);
        representative_pixel2 = sp_idx2(representative_pixel2_idx);
        
        connect_representative_pixel(i,j) = representative_pixel1;
        connect_representative_pixel(j,i) = representative_pixel2;
        
        %             distance(i,j_idx2) = sum(sum(distance_temp,1),2);
    end
    distance(i,i) = 0;
end
% Complete upper triangular part
distance = distance + distance';

N=sp_num;

% Number of nearest neighbors
knn_param = 10;

% Calculating distances of k-nearest neighbors
knn_distance = zeros(N,1);
nn_distance = zeros(N,1);
for i = 1:N
    % sort all possible neighbors according to distance
    temp = sort(distance(i,:), 'descend');
    % select k-th neighbor: knn_param+1, as the node itself is considered
    knn_distance(i) = temp(knn_param + 1);
    nn_distance(i) = temp(1 + 1);
end

% sparsification matrix
nodes_to_retain = true(N);
nodes_to_nn = true(N);
for i = 1:N
    nodes_to_retain(i, distance(i,:) < knn_distance(i) ) = false;
    nodes_to_retain(i,i) = false; % diagonal should be zero
    nodes_to_nn(i, distance(i,:) < nn_distance(i) ) = false;
    nodes_to_nn(i,i) = false; % diagonal should be zero
end
nodes_to_retain( nodes_to_retain ~= nodes_to_retain' ) = false;
nodes_to_nn( nodes_to_nn ~= nodes_to_nn' ) = true;

distance2 = zeros(sp_num);
distance2(nodes_to_retain)=distance(nodes_to_retain);
