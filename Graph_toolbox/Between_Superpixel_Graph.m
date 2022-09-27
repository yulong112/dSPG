
function [WeightMat,distance_contracted1,connect_representative_pixel] = ...
    Between_Superpixel_Graph(HyperCube, distance, superpix_img,thresh_weight, sigma_l, kernel_h, weight_beta)
     
%
%  Between Superpixel Graph Construction
%  
%  see paper: 
%
%  Yu, L.; Li, J.; et al. "dSPG: A New Discriminant Superpixel Graph Regularizer for Semi-Supervised Image Classification"
%
% -- Input Parameters ----------------
%
% HyperCube     ->   data cube
% distance      ->   precomputed similarity weight matrix
% superpix_img  ->   superpixel segmentation results
% thresh_weight ->   threshhold value for select spatially adjacent connections between superpixel
% sigma_l       ->   parameter for SPG
% kernel_h      ->   parameter for SPG
% weight_beta   ->   parameter for SPG
% 
% Author: Long Yu, September. 2022



if ~exist('thresh_weight')    thresh_weight=0.95;end
if ~exist('kernel_h','var')
    kernel_h = 15;
end
if ~exist('weight_beta','var')
    weight_beta = 0.9;
    %     weight_beta = 0.5;
end

sp_num=max(superpix_img);
[Height,Width,Bands] = size(HyperCube);
N = Width*Height; 

%% within-superpixel SNG part
% For between-superpixel graph, additional SNG connections are added in superpixels
% to ensure that a between-superpixel graph can also be used for
% classification.
nodes_to_retain = false(N);        %% spatially adjacent within superpixel
nodes_to_beSP_SNG = false(N);      %% spatially adjacent between superpixel

for i=1:sp_num
    sp_idx1 = find(superpix_img==i);
    
    record_all=[];
    for idxii = 1:length(sp_idx1)
        pixel_idx = sp_idx1(idxii);
        
        [TempYCoord, TempXCoord] = f1DTo2DCoord([Height, Width],pixel_idx);           %计算该训练样本的坐标
        
        %% pixels in spatial neighborhoods
        candidateXCoord=[TempXCoord;TempXCoord;TempXCoord-1;TempXCoord+1];
        candidateYCoord=[TempYCoord-1;TempYCoord+1;TempYCoord;TempYCoord];
        borderFilter=(candidateXCoord>=1)&(candidateXCoord<=Width)&(candidateYCoord>=1)&(candidateYCoord<=Height);
        candidateXCoord=candidateXCoord(borderFilter);
        candidateYCoord=candidateYCoord(borderFilter);
        
        candidateLocSet_tmp=(candidateXCoord-1)*Height+candidateYCoord;
        
        %% superpixel constrain
        borderFilter=ismember(candidateLocSet_tmp,sp_idx1);
        connectLocSet=candidateLocSet_tmp(borderFilter);
        
        nodes_to_retain( pixel_idx, connectLocSet ) = true;
        nodes_to_retain( pixel_idx, pixel_idx ) = false; % diagonal should be zero
        
        outsideSuperpixelLocSet_tmp = candidateLocSet_tmp(~borderFilter);
        WeightMat_tmp=distance(pixel_idx,outsideSuperpixelLocSet_tmp);
        record_ii=[repmat(pixel_idx,length(outsideSuperpixelLocSet_tmp),1) outsideSuperpixelLocSet_tmp WeightMat_tmp'];
        record_all=[ record_all;  record_ii ];
        
    end
    %% select the part of spatially adjacent between superpixel
    connect_superpixel_idxii=find(record_all(:,3)>thresh_weight);
    %         connect_superpixel_idxii=find(record_all(:,3)>0.9);
    record_idxii=record_all(connect_superpixel_idxii,:);
    nodes_to_beSP_SNG( record_idxii(:,1)+ (record_idxii(:,2)-1)*N ) = true;
    
end
nodes_to_retain( nodes_to_retain ~= nodes_to_retain' ) = true;
nodes_to_beSP_SNG( nodes_to_beSP_SNG ~= nodes_to_beSP_SNG' ) = true;

%% connections within superpixel
nodes_to_retain_within=nodes_to_retain;
nodes_to_retain_within=nodes_to_retain_within|nodes_to_beSP_SNG;

%% calculating most representative pixels
[nodes_to_super1,distance_contracted1,connect_representative_pixel] = SuperPixel_Connect_Corr(HyperCube, distance, superpix_img);

%% combine SPG with spectral adjacent connections between superpixels
[distance_contracted,nodes_to_super] = Between_Superpixel_Graph_light(HyperCube,superpix_img,nodes_to_super1,distance_contracted1, 0, sigma_l, kernel_h, weight_beta);

for i = 1:sp_num
    pixel_set_i = find(superpix_img==i);
    similar_set=find(nodes_to_super(i,:)>0);
    for j=1:length(similar_set)
        similar_superpixel=similar_set(j);
        pixel_set_j = find(superpix_img==similar_superpixel);
        distance(pixel_set_i,pixel_set_j) = distance_contracted(i,similar_superpixel);
        distance(pixel_set_j,pixel_set_i) = distance_contracted(i,similar_superpixel);
    end
end

for i = 1:sp_num
    pixel_set_i = find(superpix_img==i);
    similar_set=find(nodes_to_super1(i,:)>0);
    for j=1:length(similar_set)
        similar_superpixel=similar_set(j);
        pixel_set_j = find(superpix_img==similar_superpixel);
        distance(pixel_set_i,pixel_set_j) = distance_contracted1(i,similar_superpixel);
        distance(pixel_set_j,pixel_set_i) = distance_contracted1(i,similar_superpixel);
    end
end

nodes_to_beSP_SNG = Connect_pixel_btwsp3(nodes_to_super,connect_representative_pixel,superpix_img, Height, Width);
%     nodes_to_beSP_SNG = Connect_pixel_btwsp2(nodes_to_super,superpix_img, distance, Height, Width, knn_distance2, param2);
nodes_to_beSP_SNG=sparse(nodes_to_beSP_SNG);


nodes_to_retain_within=nodes_to_retain_within|nodes_to_beSP_SNG;
nodes_to_retain_within = sparse(nodes_to_retain_within);
%     nodes_to_retain = nodes_to_retain_within|nodes_to_add;
nodes_to_retain = nodes_to_retain_within;
nodes_to_retain( nodes_to_retain ~= nodes_to_retain' ) = true;

%% sparse
nodes_to_retain = sparse(nodes_to_retain);
nodes_to_beSP_SNG = sparse(nodes_to_beSP_SNG);

WeightMat=zeros(N);
WeightMat(nodes_to_retain) = distance(nodes_to_retain);
WeightMat = sparse(WeightMat);

WeightMat2 = zeros(N);
WeightMat2(nodes_to_beSP_SNG)=distance(nodes_to_beSP_SNG);
WeightMat2 = sparse(WeightMat2);


