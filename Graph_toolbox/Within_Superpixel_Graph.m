
function [WeightMat,WeightMat2] = ...
    Within_Superpixel_Graph(HyperCube,distance, superpix_img,...
    distance_contracted1,connect_representative_pixel,...
    thresh_weight, sigma_l,spatial_hop_n,knn_param_sp0, kernel_h, weight_beta)

%
%  Within Superpixel Graph Construction
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
% distance_contracted1  ->   similarity weight matrix of SPG
% connect_representative_pixel  ->   connected representative pixel pairs of SPG
% thresh_weight ->   threshhold value for select spatially adjacent connections between superpixel
% sigma_l       ->   parameter for SPG
% spatial_hop_n       ->   parameter for within-superpixel graph, empirically
% knn_param_sp0       ->   parameter for within-superpixel graph, empirically
% kernel_h      ->   parameter for SPG
% weight_beta   ->   parameter for SPG
%
% Author: Long Yu, September. 2022


% if ~exist('thresh_weight')    thresh_weight=0.82;end
if ~exist('thresh_weight')    thresh_weight=0.95;end
if ~exist('spatial_hop_n')    spatial_hop_n=3;end
if ~exist('knn_param_sp0')    knn_param_sp0=20;end
if ~exist('kernel_h','var')
    kernel_h = 15;
end
if ~exist('weight_beta','var')
    weight_beta = 0.9;
end

sp_num=max(superpix_img);
[Height,Width,Bands] = size(HyperCube);
ALL_data = reshape(HyperCube, Width*Height, Bands); %3D变 2D
N = Width*Height; % k = 256, l = 1100, m = 10

knn_distance2 = zeros(N,1);
%% spatially adjacent within superpixel
nodes_to_retain = false(N);
%% spectrally adjacent within superpixel
nodes_to_knn_withinsp = false(N);
%% spatially adjacent between superpixel
nodes_to_beSP_SNG = false(N);

for i=1:sp_num
    sp_idx1 = find(superpix_img==i);
    sp_num_i = length(sp_idx1);
    
    knn_param_sp=knn_param_sp0;
    if sp_num_i-1<knn_param_sp
        knn_param_sp=sp_num_i-1;
    end
    knn_distance_sp = zeros(sp_num_i,1);
    
    record_all=[];
    for idxii = 1:length(sp_idx1)
        pixel_idx = sp_idx1(idxii);
        temp = sort(distance(pixel_idx,sp_idx1), 'descend');
        knn_distance_sp(idxii) = temp(knn_param_sp + 1);
        
        %% Spa-NN
        [TempYCoord, TempXCoord] = f1DTo2DCoord([Height, Width],pixel_idx);           %计算该训练样本的坐标
        
        %% filter
        candidateXCoord=(TempXCoord-1):(TempXCoord+1);
        candidateYCoord=(TempYCoord-1):(TempYCoord+1);
        
        borderFilter=(candidateXCoord>=1)&(candidateXCoord<=Width);
        candidateXCoord=candidateXCoord(borderFilter);
        borderFilter=(candidateYCoord>=1)&(candidateYCoord<=Height);
        candidateYCoord=candidateYCoord(borderFilter);
        
        candidateXCoord1=repmat(candidateXCoord,length(candidateYCoord),1);
        candidateYCoord1=repmat(candidateYCoord',1,length(candidateXCoord));
        
        candidateLocSet_1hop=(candidateXCoord1-1)*Height+candidateYCoord1;
        
        candidateLocSet_1hop=candidateLocSet_1hop(:);
        
        %% 
        candidateXCoord=(TempXCoord-spatial_hop_n):(TempXCoord+spatial_hop_n);
        candidateYCoord=(TempYCoord-spatial_hop_n):(TempYCoord+spatial_hop_n);
        
        borderFilter=(candidateXCoord>=1)&(candidateXCoord<=Width);
        candidateXCoord=candidateXCoord(borderFilter);
        borderFilter=(candidateYCoord>=1)&(candidateYCoord<=Height);
        candidateYCoord=candidateYCoord(borderFilter);
        
        candidateXCoord1=repmat(candidateXCoord,length(candidateYCoord),1);
        candidateYCoord1=repmat(candidateYCoord',1,length(candidateXCoord));
        
        candidateLocSet_tmp=(candidateXCoord1-1)*Height+candidateYCoord1;
        
        candidateLocSet_tmp=candidateLocSet_tmp(:);
        
        %% superpixel constrain
        borderFilter=ismember(candidateLocSet_tmp,sp_idx1);
        connectLocSet=candidateLocSet_tmp(borderFilter);
        
        nodes_to_retain( pixel_idx, connectLocSet ) = true;
        nodes_to_retain( pixel_idx, pixel_idx ) = false; % diagonal should be zero
        
        %% record
        outsideSuperpixelLocSet_tmp = candidateLocSet_tmp(~borderFilter);
        
        if ~isempty(outsideSuperpixelLocSet_tmp)
            %% filter constrain
            %                 borderFilter=ismember(outsideSuperpixelLocSet_tmp,candidateLocSet_1hop);
            %                 outsideSuperpixelLocSet_tmp = outsideSuperpixelLocSet_tmp(borderFilter);
            
            if ~isempty(outsideSuperpixelLocSet_tmp)
                WeightMat_tmp=distance(pixel_idx,outsideSuperpixelLocSet_tmp);
                record_ii=[repmat(pixel_idx,length(outsideSuperpixelLocSet_tmp),1) outsideSuperpixelLocSet_tmp WeightMat_tmp'];
                
                %% all connections
                record_all=[ record_all;  record_ii ];
            end
        end
    end
    
    %% 
    connect_superpixel_idxii=find(record_all(:,3)>thresh_weight);
    record_idxii=record_all(connect_superpixel_idxii,:);
    nodes_to_beSP_SNG( record_idxii(:,1)+ (record_idxii(:,2)-1)*N ) = true;
    
    %% spectrally adjacent within superpixel
    for idxii = 1:length(sp_idx1)
        pixel_idx = sp_idx1(idxii);
        nodes_to_knn_withinsp(pixel_idx, sp_idx1( distance(pixel_idx,sp_idx1) > knn_distance_sp(idxii) ) ) = true;
        
        nodes_to_knn_withinsp(pixel_idx,pixel_idx) = false; % diagonal should be zero
    end
    
end

nodes_to_retain( nodes_to_retain ~= nodes_to_retain' ) = true;
nodes_to_beSP_SNG( nodes_to_beSP_SNG ~= nodes_to_beSP_SNG' ) = true;
nodes_to_knn_withinsp( nodes_to_knn_withinsp ~= nodes_to_knn_withinsp' ) = true;
nodes_to_beSP_SNG=sparse(nodes_to_beSP_SNG);

nodes_to_retain_within=nodes_to_retain | nodes_to_knn_withinsp;

%% connection between superpixels, ==SPG
[distance_contracted,nodes_to_super] = SPG_connections_RBF_ACS(HyperCube, distance_contracted1, superpix_img, sigma_l, kernel_h, weight_beta);

nodes_to_beSP_SNG2 = Connect_pixel_btwsp3(nodes_to_super,connect_representative_pixel,superpix_img, Height, Width);
nodes_to_beSP_SNG2=sparse(nodes_to_beSP_SNG2);

%% 
nodes_to_retain_within=nodes_to_retain_within|nodes_to_beSP_SNG;

%% sparse
nodes_to_retain_within = sparse(nodes_to_retain_within);
nodes_to_beSP_SNG = sparse(nodes_to_beSP_SNG);

WeightMat=zeros(N);
WeightMat(nodes_to_retain_within) = distance(nodes_to_retain_within);
%     WeightMat(nodes_to_beSP_SNG) = 0.1*WeightMat(nodes_to_beSP_SNG);
WeightMat = sparse(WeightMat);

WeightMat2 = zeros(N);
WeightMat2(nodes_to_beSP_SNG)=distance(nodes_to_beSP_SNG);
WeightMat2 = sparse(WeightMat2);

