%%     demo for dSPG algorithm
%
%%    dSPG: A New Discriminant Superpixel Graph Regularizer for Semi-Supervised Image Classification
%  
% Author: Long Yu
%         Jun Li,   Member, IEEE, , Fellow, IEEE, Lizhe Wang, Fellow, IEEE, Zhonghui Tang, Li Zhuo, and Yuchen Yuan
%         Lin He,
%         Antonio Plaza
%         Lizhe Wang, et al.
%         Time: September. 2022

clc
close all
clear all
addpath(genpath('./Graph_toolbox'))
addpath(genpath('./SuperPixel_toolbox'));
addpath(genpath('./self_toolbox'));

%% input param
TrainSamPerOrNum = [3,5,7,10,15]; %% number of labeled pixels per class
% TrainSamPerOrNum = [5];
alpha_ratio_all=[0.1]; %% fusion weight of between-superpixel graph and within-superpixel graph

superpix_allnum_vec=[1200]; %% preset superpixel num
thresh_weight=0.98 %% empirical parameter, for connection between superpixel in spatial neighborhood
sigma_l=180;  %% empirical parameter, for connection between superpixel

%% IndianaPine dataset
load IndianaPine
load IndianaPine_LabelMap
HyperCube= IndianaPine;
TruthMap = IndianaPine_LabelMap;
clear IndianaPine IndianaPine_LabelMap
[Height Width] = size(TruthMap);
YStart = 1; % left-up point
XStart = 1;
YEnd = 145; % right-down point
XEnd = 145;
DelBand = [1:5, 101:111,148:166,216:220];
DelClass = [1 2 8 10 17]; % use 12 classes for experiments
TrainPercent = 0.01;
PreliminarySteps1

[Height Width] = size(TruthMap);
no_classes=length(TotalSamNumAClass)

para_method='CORR';
distance = Calculate_Similarity(HyperCube, para_method);

AccRate_all={};
OA_std_all={};
AverageAcc_all={};
AA_std_all={};
IndiAccRate_all={};
kappaCoef_all={};
zscore_ESD_all={};
kappaCoef_std_all={};

for iter_sp_num=1:size(superpix_allnum_vec,1)
    
    superpix_num=superpix_allnum_vec(iter_sp_num,:)
    
    %% superpixel segmentation
    superpix_img = superpixel_cut(HyperCube,superpix_num);
    
    %% calculate the between-superpixel graph
    [WeightMat_between,distance_contracted1,connect_representative_pixel] = ...
        Between_Superpixel_Graph(HyperCube, distance,superpix_img,thresh_weight, sigma_l);

    %% calculate the within-superpixel graph
    [WeightMat_within,WeightMat2] = ...
        Within_Superpixel_Graph(HyperCube,distance, superpix_img,...
        distance_contracted1,connect_representative_pixel,thresh_weight, sigma_l);
    
    %% record the classification results
    AccRate_alpha=[];
    OA_std_alpha={};
    AverageAcc_alpha={};
    AA_std_alpha={};
    IndiAccRate_alpha={};
    kappaCoef_alpha={};
    zscore_ESD_alpha={};
    kappaCoef_std_alpha={};
    
    for alpha_ratio_i=1:length(alpha_ratio_all)
        
        alpha_ratio=alpha_ratio_all(alpha_ratio_i)
        S_W_no_normalize = alpha_ratio*WeightMat_between+(1-alpha_ratio)*WeightMat_within;
        
        N=size(S_W_no_normalize,1);
        d = sum(S_W_no_normalize,2);
        Dinv = spdiags(d,0,N,N);
        Laplace_Mat_no_normalize=Dinv-S_W_no_normalize;
        
        IndiAccRate=[];
        AA_std=[];
        kappaCoef=[];
        kappaCoef_std=[];
        zscore_ESD=[];
        OA_std=[];
        AccRate=[];
        
        FlagOfPerOrNum = 1; %the sign of setting the number of the training sample with the percent or sample number. 0 for Percent; 1 for Number;
        for PerNumNo = 1: length(TrainSamPerOrNum)
            TrainPerNum = TrainSamPerOrNum(PerNumNo);
            if FlagOfPerOrNum == 0
                NumClassSam = ceil(TrainPerNum*TotalSamNumAClass);
            else
                NumClassSam = ceil(TrainPerNum*ones(length(TotalSamNumAClass),1));
            end
            
            NeighborMask = [1];
            
            AccRate_RandNo = [];
            IndiAccRate_RandNo = [];
            kappaCoef_RandNo = [];
            zscore_ESD_RandNo = [];
            RandNo = 1:10;
            for RandNo_i=RandNo
                
                [NormTrainSam, NormTestSam, TrainLabels, TestLabels, TrainSamLoc, TestSamLoc] = fSetTrainTestSam_Neighbor_V2...
                    (HyperCube, TruthMap, SelClassNo, NeighborMask, RandNo_i, NumClassSam);
                
                NormTrainSam = NormTrainSam{1};
                NormTestSam = NormTestSam{1};
                
                %% encode the prior labels into initial_Labels
                initial_Labels = zeros(Height*Width,no_classes);
                for select_i=1:length(TrainSamLoc)
                    initial_Labels(TrainSamLoc(select_i),TrainLabels(select_i)) = 1;
                end
                initial_Labels=sparse(initial_Labels);
                
                %% graph-based classifier
                one_accuracy=[];
                para_alpha=[0.5,0.6,0.7,0.8,0.9,0.93,0.95,0.97,0.99];
                for iter_alpha=1:length(para_alpha)
                    PredLabels2 = SemiSupervised_Graph_Classifier(Laplace_Mat_no_normalize, initial_Labels, TestSamLoc, para_alpha(iter_alpha));
                    one_accuracy(1,iter_alpha) = fAccuracy(PredLabels2,TestLabels,0);
                end
                [tmp_max_OA,maxOA_idx]=max(one_accuracy,[],2);
                AccRate_RandNo = [AccRate_RandNo tmp_max_OA];
                
                [PredLabels2,PredLabels_ALL] = ...
                    SemiSupervised_Graph_Classifier(Laplace_Mat_no_normalize,...
                    initial_Labels, TestSamLoc, para_alpha(maxOA_idx(1)));
                
                IndiAccRate_RandNo = [IndiAccRate_RandNo fAccuracy(PredLabels2,TestLabels,1)];
                [kappaCoef_tmp, zscore_ESD_tmp] = fKappaCoef(PredLabels2,TestLabels);
                kappaCoef_RandNo = [kappaCoef_RandNo kappaCoef_tmp];
                zscore_ESD_RandNo = [zscore_ESD_RandNo zscore_ESD_tmp];
            end
            AccRate(:,PerNumNo) = mean(AccRate_RandNo,2)
            OA_std(:,PerNumNo) = std(AccRate_RandNo,0,2)
            IndiAccRate(:,PerNumNo) = mean(IndiAccRate_RandNo,2)
            %                     AccRate(:,PerNumNo) = max(one_accuracy,[],2)
            AA_std(:,PerNumNo) = std(mean(IndiAccRate_RandNo,1))
            
            kappaCoef(:,PerNumNo) = mean(kappaCoef_RandNo,2)
            kappaCoef_std(:,PerNumNo) = std(kappaCoef_RandNo,0,2)
            zscore_ESD(:,PerNumNo) = mean(zscore_ESD_RandNo,2);
        end
        AccRate_alpha=[AccRate_alpha;AccRate];
        OA_std_alpha{alpha_ratio_i}=OA_std;
        AverageAcc=mean(IndiAccRate,1)
        AverageAcc_alpha{alpha_ratio_i} = AverageAcc;
        AA_std_alpha{alpha_ratio_i}=AA_std;
        IndiAccRate_alpha{alpha_ratio_i} = IndiAccRate;
        kappaCoef_alpha{alpha_ratio_i} = kappaCoef;
        kappaCoef_std_alpha{alpha_ratio_i} = kappaCoef_std;
        zscore_ESD_alpha{alpha_ratio_i} = zscore_ESD;
    end
    AccRate_all{iter_sp_num}=AccRate_alpha;
    OA_std_all{iter_sp_num}=OA_std_alpha;
    AverageAcc_all{iter_sp_num} = AverageAcc_alpha;
    AA_std_all{iter_sp_num}=AA_std_alpha;
    IndiAccRate_all{iter_sp_num} = IndiAccRate_alpha;
    kappaCoef_all{iter_sp_num} = kappaCoef_alpha;
    kappaCoef_std_all{iter_sp_num} = kappaCoef_std_alpha;
    zscore_ESD_all{iter_sp_num} = zscore_ESD_alpha;
end
