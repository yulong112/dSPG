function [PrimaryComponent, PCDirection, PCVariance]=...
    fPCA_FinalVersionV2(Data,Para, VarOrNum, Style, PCOrPCDirection)
%VarOrNum：判断Para参数是主成分方差能量比，还是主成分个数
%                   VarOrNum=1 Para为主成分方差能量比，VarOrNum=0 Para为主成分个数
%Style 0 normal PCA
%        1 regularized PCA
%PCOrPCDirection-- 0 or does not exist:  both of PC and PC direction
%                                  1: PC direction only
%

DimNum = ndims(Data);
if DimNum == 3
    [Height Width Dim] = size(Data);
    SamNum = Height*Width;
    Data = reshape(Data, SamNum, Dim);
else
    [SamNum,Dim] = size(Data);%使用的源数据
end

CoVar = cov(Data);
[V,D] = eig(CoVar) ;

% if rank(CoVar) == Dim
%     [V,D] = eig(CoVar) ;
% else
%     [V,D] = eig(CoVar+0.00001*eye(Dim,Dim)) ; % to avoid the singularity of matrix
% end

EigenValueTemp = diag(D);
% for i=1:Dim
%     EigenValueTemp(i)=D(i,i);%将特征值提取到一向量中
% end

[EigenValue,Index]= sort(EigenValueTemp,'descend');%特征值由小到大排序
EigenMatrix=V(:,Index);

if  VarOrNum==1
    EigenValueSum=sum(EigenValue);
    TempEigenValueSum = 0;
    EigenVectorNum = 0;
    for i=1:Dim
        Temp=EigenValue(i);
        TempEigenValueSum = TempEigenValueSum+Temp;
        EigenValueWeigeht=TempEigenValueSum/EigenValueSum;
        if EigenValueWeigeht>Para
            break;
        end
        EigenVectorNum=EigenVectorNum+1;
    end
end

if  VarOrNum==0
    EigenVectorNum=Para;
end

if Style == 1
    for i = 1:size(EigenMatrix,2)
        TempNorm(1,i) = norm(EigenMatrix(:,i));
    end
    EigenMatrix = EigenMatrix./repmat(TempNorm,size(EigenMatrix,1),1);
end

PCDirection = EigenMatrix(:,1:EigenVectorNum);
PCVariance = EigenValue(1:EigenVectorNum);
PrimaryComponent = [];
if exist('PCOrPCDirection')
    if PCOrPCDirection==0
        PrimaryComponent = Data*PCDirection;
    end
else
    PrimaryComponent = Data*PCDirection;
end

if DimNum == 3
    PrimaryComponent = reshape(PrimaryComponent, Height, Width, EigenVectorNum);
else
    
end
