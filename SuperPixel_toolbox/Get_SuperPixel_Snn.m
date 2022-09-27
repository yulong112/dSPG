
function snn_sp_cell = Get_SuperPixel_Snn(HyperCube, superpix_img)

[Height,Width,Bands] = size(HyperCube);
sp_num=max(superpix_img);
snn_sp_cell=cell(1,double(sp_num));
for i=1:sp_num
    sp_idx1=find(superpix_img==i);
    PointSet = [];
    for j=1:length(sp_idx1)
%         [TempYCoord, TempXCoord] = f1DTo2DCoord([Height, Width],sp_idx1(j));           %�����ѵ������������
        PointLocTemp=[];
        %�����������ͼ���±߽磬ͼ���ϱ߽磬ͼ���м�����λ��
        if mod(sp_idx1(j),Height)==0
            PointLocTemp = [sp_idx1(j)-Height,sp_idx1(j)+Height,sp_idx1(j)-1];
        elseif mod(sp_idx1(j)-1,Height)==0
            PointLocTemp = [sp_idx1(j)-Height,sp_idx1(j)+Height,sp_idx1(j)+1];
        else
            PointLocTemp = [sp_idx1(j)-Height,sp_idx1(j)+Height,sp_idx1(j)-1,sp_idx1(j)+1];
        end
        %�ٽ���������<1��>Height*Width��ȥ��
        delete_idx=find(PointLocTemp<1|PointLocTemp>Height*Width);
        if ~isempty(delete_idx)
            PointLocTemp(delete_idx)=[];
        end
        
        sp_class=superpix_img(PointLocTemp);
        sp_class=unique(sp_class);
        %��0Ԫ��ȥ��
        zero_idx=find(sp_class==0);
        if ~isempty(zero_idx)
            sp_class(zero_idx)=[];
        end
        %�����ڵ���i��Ԫ��ȥ��
        upper_idx=find(sp_class>=i);
        if ~isempty(upper_idx)
            sp_class(upper_idx)=[];
        end
        PointSet=[PointSet;sp_class];
    end
    PointSet=unique(PointSet);
    %��0Ԫ��ȥ��
    zero_idx=find(PointSet==0);
    if ~isempty(zero_idx)
        PointSet(zero_idx)=[];
    end
    %�����ڵ���i��Ԫ��ȥ��
    upper_idx=find(PointSet>=i);
    if ~isempty(upper_idx)
        PointSet(upper_idx)=[];
    end
    snn_sp_cell{i}=PointSet;
end