function [newX,T,meantrain] = pca_row(X,select_k,CRate)
%ÿ����һ������
%newX  ��ά����¾���
%T �任����
%meanValue  Xÿ�о�ֵ���ɵľ������ڽ���ά��ľ���newX�ָ���X
%CRate ������
%�������Ļ���������
meantrain=mean(X);
meanValue=ones(size(X,1),1)*meantrain;
X=X-meanValue;%ÿ��ά�ȼ�ȥ��ά�ȵľ�ֵ
% X=X - repmat(mean(X,2),1,size(X,2));%ÿ��ά�ȼ�ȥ��ά�ȵľ�ֵ
% C=X'*X/(size(X,1)-1);%����Э�������
% C=X*X'/(size(X,2)-1);%����Э�������
if  select_k==0
    [U,S,V] = svd(X);
else
    [U,S,V] = svds(X,select_k);
end
%������������������ֵ
% [V,D]=eig(C);
%��������������������
% [dummy,order]=sort(diag(-D));
% V=V(:,order);%������������������ֵ��С���н�������
% d=diag(D);%������ֵȡ��������һ��������
% newd=d(order);%������ֵ���ɵ�����������������
if select_k==0
    %ȡǰn���������������ɱ任����
    newd=diag(S);%������ֵȡ��������һ��������
    sumd=sum(newd);%����ֵ֮��
    for j=1:length(newd)
        i=sum(newd(1:j,1))/sumd;%���㹱���ʣ�������=ǰn������ֵ֮��/������ֵ֮��
        if i>CRate%�������ʴ���95%ʱѭ������,������ȡ���ٸ�����ֵ
            cols=j;
            break;
        end
    end
else
    cols=select_k;
end

T=V(:,1:cols);%ȡǰcols���������������ɱ任����T
newX=X*T;%�ñ任����T��X���н�ά
end