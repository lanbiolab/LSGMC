function [acc,nmi,F,precision,AR,Purity,Recall] = AllMeasure(label,true)
% label -- Clustering results column vector
% true -- Real label column vector
type = unique(label);
newlabel = label;
for i=1:length(type)
    idx = abs(label-type(i))<0.1;
    truelabel = true(idx);
    typenum = unique(truelabel);
    nowmax = 0;
    Determine = 0;
    for j=1:length(typenum)
        temp2 = sum(abs(truelabel-typenum(j))<0.1);
        if temp2>nowmax
            Determine = typenum(j);
            nowmax = temp2;
        end
    end
    newlabel(idx) = Determine;
end
%% acc nmi f Purity Recall precision
acc = sum(abs(newlabel-true)<0.1)/length(true);
[~,nmi,~]=compute_nmi(true,newlabel);
N = length(true);
numT = 0;
numH = 0;
numI = 0;
for n=1:N
    Tn = (true(n+1:end))==true(n);
    Hn = (newlabel(n+1:end))==newlabel(n);
    numT = numT + sum(Tn);
    numH = numH + sum(Hn);
    numI = numI + sum(Tn .* Hn);
end
precision = 1;
Recall = 1;
if numH > 0
    precision = numI / numH;
end
if numT > 0
    Recall = numI / numT;
end
if (precision+Recall) == 0
    F = 0;
else
    F = 2 *  precision * Recall / (precision + Recall);
end
correnum = 0;
for ci = 1:length(type)
    incluster = true(label == type(ci));
    inclunub = hist(incluster,1:max(incluster));
    if isempty(inclunub)
        inclunub=0;
    end
    correnum = correnum + max(inclunub);
end
Purity = correnum/N;
%% ARI
matrix=zeros(max(label),max(true));
for i = 1:N
   matrix(label(i),true(i))=matrix(label(i),true(i))+1;
end
n=sum(sum(matrix));
nis=sum(sum(matrix,2).^2);		%sum of squares of sums of rows
njs=sum(sum(matrix,1).^2);		%sum of squares of sums of columns

t1=nchoosek(n,2);		%total number of pairs of entities
t2=sum(sum(matrix.^2));	%sum over rows & columnns of nij^2
t3=.5*(nis+njs);
%Expected index (for adjustment)
nc=(n*(n^2+1)-(n+1)*nis-(n+1)*njs+2*(nis*njs)/n)/(2*(n-1));

A=t1+t2-t3;		%no. agreements
if t1==nc
   AR=0;			%avoid division by zero; if k=1, define Rand = 0
else
   AR=(A-nc)/(t1-nc);		%adjusted Rand - Hubert & Arabie 1985
end
