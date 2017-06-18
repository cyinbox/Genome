function [ UPGMAtree ] = phylogeneticTree(data)
% Computate the UPGMA phylogenetic tree
% Inputs: Gene data cell containing sequence and its header.
% Output: UPGMA tree
%
% Citation
% Yin, C., & Yau, S. S. T. (2015). An improved model for whole genome phylogenetic analysis by Fourier transform. Journal of Theoretical Biology.
%
% Changchuan Yin
% University of Illinois at Chicago
% Email: cyinbox@gmail.com
% Last update 06/18/2015
  
L=zeros(1,length(data));

for k = 1:length(data)
    name=data{k,2};
    seq=getgenbank(data{k,1},'sequenceonly','true');
    L(k)=length(seq);
    Genes(k).Header   = name;
    Genes(k).Sequence = seq;
end
maxL=max(L)
N=length(Genes);
distM=zeros(N,N);

for i=1:N
    seq1=Genes(i).Sequence;
    for j=i:N
        seq2=Genes(j).Sequence;
        
        option=get2DMappingOption(seq1,seq2);
        distM(i,j)= getDistanceFFTSpace_ES2D(seq1,seq2,maxL,option)
           
    end
end

maxDist=max(max(distM));
dM=distM/maxDist;
   
dz=[];
for k=2:N
  dm=dM(k-1,k:N)';
  dz=[dz;dm];
end

UPGMAtree = seqlinkage(dz,'UPGMA',Genes);
hTree = plot(UPGMAtree,'orient','left');
xlabel('Similarity distance', 'FontName', 'AvantGarde','FontSize', 10,'FontWeight','bold')
title('UPGMA phylogenetic tree', 'FontName', 'AvantGarde','FontSize', 10,'FontWeight','bold')

end

