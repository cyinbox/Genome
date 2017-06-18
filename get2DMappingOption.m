function [option] = get2DMappingOption(seq1,seq2)
% Function to select which 2D numerical mapping will be used in similarity analysis based on sequence nucleotide compositions.
% Inputs: DNA sequences 
% Output: Choice indicator for 2D numerical mapping
%
% Citation
% Yin, C., & Yau, S. S. T. (2015). An improved model for whole genome phylogenetic analysis by Fourier transform. Journal of Theoretical Biology.
%
% Changchuan Yin
% University of Illinois at Chicago
% Email: cyinbox@gmail.com
% Last update 06/18/2015

Bases= basecount(seq1);
AV(1)=Bases.A;
TV(1)=Bases.T;
CV(1)=Bases.C;
GV(1)=Bases.G;
NV(1)=Bases.A+Bases.T+Bases.C+Bases.G;
    
Bases2= basecount(seq2);
AV(2)=Bases2.A;
TV(2)=Bases2.T;
CV(2)=Bases2.C;
GV(2)=Bases2.G;
NV(2)=Bases2.A+Bases2.T+Bases2.C+Bases2.G;

cA= sum(AV);
cT=sum(TV);
cC=sum(CV);
cG=sum(GV);
W=cA+cT+cG+cC;

%cAT=abs(0.5-(cA+cT)/W)
cAC=abs(0.5-(cA+cC)/W);
cAG=abs(0.5-(cA+cG)/W);

%cATCG(1)=cAT; %Do not consder AT
cATCG(1)=cAC;
cATCG(2)=cAG;
x=min(cATCG);

if x==cAC
    option='AC';
else 
    option='AG';
end

end

