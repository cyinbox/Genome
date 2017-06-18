function [distance] = getDistFFTGeneSpace_ES2D(seq1,seq2,N)
% Function to select which 2D numerical mapping will be used in similarity analysis based on sequence nucleotide compositions.
% Inputs: seqences seq1 and seq2, and length N that the sequences will be even scaled to 
% Output: Euclidean distance of the two Fourier power spectra of the DNA sequences using 2D numerical representations
%
% Citation
% Yin, C., & Yau, S. S. T. (2015). An improved model for whole genome phylogenetic analysis by Fourier transform. Journal of Theoretical Biology.

% Changchuan Yin
% University of Illinois at Chicago
% Email: cyinbox@gmail.com
% Last update 06/18/2015

PS1=SE_DFTDNA2D(seq1);
PS2=SE_DFTDNA2D(seq2);

%Remove first DC in DFT power spectrum
PS2(1)=[];

PS1=evenScaleVector(PS1,N);
PS2=evenScaleVector(PS2,N);

%Compute distance
distance=getDistanceSpace(PS1,PS2);

end