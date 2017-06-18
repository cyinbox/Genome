function [distance] = getDistanceFFTSpace_ES2D(seq1,seq2,N,option)
% Function to compute the Euclidean distance of two DNA sequences.
% Inputs: seqences seq1 and seq2, and length N that the sequences will be
% even scaled to, and option for which 2D will be used.
% Output: Euclidean distance of the two Fourier power spectra of the DNA
% sequences using 2D numerical representations.
%
% Citation
% Yin, C., & Yau, S. S. T. (2015). An improved model for whole genome phylogenetic analysis by Fourier transform. Journal of Theoretical Biology.
%
% Changchuan Yin
% University of Illinois at Chicago
% Email: cyinbox@gmail.com
% Last update 06/18/2015
%

  PS1=SE_DFTDNA2D_R123(seq1,option); % 1,0,-1 mapping
  PS2=SE_DFTDNA2D_R123(seq2,option);
 
  PS1(1)=[];
  PS2(1)=[];

  PS1=evenScaleVector(PS1,N);
  PS2=evenScaleVector(PS2,N);
  distance=getDistanceSpace(PS1,PS2);

end