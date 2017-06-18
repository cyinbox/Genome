% Function to select which 2D numerical mapping will be used in similarity analysis based on sequence nucleotide compositions.
% Inputs: Fourier power spectra, PS1 and PS2 
% Output: Euclidean distance of the two Fourier power spectra
%
% Changchuan Yin
% University of Illinois at Chicago
% Email: cyinbox@gmail.com
% Last update 06/18/2015
%
% Citation
% Yin, C., & Yau, S. S. T. (2015). An improved model for whole genome phylogenetic analysis by Fourier transform. Journal of Theoretical Biology.
function [ distance ] = getDistanceSpace( PS1,PS2 )
   distance = sqrt(sum(abs(PS1-PS2).^2));
end

