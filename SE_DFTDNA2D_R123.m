
function [I] = SE_DFTDNA2D_R123(seq,option )
% Function to compute FFT of a DNA seuence using different 2D numerical representation.
% Inputs: seqence and option for which 2D will be used.
% Output: Fourier power spectrum of the DNA sequence
%
% Citation
% Yin, C., & Yau, S. S. T. (2015). An improved model for whole genome phylogenetic analysis by Fourier transform. Journal of Theoretical Biology.
%
% Changchuan Yin
% University of Illinois at Chicago
% Email: cyinbox@gmail.com
% Last update 06/18/2015

seq = upper(seq);
n = length(seq);
Y = zeros(2,n);

%R3 method, same leve of accuracy as MSA: This is used in the paper
if strcmp(option,'AT') %R1
    a=[0 -1]';
    t=[0 1]';
    c=[-1 0]';
    g=[1 0]';
elseif strcmp(option,'AC') %R2
    a=[0 -1]';
    t=[-1 0]';
    c=[0 1]';
    g=[1 0]';
elseif strcmp(option,'AG') %R3
    a=[0 -1]';
    t=[-1 0]';
    c=[1 0]';
    g=[0 1]';
end
 for j = 1:n
   if strcmp(seq(j),'A') 
      Y(:,j) = a;
   elseif strcmp(seq(j),'C')
      Y(:,j) = c;
   elseif strcmp(seq(j),'G')
      Y(:,j) = g;
   elseif strcmp(seq(j),'T')
      Y(:,j) = t;
   end     
 end
 
FFT1 = fft(Y(1,:),n);
FFT2 = fft(Y(2,:),n); 

%Periodiogram 2:stable method based on different base as zero
I1 = abs(FFT1).^2;
I2 = abs(FFT2).^2;
I = I1+I2;

end





