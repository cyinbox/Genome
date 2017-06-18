function [ I] = SE_DFTDNA2D( seq )
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

seq = upper(seq);
n = length(seq);
Y = zeros(2,n);

%implementation 1:
%{
for t = 1:n
   if strcmp(seq(t),'A') 
      Y(:,t) = [0 1]';
   elseif strcmp(seq(t),'C')
      Y(:,t) = [1 1]';
   elseif strcmp(seq(t),'G')
      Y(:,t) = [1 0]';
   elseif strcmp(seq(t),'T')
      Y(:,t) = [0 0]';
   end     
end
%}
%implementation 2: used in 2014 Nov paper
% reference file: SE_DFTDNA2D_RC_Research.m 
% Test file: DFTS_primate_ES2D_KeyDev_11062014

%R3 method, same leve of accuracy as MSA: This is used in the paper
%{
for j = 1:n
   if strcmp(seq(j),'A') 
      Y(:,j) = [0 -1]';
   elseif strcmp(seq(j),'C')
      Y(:,j) = [1 0]';
   elseif strcmp(seq(j),'G')
      Y(:,j) = [0 1]';
   elseif strcmp(seq(j),'T')
      Y(:,j) = [-1 0]';
   end     
end
%}
%R1 method:not as good as R3 and MSA, same as 4D
%{
for j = 1:n
   if strcmp(seq(j),'A') 
      Y(:,j) = [0 -1]';
   elseif strcmp(seq(j),'C')
      Y(:,j) = [1 0]';
   elseif strcmp(seq(j),'G')
      Y(:,j) = [-1 0]';
   elseif strcmp(seq(j),'T')
      Y(:,j) = [0 1]';
   end     
end
%}
%R2 method, similar as R1,similar as 4D
for j = 1:n
   if strcmp(seq(j),'A') 
      Y(:,j) = [0 -1]';
   elseif strcmp(seq(j),'C')
      Y(:,j) = [0 1]';
   elseif strcmp(seq(j),'G')
      Y(:,j) = [1 0]';
   elseif strcmp(seq(j),'T')
      Y(:,j) = [-1 0]';
   end     
end


%mY = mean(Y);
%S = cov(Y);
%sz=size(S)
%Method 1: Calculate the DFT from each row of Y

FFT1 = fft(Y(1,:),n); 
FFT2 = fft(Y(2,:),n); 

%FFTY=FFTA+FFTC+FFTG;

%half=floor(n/2);
%Periodiogram 1
%I = sqrt(n)*abs(FFTY).^2;
%I(10)
%figure
%plot(I(2:half))


%Periodiogram 2:stable method based on different base as zero
I1 = abs(FFT1).^2;
I2 = abs(FFT2).^2;
I = I1+I2;

end





