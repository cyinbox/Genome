% Computate the UPGMA phylogenetic tree of the primate NADH dehydrogenase subunit 4 genes.
%
% Citation
% Yin, C., & Yau, S. S. T. (2015). An improved model for whole genome phylogenetic analysis by Fourier transform. Journal of Theoretical Biology.
%
% Changchuan Yin
% University of Illinois at Chicago
% Email: cyinbox@gmail.com
% Last update 06/18/2015

clear

% Macaca fascicularis NADH dehydrogenase subunit 4 gene, partial cds; tRNA-His, 
% tRNA-Ser, and tRNA-Leu genes, complete sequence; 
% and NADH dehydrogenase subunit 5 gene, partial cds; mitochondrial genes for mitochondrial products
data={'M22653' 'Macaca fascicular';
      'M22651' 'Macaca fuscata';
      'M22650' 'Macaca mulatta';
      'M22655' 'Saimiri sciureus';
      'M22654' 'Macaca sylvanus';
      'V00672' 'Chimpanzee';
      'M22657' 'Lemur catta';
      'V00658' 'Gorilla';
      'V00659' 'Hylobates';
      'V00675' 'Orangutan';  
      'M22656' 'Tarsisus syrichta';
      'L00016' 'Human';
      };
  
 UPGMAtree = phylogeneticTree(data);


   