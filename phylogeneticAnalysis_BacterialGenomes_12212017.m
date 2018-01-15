% Computate the UPGMA phylogenetic tree of bacterial genomes using 2D
% representations
%
% Citation
% Yin, C., & Yau, S. S. T. (2015). An improved model for whole genome phylogenetic analysis by Fourier transform. Journal of Theoretical Biology, 382, 99-110.
%
% Changchuan Yin
% University of Illinois at Chicago
% Email: cyin1@uic.edu
% Last update 12/21/2017
clear
% The input is a fasta file containing genome sequences to be compared.
% Note: the following fasta file, 'giantVirusGenome.fasta', contains the
% DNA sequences of the following bacterium genBank access ids and names
%{
data={'AE005674':'Shigell/2astr301',  
'CP001383':'Shigell/2002017',
'CP001671':'Escherichia/ABU83972',
'CP000468':'Escherichia/APECO1',
'AE014075':'Escherichia/CFT073',
'AE005174':'Escherichia/O157H7EDL933',
'NC_002695':'Escherichia/O157H7Sakai',
'HG738867':'Escherichia/K-12MC4100',
'AL513382':'Salmonella/CT18',
'CP002614':'Salmonella/UK-1',
'AE009952':'Yersinia/KIM',
'AL590842':'Yersinia/CO92',
'CP001593':'Yersinia/Z176003',
'CP001585':'Yersinia/D106004',
'AM295250':'Staphylococcus/carTM',
'AE015929':'Staphylococcus/epiATCC',
'AP006716':'Staphylococcus/haeJCSC',
'CP001837':'Staphylococcus/lugHKU',
'CP001151':'Rhodobacter/KD131',
'CP000578':'Rhodobacter/ATCC',
'AM260480':'Ralstonia/H16',
'CP000091':'Ralstonia/JMP',
'CP001598':'Bacillus/A0248',
'AE016879':'Bacillus/Ames',
'CP001215':'Bacillus/CDC',
'AE017225':'Bacillus/Stern',
'CP000048':'Borrelia/hermsii',
'CP000976':'Borrelia/dut',
'CP000049':'Borrelia/tur',
'CP000993':'Borrelia/recu',
'CP000246':'Clostridium/perATCC',
'CP000312':'Clostridium/perSM101',
'BA000016.3':'Clostridium/per13',
'CP000527':'Desulfovibrio/vulDP4',
'CP002297':'Desulfovibrio/vulRCH1',
'AE017285':'Desulfovibrio/vulHild',
'NC_002754':'Sulfolobus/solP2', 
'NC_003106':'Sulfolobus/tok7', 
'NC_002578':'Thermoplasma/DSM1728',
'NC_002689':'Thermoplasma/volGSS1', 
}
%}
%------------------------------------------------------------------------
% 1. Get DNA sequences
NCBISeqs = fastaread('bacterialGenomes.fasta');
m=length(NCBISeqs);
L=zeros(1,m);

for k = 1:length(NCBISeqs)
    header=NCBISeqs(k).Header;
    C=textscan(header, '%s', 'delimiter', '|');    
    X=C{1}; 
    id=X(1);          %Genbank access number of the DNA sequence
    name=char(X(2))   %Short name that will be on phylogenetic tree
    seq=NCBISeqs(k).Sequence;
    L(k)=length(seq);
    Genes(k).Header   = name;
    Genes(k).Sequence = seq;
end

%-------------------------------------------------------------------------
% 2.  Compute the distance matrix of the sequences
tic
maxL=max(L)
N=length(Genes);
distM=zeros(N,N);

for i=1:N
    seq1=Genes(i).Sequence;
    for j=i:N
        seq2=Genes(j).Sequence;
        option=get2DMappingOptionSeq(seq1,seq2);
        distM(i,j)= getDistanceFFTSpace_ES2D(seq1,seq2,maxL,option);
    end
end
maxDist=max(max(distM));
dM=distM/maxDist;
   
dz=[];
for k=2:N
  dm=dM(k-1,k:N)';
  dz=[dz;dm];
end

%-------------------------------------------------------------------------
%3.  Plot UPGMA tree using the distance matrix 
UPGMAtree = seqlinkage(dz,'UPGMA',Genes);
hTree = plot(UPGMAtree,'orient','left');
xlabel('Disimilarity distance', 'FontName', 'AvantGarde','FontSize', 10,'FontWeight','bold')
%title('UPGMA phylogenetic tree of bacterial genomes', 'FontName', 'AvantGarde','FontSize', 10,'FontWeight','bold')

% Set PaperPositionMode to auto so that the exported figure looks like it does on the screen.
set(gcf, 'PaperPositionMode', 'auto');
%print -depsc2 DFTS_ES2DTree_Bacteria_12212017_fig1.eps

wtime = toc

      

