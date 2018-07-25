% Plotting the phylogenetic tree by MATLAB using pair-wise distances.
% 
% Changchuan Yin, 07/22/2018
% University of Illinois at Chicago

% 1. Get data from data file
clear

fileName='distances.txt';
% Sample content in data file: pair-wise distances.
%{
gene A, gen B, 0.18
gene A, gene C, 0.12
gene A, gene D, 0.15
gene B, gene C, 0.55
gene B, gene D, 0.14
gene C, gene D, 0.92
%}

fileID = fopen(fileName);
data = textscan(fileID,'%s %s %f','Delimiter',',')
n = length(data{1});        % All names in column 1

% 2. Get the names in the order of distances and unique.
names = unique(data{1});
last=length(names)
nameLast = data{2}(n);      % And add the last element in column 2.
names(last+1) = nameLast;

% 3. Populate the gene vector.
for k = 1:length(names)
   name = names(k);
   Genes(k).Header = char(name);
end

% 4. Get the distances.
dists = data{3};

% 5. Plotting the phylogenetic tree.
UPGMAtree = seqlinkage(dists,'UPGMA',Genes);
plot(UPGMAtree,'orient','left');
xlabel('Dissimilarity distance', 'FontSize', 10,'FontWeight','bold')
title('The phylogenetic tree of genes A, B, C, and D', 'FontSize',10,'FontWeight','bold');

fclose(fileID);
