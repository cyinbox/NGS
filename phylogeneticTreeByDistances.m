% Plotting the phylogenetic tree by MATLAB using pair-wise distances.
% 
% Changchuan Yin, 07/22/2018
% Univ. of Illinois at Chicago

% 1. Get data from data file
fileName='distances.txt';
% Sample content in data file: pair-wise distances.
%{
A B 0.18
A C 0.12
A D 0.15
B C 0.55
B D 0.14
C D 0.92
%}
fileID = fopen(fileName);
data = textscan(fileID,'%s %s %f');
n = length(data{1});        % All names in column 1

% 2. Get the names in the order of distances and unique.
names = unique(data{1});
last=length(names)
nameLast = data{2}(n);      % And add the last element in column 2.
names(last+1) = nameLast;
%names = ['A','B','C','D']; % All names (column 1 and 2) in order.

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
title('The phylogenetic tree of genomes A, B, C, and D', 'FontSize',10,'FontWeight','bold');

fclose(fileID);