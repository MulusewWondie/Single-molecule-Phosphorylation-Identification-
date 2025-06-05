clc; clear; close all;

% This script clusters Raman spectra from A, B, and contamination reference C.
% KMeans (k=3) is applied to identify pure-A, pure-B, and contaminated groups.
% The cluster overlapping with spectra from C is labeled as contaminated.
% Spectra from A and B are labeled pure if they do not belong to the contaminated cluster.
% t-SNE is used for 2D visualization of all spectra and pure spectra only.
% Extracted spectra are saved into 'pure-A', 'pure-B', and 'contaminated' folders.


% Folder paths
folderA = ""; %  Directory of Target molecule 1
folderB = ""; % Directory of Target molecule 2
folderC = ""; % Directory of Citrate

% Clustering of spectra into pure-A, pure-B, and contaminated; t-SNE plotting.

% Create output folders
mkdir('pure-A');      % folder for pure-A spectra
mkdir('pure-B');      % folder for pure-B spectra
mkdir('contaminated');% folder for contaminated spectra
mkdir('figures');     % folder for figures

% Load spectra from folders
[wavenumberA, spectraA] = loadSpectra(folderA);
[wavenumberB, spectraB] = loadSpectra(folderB);
[wavenumberC, spectraC] = loadSpectra(folderC);

% Get file lists for copying later
filesA = dir(fullfile(folderA, '*.txt'));
filesB = dir(fullfile(folderB, '*.txt'));
filesC = dir(fullfile(folderC, '*.txt'));

% Combine all spectra into one matrix
allSpectra = [spectraA; spectraB; spectraC];
numClusters = 3;  % Three clusters: pure-A, pure-B, contaminated

% Perform KMeans clustering
[idx, ~] = kmeans(allSpectra, numClusters);

% Split cluster indices for each set
nA = size(spectraA, 1);
nB = size(spectraB, 1);
idxA = idx(1:nA);
idxB = idx(nA+1:nA+nB);
idxC = idx(nA+nB+1:end);

% Identify contaminated cluster as the majority cluster in C
contaminatedCluster = mode(idxC);

% Determine pure vs contaminated indices
pureA_indices         = find(idxA ~= contaminatedCluster);
contaminatedA_indices = find(idxA == contaminatedCluster);
pureB_indices         = find(idxB ~= contaminatedCluster);
contaminatedB_indices = find(idxB == contaminatedCluster);

% ================================
% t-SNE Visualization of all data
% ================================
tsneResult = tsne(allSpectra);
figure(1);
gscatter(tsneResult(:,1), tsneResult(:,2), idx, 'rgb', 'ox+');
xlabel('t-SNE 1', 'FontSize',14);
ylabel('t-SNE 2', 'FontSize',14);
title('KMeans Clustering of All Spectra');
legend('Ser','pSer','Citrate');
grid off;
saveas(gcf, fullfile('figures', 'kmeans_clusters_tsne.png'));

% ============================================
% t-SNE Visualization of pure-A vs pure-B only
% ============================================
pureSpectra = [spectraA(pureA_indices, :); spectraB(pureB_indices, :)];
pureLabels  = [ones(numel(pureA_indices),1); 2*ones(numel(pureB_indices),1)];
tsneResultPure = tsne(pureSpectra);
figure(2);
gscatter(tsneResultPure(:,1), tsneResultPure(:,2), pureLabels, 'rg', 'ox');
xlabel('t-SNE 1','FontSize',14);
ylabel('t-SNE 2','FontSize',14);
title('t-SNE of Citrate-free Ser Vs pSer ');
legend('Ser','pSer','FontSize',14);
grid off;
saveas(gcf, fullfile('figures', 'pure_spectra_tsne.png'));

% ================================
% Copy clustered files into output folders
% ================================
for i = 1:numel(pureA_indices)
    copyfile(fullfile(folderA, filesA(pureA_indices(i)).name), fullfile('pure-A', filesA(pureA_indices(i)).name));
end
for i = 1:numel(contaminatedA_indices)
    copyfile(fullfile(folderA, filesA(contaminatedA_indices(i)).name), fullfile('contaminated', filesA(contaminatedA_indices(i)).name));
end
for i = 1:numel(pureB_indices)
    copyfile(fullfile(folderB, filesB(pureB_indices(i)).name), fullfile('pure-B', filesB(pureB_indices(i)).name));
end
for i = 1:numel(contaminatedB_indices)
    copyfile(fullfile(folderB, filesB(contaminatedB_indices(i)).name), fullfile('contaminated', filesB(contaminatedB_indices(i)).name));
end

disp('Clustering complete. Results saved in pure-A, pure-B, contaminated, and figures folders.');

% ========== Helper Function ==========
function [wavenumber, spectra] = loadSpectra(folderPath)
    files   = dir(fullfile(folderPath, '*.txt'));
    nFiles  = numel(files);
    first   = readmatrix(fullfile(folderPath, files(1).name));
    wavenumber = first(:,1);
    spectra = zeros(nFiles, numel(wavenumber));
    for j = 1:nFiles
        data = readmatrix(fullfile(folderPath, files(j).name));
        spectra(j,:) = data(:,2);
    end
end

