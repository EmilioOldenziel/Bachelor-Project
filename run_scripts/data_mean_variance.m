% read the data from the files
[spectra, labels, wavelengths] = loadDataset('sugar_neoVNIR1600.csv');
[spectra_val, labels_val, wavelengths_val] = loadDataset('sugar_neoVNIR1600_val.csv');

%reduce from 9 to 3 classes #problem 1
problem1labels = labels
for i=1:length(labels)
    if problem1labels(i) < 5
        problem1labels(i) = 1;
    elseif problem1labels(i) > 6
        problem1labels(i) = 3;
    else
        problem1labels(i) = 2;
    end
end

c1 = spectra(problem1labels == 1, :);
c2 = spectra(problem1labels == 2, :);
c3 = spectra(problem1labels == 3, :);

errorbar(wavelengths, mean(c1), var(c1));
hold on
errorbar(wavelengths, mean(c2), var(c2));
hold on
errorbar(wavelengths, mean(c3), var(c3));