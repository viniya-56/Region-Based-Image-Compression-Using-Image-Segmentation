clc
clear
close all;

foregroundQuality = 90;
backgroundQuality = 35;
overallQuality = 100;
completeCompressionQuality = 35;

disp('Available image choices:');
disp('1. Desert Image (bmp)');
disp('2. Complex Image (bmp)');
 
image_choice = input("Enter the image choice (digits 1-2): ");
switch image_choice
    case 1
        img = 'desert.bmp';
    case 2
        img = 'complex_img.bmp';
    otherwise
        error('Error. Select from the available choices.');
end

image = imread(img);
grayImage = rgb2gray(image);
 
figure;
subplot(1,2,1), imshow(img), title("Original Image");
subplot(1,2,2), imshow(grayImage), title("Grayscale Image");

% thresholding using Otsu method

thresholdLevel = graythresh(grayImage);  
binaryImage = imbinarize(grayImage, thresholdLevel); % Convert to binary
 
otsuThreshold = thresholdLevel * 255;
fprintf('\nOtsu''s Threshold Value: %.2f\n', otsuThreshold);
 
figure; % Plotting the histogram
subplot(1,2,1), imhist(grayImage); hold on;          
xline(otsuThreshold, 'r', 'LineWidth', 2); hold off; 
title('Histogram of Grayscale Image');

% Clean the image
cleanedImage = imopen(binaryImage, strel('disk', 5));
cleanedImage = imfill(cleanedImage, 'holes'); subplot(1,2,2)
imshow(~cleanedImage), title("Cleaned mask after thresholding");

% Label connected components and extract the largest one (foreground)
labeledImage = bwlabel(cleanedImage);  
measurements = regionprops(labeledImage, 'Area');
allAreas = [measurements.Area];
[~, largestBlobIndex] = max(allAreas);
objectSegment = ismember(labeledImage, largestBlobIndex);
objectMask = uint8(objectSegment); % Convert logical mask to uint8

% Separate foreground object (ROI) from the background
background = bsxfun(@times, image, cat(3, objectMask, objectMask, objectMask));
foreground = image - background;

figure;
subplot(1, 2, 1), imshow(foreground), title('Foreground (High Quality)');
subplot(1, 2, 2), imshow(background), title('Background (Low Quality)');

outputDir = 'Thresholding JPG Images';
if ~exist(outputDir, 'dir')
    mkdir(outputDir); % Create the directory
end
 
% Compress the foreground with higher quality (lower compression)
imwrite(foreground, fullfile(outputDir, 'foreground_high_quality.jpg'), 'Quality', foregroundQuality); % Higher quality
 
% Compress the background with lower quality (higher compression)
imwrite(background, fullfile(outputDir, 'background_low_quality.jpg'), 'Quality', backgroundQuality); % Lower quality
 
highQualityForeground = imread(fullfile(outputDir, 'foreground_high_quality.jpg'));
lowQualityBackground = imread(fullfile(outputDir, 'background_low_quality.jpg'));
 
% Combine the high-quality foreground and low-quality background
combinedImage = highQualityForeground + lowQualityBackground;
 
% Save the image
imwrite(combinedImage, fullfile(outputDir, 'combined_selective_compression.jpg'), 'quality', overallQuality);
 
% Complete compression
imwrite(image, 'complete_compressed_img.jpg', 'jpg', 'Quality', completeCompressionQuality);
 
orig = imread(img);
complete_compressed = imread('complete_compressed_img.jpg');
selective_compressed = imread(fullfile(outputDir, 'combined_selective_compression.jpg'));

figure;
subplot(1, 3, 1), imshow(image), title('Original Image');
subplot(1, 3, 2), imshow(combinedImage), title('Selectively Compressed');
subplot(1, 3, 3), imshow(complete_compressed), title('Complete Compressed');
 
% Calculating performance metrics
mse1 = immse(orig, complete_compressed);
mse2 = immse(orig, selective_compressed);
fprintf('\nMSE Selective Compression: %f\n', mse2);
fprintf('MSE Complete Compression: %f\n', mse1);
 
% SNR calculation for the foreground and foreground
snrForeground = calculate_snr(foreground, highQualityForeground);
fprintf('\nSNR Foreground: %f\n', snrForeground);
snrBackground = calculate_snr(background, lowQualityBackground);
fprintf('SNR Background: %f\n', snrBackground);

snr1 = calculate_snr(orig, complete_compressed);
snr2 = calculate_snr(orig, selective_compressed);
fprintf('SNR Selective Compression: %f\n', snr2);
fprintf('SNR Complete Compression: %f\n', snr1);
 
% Compare file sizes
fileInfoOriginal = dir(img);
fileInfoForeground = dir(fullfile(outputDir, 'foreground_high_quality.jpg'));
fileInfoBackground = dir(fullfile(outputDir, 'background_low_quality.jpg'));
fileInfoCombined = dir(fullfile(outputDir, 'combined_selective_compression.jpg'));
fileInfoComplete = dir('complete_compressed_img.jpg');
 
fprintf('\nOriginal Image Size: %.2f MB\n', fileInfoOriginal.bytes / (1024 * 1024));
fprintf('Foreground (High Quality) Image Size: %.2f MB\n', fileInfoForeground.bytes / (1024 * 1024));
fprintf('Background (Low Quality) Image Size: %.2f MB\n', fileInfoBackground.bytes / (1024 * 1024));
fprintf('Selectively Compressed Image Size (JPG): %.2f MB\n', fileInfoCombined.bytes / (1024 * 1024));
fprintf('Complete Compressed Image Size (JPG): %.2f MB\n', fileInfoComplete.bytes / (1024 * 1024));

uncompressedSize = fileInfoOriginal.bytes;
selectiveCompressedSize = fileInfoCombined.bytes;
completeCompressedSize = fileInfoComplete.bytes;
selectiveCompressionRatio = uncompressedSize / selectiveCompressedSize;
completeCompressionRatio = uncompressedSize / completeCompressedSize;
fprintf('\nCompression Ratio (selective): %.2f\n', selectiveCompressionRatio);
fprintf('Compression Ratio (complete): %.2f\n', completeCompressionRatio);

function snrValue = calculate_snr(originalImage, noisyImage)
    % Convert images to double precision to avoid overflow
    originalImage = double(originalImage);
    noisyImage = double(noisyImage);
 
    % Compute the signal power (sum of squares of the original image)
    signalPower = sum(originalImage(:).^2);
 
    % Compute the noise (difference between the original and noisy/compressed image)
    noise = originalImage - noisyImage;
 
    % Compute the noise power (sum of squares of the noise)
    noisePower = sum(noise(:).^2);
 
    % Calculate the SNR (in dB)
    snrValue = 10 * log10(signalPower / noisePower);
end
