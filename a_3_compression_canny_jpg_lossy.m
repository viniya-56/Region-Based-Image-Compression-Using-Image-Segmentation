clc;
clear all;
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
       image = 'desert.bmp';
   case 2
       image = 'complex_img.bmp';
   otherwise
       error('Error. Select from the available choices.');
end
outputDir = 'Canny JPG Images';

img = imread(image);
img_gray = rgb2gray(img);

figure;
subplot(1,2,1), imshow(img), title("Original Image");
subplot(1,2,2), imshow(img_gray), title("Grayscale Image");
figure;
subplot(1,2,1), imhist(img_gray), title('Histogram of Gray Image');

% Edge detection using Canny to find the ROI
edges = edge(img_gray, "canny");
dilatedEdges = imdilate(edges, strel('disk', 5));   % Dilate the edges
roi_mask = imfill(dilatedEdges, 'holes');   % Fill the holes
subplot(1,2,2), imshow(roi_mask), title('ROI mask using Canny Edge Detection');

figure;
subplot(2,2,1), imshow(img), title("Original Image");
subplot(2,2,2), imshow(edges), title("Canny Edges");
subplot(2,2,3), imshow(dilatedEdges), title("Dilated Edges");
subplot(2,2,4), imshow(roi_mask), title("ROI Mask (Filled)");

% Separate the foreground (ROI) and background
mask = roi_mask; background_mask = ~mask;
img = im2uint8(img);

% Extract foreground using logical indexing
foreground = img;
foreground(~cat(3, mask, mask, mask)) = 0; % Set non-ROI pixels to 0 (black)

% Extract background
background = img;
background(~cat(3, background_mask, background_mask, background_mask)) = 0; % Set non-ROI pixels to 0 (black)

% Compressing the foreground with higher quality (lower compression)
imwrite(foreground, fullfile(outputDir, 'foreground_high_quality.jpg'), 'Quality', foregroundQuality); % Higher quality

% Compressing the background with lower quality (higher compression)
imwrite(background, fullfile(outputDir, 'background_low_quality.jpg'), 'Quality', backgroundQuality); % Lower quality
highQualityForeground = imread(fullfile(outputDir, 'foreground_high_quality.jpg'));
lowQualityBackground = imread(fullfile(outputDir, 'background_low_quality.jpg'));

% Combine the foreground and background for the final result
combinedImage = foreground + background;
imwrite(combinedImage, fullfile(outputDir, 'combined_selective_compression.jpg'), 'Quality', overallQuality);

% Complete compression
imwrite(img, fullfile(outputDir, 'complete_compressed_img.jpg'), 'jpg', 'Quality', completeCompressionQuality);

orig = imread(image);
complete_compressed = imread(fullfile(outputDir, 'complete_compressed_img.jpg'));
selective_compressed = imread(fullfile(outputDir, 'combined_selective_compression.jpg'));

figure;
subplot(1, 3, 1), imshow(img), title('Original Image');
subplot(1, 3, 2), imshow(combinedImage), title('Final Compressed Image (jpg, canny)');
subplot(1, 3, 3), imshow(complete_compressed), title('Complete Compressed Image (jpg)');

% Calculating performance metrics
mse1 = immse(orig, complete_compressed); 
mse2 = immse(orig, selective_compressed); 
fprintf('\nMSE Selective Compression: %f\n', mse2); 
fprintf('MSE Complete Compression: %f\n', mse1); 
 
snr1 = calculate_snr(orig, complete_compressed); 
snr2 = calculate_snr(orig, selective_compressed); 
fprintf('\nSNR Selective Compression: %f\n', snr2); 
fprintf('SNR Complete Compression: %f\n', snr1); 
 
% Comparing file sizes 
fileInfoOriginal = dir(image);
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
fprintf('Compression Ratio (selective): %.2f\n', selectiveCompressionRatio);
fprintf('Compression Ratio (complete): %.2f\n', completeCompressionRatio);


function snrValue = calculate_snr(originalImage, noisyImage) 
    % Convert images to double precision to avoid overflow 
    originalImage = double(originalImage); 
    noisyImage = double(noisyImage); 
 
    % Compute the signal power (sum of squares of the original image pixels) 
    signalPower = sum(originalImage(:).^2); 
 
    % Compute the noise (difference between the original and noisy/compressed image) 
    noise = originalImage - noisyImage; 
 
    % Compute the noise power (sum of squares of the noise) 
    noisePower = sum(noise(:).^2); 
 
    % Calculate the SNR (in dB) 
    snrValue = 10 * log10(signalPower / noisePower); 
end 