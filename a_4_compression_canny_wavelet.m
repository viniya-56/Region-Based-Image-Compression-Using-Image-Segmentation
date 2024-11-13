clc; 
clear all; 
close all; 

% foregroundRatio is lossless, when unspecified
backgroundRatio = 90;
completeCompressionRatio = 90;

disp('Available image choices:');
disp('1. Desert Image (bmp)');
disp('2. Horse Image (jpg)');

image_choice = input("Enter the image choice (digits 1-2): ");
switch image_choice
    case 1
        image = 'desert.bmp';
    case 2
        image = 'complex_img.bmp';
    otherwise
        error('Error. Select from the available choices.');
end

outputDir = 'Canny Wavelet Images';
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

img = imread(image); 
img_gray = rgb2gray(img); 

figure;
subplot(1,2,1), imshow(img), title("Original Image");
subplot(1,2,2), imshow(img_gray), title("Grayscale Image");

% histogram
figure;
subplot(1,2,1), imhist(img_gray), title('Histogram of Grayscale Image');

% Edge detection using Canny to find the ROI (Region of Interest) 
edges = edge(img_gray, "canny"); 
dilatedEdges = imdilate(edges, strel('disk', 5)); % Dilate the edges 
roi_mask = imfill(dilatedEdges, 'holes'); % Fill the holes 

subplot(1,2,2), imshow(roi_mask), title('ROI mask using Canny Edge Detection');

figure;
subplot(2,2,1), imshow(img), title("Original Image"); 
subplot(2,2,2), imshow(edges), title("Canny Edges");
subplot(2,2,3), imshow(dilatedEdges), title("Dilated Edges");
subplot(2,2,4), imshow(roi_mask), title("ROI Mask (Filled)");


% Separate the foreground (ROI) and background 
mask = roi_mask; 
background_mask = ~mask; 

% Ensure img is in uint8 format (common for images) 
img = im2uint8(img);

% Extract foreground using logical indexing 
foreground = img; 
foreground(~cat(3, mask, mask, mask)) = 0; % Set non-ROI pixels to 0 (black) 

% Extract background 
background = img; 
background(~cat(3, background_mask, background_mask, background_mask)) = 0; % Set non-ROI pixels to 0 (black) 

figure; 
subplot(1, 2, 1), imshow(foreground), title('Foreground (High Quality)'); 
subplot(1, 2, 2), imshow(background), title('Background (Low Quality)'); 

% Combine the foreground and background for the final result
final_combined = foreground + background; 

% Save images in JP2 formats with lossy compression
imwrite(foreground, fullfile(outputDir, 'foreground_lossless.jp2')); 
imwrite(background, fullfile(outputDir, 'background_lossy.jp2'), 'CompressionRatio', backgroundRatio); % 90 times smaller

% Save the final result in JP2 formats
imwrite(final_combined, fullfile(outputDir, 'Final_Combined_Image.jp2')); 

% Complete compression
imwrite(img, 'complete_compressed_img.jp2', 'CompressionRatio', completeCompressionRatio);

orig = imread(image);
complete_compressed = imread('complete_compressed_img.jp2');
selective_compressed = imread(fullfile(outputDir, 'Final_Combined_Image.jp2'));

figure; 
subplot(1,3,1), imshow(img), title("Original Image"); 
subplot(1,3,2), imshow(final_combined), title("Selectively Compressed Image");
subplot(1,3,3), imshow(complete_compressed), title("Completely Compressed Image");

% Calculating performance metrics
mse1 = immse(orig, complete_compressed); 
mse2 = immse(orig, selective_compressed); 
fprintf('\nMSE Selective Compression: %f\n', mse2); 
fprintf('MSE Complete Compression: %f\n', mse1); 

snr1 = calculate_snr(orig, complete_compressed); 
snr2 = calculate_snr(orig, selective_compressed); 
fprintf('\nSNR Selective Compression: %f\n', snr2); 
fprintf('SNR Complete Compression: %f\n', snr1); 

fileInfoOriginal = dir(image); 
fileInfoForeground = dir(fullfile(outputDir, 'foreground_lossless.jp2')); 
fileInfoBackground = dir(fullfile(outputDir, 'background_lossy.jp2')); 
fileInfoCombinedJP2 = dir(fullfile(outputDir, 'Final_Combined_Image.jp2')); 
fileInfoCompleteCompressedJP2 = dir('complete_compressed_img.jp2'); 

fprintf('\nOriginal Image Size: %.2f MB\n', fileInfoOriginal.bytes / (1024 * 1024));
fprintf('Foreground (High Quality) Image Size: %.2f MB\n', fileInfoForeground.bytes / (1024 * 1024));
fprintf('Background (Low Quality) Image Size: %.2f MB\n', fileInfoBackground.bytes / (1024 * 1024));
fprintf('Selectively Compressed Image Size (JP2): %.2f MB\n', fileInfoCombinedJP2.bytes / (1024 * 1024));
fprintf('Complete Compressed Image Size (JP2): %.2f MB\n', fileInfoCompleteCompressedJP2.bytes / (1024 * 1024));

uncompressedSize = fileInfoOriginal.bytes;
selectiveCompressedSize = fileInfoCombinedJP2.bytes;
completeCompressedSize = fileInfoCompleteCompressedJP2.bytes;
% Calculating compression ratio
selectiveCompressionRatio = uncompressedSize / selectiveCompressedSize;
completeCompressionRatio = uncompressedSize / completeCompressedSize;
% Display the result
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

