clc
clear all
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

img = imread(image);

redChannel = img(:, :, 1); % Separate the RGB channels
greenChannel = img(:, :, 2);
blueChannel = img(:, :, 3);

% Otsus thresholding on each RGB channel independently
redThreshold = graythresh(redChannel);  
greenThreshold = graythresh(greenChannel);  
blueThreshold = graythresh(blueChannel);  

figure; % Plot the histogram
subplot(2,2,1), imhist(redChannel), hold on;      
xline(redThreshold*255, 'r', 'LineWidth', 2), hold off;
title('Histogram of Red');
subplot(2,2,2), imhist(greenChannel), hold on;          
xline(greenThreshold*255, 'g', 'LineWidth', 2), hold off;
title('Histogram of Green');
subplot(2,2,3), imhist(blueChannel), hold on;          
xline(blueThreshold*255, 'b', 'LineWidth', 2), hold off;
title('Histogram of Blue');
subplot(2,2,4), imhist(rgb2gray(img)), hold on;          
xline(blueThreshold*255, 'b', 'LineWidth', 2);
xline(greenThreshold*255, 'g', 'LineWidth', 2);
xline(redThreshold*255, 'r', 'LineWidth', 2); hold off;
title('Histogram Grayscale with rgb thresholds');

% Create binary masks for each channel based on the thresholds
redMask = imbinarize(redChannel, redThreshold);   % Binary mask for red
greenMask = imbinarize(greenChannel, greenThreshold);  % For green
blueMask = imbinarize(blueChannel, blueThreshold);   % Binary mask for blue

% Clean the masks using morphological operations
redMask = imopen(redMask, strel('disk', 5));  % Open to remove noise
redMask = imfill(redMask, 'holes');           % Fill holes
greenMask = imopen(greenMask, strel('disk', 5));  % Open to remove noise
greenMask = imfill(greenMask, 'holes');           % Fill holes
blueMask = imopen(blueMask, strel('disk', 5));  % Open to remove noise
blueMask = imfill(blueMask, 'holes');           % Fill holes

% Complement masks for the foreground
redCompMask = imcomplement(redMask);
greenCompMask = imcomplement(greenMask);
blueCompMask = imcomplement(blueMask);

% Combine the masks (for representation only)
combinedMask = redCompMask | greenCompMask | blueCompMask;

% Display the individual RGB masks and the concatenated mask
figure;
subplot(2, 2, 1), imshow(redCompMask), title('Red Mask');
subplot(2, 2, 2), imshow(greenCompMask), title('Green Mask');
subplot(2, 2, 3), imshow(blueCompMask), title('Blue Mask');
subplot(2, 2, 4), imshow(combinedMask), title('Overlapping Mask (for Representation)');

% Apply the respective masks to each color channel
foreground_red = uint8(redCompMask) .* redChannel; 
foreground_green = uint8(greenCompMask) .* greenChannel;  
foreground_blue = uint8(blueCompMask) .* blueChannel;

% Combine the masked RGB channels to get the foreground
foreground = cat(3, foreground_red, foreground_green, foreground_blue);

% Apply complement mask for the background
background_red = uint8(redMask) .* redChannel; % Background red channel
background_green = uint8(greenMask) .* greenChannel; % Background green channel
background_blue = uint8(blueMask) .* blueChannel; % Background blue channel

% Combine the masked RGB channels to get the background
background = cat(3, background_red, background_green, background_blue);

% Save the foreground with high-quality compression
imwrite(foreground, 'foreground_high_quality.jpg', 'jpg', 'Quality', foregroundQuality);

% Save the background with low-quality compression
imwrite(background, 'background_low_quality.jpg', 'jpg', 'Quality', backgroundQuality);

% Read the compressed images back
highQualityForeground = imread('foreground_high_quality.jpg');
lowQualityBackground = imread('background_low_quality.jpg');

% Combine the high-quality foreground and low-quality background
final_combined_image = highQualityForeground + lowQualityBackground;

% Save the final combined image
imwrite(final_combined_image, 'Final_combined_Image.jpg', 'jpg', 'Quality', overallQuality);

% Complete compression
imwrite(img, 'complete_compressed_img.jpg', 'jpg', 'Quality', completeCompressionQuality);

% Display the results
figure;
subplot(1,2,1); imshow(highQualityForeground); title('Foreground (High Quality)');
subplot(1,2,2); imshow(lowQualityBackground); title('Background (Low Quality)');

% Complete compression
imwrite(img, 'complete_compressed_img.jpg', 'jpg', 'Quality', completeCompressionQuality);
orig = imread(image);
complete_compressed = imread('complete_compressed_img.jpg');
selective_compressed = imread('Final_combined_Image.jpg');

% Display the final results
figure;
subplot(1,3,1); imshow(img); title('Original Image');
subplot(1,3,2); imshow(final_combined_image); title('Final Selective Compressed Image (jpg, Unified Mask)');
subplot(1,3,3); imshow(complete_compressed); title('Final Complete Compressed Image (jpg)');

% Calculating performance metrics
mse1 = immse(orig, complete_compressed);
mse2 = immse(orig, selective_compressed);
fprintf('MSE Complete Compression: %f\n', mse1);
fprintf('MSE Selective Compression: %f\n', mse2);

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

snr1 = calculate_snr(orig, complete_compressed);
snr2 = calculate_snr(orig, selective_compressed);
fprintf('SNR Complete Compression: %f\n', snr1);
fprintf('SNR Selective Compression: %f\n', snr2);

% Compare file sizes
fileInfoOriginal = dir(image);
fileInfoForeground = dir('foreground_high_quality.jpg');
fileInfoBackground = dir('background_low_quality.jpg');
fileInfoSelective = dir('Final_combined_Image.jpg');
fileInfoComplete = dir('complete_compressed_img.jpg');

fprintf('Original Image Size: %.2f MB\n', fileInfoOriginal.bytes / (1024 * 1024));
fprintf('Foreground (High Quality) Image Size: %.2f MB\n', fileInfoForeground.bytes / (1024 * 1024));
fprintf('Background (Low Quality) Image Size: %.2f MB\n', fileInfoBackground.bytes / (1024 * 1024));
fprintf('Selectively Compressed Image Size (JPG): %.2f MB\n', fileInfoSelective.bytes / (1024 * 1024));
fprintf('Completely Compressed Image Size (JPG): %.2f MB\n', fileInfoComplete.bytes / (1024 * 1024));

uncompressedSize = fileInfoOriginal.bytes;
selectiveCompressedSize = fileInfoSelective.bytes;
completeCompressedSize = fileInfoComplete.bytes;

% Calculate compression ratio
selectiveCompressionRatio = uncompressedSize / selectiveCompressedSize;
completeCompressionRatio = uncompressedSize / completeCompressedSize;

% Display the result
fprintf('Compression Ratio (selective): %.2f\n', selectiveCompressionRatio);
fprintf('Compression Ratio (complete): %.2f\n', completeCompressionRatio);