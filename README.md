# Region-Based Image Compression Using Image Segmentation

### Overview
This project implements a region-based image compression system utilizing image segmentation techniques such as thresholding and Canny edge detection. The primary aim is to enhance image compression by preserving high-quality details in specified Regions of Interest (ROI),  while applying greater compression to less important areas, thus optimizing storage and bandwidth.

### Features
- Segmentation Techniques: Implements thresholding and Canny edge detection to accurately identify ROIs in images.
- Compression Methods: Utilizes JPEG </h4>(Discrete Cosine Transform) and wavelet-based compression algorithms to achieve efficient data storage and preservation of critical image details.
- Applications: Suitable for use cases such as medical imaging, satellite imagery, surveillance systems, and more.

### Technologies Used
- MATLAB: For image processing, segmentation, and compression.
- Image Segmentation Techniques: Thresholding, Canny edge detection.
- Compression Algorithms: JPEG (DCT) and wavelet-based compression.

### Installation and Usage
1. Clone the repository:
   ```bash
   git clone https://github.com/viniya-56/Region-Based-Image-Compression-Using-Image-Segmentation.git

Open the project files in MATLAB.
Load an image into the MATLAB environment.
Run the main script to perform region-based image compression.
Apply selected compression techniques to the image, preserving ROI quality.