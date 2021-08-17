# ROI-analysis
ROI analysis package for ISS and INO image cubes
![image](https://user-images.githubusercontent.com/84033812/129761866-83b31096-4b1d-437f-8a8f-1ffb65cb51c3.png)

File descriptions:
- ROI_GUI.py
  - pyqt5 gui functions/classes
-ROI_functions_images
  - fluorescence_functions.py
    - image file handling, data analysis functions
  - bad_path.png, none.png
    - background images
- requirements.txt
  - required packages
  
# Requirements:
- python (>3.8)
- pyqt5
- matplotlib
- numpy
- qtrangeslider
- tifffile
- imagecodecs
- scipy
- scikit-image
- opencv-python
- pandas

# Running the package
To install the requirements: ensure that you already have python and pip installed, then run

    pip install -r requirements.txt

Next, ensure that `ROI_GUI.py` and the `ROI_functions_images` folder are both in the same working directory. 
Running `ROI_GUI.py` will open the GUI and allow access to the ROI analysis functionalities.
