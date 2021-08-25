# ROI-analysis
ROI analysis package for ISS and INO image cubes
![image](https://user-images.githubusercontent.com/84033812/129761866-83b31096-4b1d-437f-8a8f-1ffb65cb51c3.png)

File descriptions:
- ROI_GUI.py
  - pyqt5 gui functions/classes
- ROI_functions_images
  - fluorescence_functions.py
    - image file handling, data analysis functions
  - bad_path.png, none.png
    - background images
- requirements.txt
  - required packages
- filter_bin_data.ipynb
  - jupyter notebook for data filtering, binning and curve fitting
  
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

    > pip install -r requirements.txt

Next, ensure that `ROI_GUI.py` and the `ROI_functions_images` folder are both downloaded and extracted directly into the desired python environment. 
Running `ROI_GUI.py` will open the GUI and allow access to the ROI analysis functionalities.

If using Anaconda, it is best to run the script from the Anaconda command prompt with

    > python ROI_GUI.py
    
This ensures that any error messages that may occur are accessible if the GUI crashes.

The binning script can be run from the Anaconda command prompt with

    > jupyter notebook filter_bin_data.ipynb

To close the notebook, type `Ctrl+c` in the Anaconda command prompt.
