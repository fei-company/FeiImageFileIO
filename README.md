$ cd <path>
$ ./FeiSer2Img.py --help
ser2img.py: convert .ser (TIA format) to .img
options:
-i, --input 	<inputfile>
-o, --output 	<outputfile>
-x, --cx 	<center x(px)> (required for images with a beam stopper, otional for images without a beam stopper)
-y, --cy 	<center y(px)> (required for images with a beam stopper, otional for images without a beam stopper)
-k,- -HT 	<voltage (kV)>
-c, --camera 	<camera length (m)>
-a,--osc 	<rotation angle per frame (deg)>
-s, --sigma 	<sigma of Gaussian filter> (default 2.0)
-g, --gain 	<multiply the data by the provided gain> (default 2.0)
-b, --bias_correction	 (Optional processing) Bias correction, sometimes the bias drifted during the exposure and you observe strips in the image, this may help
-d, --beam_centering	 (Optional processing) align the beam center, sometimes the beam center drifted during the exposure, if it is not blocked by the beam stopper and not saturated, this simple routine can fix it
-l, --lowpass_filtering	(Optional processing) low-pass filtering, may help in spots finding and indexing)

$ ./FeiSer2Img.py -i ~/Projects/MED/GV/gv23_1.ser -o ~/Projects/MED/GV/gv23/gv23 -k 200 -c 4.3 -a 0.5 
Input .ser file: /home/lyu/Projects/MED/GV/gv23_1.ser
	Number of images: 50
Rounding to unsigned int (gain=2.0)
	Applying pixel value offset to all images of 132.74
	Scaling all images by  2.00
Finding Beam center: 979.94, 1027.75
Saving to /home/lyu/Projects/MED/GV/gv23/gv23###.img

$ ./FeiSer2Img.py -i ~/Projects/MED/GV/gv23_1.ser -o ~/Projects/MED/GV/gv23/gv23 -k 200 -c 4.3 -a 0.5 -d
Input .ser file: /home/lyu/Projects/MED/GV/gv23_1.ser
	Number of images: 50
Align beam center
Rounding to unsigned int (gain=2.0)
	Applying pixel value offset to all images of 132.74
	Scaling all images by  2.00
Finding Beam center: 985.45, 1028.13

$ ./FeiSer2Img.py -i ~/Projects/MED/GV/gv23_1.ser -o ~/Projects/MED/GV/gv23/gv23 -k 200 -c 4.3 -a 0.5 -d -l
Input .ser file: /home/lyu/Projects/MED/GV/gv23_1.ser
	Number of images: 50
Align beam center
Low-pass filter (Gaussian kernel, sigma=2.0)
Rounding to unsigned int (gain=2.0)
	Applying pixel value offset to all images of 41.08
	Scaling all images by  2.00
Finding Beam center: 985.45, 1028.13
Saving to /home/lyu/Projects/MED/GV/gv23/gv23###.img
