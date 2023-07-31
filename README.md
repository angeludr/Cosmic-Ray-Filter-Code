# Cosmic Ray Filter Code

clean_cr.py is for the use of the cleaning out cosmic rays from dark calibration images.

## Download ```clean_cr.py```

## In terminal ```open clean_cr.py```

In the ```MAIN``` of the code
- Input the fits file you would like to correct for object ```infile```
- Input the desired name of the outfile for object ```outfile```

## Running filter
Once ```infile``` and ```outfile``` names are inserted

Start correcting by executing
        ```python clean_cr.py``` in terminal

## As code runs
Code algorithm:
- Open the fits files and obtain detectors dimensions
    - Length, width, image frames
- Iterate through each pixel within the detector
- <sub> Compute the difference between image frames <sub>
- <sub> Calculate the robust standard deviation <sub>
- If a difference value is more than the calculated robust standard deviation, it will be flagged for having a cosmic ray
- If a pixel is flagged for having a cosmic ray:
-   <sub> Line fitting via linear regression is applied before and after the cosmic ray jump <sub>
-   <sub> An offset is calculated from the before and after fitted slopes and the offset is applied to the affected section of the pixel <sub>
-   <sub> Corrected pixel data gets saved into a corrected output file at the same coordinates within detector <sub>
- If a pixel is not flagged, it will be saved regardless into a corrected output file at the same coordinates within detector

## In terminal
Once the code is executed and begins correcting, the following will be printed:
- <sub> The pixel coordinates of the pixel that has been flagged for having a cosmic ray <sub>
- <sub> The image frames within the pixel where the cosmic ray was detected <sub>
- <sub> The amplitude of the jump of the cosmic ray <sub>

## Output
- <sub> A fits file of the corrected pixel <sub>
- <sub> A npz file containing the information printed in terminal + more _(coming soon)_ <sub>
