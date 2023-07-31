# Cosmic Ray Filter Code

clean_cr.py is for the use of the cleaning out cosmic rays from dark calibration images.

## Download ```clean_cr.py```

## In terminal ```open clean_cr.py```

In the ```MAIN``` of the code
- Replace ```in_filename``` with the name of the fits file that you would like to correct
- Replace ```out_filename``` with the name of the desired name of file once corrections are done

Save inputs with ```cmd + s``` (_MACOS_) or ```cntrl + s``` (_WINDOWS_)

## Execute code:
```python clean_cr.py``` in terminal

## Upon execution:
Code algorithm:
- Open the fits files and obtain detectors dimensions
    - Length, width, image frames
- Iterate through each pixel within the detector
    - Compute the difference between image frames
    - Calculate the robust standard deviation
        - _If a difference value is more than the calculated robust standard deviation, the pixel will be flagged for having a cosmic ray_
- **If a pixel is flagged for having a cosmic ray**:
    - Line fitting via linear regression is applied before and after the cosmic ray jump
    - An offset is calculated from the before and after fitted slopes and the offset is applied to the affected section of the pixel
    - Corrected pixel data gets saved into a corrected output file at the same coordinates within detector
- **If a pixel is not flagged, it will be saved regardless into a corrected output file at the same coordinates within detector**

## In terminal
Once the code is executed and begins correcting, the following will be printed:
- The pixel coordinates of the pixel that has been flagged for having a cosmic ray
- The image frames within the pixel where the cosmic ray was detected
- The amplitude of the jump of the cosmic ray
- Once code is finished, computed run time will print as _minutes_.

## Output
- A fits file of the corrected pixel 
- A npz file containing the information printed in terminal + more (_coming soon_)
