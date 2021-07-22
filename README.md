# LGN

The prerequisites for the image processing are installed IgorPro and ITK-SNAP software.

1. Open command line window (macOS), go to the folder with NIFTI MRI images, set a path for c3d plugin in ITK-SNAP software (e.g., export PATH="/Applications/ITK-SNAP.app/Contents/bin/":$PATH). Copy shell script 3D_to_1D to this folder, and grant it executive privilege with command chmod.
2. Open MRI image in ITK-SNAP. Identify the center of LGN (in a sagittal slice). Run the script 3D_to_1D with five arguments. The 1st, 2nd, and 3rd arguments are the x-, y-, and z-coordinates of the LGN center. The 4th argument is the name of the image file without extension ".nii.gz". The 5th argument is either R (right LGN), or L (left LGN). For example: ./3D_to_1D  120 120 120 MRimage R
The outputs of the shell script are the text file MRimage_R.txt and 22x22x22 excise of LGN containing image.
3. Open IgorPro experiment EdgeEnhancement.pxp in IgorPro. A permanent link for downloading IgorPro experiment EdgeEnhancement.pxp is here:
https://drive.google.com/drive/folders/1zTe98Gw6bCtnQUJpOLICDqUQEz1Rieec?usp=sharing
In menu Data => Load Waves => Load Delimited Text   load file MRimage_R.txt and save it as MRimage_R
4. In IgorPro command line, run Ar2Im(MRimage_R). The output will be the 3D 22x22x22 image with the name of MRimage_R_Im. Type this name into the table ImName. The table ImName may contain as many 3D images as you want to run at once in the next step.
5. In IgorPro command line, run Batch(SL,Nmeans,Niterations), where the first argument SL is either 0 (Shorter processing time, preferred) or 1 (Longer processing time). The second argument Nmeans is 1 for shapes with vertices that can be represented by 3 crossing planes (e.g., cube or triangular pyramid), and 4 to 20 for shapes with vertices that cannot be represented by 3 crossing planes (e.g., sphere, or LGN). The third argument is the number of consecutive filtrations (6 for images with moderate noise level, 18 to 24 for images with high noise level). The processing takes ~50 min/iteration/image. For example, to reproduce the image Cube_Im_1i in the root directory of IgorPro data browser, execute Batch(0,1,1) ("1i" at the end of the name Cube_Im_1i stands for "1 iteration").
6. When the previous step is finished, run Im2Ar(ImageName) in IgorPro command line, where ImageName is the name of the edge-enhanced image. Save the output in the folder with original MR image.
7. In the macOS command line, run c3d x.nii.gz -region 0x0x0vox 22x22x22vox -resample 200x200x200% -scale 0 -landmarks-to-spheres xx.txt 0.2 -o xxx.nii.gz
where x.nii.gz is the image resulted from execution of 3D_to_1D in the step 2,
xx.txt is the name of the text image resulted from Im2Ar(ImageName) in IgorPro command line in the previous step,
xxx.nii.gz is the desired name of 44x44x44 LGN containing image representing edge-enhanced and upsampled version of the original LGN-containing image.
