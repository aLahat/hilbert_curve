Hilbert curve is a space filling fractal which has some intrinsic clustering properties [16,17], meaning that regions positioned closer to each other are very to retain localicity in any N-dimensional iteration of the Hilbert curve. 

This property can be exploited to plot linear information, and by using multiple RGB channels we can visually observe up to three large datasets simultaneously. 
Code for the generation of Hilbert curves was forked from Galtay [16]. The original code was capable of producing N-dimensional Hilbert shapes. Here, we are only interested in the use of 2-d Hilbert shapes. 
Results
In order to generate a Hilbert plot we need to import a couple of modules; a plotting module such as matplotlib and the HilbertColor function from the Hilbert package.
```
import matplotlib.pyplot as plt
from hilbert import HilberColor

plt.imshow(         # plots an image based on a 2d or 3d (RGB) matrix.
    HilberColor(    # turns a linear numerical list into a 3d (RGB) matrix.
        range(4**5) # generates a list of numbers ranging from 0 to 4^5
    )
)
plt.axis('off') # removes the axis annotations
plt.show()      # Displays plot
```
This results in a Hilbert fractal being drawn containing a grayscale from white to black.

For a more biologically relevant Hilbert plot, data from epigenetic markers in bigwig files can be read and plotted:
```
import matplotlib.pyplot as plt
import pyBigWig # needed to read bigwig files
from hilbert import HilberColor

N = 4**6	# this will be the number of bins. 
		# If the number of bins is a 4x, it will utilize all the available canvas.

CHROMS = 	[
	'chr1','chr2','chr3','chr4','chr5',
	'chr6','chr7','chr8','chr9','chr10',
	'chr11','chr12','chr13','chr14','chr15',
	'chr16','chr17','chr18','chr19','chr20',
	'chr21','chr22','chrX','chrY'
	]
f, axarr = plt.subplots(len(CHROMS),4,
                    	figsize = (2*4,len(CHROMS)*2),
                    	)	# generates a figures with four subfigures for each chromosome
				# three for showing individual channels and one for showing the 
				# merged channels
repliTime	= pyBigWig.open('./ENCFF000KUF.bigWig')
H3K27me3	= pyBigWig.open('./ENCFF832TSN.bigWig')
H3K4me3 	= pyBigWig.open('./ENCFF372YWG.bigWig')

for axes,CHR in zip(axarr,CHROMS):
	A,B,C,D=axes
	[ax.get_xaxis().set_ticks([]) for ax in [A,B,C,D]] # removes x ticks
	[ax.get_yaxis().set_ticks([]) for ax in [A,B,C,D]] # removes y ticks
    
	A.set_ylabel(CHR)	# sets the chromosome name on the left of the first row
	A.imshow(HilberColor( repliTime.stats(CHR,nBins = N) ,[],[]))	# Uses the red channel for RepliChip
										# Empty lists force the RGB channels
										# otherwise it would interpret a grayscale figure
	B.imshow(HilberColor( [],H3K27me3.stats(CHR,nBins = N),[] ))	# Green channel for H3K27me3
	C.imshow(HilberColor( [],[],H3K4me3.stats(CHR,nBins = N) ))	# Blue channel for H3K4me3
	D.imshow(HilberColor( repliTime.stats(CHR,nBins = N),
                      	H3K27me3.stats(CHR,nBins = N),
                      	H3K4me3.stats(CHR,nBins = N) ))	# Uses all the channels this time
	if CHR =='chr1': # only writes titles on the top most subfigures
    		A.set_title('RepliChip')
    		B.set_title('H3K27me3')
    		C.set_title('H3K4me3')
    		D.set_title('Combined')
plt.subplots_adjust(wspace=0.05, hspace=0.05)
plt.show()
```

![cropped example](/hilbert_cropped.png)

As Hilbert curves naturally cluster nearby regions, we can visualize distinct early/late regions in the chromosome. 
In this example we can observe that strong H3K27me3 (enhancer tag) and H3K4me (promoter tag) do not overlap (promoter and enhancer are distant to each other), but are mostly located in replicative early regions.

