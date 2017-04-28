// start with an open image

roiManager("Reset");
setTool("rectangle");
waitForUser("Select an ROI","parse the video and select the region to be analyzed (such as a protrusion), \n Note the first and last frames of interest \n, Click OK when the Roi is drawn on the video");
fullimageID=getImageID();
title = "Sub time interval to process";
Dialog.create(title);
Dialog.addNumber("First Frame:", 1);
Dialog.addNumber("Last frame:", nSlices);
Dialog.show();
 
ff = Dialog.getNumber();
lf = Dialog.getNumber();;
 
// Create the crop, smooth it
run("Duplicate...", "duplicate range="+ff+"-"+lf);
cropID=getImageID();
selectImage(fullimageID);
close();
selectImage(cropID);
run("Duplicate...", "duplicate");
cropMaskID=getImageID();
run("Smooth", "stack");
nbexpectedRois=nSlices;
//let the user select the Threshold

setAutoThreshold("Default dark");
run("Threshold...");
waitForUser("Select an ROI","Threshold the protrusions. \n Holes and some small additional pixels are acceptable at this stage, mask will be further cleaned\n, Click OK when Threshold have been applied");
selectImage(cropMaskID);
run("Analyze Particles...", "size=100-Infinity include add stack");

selectImage(cropID);
nbrois=roiManager("Count");
if (nbrois!=nbexpectedRois){
	rois = newArray();
	slices = newArray();
	for (i=0;i<nbrois;i++){
		roiManager("Select",i);
		rois=Array.concat(rois,i);
		slices=Array.concat(slices,getSliceNumber());
		
	}	
}
if (nbrois<nbexpectedRois){
	missingSlices=createmissingSlicesArray(slices,nbexpectedRois);
	Array.show(missingSlices);
	selectImage(cropMaskID);
	
	roiManager("Show None");
	run("Select None");
	roiManager("Reset");
	waitForUser("Not enough Rois, check your segmentation, missing frames are listed. \n"+"Edit the mask \n click OK when done");
	selectImage(cropMaskID);
	run("Analyze Particles...", "size=100-Infinity include add stack");
	
}
if (nbrois>nbexpectedRois){
	duplicatedSlices=createduplicatedSlicesArray(slices,nbexpectedRois);
	Array.show(duplicatedSlices);
	selectImage(cropMaskID);
	
	roiManager("Show None");
	run("Select None");
	roiManager("Reset");
	waitForUser("Too much Rois, check your segmentation, frames with problem are listed. \n"+"Edit the mask \n click OK when done");
	selectImage(cropMaskID);
	run("Analyze Particles...", "size=100-Infinity include add stack");
}

selectImage(cropMaskID);
close();
selectImage(cropID);
roiManager("Show All");

waitForUser("You can now launch Recrutment Edge Analysis");

function createmissingSlicesArray(slices,nbexpectedRois) {
	count=newArray();
	 output = newArray();
	for (s=1;s<=nbexpectedRois;s++)
		count=Array.concat(count,0);	
	for (s=1;s<nbexpectedRois;s++){
		
		for (i=0;i<slices.length; i++){
			if (slices[i]==s)
				count[s-1]=count[s-1]+1;
		}
	}
	for (s=1;s<nbexpectedRois;s++){
		if (count[s-1]<1)
			output=Array.concat(output,s);
	}
     
      return output;
   }
  // in a future version, could correct it automatically based on area (keeping the biggest for each slide) 
function createduplicatedSlicesArray(slices,nbexpectedRois) {
    count=newArray();
	 output = newArray();
	for (s=1;s<=nbexpectedRois;s++)
		count=Array.concat(count,0);	
		
	for (s=1;s<nbexpectedRois;s++){
		
		for (i=0;i<slices.length; i++){
			if (slices[i]==s)
				count[s-1]=count[s-1]+1;
		}
	}
	for (s=1;s<nbexpectedRois;s++){
		if (count[s-1]>1){
			
			
			output=Array.concat(output,s);
		}
	}
     
      return output;
   }