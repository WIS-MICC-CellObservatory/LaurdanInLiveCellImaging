/* ImageJ macro for GP images analysis
   Quantitative Imaging of Membrane Lipid Order in Cells and Organisms
   Dylan M. Owen, Carles Rentero, Astrid Magenau, Ahmed Abu-Siniyeh and Katharina Gaus
   Nature Protocols 2011
   version July 2011   

   Following Nature protocols paper: “Quantitative imaging of membrane lipid order in cells and organisms” , Owen et al, 2011. doi:10.1038/nprot.2011.419
   https://www.nature.com/articles/nprot.2011.419
   Adapting ImageJ macro of the Nature protocol (NP) paper, which is available at: https://media.nature.com/original/nature-assets/nprot/journal/v7/n1/extref/nprot.2011.419-S1.txt  
*/

/* This macro calculates GP images from Laurdan and di-4-ANEPPDHQ ratiometric images in
   bactch mode (whole chosen folder) obtained using a Leica microscope. The generation
   of HSB images of these GP images has been also implemented.
*/

/* 
 *  The original macro was minimally modified and tuned 
 *  by Ofra Golani, March 2019
 * 
 * Modifications: 
 * ==============
 * 	 
 * - GP calculation corrected (the original GP – which I left as is - gave “all zero results”), 
 * 	 the corrected calculation is stored in GP1 folder. Basically GP is a masked version of the PreGP image, 
 * 	 the correction is in the implementation of masking: instead of multiplying the PreGP image by a (1bit) mask and converting to 8 bit, 
 * 	 the mask is used to set the out-of-mask pixels to NaN to ensure they are discarded from the histogram calculation. 
 * 	 Setting the out-of-mask values to zero is not correct as zero is a valid value for PreGP/GP
 * 	 
 * - When GP channles are used for masking, Masking is always based on original (non-scaled) disordered file : (Ord+Disord)>2*Th
 * 	 
 * - User can select regions of interests (ROIs) of different types (eg cells / membrane), 
 * 	 ROIs are saved and histograms of GP values within those regions are calculated. 
 * 	 
 * - Results folder name is fixed and does not include date/time. 
 * 	 This way, one can use existing ROIs for repeated calculation.
 * 	 
 * - Histograms are saved for ALL (in GP1\Histograms) and for each ROI type in ROI subfolder
 * 
 * - Number of Bins in the Histograms can be tuned using DefaultNBins
 * - histograms are also saved as jpg graphs
 * 
 * - Pixel values are saved for each ROI in a table (default is one table for all the ROIs of an image), 
 * 	 for ordered channel, disordered channel, PreGP and GP1
 * 
 * See accompny documentation file for further instructions 
 */

// =============== Tunable Parameters ===================
var fileExtention = ".tif";

// Parameters for the First Dialog that controls GP calculation
CHANNEL1 = newArray("none","ch00","ch01","ch02","ch03","ch04","ch05");
//CHANNEL1 = newArray("none","00","01","02","03","04","05");
//CHANNEL1 = newArray("none","c001","c002","c003");

// Used for selecting ordered / disordered channels in the GUI
CHANNEL2 = newArray("ch00","ch01","ch02","ch03","ch04","ch05");
//CHANNEL2 = newArray("c001","c002","c003");
var QUESTION = newArray("Yes","No"); 
DefaultOrderedChannel = "ch00";
DefaultDisOrderedChannel = "ch01";
DefaultLowerThresholdGPMask = 15;
DefaultLowerThresholdIFMask = 50;
DefaultGFactor = 1;
DefaultUseGFactor = "Yes";
DefaultNBins = 40; //256; // number of Bins in the histograms
DefaultNormHistMax = 10; //2.5;  // yMax for the histogram plot
UseDateTimeInResultsFolderName = 0;
var DeafultSavePerROIPixelValuesInTable = true;
var DeafultValuesOfRoisInSeparateFiles = false; 	// 0 - values for all ROIs ofan image are saved in one file, 1 - values for each ROI are saved into independent file
var DeleteOldFilesFromResultsFolder = 1;

// Parameters for the 3rd Dialog, used to control ROI based analysis (in the function ROIanalysis)
var ROI_SELECTION_IMAGE = newArray("Order channel","Disorder channel","Order+Disorder","IF channel","HSB Image");
var ROI_TYPES = newArray("Membrane","Cell","Other", "none");
var ROI_COLORS = newArray("red","cyan","magenta");
//var QUESTION = newArray("Yes","No");
var DefaultNRoiTypes = 2;
var DrawRoisStartMessage = "The script will prompt you to select ROIs of selected type.\n It will prompt you to click OK after each ROI, and ask if you want to select another ROI of the same type..";
var DefaultDisplayImage = "Order+Disorder";

var RoiLineWidth = 2;
// =============== End of Tunable Parameters ===================

var year;
var months;
var dayOfMonths;
var hours;
var minutes;
var dir;
var results_Dir;
var ordered_images_Dir;
var disordered_images_Dir;
var GP_images_Dir;
var x;
var y;
var nBins;
var NormHistMax;  // yMax for the histogram plot
//var Folder;
//var s;
//var histoDir;
var chA;
var chB;

print("\\Clear");
run("Options...", "iterations=1 count=1 black do=Nothing");

// Initialization
requires("1.44d");
closeAllImages();
roiManager("Reset");

// Select images folder
dir = getDirectory("Choose a Directory of images to process");

/*// Inicialize choice variables
//CHANNEL1 = newArray("none","ch00","ch01","ch02","ch03","ch04","ch05");
//CHANNEL1 = newArray("none","00","01","02","03","04","05");
CHANNEL1 = newArray("none","c001","c002","c003");

// Used for selecting ordered / disordered channels in the GUI
//CHANNEL2 = newArray("ch00","ch01","ch02","ch03","ch04","ch05");
CHANNEL2 = newArray("c001","c002","c003");
QUESTION = newArray("Yes","No");*/
Lut_Dir = getDirectory("luts");
lut = getFileList(Lut_Dir);

// Choose image channels and threshold value
Dialog.create("GP analysis parameters");
//Dialog.addChoice("Acquisition ordered channel:  ", CHANNEL2, "ch00");
Dialog.addChoice("Acquisition ordered channel:  ", CHANNEL2, DefaultOrderedChannel);
//Dialog.addChoice("Acquisition disordered channel:  ", CHANNEL2, "ch01");
Dialog.addChoice("Acquisition disordered channel:  ", CHANNEL2, DefaultDisOrderedChannel);
//Dialog.addNumber("Lower Threshold Value for GP the mask:  ", 15);
Dialog.addNumber("Lower Threshold Value for GP the mask:  ", DefaultLowerThresholdGPMask);
Dialog.addChoice("Scale color for GP images:  ", lut, "Rainbow RGB.lut");
Dialog.addMessage("\n");
Dialog.addChoice("Immunofluorescence channel:  ", CHANNEL1, "none");
//Dialog.addNumber("Lower Threshold Value for the IF mask:  ", 50);
Dialog.addNumber("Lower Threshold Value for the IF mask:  ", DefaultLowerThresholdIFMask);
Dialog.addMessage("\n");
//Dialog.addNumber("G factor (1 if unknown):  ", 1);
Dialog.addNumber("G factor (1 if unknown):  ", DefaultGFactor);
//Dialog.addChoice("Do you want to use G factor for GP image calculation?",QUESTION, "No");
Dialog.addChoice("Do you want to use G factor for GP image calculation?",QUESTION, DefaultUseGFactor);
Dialog.addMessage("\n");
Dialog.addChoice("Do you want to generate HSB images?",QUESTION, "Yes");
Dialog.addMessage("\n");
Dialog.addNumber("Number of Histogram Bins:  ", DefaultNBins);
Dialog.addNumber("Maximum Hist plot Y Axis:  ", DefaultNormHistMax);
Dialog.addMessage("\n");
Dialog.show();
//waitForUser; // commented by OG 
// Feeding variables from dialog choices
chA = Dialog.getChoice();
chB = Dialog.getChoice();
t = Dialog.getNumber();
lut = Dialog.getChoice();
chC = Dialog.getChoice();
tc = Dialog.getNumber();
Gf = Dialog.getNumber();
Ques = Dialog.getChoice();
HSB = Dialog.getChoice();
nBins = Dialog.getNumber();
NormHistMax = Dialog.getNumber();

lut1= substring(lut,0,lengthOf(lut)-4);
suffixLen = lengthOf(chA)+lengthOf(fileExtention)+1;
time0 = getTime();
setBatchMode(true); 
//setBatchMode(false);

// Folder management
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
MonthNames = newArray("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
DayNames = newArray("Sun", "Mon","Tue","Wed","Thu","Fri","Sat");
if (hour<10) {hours = "0"+hour;}
else {hours=hour;}
if (minute<10) {minutes = "0"+minute;}
else {minutes=minute;}
if (month<10) {months = "0"+(month+1);}
else {months=month+1;}
if (dayOfMonth<10) {dayOfMonths = "0"+dayOfMonth;}
else {dayOfMonths=dayOfMonth;}

var Folder = 0;
if (UseDateTimeInResultsFolderName)
	results_Dir = dir + "Results "+year+"-"+months+"-"+dayOfMonths+" "+hours+"h"+minutes+ File.separator;
else
	results_Dir = dir + "Results"+ File.separator;

File.makeDirectory(results_Dir);

ordered_images_Dir = results_Dir + "Ordered Images" + File.separator;
Folder = ordered_images_Dir;
newFolder();
disordered_images_Dir = results_Dir + "Disordered Images" + File.separator;
Folder = disordered_images_Dir;
newFolder();
GP_images_Dir = results_Dir + "GP images" + File.separator;
Folder = GP_images_Dir;
newFolder();
histogramGP_Dir = GP_images_Dir + "Histograms" + File.separator;
Folder = histogramGP_Dir;
newFolder();

// GP1 added by OG - to Fix for GP calculation - convert low intensity pixels to NaN
GP1_images_Dir = results_Dir + "GP1 images" + File.separator;
Folder = GP1_images_Dir;
newFolder();
histogramGP1_Dir = GP1_images_Dir + "Histograms" + File.separator;
Folder = histogramGP1_Dir;
newFolder();

rawGP_images_Dir = results_Dir + "raw GP images" + File.separator;
Folder = rawGP_images_Dir;
newFolder();

IF_images_Dir = results_Dir + "Immunofluorescence Images" + File.separator;
GP_IF_images_Dir = results_Dir + "GP-IF images" + File.separator;
histogramIF_Dir = GP_IF_images_Dir + "Histograms" + File.separator;
HSB_Dir = results_Dir + "HSB images" + File.separator;

// Open and save images in 32bit format 
listDir = getFileList(dir);
var s = 0;
for (i = 0; i < listDir.length; i++) {
	//print(i, listDir[i], chA+".tif", chB+".tif", endsWith(listDir[i], chA+".tif"), endsWith(listDir[i], chB+".tif"));	
	print("OG: i, listDir[i]: ", i, listDir[i]);
	if (endsWith(listDir[i], chA+fileExtention)) {
		open(dir + listDir[i]);
		prepareImage();
		//print("A: s=",s," substring=",substring(s,0,lengthOf(s)-suffixLen));
		saveAs("tiff", ordered_images_Dir + substring(s,0,lengthOf(s)-suffixLen) + "_chA_32bits");
		close(); }
	if (endsWith(listDir[i], chB+fileExtention)) {
		open(dir + listDir[i]);
		prepareImage();
		Gf1=Gf;
		//print("B: s=",s," substring=",substring(s,0,lengthOf(s)-suffixLen));
		print("Gf1=",Gf1);		
		if (Ques=="Yes") {
			run("Multiply...","value=" + Gf);
			Gf1=1; }
		saveAs("tiff", disordered_images_Dir + substring(s,0,lengthOf(s)-suffixLen) + "_chB_32bits");
		close(); }
}

// OG: 
print("Gf1=",Gf1);

// GP and GPc Arrays 
/*GP=newArray(256);
for (i=0; i<256; i++) {
	GP[i]=((i-127)/127); }

GPc=newArray(256);
for (i=0; i<256; i++) {
	GPc[i]=-(1+GP[i]+Gf1*GP[i]-Gf1)/(-1-GP[i]+Gf1*GP[i]-Gf1); }
*/
half_nBins = nBins/2 -1;
GP=newArray(nBins);
for (i=0; i<nBins; i++) {
	GP[i]=((i-half_nBins)/half_nBins); }

GPc=newArray(nBins);
for (i=0; i<nBins; i++) {
	GPc[i]=-(1+GP[i]+Gf1*GP[i]-Gf1)/(-1-GP[i]+Gf1*GP[i]-Gf1); }

//Array.show("GP VAL", GP, GPc);
//waitForUser;

// GP analysis
listOrd = getFileList(ordered_images_Dir);
listDisord = getFileList(disordered_images_Dir);

var histoDir=0;
for (h = 0, j = h; h < listDisord.length; h++, j++) 
{
	Name=newArray(listOrd.length);
	open(ordered_images_Dir+listOrd[h]);
	name = getTitle;
	Name[j] = substring(name,0,lengthOf(name)-15);
	rename("Image_1a.tif");
	run("Duplicate...","title=Image_1b.tif");

	// Open also original disorder file for masking 
	name_disordered = Name[j] + "_"+chB+fileExtention; 
	//print (ordered_images_Dir+listOrd[j], dir+name_disordered);
	open(dir+name_disordered);
	rename("Image_2a_for_masking.tif");
	open(disordered_images_Dir+listDisord[j]);
	rename("Image_2a.tif");
	run("Duplicate...","title=Image_2b.tif");

	imageCalculator("Substract create 32-bit", "Image_1a.tif", "Image_2a.tif");
	rename("Image_Subs.tif");
	imageCalculator("Add create 32-bit", "Image_1b.tif", "Image_2b.tif");
	rename("Image_Add.tif");
	imageCalculator("Divide create 32-bit", "Image_Subs.tif", "Image_Add.tif");
	saveAs("tiff", rawGP_images_Dir + Name[j] + "_preGP");
	rename("Image_preGP.tif");
	setMinAndMax(-1.0000, 1.0000);
	call("ij.ImagePlus.setDefault16bitRange", 0);

	imageCalculator("Add create 32-bit", "Image_1b.tif", "Image_2a_for_masking.tif");
	rename("Image_Add_forMasking.tif");
	//selectImage("Image_Add.tif");
	selectImage("Image_Add_forMasking.tif");
	setThreshold(t*2,510);
	run("Convert to Mask");
	run("Duplicate...","title=Image_Mask.tif"); // OG
	//run("Invert");
	//selectImage("Image_Add.tif"); // OG
	selectImage("Image_Add_forMasking.tif");
	run("Subtract...", "value=254");
	rename("Image_1bit.tif");
	imageCalculator("Multiply create", "Image_1bit.tif", "Image_preGP.tif"); // OG: original GP in the macro  - but does not really work
	run(lut1);
	saveAs("tiff", GP_images_Dir + Name[j] + "_GP");

// Histogram generation
	//histoDir=histogramGP_Dir;
	//HistogramGeneration();

	// OG: fixed to exclude the low intensity values
	//imageCalculator("Multiply create 32-bit", "Image_1bit.tif", "Image_preGP.tif"); 
	selectImage("Image_Mask.tif"); // OG
	run("Create Selection");
	selectImage("Image_preGP.tif"); // OG
	run("Duplicate...","title=Image_GP1.tif");
	run("Restore Selection");
	run("Set...", "value=NaN");
	run("Select None");	
	run(lut1);
	saveAs("tiff", GP1_images_Dir + Name[j] + "_GP");
	rename("Image_GP1.tif");

// Histogram generation
	histoDir=histogramGP1_Dir;
	HistogramGeneration("All");
	run("Close All");
}

// If there is not fluorescent immunostaining
if (chC == "none") {
// HSB image generation
	if (HSB=="Yes") {
		HSBgeneration();
	}
// Print information
	listGP = getFileList(GP_images_Dir);
	time1 = getTime();
	TOTALtime = (time1-time0)/1000;
	printInfo();
}

// If there is fluorescent immunostaining
else {
	print("OG: before fluorescent immunostaining part ");
// Folder management
	temp_Dir = getDirectory("temp");

	temp_IFmask_Dir = temp_Dir + "IF mask" + File.separator;
	Folder = temp_IFmask_Dir;
	newFolder();

	File.makeDirectory(IF_images_Dir);
	File.makeDirectory(GP_IF_images_Dir);
	File.makeDirectory(histogramIF_Dir);

// Create a 1-bit mask from Immunofluorescence image and save for each image
	listDir = getFileList(dir);
	Name=newArray(listDir.length);
	for (i = 0; i < listDir.length; i++) {
		if (endsWith(listDir[i], chC + fileExtention)) {
			open(dir + listDir[i]);
			name = getTitle;
			Name[i]=substring(name,0,lengthOf(name)-suffixLen);
			run("32-bit");
			saveAs("tiff",IF_images_Dir + Name[i] + "_IF_32bits");
			run("8-bit");
			setThreshold(tc,255);
			run("Convert to Mask");
			run("Divide...","value=255");
			saveAs("tiff", temp_IFmask_Dir + Name[i] + "_GP-IF_1bit");
			close();
		}
	}

// GP-IF image
	listTempIF = getFileList(temp_IFmask_Dir);
//	listL = getFileList(GP_images_Dir);
	Name=newArray(listTempIF.length);
	for (h = 0, j = h; h < listTempIF.length; h++, j++) {
		open(temp_IFmask_Dir+listTempIF[h]);
		k = getTitle;
		Name[j]=substring(k,0,lengthOf(k)-15);
		open(GP_images_Dir + Name[j] + "_GP.tif");
		l = getTitle;
		run("8-bit");
		imageCalculator("Multiply create", k, l);
		run(lut1);
		saveAs("tiff", GP_IF_images_Dir + Name[j] + "_GP-IF");

// Histogram and Normal Distribution
		histoDir=histogramIF_Dir;
		HistogramGeneration("All");
	}

// Folder management
	listDel = getFileList(temp_IFmask_Dir);
	for (i = 0; i < listDel.length; i++) {
		File.delete(temp_IFmask_Dir+listDel[i]); }
	File.delete(temp_IFmask_Dir);

// HSB image generation
	if (HSB=="Yes") {
		HSBgeneration();
	}

// Print information
	listGP = getFileList(GP_images_Dir);
	time2 = getTime();
	TOTALtime = (time2-time0)/1000;
	printInfo();
} // OG: end of "If there is fluorescent immunostaining"

ROIanalysis(dir, listDir);

// Cleanup
run("Close All");
roiManager("Reset");
time3 = getTime();
TOTALtime = (time3-time0)/1000;
print("Execution time: " + d2s(TOTALtime,2) + " seg.");
print("Analysis Done !");


///////////////// FUNCTIONS ////////////////////

function closeAllImages() {				// This function closes all images
	while (nImages>0) {
		selectImage(nImages);
		close(); }
}

function newFolder() {					// This function creates a folder, removing any existing file in a folder with the same name
	File.makeDirectory(Folder);
	if (DeleteOldFilesFromResultsFolder)
	{
		listFolder = getFileList(Folder);
		for (i = 0; i < listFolder.length; i++) {
			File.delete(Folder+listFolder[i]); }
	}
}

function prepareImage () {				// This funcion prepares each image for the analysis
	s=getTitle;
	run("8-bit");
	run("Grays");
	run("32-bit");
	return s;
}

function HistogramGeneration (type) {		// This funcion obtains the intensity frequency histogram of the images, makes it smoother,
										// calculates the normal average distribution and also includes the GP value (and GP value
										// corrected by the Gfactor) corresponding to each intensity. An MS Excel file is generated
										// with all this data
/*	Int=newArray(256);
	Cou=newArray(256);
	Smo=newArray(256);
	NAvDist=newArray(256);
	nBins = 256;*/

	Int=newArray(nBins);
	Cou=newArray(nBins);
	Smo=newArray(nBins);
	NAvDist=newArray(nBins);
	//nBins = 256;

	if ( bitDepth() == 8)
		getHistogram(values, counts, nBins);
	else
		//getHistogram(values, counts, 256, -1, 1); // OG: fix to make sure to use [-1,1] range
		getHistogram(values, counts, nBins, -1, 1); // OG: fix to make sure to use [-1,1] range
	//Array.show("Orig Hist", values, counts); // OG

	if (0) {
		close();
		while (nImages>0) {
			selectImage(nImages);
			close(); }
	}
	
	for (u=0; u<nBins; u++) {
		Int[u]=u;
		Cou[u]=counts[u];
		if (u<=1) {
			Smo[u]=0; }
		//else if (u==255) {
		else if (u==nBins-1) {
			Smo[u]=0; }
		else {
			Smo[u]=(counts[u-1]+counts[u]+counts[u+1])/3;}
	}
	Array.getStatistics(Cou,min,max,mean,stdDev);
	//Sa=(mean*256)-counts[0]-counts[255];
	Sa=(mean*nBins)-counts[0]-counts[nBins-1];
	d=File.open(histoDir + Name[j]+"_"+type+"_Histogram.xls");
	print(d, "Intensity	Counts	Smooth	Norm Av Dist	GP	GP corrected");
	//for (k=0; k<256; k++) {
	for (k=0; k<nBins; k++) {
		NAvDist[k]=100*Smo[k]/Sa;
		print(d, Int[k]+"	"+Cou[k]+"	"+Smo[k]+"	"+NAvDist[k]+"	"+GP[k]+"	"+GPc[k]);
	}
	File.close(d);

	// Plot Histogram
	Plot.create("Normalized smoothed Histogram "+type, "GP", "NAvDist", GPc, NAvDist);
	//Plot.setLimits(xMin, xMax, yMin, yMax);
	Plot.setLimits(-1, 1, 0, NormHistMax);
	Plot.setFontSize(0.0);
	Plot.setAxisLabelSize(12.0, "plain");
	Plot.setXYLabels("GP", "NAvDist");
	Plot.setFormatFlags("1100111111");	
	Plot.show();
	selectWindow("Normalized smoothed Histogram "+type);
	saveAs("Jpeg", histoDir + Name[j]+"_"+type+"_Histogram.jpg");
	//Array.show("Full Table", Int, Cou, Smo, NAvDist, GP, GPc); // OG	
}



function HSBgeneration() {				// This function generates the HSB image of the GP images as explained in the paper
	closeAllImages();
	setBatchMode(false);

// Select images folder
	listRAW = getFileList(rawGP_images_Dir);
	
	BRIGHTNESS = newArray("Order channel","Disorder channel","IF channel");
	BRIGHTNESSnoIF = newArray("Order channel","Disorder channel");

	HSB_Dir = results_Dir + "HSB images" + File.separator;
	Folder = HSB_Dir;
	newFolder();

	Lut_Dir = getDirectory("luts");
	lut = getFileList(Lut_Dir);

	Dialog.create("Choose images and LUT");
	Dialog.addChoice("Select Hue (GP) image: ", listRAW, "none");
	if (chC == "none") {
		Dialog.addChoice("Brightness folder: ", BRIGHTNESSnoIF, "Order channel");}
	else {
		Dialog.addChoice("Brightness folder: ", BRIGHTNESS, "Order channel");}
	Dialog.addMessage("\n");
	Dialog.addChoice("Select color LUT: ",lut, "Rainbow RGB.lut");
	Dialog.addMessage("\n");
	Dialog.addChoice("Do you want to convert the whole folder into HSB images?", QUESTION, "Yes");
	Dialog.show();

// Feeding variables from dialog choices
	H = Dialog.getChoice();
	BRIGHT = Dialog.getChoice();
	Lut = Dialog.getChoice();
	WholeDir=Dialog.getChoice();

	if (BRIGHT=="Order channel") {
		mark = "_chA_32bits.tif";
		brightness_Dir = results_Dir + "Ordered Images" + File.separator;
	}
	else if (BRIGHT=="Disorder channel") {
		mark = "_chB_32bits.tif";
		brightness_Dir = results_Dir + "Disordered Images" + File.separator;
	}
	else {
		mark = "_IF_32bits.tif";
		brightness_Dir = results_Dir + "Immunofluorescence Images" + File.separator;
	}

	run("Set Measurements...", "  min limit display redirect=None decimal=5");
	index2=indexOf(Lut,".lut");
	L=substring(Lut,0,index2);

	open(rawGP_images_Dir+H);
	name=getTitle();
	Name=substring(name,0,lengthOf(name)-10);
	run(L);
	rename("Hue");

	run("Brightness/Contrast...");
	waitForUser("set min & max","set min & max");
	getMinAndMax(min,max);
	time0 = getTime();

	open(brightness_Dir+Name+mark);
	Bri=getTitle();
	run("Enhance Contrast", "saturated=0.5 normalize");

	selectWindow("Hue");
	run("RGB Color");
	run("Split Channels");

	imageCalculator("Multiply create 32-bit", Bri,"Hue (red)");
	rename("bR");
	run("8-bit");

	imageCalculator("Multiply create 32-bit", Bri,"Hue (green)");
	rename("bG");
	run("8-bit");

	imageCalculator("Multiply create 32-bit", Bri,"Hue (blue)");
	rename("bB");
	run("8-bit");

	run("Merge Channels...", "red=bR green=bG blue=bB gray=*None*");
	saveAs("tiff", HSB_Dir + Name + "_HSB");

	closeAllImages();

// HSB whole folder processing
	if (WholeDir == "Yes") {
		for (j = 0; j < listRAW.length; j++) {
			open(rawGP_images_Dir+listRAW[j]);
			name1=getTitle();
			Name1=substring(name1,0,lengthOf(name1)-10);
			rename("Hue");
			run(L);
			setMinAndMax(min,max);

			open(brightness_Dir+Name1+mark);
			run("Enhance Contrast", "saturated=0.5 normalize");
			Bri=getTitle;

			selectWindow("Hue");
			run("RGB Color");
			run("Split Channels");

			imageCalculator("Multiply create 32-bit", Bri,"Hue (red)");
			rename("bR");
			run("8-bit");

			imageCalculator("Multiply create 32-bit", Bri,"Hue (green)");
			rename("bG");
			run("8-bit");

			imageCalculator("Multiply create 32-bit", Bri,"Hue (blue)");
			rename("bB");
			run("8-bit");

			run("Merge Channels...", "red=bR green=bG blue=bB gray=*None*");
			saveAs("tiff", HSB_Dir + Name1 + "_HSB");

			closeAllImages();
		}
	}

// Make HSB LUT bar
	newImage("Hue", "8-bit Ramp", 256, 256, 1);
	run(L);
	run("Duplicate...", "title=Brightness");
	run("Rotate 90 Degrees Left");
	run("32-bit");
	selectWindow("Hue");
	run("RGB Color");
	run("Split Channels");
	imageCalculator("Multiply create 32-bit", "Brightness","Hue (red)");
	rename("bR");
	run("8-bit");
	imageCalculator("Multiply create 32-bit", "Brightness","Hue (green)");
	rename("bG");
	run("8-bit");
	imageCalculator("Multiply create 32-bit", "Brightness","Hue (blue)");
	rename("bB");
	run("8-bit");
	run("Merge Channels...", "red=bR green=bG blue=bB gray=*None*");
	rename("HSB LUT");
	selectWindow("Hue (red)");
	close();
	selectWindow("Hue (green)");
	close();
	selectWindow("Hue (blue)");
	close();
	selectWindow("Brightness");
	close();

	selectWindow("HSB LUT");
	run("Rotate 90 Degrees Left");
	run("Size...", "width=32 height=256 interpolation=None");


// Copy LUT bar to the image
	run("Copy");
	newImage("Panel LUT", "RGB White", 70, 264, 1);
	setBatchMode(false);
	makeRectangle(4,4,32,256);
	run("Paste");
	run("Colors...", "foreground=black background=black selection=yellow");
	run("Line Width...", "line="+2);
	run("Draw");
	run("Select None");
	setFont("Arial", 12);
	setColor(0, 0, 0);
	drawString(d2s(max,2),39,15);
	drawString(d2s(min,2),39,264);
	run("Select None");

	if (WholeDir == "Yes") {
		saveAs("tiff", HSB_Dir + "Lut bar");
		close();
	}
	else {
		saveAs("tiff", HSB_Dir + Name + "_Lut bar");
		open(HSB_Dir + Name + "_HSB.tif");
	}
	closeAllImages();
}

function printInfo () {					// This function prints and saves a summary of the results of the macro
	//waitForUser("Before print info, look at the Log window"); // added by OG
	print("\\Clear");
	print("GP images analysis macro");
	print("   Quantitative Imaging of Membrane Lipid Order in Cells and Organisms");
	print("   Dylan M. Owen, Carles Rentero, Astrid Magenau, Ahmed Abu-Siniyeh and Katharina Gaus");
	print("   Nature Protocols 2011");
	print("   version July 2011");
	print("\n");
	print("Ordered channel: " + chA);
	print("Disordered channel: " + chB);
	print("   Lower threshold value: " + t);
	if (chC != "none") {
		print("Channel IF: " + chC);
		print("   Lower threshold value IF mask: " + tc); }
	print("");
	print("GP images (" + (listGP.length-1) + ") saved at:");
	print("  " + GP_images_Dir);
	print("G factor: " + Gf);
	print("G factor Used: " + Ques);
	print("nBins: " + nBins);
	print("NormHistMax: " + NormHistMax);
	
	print("Scale color for GP images: " + lut1);
	print("");
	if (chC != "none") {
		listGP2 = getFileList(GP_IF_images_Dir);
		print("GP-IF images (" + (listGP2.length-1) + ") saved at:");
		print("  " + GP_IF_images_Dir);
		print(""); }
	if (HSB=="Yes") {
		print("HSB images saved at:");
		print("  " + HSB_Dir);
		print("");
	}
	print("");
	print("Execution time: " + d2s(TOTALtime,2) + " seg.");
	print("");
	print("Analysis date: "+DayNames[dayOfWeek]+", "+dayOfMonth+" "+MonthNames[month]+" "+year+" - "+hours+":"+minutes);
	print("ImageJ version " + getVersion());

	selectWindow("Log");
	saveAs("Text", results_Dir + "Log GP images generation "+year+"-"+months+"-"+dayOfMonths+" "+hours+"h"+minutes);

	setBatchMode("exit and display");	
}

// OG: ROI Analysis 
// for each image - prompt the user for selecting the diferent types of ROIs - number of types is set in the first dialog 
// use the selected type of image for the ROI selection - HSB / ordered / disordered  
// save the ROIs in the Roi Manager 
// calc hist per ROI / ROI Type
// for each ROI - Save the info: pixels (x,y), Ch1 Intensity, Ch2 Intensity, GP Intensity
// Curve fitting ? 
// It is based on Histogramgeneration function
function ROIanalysis(dir, listDir)
{
	closeAllImages();
	setBatchMode(false);

	ROI_Dir = results_Dir + "ROI images" + File.separator;
	Folder = ROI_Dir;
	oldFlagVal = DeleteOldFilesFromResultsFolder;
	DeleteOldFilesFromResultsFolder = 0;
	newFolder();
	DeleteOldFilesFromResultsFolder = oldFlagVal;

	Dialog.create("Select ROIs");
	Dialog.addChoice("Select Display Image: ", ROI_SELECTION_IMAGE, DefaultDisplayImage);
	Dialog.addChoice("Look for existing ROIs? ", QUESTION, "No");
	Dialog.addNumber("Number of ROI Types:  ", DefaultNRoiTypes);
	Dialog.addChoice("Select 1st ROI Type: ",ROI_TYPES, "Membrane");
	Dialog.addChoice("Select 2nd ROI Type: ",ROI_TYPES, "Cell");
	Dialog.addChoice("Select 3rd ROI Type: ",ROI_TYPES, "none");
	Dialog.addCheckbox("Save Per ROI Pixel Values In Table", DeafultSavePerROIPixelValuesInTable);
	Dialog.addCheckbox("Values Of Rois In Separate Files", DeafultValuesOfRoisInSeparateFiles);
	Dialog.addMessage("\n");
	Dialog.show();

	// Feeding variables from dialog choices
	D = Dialog.getChoice();
	Ques1 = Dialog.getChoice();	
	nT = Dialog.getNumber();
	type1 = Dialog.getChoice();
	type2 = Dialog.getChoice();
	type3 = Dialog.getChoice();
	SavePerROIPixelValuesInTable = Dialog.getCheckbox();
	ValuesOfRoisInSeparateFiles = Dialog.getCheckbox();

	run("Set Measurements...", "  min limit display redirect=None decimal=5");
	
	listOrd = getFileList(ordered_images_Dir);
	//listDisord = getFileList(disordered_images_Dir);
	/*listGP = getFileList(GP_images_Dir);
	listGP1 = getFileList(GP1_images_Dir);
	listPreGP = getFileList(rawGP_images_Dir);
	if (File.exists(IF_images_Dir))
		listIF = getFileList(IF_images_Dir); */

	//Name=newArray(listOrd.length);
	Name=newArray(listOrd.length); 

	//listDir = getFileList(Disp_Dir);
	//listDir = listOrd;
	s = 0;
	for (j = 0; j < listOrd.length; j++) 
	{
		//ord_filename = dir + listDir[i], chA+fileExtention; // same in original folder and in the ordered file folder
		open(ordered_images_Dir+listOrd[j]);
		name = getTitle;
		name = replace(name,"_chA_32bits.tif",""); 
		name_disordered = name + "_"+chB+fileExtention; 
		name_HSB = name + "_HSB.tif"; 
		// just to save the file name for HistogramGeneration function
		Name[j] = name;
		rename("OrdIm");

		print (ordered_images_Dir+listOrd[j], dir+name_disordered);
		open(dir+name_disordered);
		//open(disordered_images_Dir+listDisord[j]);
		rename("DisOrdIm");
		open(GP_images_Dir+name+"_GP.tif");
		rename("GPIm");
		open(GP1_images_Dir+name+"_GP.tif");
		rename("GP1Im");
		//open(rawGP_images_Dir+listPreGP[j]);
		open(rawGP_images_Dir+name+"_preGP.tif");
		rename("PreGPIm");
		open(HSB_Dir+name_HSB);
		rename("HSBIm");
		if (File.exists(IF_images_Dir)) {
			open(IF_images_Dir+listIF[j]);
			rename("IFIm");
		}
		if (D=="Order channel") 		{ dispIm = "OrdIm"; }
		else if (D=="Disorder channel") { dispIm = "DisOrdIm"; }
		else if (D=="Order+Disorder")	{ 
			imageCalculator("Add create 32-bit", "OrdIm", "DisOrdIm");
			rename("AddIm");
			dispIm = "AddIm"; }
		else if (D=="IF channel") 		{ dispIm = "IFIm"; }
		else 							{ dispIm = "HSBIm"; }

		// select ROIs
		roiManager("Reset");
		RoiFound = 0;
		if (Ques1=="Yes") { // look for exsiting ROIs
			// Open contact ROIs 
			if (File.exists(ROI_Dir+name+"_ROIs.zip")) 
			{
				//print("Opening ROI file");
				roiManager("Reset");
				roiManager("Open", ROI_Dir+name+"_ROIs.zip");		
			
				// loop over ROIS and count different types 
				nRois = roiManager("count");
				lastType1 = 0;
				nType2 = 0;
				nType3 = 0;
				for (n = 0; n < nRois; n++) 
				{
					// it is assumed - but not checked that all type1 ROIs come first follows by all type2 ROIs followed by type3
					roiManager("Select", n);
					roiName=call("ij.plugin.frame.RoiManager.getName", n);
					if (startsWith(roiName, type1)) lastType1 = lastType1 + 1;
					else if (startsWith(roiName, type2)) nType2 = nType2 + 1;
					else if (startsWith(roiName, type3)) nType3 = nType3 + 1;
				}	
				lastType2 = nType2 + lastType1;			
				lastType3 = nType3 + lastType2;			
				RoiFound = 1;
			}
			else {
				print("ROI File does not exist: ",ROI_Dir+name+"_ROIs.zip");
			}
		} 
		if (RoiFound == 0) {
			// Display intro message
			showMessage("Draw ROIs",DrawRoisStartMessage);
			getROIs(dispIm, type1, ROI_COLORS[0]);
			lastType1 = roiManager("count");
			if (nT > 1) {
				getROIs(dispIm, type2, ROI_COLORS[1]);
				lastType2 = roiManager("count");
			}
			if (nT > 2) {
				getROIs(dispIm, type3, ROI_COLORS[2]);
				lastType3 = roiManager("count");
			}
			// Save the ROIs
			roiManager("Save", ROI_Dir+name+"_ROIs.zip");
		}			
		// Save flattened image with overlay of selected contacts
		//selectImage(name1);
		selectImage(dispIm);
		roiManager("show all with labels");
		//roiManager("show all without labels");
		run("Flatten");
		saveAs("Tiff", ROI_Dir+name+"_ROIsOverlay.tif");
		close();
		
		setTool("hand");		

		// Measure ROIs
		// calculate histograms for each category ?
		type1_indexes = newArray(lastType1);
		for (n=0; n < lastType1; n++) type1_indexes[n] = n;
		HistogramGenerationForROIs(type1, type1_indexes);

		if (nT > 1) {
			nType2 = lastType2-lastType1;
			type2_indexes = newArray(nType2);
			for (n=0; n < nType2; n++) type2_indexes[n] = n + lastType1;
			HistogramGenerationForROIs(type2, type2_indexes);
		}
					
		if (nT > 2) {
			nType3 = lastType3-lastType2;
			type3_indexes = newArray(nType3);
			for (n=0; n < nType3; n++) type3_indexes[n] = n + lastType2;
			HistogramGenerationForROIs(type3, type3_indexes);
		}

		if (SavePerROIPixelValuesInTable)
		{
			// get values of ord/disord preGP ([-1,1]) values for each Roi - print them into one text file for each iamge
			run("Clear Results");
			setBatchMode(true);
			roiValuesBaseName = ROI_Dir+name+"_RoiVal";
			for (n = 0; n < roiManager("count"); n++)
			{
				roiName=call("ij.plugin.frame.RoiManager.getName", n);
				//print("Save Per Pixel Values:", name," ROI #",n, roiName, "Getting values ...");
				ord = Get_xyVal_Column(n, "OrdIm");			
				disord = Get_xyVal_Column(n, "DisOrdIm");
				pregp = Get_xyVal_Column(n, "PreGPIm");
				gp = Get_xyVal_Column(n, "GPIm");
				gp1 = Get_xyVal_Column(n, "GP1Im");
				//run("Clear Results");
				np = Array.getSequence(x.length); // x&y array are filled by Get_xyVal_Column function
				/*namearr = newArray(roiName);	
				for (m = 1; m < x.length; m++) namearr = Array.concat(namearr,roiName);*/
				namearr = newArray(x.length);	
				for (m = 0; m < x.length; m++) namearr[m] = roiName;

				//Array.show(roiName+"_Table",np, namearr, x, y, ord, disord, pregp, gp, gp1);
				//waitForUser("Check table for Roi #",n);
				if (ValuesOfRoisInSeparateFiles == 1)
				{			
					//print("Save Per Pixel Values:", name," ROI #",n, roiName, "Save ROI info to file ...");
					Array.show(roiName+"_Table",np, namearr, x, y, ord, disord, pregp, gp, gp1);
					selectWindow(roiName+"_Table");	
					// save the table to file , reset table
					saveAs("Results", roiValuesBaseName+"_"+roiName+".xls");
					run("Clear Results");
				} else 
				{
					if (n == 0) 
					{
						xA = x;
						yA = y;
						npA = np;
						namearrA = namearr;
						ordA = ord;
						disordA = disord;
						pregpA = pregp;
						gpA = gp;
						gp1A = gp1;
					} else 
					{
						print("Save Per Pixel Values:", name," ROI #",n, roiName, "Concat arrays ...");
						xA = Array.concat(xA,x);
						yA = Array.concat(yA,y);
						npA = Array.concat(npA,np);
						namearrA = Array.concat(namearrA,namearr);
						ordA = Array.concat(ordA,ord);
						disordA = Array.concat(disordA,disord);
						pregpA = Array.concat(pregpA,pregp);
						gpA = Array.concat(gpA,gp);
						gp1A = Array.concat(gp1A,gp1);
					}
				}
			}
			if (ValuesOfRoisInSeparateFiles == 0)
			{
				Array.show("Results",npA, namearrA, xA, yA, ordA, disordA, pregpA, gpA, gp1A);
				//waitForUser("Check table for All Rois");
				// save the table to file, reset table
				selectWindow("Results");
				saveAs("Results", roiValuesBaseName+".xls");
				run("Clear Results");
			}
			setBatchMode(false);
		}	
		closeAllImages();
	} // for
	//print("End of ROIanalysis");
} // end of ROIanlysis


//function Add_xyVal_Column(x, y, n, roiName, Im, ColName, FirstRow)
function Get_xyVal_Column(n, Im)
{
	selectImage(Im);
	roiManager("Show None");
	roiManager("select", n);
	run("Save XY Coordinates...", "save=["+ROI_Dir+"temp.csv]");
	if (isOpen("temp.csv")) 
	{
		selectWindow("temp.csv");
		run("Close");
	}
	//if (isOpen("Results")) 
	//	run("Clear Results");
	open(ROI_Dir+"temp.csv");
	//selectWindow("Results");
	x = Table.getColumn("X","temp.csv");
	y = Table.getColumn("Y","temp.csv");
	val = Table.getColumn("Value","temp.csv");
	//Array.show("Tmp",x, y, val);
	/*print("In Get_xyVal, length=", x.length, y.length, val.length);
	Array.show("temp",x,y,val);
	waitForUser("Check temp table");*/
	return val;
}


function HistogramGenerationForROIs(type, indexes)
{
	if(isOpen("Mask")) 
	{
		selectWindow("Mask");
		close(); 
	}
	selectImage("GP1Im");
	imName = "Image_GP_"+type+".tif";
	run("Duplicate...","title="+imName); 
	roiManager("select", indexes);
	if (lengthOf(indexes) > 1)
		roiManager("Combine");
	run("Create Mask");
	selectWindow("Mask");
	run("Create Selection");
	selectWindow(imName);
	run("Restore Selection");
	run("Set...", "value=NaN");
	run("Select None");
	// JAKUB
	saveAs("tiff", GP_images_Dir + Name[j] + "_thresholdGP");
	//histoDir=histogramGP1_Dir;
	histoDir=ROI_Dir;
	HistogramGeneration(type);
	
}


function getROIs(Im, type, color)
{
	getAnotherROI = 1;
	RoiNum = 1;

	// Prompt the user to draw the pairs of contact lines 
	while(getAnotherROI == 1)
	{
		SelectOneRoi(Im, type, color, RoiNum);
		getAnotherROI = getBoolean("Do you want to select another ROI for  "+type+" ?");
		if (getAnotherROI == 1) 
			RoiNum = RoiNum + 1;
	}	
}

// prompt the user to draw a ROI, add it to Roi manager, rename it and change its color
function SelectOneRoi(Im, type, color, RoiNum)
{
	roiName = type + "_" + RoiNum;
	
//	mitoSideName = "M"+	MitoNum + "_C" + ContactNum + "_M";
//	erSideName   = "M"+	MitoNum + "_C" + ContactNum + "_E";

	selectImage(Im);
	run("Select None");
	roiManager("show all without labels");
	//setTool("freeline");
	setTool("polygon");
	waitForUser("Please draw roi for " + type + " "+ RoiNum +" , Click OK when done");
	roiManager("Add");
	roiManager("Select", roiManager("Count")-1);
	roiManager("Rename", roiName);	
	roiManager("Set Color", color);
	roiManager("Set Line Width", RoiLineWidth);
}

