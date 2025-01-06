
path_original = "C:/Users/satheps/PycharmProjects/Results/2022/May27/Stacks_hotl/";
path_segmented = "C:/Users/satheps/PycharmProjects/Results/2022/May27/Hotl_phases/phase3_seg/";
path_combined = "C:/Users/satheps/PycharmProjects/Results/2022/May27/Hotl_phases/phase3_overlay_june3/";


subdirlist_original  = getFileList(path_original);
subdirlist_segmented = getFileList(path_segmented);
//subdirlist_combined = getFileList(path_combined);

Array.print(subdirlist_original);
Array.print(subdirlist_segmented);
if (subdirlist_original.length == subdirlist_segmented.length){
	//Do nothing ?
}
else {
	if (subdirlist_segmented.length < subdirlist_original.length){
		subdirlist_original = subdirlist_segmented;
		print("Original set to be equal to segmented");
		Array.print(subdirlist_original);
		Array.print(subdirlist_segmented);
		// Temporary. ideally check if all elements of segmented are also in original
	}
	else{
		exit("original and segmented subdirectories dont match in size");
	}
}


for (i=0; i<subdirlist_original.length; i++) {
	subdirfilelist_original = getFileList(path_original+ File.separator +subdirlist_original[i]);
	subdirfilelist_segmented = getFileList(path_segmented+ File.separator +subdirlist_segmented[i]);
	if (!File.exists(path_combined + File.separator +subdirlist_original[i])){
		File.makeDirectory(path_combined+File.separator +subdirlist_original[i]);
	}
		
	
	for (og=0;og<subdirfilelist_original.length;og++){
		
		basename = split(subdirfilelist_original[og],".");//[0];
		basename = basename[0];
		print("Starting" + basename);
		
		for (sg=0;sg<subdirfilelist_segmented.length;sg++){
			
			if (subdirfilelist_segmented[sg].contains(basename)){
				print("Working on" + subdirfilelist_segmented[sg]);
				//original
				origfilepath = path_original + File.separator +subdirlist_original[i] + File.separator + subdirfilelist_original[og];
				segfilepath = path_segmented + File.separator +subdirlist_segmented[i] + File.separator + subdirfilelist_segmented[sg];
				run("Bio-Formats Importer", "open="+origfilepath+" autoscale color_mode=Colorized rois_import=[ROI manager] split_channels view=Hyperstack stack_order=XYCZT");
				close("["+subdirfilelist_original[og]+" - C=1]"); //try simple pattern
				close("["+subdirfilelist_original[og]+" - C=2]");
				close("["+subdirfilelist_original[og]+" - C=3]");
				selectWindow(""+subdirfilelist_original[og]+" - C=0");
				run("Enhance Contrast...", "saturated=0.3 normalize process_all use");
				//open(path_original + File.separator +subdirlist_original[i] + File.separator + subdirfilelist_original[og]);
				run("8-bit");
//				run("RGB Color");
				//segmented
				run("Bio-Formats Importer", "open="+segfilepath+" autoscale color_mode=Colorized rois_import=[ROI manager] split_channels view=Hyperstack stack_order=XYCZT");
//				open(path_segmented + File.separator +subdirlist_segmented[i] + File.separator + subdirfilelist_segmented[sg]);
				run("8-bit");
//				run("RGB Color");
				//Merge
				winlist = getList("image.titles");
				Array.print(winlist);
				//c5: cyan and c1: red. cyan + red = white. Original 
				print(subdirfilelist_original[og], subdirfilelist_segmented[sg]);
				print("Merge Channels...", "c1=["+subdirfilelist_original[og]+" - C=0] c5=["+subdirfilelist_segmented[sg]+" - C=0] create ignore");
				run("Merge Channels...", "c1=["+subdirfilelist_original[og]+" - C=0] c5=["+subdirfilelist_segmented[sg]+" - C=0] create ignore");
				selectWindow("Composite");
//				combinedfilepath = path_combined + File.separator +subdirlist_original[i] + File.separator + subdirfilelist_segmented[sg];
				combinedfilepath = path_combined +subdirlist_original[i] + subdirfilelist_segmented[sg];
				print("Saving: ", combinedfilepath);
				saveAs("Tiff", combinedfilepath);
				close("*");
			}
		}
		 
	}
}
print("Completed");
//When comparing use: run("Channels Tool...");
//Stack.setDisplayMode("color");
//Stack.setDisplayMode("color");
//Stack.setChannel(1);
//Stack.setChannel(2);