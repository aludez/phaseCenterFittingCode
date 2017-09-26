This is the code i'm going to use for fitting phase centers

I drew HEAVILY from what my dear friend, Linda Cremonesi (l.cremonesi@ucl.ac.uk) wrote to do this in Anita 3. I also used a lot of what Stephen Hoover wrote in his thesis to inform this.

WAIS runs are about 122 - 153, with some weird parts along the way, especially in runs 129-132 where the pulser was at 45 degrees.

Procedure goes as follows:
First - make trees of CorrelationSummaryAnita4's using the macros/makePulserTreeAnalysisWf.C

Second - make phase center calibration files using macros/fitLindaSplitCost.C (V and H are made separately).  this is where you actually fit things and generate the phase center numbers.
	note - this was just one of many ways i tried fitting for the phase centers, probably the most successful way.

Third - turn the separate V and H files into one phaseCenterCalibration file with macros/makePhaseCenterCalibFile.C

Fourth - you can make summary resolutions for comparison purposes using macros/makeResolutions.C, this basically just replaces all your phase center numbers with ones read in from a calibration file and generates eventSummarys you can look at

Fifth - use macros/doResolutions.C to generate the summary graphs for you. These are the actualy outputs that are worth looking at


