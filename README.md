# Diabetes-LFP-Analysis

Written by Zachary Ip (ip.zachary.t@gmail.com) 3/7/23.

This repository is expected to be used in conjunction with
[spectral-analysis-tools](https://github.com/Zachary-Ip/spectral-analysis-tools),
and [utils-toolbox](https://github.com/Zachary-Ip/utils-toolbox). Code
will not run without those repositories on path. This repository was
developed in MATLAB, and requires a license with Signal processing
toolbox among others.

This code was developed to process a dataset collected by Jialing Liu
(Jialing.liu@ucsf.edu) and co: Gratianne Rabiller
(Gratianne.rabiller@gmail.com), and Shahram Zarrabian
(szarrabian@gmail.com ). This data was collected using control db/+ and
type-2 diabetic db/db mice. Data was collected using 16-site linear
electrodes. Meta-data and information about the organization of the
animals can be found in `init_Spk_info.m` and `init_Synopsis_Dbdb.xlsx`.

Raw data begins being processed in `notebook_Data_Analysis`, and
calculated measures are compared in `notebook_Statistics`. Some values
used for the manuscript were calculated by Gratianne Rabiller and stats
for those values are calculated in `notebook_Statistics_Prism`.

Auxiliary test codes `notebook_detect_SPWR`,
`notebook_CSD_individual_animals`, and `notebook_SWR_Spectrogram` are
not used in final manuscript (as of writing 3/7/23), and were only used
during preliminary data analysis.

Final figure generation for the manuscipt (as of writing 3/7/23) are
generated in prism by Shahram Zarrabian and Jialing Liu. To facilitate
this, `notebook_Data_Analysis` saves excel files of the final measures
to share, however, `notebook_Statistics` also generates figures in
MATLAB.
