# How many neurons are sufficient for perception of cortical activity?
## Dalgleish-eLife-2020
Figure data/code and raw data import/processing code used in:

Dalgleish HWP*, Russell LE*, Packer AM*, Roth A, Gauld OM, Greenstreet F, Thompson EJ, Hausser M (2020). __How many neurons are sufficient for perception of cortical activity?__ _eLife_

Data for use with this repository for should be downloaded from:
https://doi.org/10.6084/m9.figshare.13135505
https://figshare.com/projects/How_many_neurons_are_sufficient_for_perception_of_cortical_activity_/91031

These data were collected from mice at various stages through being trained to report optogenetic activation of barrel cortex (both 1P and 2P). Final experiments tested animal behaviour and neuronal responses in mice tasked with reporting perception of two-photon optogenetic activation of varying numbers of cortical excitatory neurons by licking for sucrose rewards. See paper for more details.

Raw extracted calcium imaging data (Suite2p ROIs and traces) and photostimulation data can also be found on [Figshare](https://doi.org/10.6084/m9.figshare.13128950) along with some description of data format. Raw calcium imaging movies (~1TB) are available upon reasonable request.

All code is written by [Henry Dalgleish](https://github.com/hwpdalgleish) and [Lloyd Russell](https://github.com/llerussell), except for code in `~/utils/external`.

If anything is unclear please contact hwpdalgleish@gmail.com.

If you use these data in a paper, please cite the above Dalgleish et al. 2020 (eLife) paper and reference this Github repo:<br/>
https://github.com/alloptical/Dalgleish-eLife-2020<br/>
and this Figshare databases:<br/>
https://doi.org/10.6084/m9.figshare.13135505
https://doi.org/10.6084/m9.figshare.13128950

## Setup
* (required) Ensure that `DalgleishHausser2020_setup.m`, `DalgleishHausser2020_plots_main.m` and `DalgleishHausser2020_analysis_main.m` are in the repository's base directory (`~/DalgleishHausser2020Repo/`) as they use this to add all necessary paths. Unless you want to add paths manually then simply make sure that all repository contents (and any raw data downloaded from Figshare) remain somewhere in `~/DalgleishHausser2020Repo/`.
* (required) [Download processed data from Figshare](https://doi.org/10.6084/m9.figshare.13135505). Just ensure that all downloaded files are put in a directory in the repository's base directory.
* (optional) [Download raw data from Figshare](https://doi.org/10.6084/m9.figshare.13128950) (see _Raw data download and import_ section below for instructions).

## Plotting
`DalgleishHausser2020_plots_main.m` quick access to figure/relevant figure supplement plots. This calls plotting scripts in `~/figureScripts`, the names of which should be self-explanatory.

## Analysis
* `DalgleishHausser2020_analysis_main.m`: quick access to run analysis routines.
* `~/analysisScripts` contains:
  * `DalgleishHausser2020_importProcessing.m`: import/processing pipeline used for the majority of data in the paper (Figures 2 - 4). Note this needs the Figshare data.
  * `DalgleishHausser2020_threshold_analysis.m`: the threshold estimation procedure.
  * `DalgleishHausser2020_Fig1_analysis.m` and `DalgleishHausser2020_Fig2_3_4_analysis.m`: intermediate analysis routines used for Figure 1 and Figures 2 - 4 respectively.

## Exploring
* Run `DalgleishHausser2020_setup.m`.
* Import desired data from `~/data`:
  * `DalgleishHausser2020_All1PLearning_raw.mat` raw training data during initial 1P priming phase (primarily Figure 1 - figure supplement 3).
  * `DalgleishHausser2020_All1PLearning_analysed.mat` same data as above but processed to return variables plotted in manuscript figures.
  * `DalgleishHausser2020_All1PLearning_raw.mat` raw training data during 1P --> 2P transition sessions when 1P and 2P photostimuli are delivered together (primarily Figure 1, Figure 1 - figure supplement 5).
  * `DalgleishHausser2020_All2PLearning_analysed.mat` same data as above but processed to return variables plotted in manuscript figures.
  * `DalgleishHausser2020_imaging_raw.mat` calcium imaging/photostimulation data from number of cells psychometric curve experiments (Figures 2 - 4). Note this doesn't include traces as they are too large, but these can be downloaded from Figshare (see above).
  * `H527_20190206_165853.mat` individual Phase 2 behavioural session (1P and 2P 200 neuron photostimulation + catch trials) plotted in Figure 1d.
* Use scripts in `analysisScripts` to explore data.
* `~/utils` contains lower level analysis routines that might be useful.

## Data format
* 1P/2P training data (Figure 1 and supplements) in `DalgleishHausser2020_All1PLearning_raw.mat` and `DalgleishHausser2020_All2PLearning_raw.mat` with relevant variable:
  * `all_bhv`: a struct array where each element is the training data for a given animal. The most relevant fields are `bhv_daily` which contains the training data for each training day, `params` which contains import parameters and `powers` which contains the LED powers in mW (NB these are only present for 1P only training data where LED power is varied. For 1P/2P training data all LED powers are 0.05 mW). 
* 2P imaging, photostimulation and training data during number of cells psychometric curve sessions (Figures 2 - 4 and supplements) in `DalgleishHausser2020_imaging_raw.mat`:
  * `session`: a struct array where each element is a single animal, single behavioural session. This has relevant fields:
    * `s2p` which contains Suite2P ROIs and metadata.
    * `targets` which has all target locations for different trial-types (note there are 8: 1 catch trial and 7 stimulus trials; stimulus trials are _stim_type_ 7 and catch trials are _stim_type_ 6).
    * `cat_paq` which contains the timing information for stimuli and behavioural responses. Note that this is a struct of dimensions nStimulusChannels (these are the possible stimulus channels our behavioural hardware can use to trigger things). Stimulus channel 7 (`cat_paq(7)`) contains photostimulus trigger times (organised into 7 trial-types = [200 100 75 50 25 10 5] neurons stimulated) and stimulus channel 6 (`cat_paq(6)`) contains catch trial trigger times (1 trial-type = 0 neurons stimulated). Note these are then "flattened" into the `num_cells` variable (and all subsequent constructed variables) with the order [0 200 100 75 50 25 10 5] neurons (see below).
    * Note that if you download and import raw data from Figshare this will also include full neuropil subtracted fluorescence traces as raw F in the `f_sub` field and as dF/F (calculated with rolling baseline) in the `dff` field. Raw F was used for all STA analyses as we convert to z-scores on each trial epoch separately on the basis of the pre-trial mean and s.d.
  * `num_cells`: a vector containing the _trial-type_ order of number of cells (i.e. trial-type 2 is 200 neurons since num_cells(2) = 200 etc.). This corresponds to the columns of all of the analysis variables produced (i.e. num_cells(2) = 200 meaning that variable(:,2) corresponds to values for that variable on 200 cell stimulation trials).
  * Most other variables adhere to the following convention: dimensions = numSessions * numTrialTypes where rows = sessions (corresponding to elements of `session`) and columns = trial-types (dictated by values in `num_cells`). Variable names are somewhat logical, and can hopefully be worked out from context. Some relevant variables of this shape are:
    * `targeted_flag`: cell array where each cell contains a logical vector dictating whether corresponding ROIs were targeted in that session and on that trial type. So targeted_flag{1,2}(20) would tell you whether ROI 20 was targeted on trial type 2 (which happens to be 200 cell trials as num_cells(2) = 200) in session 1.
    * `tw_traces_zs`: cell array of z-scored fluorescence traces aligned to stimulus onset (stimulus triggered averages; STAs). Each cell contains STAs of dimensions numROIs * numTrials * numTimepoints.
    * `tw_response_zs`: cell array of post-stimulus response (z-scored fluorescence) calculated from post-stimulus epoch in above STA traces. Used for most analyses. Each cell is numROIs * numTrials.
    * `tw_response`: cell array of animals' behavioural responses (numTrials * 1).
    * `tw_rxn`: cell array of animals' reaction times (numTrials * 1). NB trials with no lick are filled with NaN.
  
## Raw data download and import
As mentioned, this is available on [Figshare](https://doi.org/10.6084/m9.figshare.13128950).
  
Download all .zip files. These contain `sessionDate_animalID` directories containing data from individual psychometric curve sessions (see Figshare description for more details). For current import/processing pipelines to run, these need to be placed in a directory structure somewhere in the DalgleishHausser2020Repo base directory:
<pre>
├── ~/DalgleishHausser2020Repo/rawDataDirectory/
    ├── 20191203_L541/
    |    └── experimentDataFile1
    |    └── experimentDataFile2
    ├── 20191009_L544/
    |    └── experimentDataFile1
    |    └── experimentDataFile2
    . etc.
</pre>

Files were zipped on a mac and should hopefully unzip fine onto a mac (i.e. myDirectory/ --> myDirectory.zip --> myDirectory/), however sometimes on Windows this can get mangled (i.e. myDirectory/ --> myDirectory.zip --> __myDirectory/__ myDirectory/). To avoid this either:
* Put all .zip files into a directory and pass to `unzipDirectories.m` (`~/DalgleishHausser2020/utils/import/`) to unzip into correct structure.
* Manually unzip and remove any added upstream directories.
