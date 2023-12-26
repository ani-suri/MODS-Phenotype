#!/bin/bash
YEAR=2015
NUM_OF_PROCESSES=180
CSN_LIST=/labs/kamaleswaranlab/MODS/EliteDataHacks/pt_lists/pt_list_2015.psv
PICKLE_PATH=/opt/scratchspace/KLAB_SAIL/Emory_Data/Pickles/Yearly_Pickles
OUTPUT_PATH=/opt/scratchspace/KLAB_SAIL/Emory_Data/Pickles/Encounter_Pickles

echo This is the csn list- $CSN_LIST
echo This is the pickle path- $PICKLE_PATH
echo This is the output path- $OUTPUT_PATH
echo This is the num processes - $NUM_OF_PROCESSES
echo The year is - $YEAR

python /home/gmatlin/MODSPhenotypes/mods/sepy/em_make_dicts_parallel.py $CSN_LIST $PICKLE_PATH $OUTPUT_PATH $NUM_OF_PROCESSES $YEAR