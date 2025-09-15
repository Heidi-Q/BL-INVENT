# BL-INVENT
Molecular generative model to design bivalent MsbA inhibitors 
## Quick start

Installation
-------------

1. Install [Conda](https://conda.io/projects/conda/en/latest/index.html)
2. Clone this Git repository
3. Open a shell, and go to the repository and create the Conda environment:
   
        $ conda env create -f blinvent.yml

4. Activate the environment:

        $ conda activate blinvent_env

5. Install in-house reinvent_scoring

        $ cd reinvent_scoring

        $ pip install reinvent_scoring-0.0.73_bq-py3-none-any.whl

6. Open another shell, and clone in-house [DockStream](https://github.com/jidushanbojue/DockStream-master) repository

   This is the docking component special for Protac-invent.

        $ conda env create -f environment.yml






## Usage
1. Edit template Json file (for example in MsbA/msba_blinvent_config.json).

   Templates can be manually edited before using. The only thing that needs modification for a standard run are the file and folder paths. Most running modes produce logs that can be monitored by tensorboard
2. python input.py config.json

## Analyse the results

1. tensorboard --logdir "progress.log"

    progress.log is the "logging_path" in template.json
