# SAMPLE-ICM
A tool to perform adaptive sampling using RADICAL ENTK

## Installation
(**Note still in progress)
- Clone the repository on both the client and server in which you intend to run SAMPLE
- execute the command `pip install .` from within the repository

## Usage
### Command-line interface
`python sample.py [model] [analysis] [param_space] [init_params] [resource_reqs] [n_samples] [n_iterations] [output]`
As seen above, SAMPLE requires certain user input files in order to execute properly. These files include:
- **model**: The location of the model executable on the server
- **analysis**: The location of the selection algorithm executable on the server. Note: ATM this needs to be a command-line which accepts as inputs the parameter space csv file, the number of sample and the output file csv. The contents of the output file CSV should contain the input parameter values for the next iteration. See radical/selection_random.py and radical/selection_adaptive.py for examples.
- **param_space**: the location of the CSV file on the server containing containing all possible parameter combinations (See example: [parameter_space.csv](https://github.com/radical-collaboration/SAMPLE-ICM/blob/radical/wrapper/parameter_space.csv))
- **init_params**: the location of the CSV file on the server containing only the initial parameter combinations (See example: [start_params.csv](https://github.com/radical-collaboration/SAMPLE-ICM/blob/radical/wrapper/start_params.csv)
- **resource_reqs**: the location of the json file containing resource requirements of the server (see: radical.pilot docs for instructions). Example [bridges](https://github.com/radical-collaboration/SAMPLE-ICM/blob/radical/radical/bridges_res.json)
- **n_samples**: The number of samples to run per iteration (note: will maybe remove/ infer from init_params length.
- **n_iterations**: Number of iterations. Currently the stop condition is a maximum number of iterations.
- **output**: Output data location of the server

## Sample command
(on Bridges)
`python sample.py /pylon5/mc3bggp/vhayot/SAMPLE-ICM/wrapper/wrapper.py /pylon5/mc3bggp/vhayot/SAMPLE-ICM/radical/selection_adaptive.py /pylon5/mc3bggp/vhayot/SAMPLE-ICM/wrapper/parameter_space.csv /pylon5/mc3bggp/vhayot/SAMPLE-ICM/wrapper/start_params.csv bridges_res.json 14 3 /pylon5/mc3bggp/vhayot/SAMPLE-ICM/wrapper`
