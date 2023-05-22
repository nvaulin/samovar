# Run the cNODE prediction

## Requirements instalation
Please, install Julia prior than running the main script. You can find the corresponding version [here](https://julialang.org/downloads/). Then add julia to your PATH by running 
```console
$  export PATH="$PATH:<julia folder>/bin"
```
where <julia folder> is the folder in the directory where you run the instalation. 

Then you have to install dependencies, it will take some time, but will make further runs faster. In order to do this, please, run the following command from the current directory (*Manifest.toml* and *Project.toml* files should be there):
```console
$  julia -e "using Pkg; Pkg.instantiate(); Pkg.precompile()"
```

## Models training
In order to run training, you have to pass a path to your dataset to **--data** flag for *cnode_train.jl*. In order to be used your data have to be in numerical format, where rows correspond 
 to unique bacteria, columns to experiments, and values of the matrix represent the relative abundance of species in each sample. Put this table to the *data* folder and add an empty directory ./hyperparameters. <br>
  
  To run the training run:
```console
$  julia --project="path_to_dir_with_manifest_file" cnode_train.jl --data=PATH_TO_DATA
```
  This script will perform training and save a model weights to the current folder in the *bson* format. You can later use them to make predictions.
  
## Predictions
  
It is important to notice that predictions can be made only with the set of bacteria that was used for the training, the model is robust and unfortunately unable to add new species during prediction.
To make prediction you have to create *.csv* file with only one column (array n*1) representing presence of bacteria species in the sample: 1 or 0. The path to this file has to be passed as **--data** argument, whereas path to weights in bson format have to be passed to **--data** parameter. If you want to specify name of the prediction file, use flag **--output**.  
Example of the full command for the prediction is below:
  
```console
$  julia --project="path_to_dir_with_manifest_file" cnode_predict.jl --weights=PATH_TO_WEIGHTS --data=PATH_TO_DATA --output=PREDICTION_PATH
```

This will output predictions-ralative abundancies to the file *PREDICTION_PATH*, where PREDICTION_PATH ends with .csv.

For more details on cNODE package and the official documentation, please, check [project github](https://github.com/michel-mata/cNODE.jl).

