# Usage example: 
# julia cnode_predict.jl --weights=PATH_TO_WEIGHTS --data=PATH_TO_DATA  --output=PATH_TO_PREDICTION

# Required imports
using Flux, DifferentialEquations, cNODE, DelimitedFiles, SharedArrays, Statistics
using Base.Iterators: product
using MLDataPattern: kfolds
using Optim, Optimisers
using Tracker
using BSON: @load
using ArgParse
using CSV

# Creating an argument parser
parser = ArgParse.ArgumentParser()
ArgParse.add_argument!(parser, "--weights", required=true, help="Path to CNODE weights file for the phenotype, .bson format")
ArgParse.add_argument!(parser, "--data", required=true, help="Path to data file to make prediction, vector in csv format with 0/1 for each species")
ArgParse.add_argument!(parser, "--output", default="prediction.csv", help="Path to output file")
args = ArgParse.parse_args(parser)

weights_path = args.weights
X = readcsv(args.data)

# Building the architecture-backbone for the predicting model
cnode = getModel(N)

# Loading weights to the empty model
@load weights_path weights
Flux.loadparams!(cnode, weights)

# Making a prediction and saving it to the csv file
pred = predict(cnode, X)
println("Check your prediction in the cnode_prediction.csv file!")
CSV.write(args.output, DataFrame(pred))

