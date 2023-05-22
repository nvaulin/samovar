################################################################################
#       0. CONFIGURATION
 ################################################################################

# 1. Make sure that:
# 	1) setup.jl is in the same directory you are running this script from!
# 	2) you have folder ./data/ with your data files
# 	3) you have folder ./hyperparameters/ for intermediate results

# 2. Setting up the initial configuration.Recommended number of workerd - 4, but you can change to more (will cause more overheads). 
begin
    num_workers = 4
    include("setup.jl")
end
# NOTE: step above is bottleneck and might take >10 minutes

parser = ArgParse.ArgumentParser()
ArgParse.add_argument!(parser, "--data", required=true, help="Path to data file")
args = ArgParse.parse_args(parser)

################################################################################
#       1. HYPERPARAMETER SEARCH and loading the data!
################################################################################

# Here we perform hyperparameters search. @everywhere ensures the implementation of this code on all workers while parallelization.

@everywhere begin
    # Define potential hyperparameter values
    LearningRates = [[0.001,0.0025],[0.001,0.005],[0.001,0.01],[0.01,0.025],[0.01,0.05],[0.01,0.1]]
    Minibatches = [1,5,10]
  
    # Parameters for run
    max_epochs = 1000
    early_stoping = 100
    report = 50
  
    # Iterator over hyperparams,params and repetitions
    inx = collect(product(enumerate(LearningRates),enumerate(Minibatches)))[:]
    "search hyperparameters..." |> println
end


i = 1
DATA = "your data" 

# Loading the data
path = args.data
Z,P = import_data(path)
N,M = size(Z)

# Exploring hyperparameters in small datasets
"training $DATA..." |> println

for it in inx
    ((j,lr),(k,mb)) = it
    mb = Minibatches[k]
    "LR: $lr MB: $mb" |> println
    Qtst = SharedArray{Float64}(M,N)
    Ptst = SharedArray{Float64}(M,N)
    LossTrain = SharedArray{Float64}(M)
    LossTest = SharedArray{Float64}(M)

    # Using Leave-one-out cross validation
    LeaveOneOut = kfolds((Z,P); k = M) |> enumerate |> collect
    @sync @distributed for l in LeaveOneOut
        (l,((ztrn,ptrn),(ztst,ptst))) = l
        "training $l..."|>println
        
        # Initializing cNODE model
        cnode = getModel(N)
        # Training cNODE
        W, loss_train, loss_val, loss_test = train_reptile(
                                                cnode, max_epochs,
                                                mb, lr,
                                                ztrn, ptrn, ztst, ptst, ztst, ptst,
                                                report, early_stoping
                                             )
        # Saving results to corresponding arrays
        Ptst[l,:] = ptst
        Qtst[l,:] = predict(cnode,ztst)
        LossTrain[l] = loss_train[end]
        LossTest[l] = loss_test[end]
        # Report
        println(l,'\t',loss_train[end],'\t',loss_test[end])
        println('#'^30)
    end

    # Saving results
    results = "./hyperparameters/"
    !ispath(results) && mkpath(results)
    writedlm(results*"real_sample_$(j)$(k).csv", Ptst, ',')
    writedlm(results*"pred_sample_$(j)$(k).csv", Qtst, ',')
    writedlm(results*"test_loss_$(j)$(k).csv",   LossTest, ',')
    writedlm(results*"train_loss_$(j)$(k).csv",  LossTrain, ',')
    
end

################################################################################
#       2. MAIN TRAINING LOOP
################################################################################

# Selecting hyperparameters
results = "./hyperparameters/"
# 1) calculating means losss for all test runs
# 2) selecting the index of the the one with minimum value (argmin)
_mean = [ mean(readdlm(results*"test_loss_$i$j.csv",',',Float64,'\n') ) for i in 1:6, j in 1:3] |> argmin
# the index of this least mean will be (index of best lr, index of best bs)
mb = Minibatches[_mean[2]]
lr = LearningRates[_mean[1]]

# Allocating variables
Qtst = SharedArray{Float64}(M,N)
Ptst = SharedArray{Float64}(M,N)
LossTrain = SharedArray{Float64}(M)
LossTest = SharedArray{Float64}(M)
  
# Runing validation
LeaveOneOut = kfolds((Z,P); k = M) |> enumerate |> collect
"training $DATA..." |> println
  
@sync @distributed for l in LeaveOneOut
    (l,((ztrn,ptrn),(ztst,ptst))) = l
    "training $l..."|>println
    
    # Get cNODE
    cnode = getModel(N)
    
    # Train cNODE
    W, loss_train, loss_val, loss_test = train_reptile(
                                            cnode, max_epochs,
                                            mb, lr,
                                            ztrn, ptrn, ztst, ptst, ztst, ptst,
                                            report, early_stoping
                                        )
    # Save values
    Ptst[l,:] = ptst
    Qtst[l,:] = predict(cnode,ztst)
    LossTrain[l] = loss_train[end]
    LossTest[l] = loss_test[end]
    
    # Save realization
    !ispath(results*"loss_epochs/") && mkpath(results*"loss_epochs/")
    writedlm(results*"loss_epochs/train$l.csv",loss_train, ',')
    writedlm(results*"loss_epochs/test$l.csv",loss_test, ',')
    writedlm(results*"loss_epochs/val$l.csv",loss_val, ',')
    
    # Report
    println(i,'\t',loss_train[end],'\t',loss_test[end])
    # Saving weights to the bson file
    weights = Tracker.data.(Flux.params(cnode));
    @save "cnode_healthy_weights.bson" weights

end

    
# Write full results
writedlm(results*"real_sample.csv",Ptst, ',')
writedlm(results*"pred_sample.csv",Qtst, ',')
writedlm(results*"test_loss.csv",LossTest, ',')
writedlm(results*"train_loss.csv",LossTrain, ',')
