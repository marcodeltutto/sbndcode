import sys

from tempfile import mkstemp
from shutil import move
from os import fdopen, remove
from bayes_opt import BayesianOptimization
from bayes_opt import UtilityFunction
import os

# Bounded region of parameter space
bdt_paramterbounds = { 'float NTrees': (5,1800),
                       'float MinNodeSize': (1,25),
                       'float Shrinkage': (0.005,0.9),
                       'float BaggedSampleFraction': (0.1,0.9),
                       'float MaxDepth': (2,6)
                      }


bdt_paramtertypes = {'float NTrees': "int",
                     'float MinNodeSize': "float", 
                     'float Shrinkage': "float",
                     'float BaggedSampleFraction': "float", 
                     'float MaxDepth': "int"
}

   
bdt_filename="./NueRecoOptimiser.cc"
equals="="
dot="."
newline="\n"
semi=";"

def black_box_function(StartFitSize, y):
    """Function with unknown internals we wish to maximize.

    This is just serving as an example, for all intents and
    purposes think of the internals of this function, i.e.: the process
    which generates its output values, as unknown.
    """
    return -StartFitSize ** 2 - (y - 1) ** 2 + 1


def replace_bdt_param(bdt_param, bdt_param_name):

    replace_bdt_param_name=bdt_param_name+equals

    if(bdt_paramtertypes[bdt_param_name] == "int"):
        bdt_param=bdt_param.round()

    #Create temp file
    fh, abs_path = mkstemp()
    with fdopen(fh,'w') as new_file:
        with open(bdt_filename) as old_file:
            for line in old_file:
                new_bdt_param_name=line
                if(replace_bdt_param_name in line):
                    equal_loc = line.find('=')
                    new_bdt_param_name=line[:equal_loc+1]+str(bdt_param) + semi + newline
                    line=line[equal_loc+1:]

                new_file.write(line.replace(line,new_bdt_param_name))

    #Remove original file
    remove(bdt_filename)
    #Move new file
    move(abs_path,bdt_filename)


def function(**fcl_parms):
    
    #Change the fcl config 
    for param in fcl_parms: 
        replace_bdt_param(fcl_parms.get(param),param)
        
    #Run the bash script required script must return the value to be maximised.
    bashCommand = "source ./runbdt.sh"
    os.system(bashCommand)

    fileHandle = open ( 'histcomp.txt',"r" )
    lineList = fileHandle.readlines()
    fileHandle.close()
    target=float(lineList[-1])

    return target


optimizer = BayesianOptimization(
    f=function,
    pbounds=bdt_paramterbounds,
    verbose=2, # verbose = 1 prints only when a maximum is observed, verbose = 0 is silent
    random_state=2,
)

#i=0;

optimizer.probe(
    params={'float NTrees': 600,
            'float MinNodeSize': 2.5,
            'float Shrinkage': 0.1,
            'float BaggedSampleFraction': 0.5,
            'float MaxDepth': 2
        }
)

optimizer.maximize(
    init_points=30,
    n_iter=1000,
    acq='ucb',
    kappa=2,
)



# utility = UtilityFunction(kind="ucb", kappa=2, xi=1e-2)


# for _ in range(1000):

#     #Get the suggested next parameters 
#     next_point_to_probe = optimizer.suggest(utility)
#     print("Next point to probe is:", next_point_to_probe)

#     #Change the bdt config 
#     for param in next_point_to_probe: 
#         replace_bdt_param(next_point_to_probe.get(param),param)

#     #Run the bash script required script must return the value to be maximised.
#     bashCommand = "source ./runbdt.sh"
#     os.system(bashCommand)

#     #process = subprocessPopen(bashCommand.split(), stdout=subprocess.PIPE)
#     #target, error = process.communicate()
#     #target = black_box_function(**next_point_to_probe)
#     fileHandle = open ( 'histcomp.txt',"r" )
#     lineList = fileHandle.readlines()
#     fileHandle.close()
#     target=float(lineList[-1])
    
#     #Add to the opttimsation processes 
#     optimizer.register(
#         params=next_point_to_probe,
#         target=target,
#     )
#     print("Found the target value to be:", target)

#     print(target, next_point_to_probe)

#     print(optimizer.max)

#     #i = i + 1
    
#     #print("this is i", i)

print(optimizer.max)

