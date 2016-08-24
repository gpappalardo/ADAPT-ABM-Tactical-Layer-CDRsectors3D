# ELSA-ABM-Tactical-Layer

The Tactical layer aims at transforming planned trajectories into "actual" trajectories by deconflicting them with a simulated air traffic controller. The standard input is a flight plan whereas the output are the realized routes, similarly to what exists in DDR data. The code of the Tactical ABM is stored in src/. In the config folder there is the config.cfg that contains the values of most of the parameters.

You can directly compile the source code in the following way:

```
cd src
LC_ALL=C gcc -O3 -c *.c && gcc *.o -o ../ElsaABM.so -lm
```

Once the binary ElsaABM.so is produced it is possible to execute the simulations:

```
./ElsaABM.so input_file output_file Config_file seed,
```

Using the files provided as examples in the tests/ folder, the user can run:

```
./ElsaABM.so tests/inputABM.dat tests/outputABM.dat config/config.cfg 4967
```

which should produce the outputABM.dat file with deconflicted trajectories. 

# Input-Output File Format

Either input file and output files have the same format. The first line must contain the number of flights (trajectories contained in the files):
```
1475\tNflight
```

the following lines are the trajectory of each aircraft. In particular the first two positions, separated by a backslash t,  represent respectively the ID, a unique int number that identify the flight, and the number of NVP in the related route. 

The following positions are the nvps of the route. They consist in 5 values, separated by a comma: the horizontal position (latitude and longitude), the Flight Level, the time in the format `2010-05-06 10:20:32', and an integer number that refers to the sector to which the navpoint belongs. The number 0 is the null sector, in fact it refers to a nvp that not belong to any sector.

It will be generate also an outputfile_Counter.dat that contains 4 columns, each line refers to an ATC operation: the first is column is a identifier of the operation (R reroute, D direct, H flight level change), the second column is the timestamps of the operation, the third column is the previous closet nvp where the operation occurs, the fourth column is ID of the flight.

# Configuration File
Most of the parameters are stored in the configuration file. Changing these values do not require to re-compile the code. 
An example of a config file can be found in config/. Each line starting with # is a comment, whereas all other lines are values. The name of the variable it refers to needs to be written after a \t#. For instance:

```
24\t#t_w\n
```

gives t_w = 24. The position of the values in the configuration file does not matter.

Note that the ABM does not perform any consistency checks on the variables. For example if one wants to increase the lookahead t\_w, one needs to fix the product t_w x t_r x t_i which represents the time-step. These three values are the most crucial ones for the simulations and should be chosen with care. The user can find their descriptions in the xxx.
\\
Most of the values of the configuration file can be changed without too much trouble. For example the most simple features are tunable with the parameters nsim (number of simulations).


# Additional files

Other files are required to run the simulations, which are called respectively temp_nvp.dat, shock_tmp.dat,  sector_capacities.dat, and bound_latlon.dat and the boundary of the sectors stored in config/boundary. Some examples are stored in the config folder config/. Their paths have to be specified in the config file, in the corresponding entries.\\

The temp_nvp.dat consists in a two columns text file with latitude and longitude of the temporary points used for rerouting. The number of temporary nvps used in the simulation is defined via a #define variable called NTMP in mSector.h. The user must be careful because the code:
* does not check if the points are within the ACC. 
* does not check the number of points matches NTMP.

The shock_tmp.dat has the same format of temp_nvp.dat and it contains the centers of the shocks used in the simulations. This also has to be created via an external script, which is not included at the moment in the repository. We aware the users that the shock code it was not sufficiently tested.

The sector_capacities.dat has two columns: the first one stores the ID of the sector, the second the related capacity. These labels need to be consistent with the ones present in the input file. Note that the user can use the '0' for navpoints which do not have to be under control. This is more specifically used for the first and the last points of the trajectories.

The file bound_latlon.dat contains the boundaries of the airspace as a list of latitudes/longitudes, like the temp_nvp.dat file. Finally in the config/boundary folder are stored the boundary of the sectors in the format 3_bound_latlon.dat where the number indicates the ID of the sector.

# Define 

Other parameters are defined in the header files. Every time you want to change these parameters you need to recompile the code. For example in mABM.h is defined N_TRY that is the number of reshuffling allowed in case of unsolved conflicts:

```
 #define N_TRY 50 1000
```

Another useful #define is DTMP_P that is the Maximum distance in meters of the selected temporary points for the rerouting operations.

Others #define do not have associated values. For example:

```
#define SINGLE_TOUCH
```

If you comment this #define the ABM can perform a multiple modification of the route of the same aircraft in the same time-step if it does not find any solution.


#Caveats

Here are some potential traps for new users concerning the code:

* The time step t_s is the one the most important parameter of the model. It is equal to t_s = t_i x t_w x t_r:
* The time increment t_i is the time resolution for the trajectories. It should be very small (around 10 seconds),
* The time window t_w represents the time horizon of the controller on which it will compute the future potential conflicts.
* The time roll t_r, which is a fraction of the time window, represents the time after which the controller updates trajectories of the flights.

In particular, a user who would like to increase the time horizon of the controller would need to keep the time-step t_s = t_i x t_w x t_r fixed by reducing the time-roll, otherwise different effects will mix up.
* The life duration of shocks is computed with respect to the time step. Hence changing t_i, t_w, or t_r will change the duration of the flights.
* The shocks appear only on 10 flight levels, not on a whole column of air. The variables shock_f_lvl_min and shock_f_lvl_max do not change this fact. Instead, they are fixing the possible interval of flight levels of apparition of the shocks.
* The starting and ending dates have to be informed in the config file at the corresponding lines. They need to be consistent with the trajectories provided. 

