# SSDFA  
This is the repository for Scatter Search for DNA Fragment Assembly tool (SS-DFA).

# Usage
java ScatterSearch input_fragment_file output_file

The output file will contain data from 10 runs. Each line contains:
RunId fitness totalOverlap ContigCount WallClockTime

# configuration
Default configuration is already hard-coded in the source code. The default (tuned) values are specified in the paper. Certain configurations can be overridden using a config.txt file. The file should be placed in the current directory (i.e. the directory from where the ScatterSearch class gets loaded).

A sample configuration file is available in the /Code/SSDFA/src folder. The sample also uses default values, but can be changed by the user as needed.

The parameters "Selection strategy for fit individuals" and "Fitness function" cannot be changed using the configuration file. We made changes in the code directly when we experimented with different selection strategies and fitness functions. We do not recommend users changing these. However, if you must, you can go into the source code and make necessary changes.

To use traditional fitness function instead of our proposed one:
- Open up the code for Utility.fitness().
- Uncomment the line: //fitness = totalOverlap;
- Comment out the line right above it: fitness=((double)totalOverlap)/contigCount;

To use tournament selection instead of truncation selection for fit individuals:
- Open up the code for ScatterSearch.getFitIndividuals()
- Uncomment the line: //return getTournamentSelectedIndividuals(population, fitnessArray, n);
- Comment out the line right above it: return getTopRankedIndividuals(population, fitnessArray, n);

To enable exploitative mutation in hill-climbing (HC):
- Open up the code for LocalSearch.HillClimbing()
- Uncomment the line: //int[] r = Utility.mutation(best, threshold, ov);
- Comment out the line right above it: int[] r = Utility.mutation(best);

