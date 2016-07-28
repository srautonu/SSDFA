
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import suffixtree.GeneralizedSuffixTree;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author pritom
 */
public class Utility {
    static final Random rng = new Random(1);

    static final int DivMeasure_HamDistance = 0;
    static final int DivMeasure_PDistance = 1;

    static int divMeasure = DivMeasure_HamDistance;

    public static int [] mutation(int individual[])
    {
        int swapIndex1;
        int swapIndex2;
        
        int length = individual.length;
        
        swapIndex1 = rng.nextInt(length);

        while(swapIndex1 == (swapIndex2 = rng.nextInt(length)));

        int[] mutatedIndividual = individual.clone();
        int temp = mutatedIndividual[swapIndex1];
        mutatedIndividual[swapIndex1] = mutatedIndividual[swapIndex2];
        mutatedIndividual[swapIndex2] = temp;

        //
        // Make the swapped fragments reverse complement with 50% probability each.
        //
        if (rng.nextBoolean())
            mutatedIndividual[swapIndex1] = (mutatedIndividual[swapIndex1] + length) % (2*length);

        if (rng.nextBoolean())
            mutatedIndividual[swapIndex2] = (mutatedIndividual[swapIndex2] + length) % (2*length);

        return mutatedIndividual;
    }
        
    public static int [] mutation(int individual[],long threshold, int[][] overlapArray)
    {
        int length = individual.length;
        boolean[] flagArray = new boolean[length];
        
        int count = 0;
        if(overlapArray[individual[0]][individual[1]] < threshold)
        {
            flagArray[0] = true;
            count++;
        }
        
        if(overlapArray[individual[length - 2]][individual[length - 1]] < threshold)
        {
            flagArray[length - 1] = true;
            count++;
        }
        
        for(int i = 1; i < length-1; i++)
        {
            if(overlapArray[individual[i-1]][individual[i]] <threshold
                && overlapArray[individual[i]][individual[i+1]] < threshold)
            {
                flagArray[i] = true;
                count++;
            }
        }
        
        List<Integer> changeList = new ArrayList<>(count);
        
        for(int i = 0; i < length; i++)
        {
            if(flagArray[i] == true)
            {
                changeList.add(individual[i]);
            }
        }
        
        Collections.shuffle(changeList, rng);
        
        int[] mutatedIndividual = individual.clone();
        
        for(int i = 0, j = 0; i < length; i++)
        {
            if(flagArray[i] == true)
            {
                mutatedIndividual[i] = changeList.get(j);
                j++;
            }
        }
        
        return mutatedIndividual;
    }
    
    public static int [][] crossover(int [] parent0, int [] parent1)
    {
        int children[][] = new int[2][];
        int length = parent0.length;
        
        
        children[0] = new int[length];
        children[1] = new int[length];
        
        boolean[] crossoverFlag = new boolean[length]; //To choose parts of a parent

        for(int i = 0; i < length; i++)
        {
            crossoverFlag[i] = rng.nextBoolean();
        }
        
        boolean[] presentFlag0 = new boolean[length];
        boolean[] presentFlag1 = new boolean[length];
        
        ArrayList<Integer> queue0 = new ArrayList<Integer>();
        ArrayList<Integer> queue1 = new ArrayList<Integer>();
        
        for(int i = 0; i < length; i++)
        {
            if(crossoverFlag[i])
            {
                children[0][i] = parent0[i];
                presentFlag0[parent0[i] % length] = true;
                queue1.add(i);
            }
            else
            {
                children[1][i] = parent1[i];
                presentFlag1[parent1[i] % length] = true;
                queue0.add(i);
            }
        }

        for (int i = 0; !queue0.isEmpty(); i++)
        {
            if(presentFlag0[parent1[i] % length] == false)
            {
                children[0][queue0.remove(0)] = parent1[i];
            }
        }

        for (int i = 0; !queue1.isEmpty(); i++)
        {
            if(presentFlag1[parent0[i] % length] == false)
            {
                children[1][queue1.remove(0)] = parent0[i];
            }
        }
        
        return children;
    }

    public static double fitness(int[] individual, int[][] overlapArray, long threshold)
    {
        double fitness = 0.0;
        long contigCount = 1;
        long totalOverlap = 0;
        
        for(int i = 0; i < individual.length - 1; i++)
        {
            int overlapCount = overlapArray[individual[i]][individual[i+1]];
            
            if(overlapCount <= threshold)
            {
                contigCount++;
            }
            else
            {
                totalOverlap += overlapCount;
            }
        }
        
        fitness=((double)totalOverlap)/contigCount;
        //fitness = totalOverlap;
        
        return fitness;
    }
    
    public static double diversity(int [][] population, int [] individual, int currentPopSize)
    {
        double dist=0;
        //double size=population.length;

        if (divMeasure == DivMeasure_PDistance) {
            for (int i = 0; i < currentPopSize; i++) {
                dist += PDistance(population[i], individual);
            }
        }
        else {
            for (int i = 0; i < currentPopSize; i++) {
                dist += hammingDistance(population[i], individual);
            }
        }

        return dist/currentPopSize;
    }
    
    public static double hammingDistance(int [] individual1, int [] individual2)
    {
        int size=individual1.length;
        double dist=0;
        
        for(int i=0; i < size; i++)
        {
            if(individual1[i] != individual2[i])
            {
                dist++;
            }
        }
        
        return dist / size;
    }

    public static double PDistance(int [] individual1, int [] individual2)
    {
        int dist = 0;

        dist += PDistanceAsymmetric(individual1, individual2);
        dist += PDistanceAsymmetric(individual2, individual1);

        return dist / 2.0;
    }

    public static double PDistanceAsymmetric(int [] individual1, int [] individual2)
    {
        int i, center, prev, next, dist;
        int size = individual1.length;

        int[] index2 = new int[2*size];
        for (i = 0; i < index2.length; i++)
        {
            //
            // Initialize each fragment's (and RC complement's) position to -1.
            // Note that we have 2n strings (regular + RC), but the individual has
            // n strings. So, half of the strings will not have a index in the
            // individual
            //
            index2[i] = -1;
        }

        for(i = 0; i < size; i++)
        {
            index2[individual2[i]] = i;
        }

        for(i = 0, dist = 0; i < size; i++)
        {
            center = index2[individual1[i]];
            prev = next = center;

            if (i != 0)
                prev = index2[individual1[i - 1]];
            if (i != size -1)
                next = index2[individual1[i + 1]];

            if (i != 0)
            {
                if (center == -1 || prev == -1)
                    dist += size - 1;
                else
                    dist += Math.abs(center - prev);
            }

            if (i != size - 1)
            {
                if (center == -1 || next == -1)
                    dist += size - 1;
                else
                    dist += Math.abs(center - next);
            }
        }

        dist -= 2*(size - 1);
        dist /= 2*(size - 1)*(size - 2);

        return dist;
    }
    
    
    public static int [] tournamentSelection(int[][] population, int tournamentSize, double[] propArray)
    {
        int populationSize = population.length;
        int index = rng.nextInt(populationSize);

        int []s = population[index];
        double propS = propArray[index];

        while(tournamentSize > 0)
        {
            index = rng.nextInt(populationSize);
            int[] r = population[index];

            double propR = propArray[index];

            if(propR > propS)
            {
                s = r;
                propS = propR;
            }

            tournamentSize--;
        }

        return s;
    }
        
    public static long contigCount(int[] individual, int[][] overlapArray, long threshold)
    {
        long contigCount = 1;
        
        for(int i = 0; i < individual.length - 1; i++)
        {
            int overlapCount = overlapArray[individual[i]][individual[i+1]];
            
            if(overlapCount <= threshold) {
                contigCount++;
            }
        }
        
        return contigCount;
    }

    public static long overlapCount(int[] individual, int[][] overlapArray, GeneralizedSuffixTree tree)
    {
        long totalOverlap = 0;
        for (int i = 0; i < individual.length - 1; i++) {
            int ov = overlapArray[individual[i]][individual[i + 1]];
            totalOverlap += ov;
        }

        if (null != tree)
        {
            //
            // validate that the pairwise overlaps as mentioned in the
            // overlap array are indeed accurate
            //
            for (int i = 0; i < individual.length - 1; i++) {
                int ov = overlapArray[individual[i]][individual[i + 1]];

                String strCurr = tree.getString(individual[i]);
                String strNext = tree.getString(individual[i + 1]);

                for (int j = strCurr.length() - ov, k = 0; j < strCurr.length(); j++, k++) {
                    if (strCurr.charAt(j) != strNext.charAt(k)) {
                        System.out.println("ERROR!!! Incorrect overlap size. Fragments: " + i + " " + j);
                    }
                }
            }
        }

        return totalOverlap;
    }

    public static int[] generateIndividual(int SIZE) {

        List<Integer> individualList = new ArrayList<>(SIZE);

        for (int i = 0; i < SIZE; i++) {
            individualList.add(i);
        }
        Collections.shuffle(individualList, rng);

        Integer[] individualI = new Integer[SIZE];
        individualList.toArray(individualI);

        int[] individual = new int[SIZE];
        for (int i = 0; i < SIZE; i++) {
            //
            // The fragment is equally likely to be from any of the strands.
            // We use (fragment index + SIZE) to represent the fragment in
            // reverse complement (RC) strand
            //
            if (rng.nextBoolean())
                individualI[i] += SIZE;

            individual[i] = individualI[i];
        }

        return individual;
    }
}
