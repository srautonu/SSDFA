/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author pritom
 */
import java.io.*;
import java.util.*;
import suffixtree.*;

public class ScatterSearch {

    public static void main(String args[]) {
        String strResultFile = args[1];

        long fileReadStartTime =System.currentTimeMillis();
        GeneralizedSuffixTree tree = ReadAndStoreFragments(args[0]);
        int[][] overlapArray = tree.allPairSuffixPrefix();

        //
        // SIZE is total number of fragments in the tree, divided by 2. Because
        // Each fragment is kept twice in the tree; once in forward direction
        // another one as RC (reverse complement)
        //
        int SIZE = tree.getStringCount() / 2;

        long threshold = 0;
        long totalLength = 0;
        long totalOverlap = 0;

        //
        // Default configuration. Can be overridden by using
        // configuration file (config.txt). The file should
        // be placed in the current directory (i.e. the directory
        // from where the ScatterSearch class gets loaded).
        //
        int n = 30;
        int m = 20;
        int nHCIter = 500;
        double thresholdWeight = 0.20;
        int totalTime = 200;
        String diversityMeasure = "PDistance";

        int popSize = n + m;

        try {
            Properties prop = new Properties();

            prop.load(new FileInputStream("config.txt"));

            n = Integer.parseInt(prop.getProperty("n"));
            m = Integer.parseInt(prop.getProperty("m"));
            popSize = n + m;
            nHCIter = Integer.parseInt(prop.getProperty("HCIteration"));
            thresholdWeight = Double.parseDouble(prop.getProperty("thresholdweight"));
            totalTime = Integer.parseInt(prop.getProperty("totaltime"));
            diversityMeasure = prop.getProperty("diversity");
            if (0 == diversityMeasure.compareTo("PDistance"))
                Utility.divMeasure = Utility.DivMeasure_PDistance;
            else
                Utility.divMeasure = Utility.DivMeasure_HamDistance;
            System.out.println("Read configurations from config.txt.");
        } catch (Exception e) {
            System.out.println("Using hard-coded default configuration.");
        }

        //
        // We set the threshold to be thresholdWeight fraction of mean fragment length
        //
        for (int i = 0; i < SIZE; i++) {
            totalLength += tree.getString(i).length();
        }
        threshold = (int)(thresholdWeight * totalLength / SIZE);
        System.out.println("Overlap threshold (Auto-tuned):  " + threshold);

        long fileReadEndTime = System.currentTimeMillis();

        System.out.println("Fragments, config read and store time: " + (fileReadEndTime - fileReadStartTime) + " ms");

        //
        // Setup 10 random numbers as seeds for 10 runs
        //
        Utility.rng.setSeed(1);
        long[] runRandSeeds = new long[10];
        for (int i = 0; i < runRandSeeds.length; i++)
            runRandSeeds[i] = Utility.rng.nextLong();

        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(strResultFile, true));

            bw.newLine();
            bw.append(new Date().toString());
            bw.newLine();
            bw.append("n = " + n);
            bw.newLine();
            bw.append("m = " + m);
            bw.newLine();
            bw.append("HCIteration = " + nHCIter);
            bw.newLine();
            bw.append("thresholdweight = " + thresholdWeight);
            bw.newLine();
            bw.append("totaltime = " + totalTime);
            bw.newLine();
            bw.append("diversity = " + diversityMeasure);
            bw.newLine();
            bw.close();
        } catch (Exception e) {

        }

        for (int runId = 1; runId <= 10; runId++)
        {
            long startTime = System.currentTimeMillis();

            System.out.println("RUN " + runId);
            System.out.println("-------------");

            // Set random number seed
            Utility.rng.setSeed(runRandSeeds[runId - 1]);

            int[][] seededPopulation = generateSeededPopulation(n, SIZE, overlapArray, nHCIter, threshold);
            System.out.println("Generated Seeded Population");

            int[][] population = generateDiversePopulation(m, SIZE, seededPopulation);
            System.out.println("Generated Diverse Population");
            double[] fitnessArray = new double[popSize];
            double[] diversityArray = new double[popSize];

            int[] BEST = null;
            double fitnessBest = 0;

            for (int i = 0; i < popSize; i++) {
                population[i] = LocalSearch.HillClimbing(overlapArray, population[i], nHCIter, threshold);
                double fitnessI = Utility.fitness(population[i], overlapArray, threshold);

                if (BEST == null || fitnessI > fitnessBest) {
                    BEST = population[i];
                    fitnessBest = fitnessI;
                }
            }

            totalOverlap = Utility.overlapCount(BEST, overlapArray, tree);
            System.out.println("0 " + fitnessBest + " " + fitnessBest + " " + totalOverlap);

            //////////////////////////////////////////////////////////////////////////////////////////////
            //MAIN LOOP STARTS FROM HERE
            /////////////////////////////////////////////////////////////////////////////////////////////
            for (int time = 1; time <= totalTime; time++) {
                double fitnessRunBest = 0;

                //Evaluating fitness and Diversity
                fitnessArray = new double[population.length];
                diversityArray = new double[population.length];

                for (int i = 0; i < population.length; i++) {
                    fitnessArray[i] = Utility.fitness(population[i], overlapArray, threshold);
                    diversityArray[i] = Utility.diversity(population, population[i], population.length);
                }

                int[][] B = getFitIndividuals(population, fitnessArray, n);
                int[][] D = getDiverseIndividuals(population, diversityArray, m);

                //Generating Parents for new Population
                for (int i = 0; i < n; i++) {
                    population[i] = B[i];
                }

                for (int i = n; i < popSize; i++) {
                    population[i] = D[i - n];
                }

                //
                // For each pair, there will be 2 children generated. So, the size of newly generated
                // set of individuals would be: popSize * (popSize - 1)
                // Plus, the existing set size is popSize.
                //
                // So, the total size of new population: popSize * popSize
                //
                int[][] newPopulation = new int[popSize * popSize][];

                for (int i = 0; i < popSize; i++) {
                    newPopulation[i] = population[i];
                }

                for (int i = 0, k = popSize; i < popSize; i++) {
                    for (int j = 0; j < i; j++) {
                        int[][] children = Utility.crossover(population[i], population[j]);

                        //
                        // Normal mutation operation
                        //
                        children[0] = Utility.mutation(children[0]);
                        children[1] = Utility.mutation(children[1]);

                        //
                        // Exploitative mutation operation
                        //
                        children[0] = Utility.mutation(children[0], threshold, overlapArray);
                        children[1] = Utility.mutation(children[1], threshold, overlapArray);

                        children[0] = LocalSearch.HillClimbing(overlapArray, children[0], nHCIter, threshold);
                        children[1] = LocalSearch.HillClimbing(overlapArray, children[1], nHCIter, threshold);

                        double fitness0 = Utility.fitness(children[0], overlapArray, threshold);
                        double fitness1 = Utility.fitness(children[1], overlapArray, threshold);

                        //
                        // best fitness in the current run
                        //
                        if (fitness0 > fitnessRunBest) {
                            fitnessRunBest = fitness0;
                        }
                        if (fitness1 > fitnessRunBest) {
                            fitnessRunBest = fitness1;
                        }

                        if (fitness0 > fitnessBest) {
                            BEST = children[0];
                            fitnessBest = fitness0;
                        }
                        if (fitness1 > fitnessBest) {
                            BEST = children[1];
                            fitnessBest = fitness1;
                        }

                        newPopulation[k++] = children[0];
                        newPopulation[k++] = children[1];
                    }
                }

                population = newPopulation;

                totalOverlap = Utility.overlapCount(BEST, overlapArray, tree);

                long nContig = Utility.contigCount(BEST, overlapArray, threshold);
                long endTime = System.currentTimeMillis();
                long duration = (endTime - startTime) +  (fileReadEndTime - fileReadStartTime);
                System.out.println(time + " " + fitnessRunBest + " " + fitnessBest + " " + totalOverlap + " " + nContig + " " + duration);
            }

            try {
                BufferedWriter bw = new BufferedWriter(new FileWriter(strResultFile, true));

                long nContig = Utility.contigCount(BEST, overlapArray, threshold);
                long endTime = System.currentTimeMillis();
                long duration = (endTime - startTime) +  (fileReadEndTime - fileReadStartTime);

                bw.append(runId + " " + fitnessBest + " " + totalOverlap + " " + nContig + " " + duration);
                bw.newLine();

                bw.close();
            } catch (Exception e) {

            }

        }
        //////////////////////////////////////////////////////////////////////////////////////////////
        //MAIN LOOP ENDS HERE
        /////////////////////////////////////////////////////////////////////////////////////////////
    }

    public static GeneralizedSuffixTree ReadAndStoreFragments(String strFragmentFile)
    {
        GeneralizedSuffixTree tree = new GeneralizedSuffixTree();
        String strFragment = "";
        int nOrigStr;

        try
        {
            BufferedReader br = new BufferedReader(new FileReader(strFragmentFile));
            for(String line; (line = br.readLine()) != null; ) {
                if (line.charAt(0) == '>')
                {
                    if (!strFragment.isEmpty())
                    {
                        tree.put(strFragment);
                        strFragment = "";
                    }
                }
                else
                {
                    strFragment += line;
                }
            }
            tree.put(strFragment);

            br.close();
        }
        catch (Exception e)
        {
            System.out.println(e);
        }

        //
        // Now add the RC strings
        //
        nOrigStr = tree.getStringCount();
        for (int i = 0; i < nOrigStr; i++)
        {
            String strOrg = tree.getString(i);
            StringBuilder strRC = (new StringBuilder(strOrg)).reverse();
            for (int j = 0; j < strRC.length(); j++)
            {
                char rc = '\0';
                switch (strRC.charAt(j))
                {
                    case 'A':
                        rc = 'T';
                        break;
                    case 'T':
                        rc = 'A';
                        break;
                    case 'C':
                        rc = 'G';
                        break;
                    case 'G':
                        rc = 'C';
                        break;
                }

                strRC.setCharAt(j, rc);
            }

            tree.put(strRC.toString());
        }
        return tree;
    }

    public static int[][] generatePopulation(int popoulationSize, int individualSize) {

        int[][] population = new int[popoulationSize][];

        for (int i = 0; i < popoulationSize; i++) {
            population[i] = Utility.generateIndividual(individualSize);
        }

        return population;

    }

    public static int[][] generateSeededPopulation(int popoulationSize, int individualSize, int[][] ov, long maxIter, long threshold) {

        int[][] population = new int[popoulationSize][];

        for (int i = 0; i < popoulationSize; i++) {
            population[i] = LocalSearch.HillClimbing(ov, Utility.generateIndividual(individualSize), maxIter, threshold);
        }

        return population;

    }

    public static int[][] generateDiversePopulation(int diversePopulationSize, int individualSize, int[][] seed) {

        double maxDiversity = 0;
        int seedSize = seed.length;

        for (int i = 0; i < seed.length; i++) {
            double diversity = Utility.diversity(seed, seed[i], seedSize);

            if (diversity > maxDiversity) {
                maxDiversity = diversity;
            }
        }

        maxDiversity /= 2;

        int[][] population = new int[diversePopulationSize + seed.length][];

        for (int i = 0; i < seedSize; i++) {
            population[i] = seed[i];
        }

        for (int i = 0; i < diversePopulationSize; i++) {

            double diversity = 0;
            int[] diverseIndividual;

            do {
                diverseIndividual = Utility.generateIndividual(individualSize);
                diversity = Utility.diversity(population, diverseIndividual, seedSize + i);
            } while (diversity < maxDiversity);

            population[seedSize + i] = diverseIndividual;
        }

        return population;
    }

    private static int[][] getTournamentSelectedIndividuals(int [][] population, final double[] prop, int n)
    {
        int[][] selectedPopulation = new int[n][];
        for(int i = 0; i < n; i++)
        {
            selectedPopulation[i] = Utility.tournamentSelection(population, 7, prop);
        }

        return selectedPopulation;
    }

    private static int[][] getTopRankedIndividuals(int[][] population, final double[] prop, int n) {
        //
        // Inner class for index-sorting prop array in descending order
        // the indices array (declaration follows) holds the indices.
        //
        Comparator<Integer> comparator = new Comparator<Integer>() {
            public int compare(Integer arg0, Integer arg1) {

                Double a = new Double(prop[arg0]);
                Double b = new Double(prop[arg1]);

                int cmp = b.compareTo(a);
                if (cmp == 0)
                    return arg0.compareTo(arg1);

                return cmp;
            }
        };

        int[][] selectedPopulation = new int[n][];
        Integer[] indices = new Integer[prop.length];

        for (int i = 0; i < indices.length; i++) {
            indices[i] = i;
        }
        Arrays.sort(indices, comparator);

        for (int i = 0; i < n; i++) {
            selectedPopulation[i] = population[indices[i]];
        }

        return selectedPopulation;
    }

    public static int[][] getFitIndividuals(int[][] population, double[] fitnessArray, int n)
    {

        return getTopRankedIndividuals(population, fitnessArray, n);
        //return getTournamentSelectedIndividuals(population, fitnessArray, n);
    }

    public static int[][] getDiverseIndividuals(int[][] population, final double[] diversityArray, int m)
    {
        return getTopRankedIndividuals(population, diversityArray, m);
    }
}
