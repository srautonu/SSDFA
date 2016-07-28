/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author pritom
 */
public class LocalSearch {

    public static int[] HillClimbing(int[][] ov, int[] individual, long nNoChangeLimit, long threshold)
    {
        int nNoChange = 0;
        int[] best = individual.clone();
        double fitnessBest = Utility.fitness(best, ov, threshold);

        while(nNoChange < nNoChangeLimit)
        {
            int[] r = Utility.mutation(best, threshold, ov);
            //int[] r = Utility.mutation(best);
            double fitness = Utility.fitness(r, ov, threshold);

            if(fitness > fitnessBest)
            {
                nNoChange = 0;
                best = r.clone();
                fitnessBest = fitness;
            }
            else
            {
                nNoChange++;
            }
        }

        return best;
    }
}
