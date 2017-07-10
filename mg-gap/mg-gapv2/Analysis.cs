using System;
using System.Data;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
namespace mg_gapv2
{
    public class Analysis
    {
        public static List<string> GetB(string filepath)
        {
            //set up list that will contain all the B values
            List<string> bList = new List<string>();
            //parameter for minimum reads
            int min_reads = 20;
            //Set up for snp specific variance
            double Var_snp_specific = 0.0;
            //set up list for holding intermediates for B standard
            List<string> accepted_snps = new List<string>();
            //Set up list for holding z score raw data
            List<double> zraw = new List<double>();


            //start processing
            using (var fileStream = File.OpenRead(filepath))
            using (var streamReader = new StreamReader(fileStream))
            {
                String line;
                while ((line = streamReader.ReadLine()) != null)
                {
                    string[] cols = line.Replace("\n", "").Split('\t');
                    if (!cols.Contains("CHROM")) //weed out header row
                    {
                        if (Convert.ToInt16(cols[0]) < 15)//if chr 1-14 only
                        {
                            int position = Convert.ToInt32(cols[1]);
                            double C_count = 0;
                            List<double> C_dat = new List<double> { };
                            double T_count = 0;
                            List<double> T_dat = new List<double> { };
                            //change this line for what to pool. C will be CA,CB,S
                            for (int i = 2; i < 6; i++)
                            {
                                string[] AD = cols[i].Split(',');
                                if (Convert.ToInt16(AD[0]) + Convert.ToInt16(AD[1]) > 0)
                                {
                                    if (i < 4) //pool of CA, CB, and S
                                    {
                                        C_count++;
                                        C_dat.Add(Convert.ToInt16(AD[0]));
                                        C_dat.Add(Convert.ToInt16(AD[1]));
                                    }
                                    else //pool of TA, TB
                                    {
                                        T_count++;
                                        T_dat.Add(Convert.ToInt16(AD[0]));
                                        T_dat.Add(Convert.ToInt16(AD[1]));
                                    }

                                }
                            }

                            //So we put the separated allele depths into the two pools of C and T for the SNP
                            if (C_count > 0 && T_count > 0)//must have something to work with
                            {
                                //reset the variables here
                                double qC_0 = 0.0;
                                double qC_1 = 0.0;
                                double qT_0 = 0.0;
                                double qT_1 = 0.0;

                                for (int i = 0; i < C_dat.Count / 2; i++)
                                {
                                    double m = C_dat[2 * i] + C_dat[2 * i + 1];
                                    qC_0 += C_dat[2 * i]; //ref base
                                    qC_1 += m; //ref + alt base
                                }

                                for (int i = 0; i < T_dat.Count / 2; i++)
                                {
                                    double m = T_dat[2 * i] + T_dat[2 * i + 1];
                                    qT_0 += T_dat[2 * i];
                                    qT_1 += m;
                                }

                                if (qT_1 >= min_reads && qC_1 >= min_reads)
                                {

                                    //qC_hat and qT_hat should be different almost all of the time
                                    double qC_hat = qC_0 / qC_1;
                                    double qT_hat = qT_0 / qT_1;

                                    if (qC_hat != qT_hat) //this catches any time qC_hat and qT_hat are the same - which shouldn't really happen and so throw them out
                                    {
                                        //moved down one logic level to prevent num_snp's from being off in final calculations
                                        double var_C = 1.0 / qC_1;
                                        double var_T = 1.0 / qT_1;
                                        Var_snp_specific += (var_C + var_T);

                                        //note double provides more precision but some of the Math class methods only take double - so there may be some small loss of precision past n^-16 
                                        double diverge = (2.0 * ((Math.Asin(Math.Pow(qT_hat, 0.5))) - (Math.Asin(Math.Pow(qC_hat, 0.5)))));
                                        if (diverge == 0)
                                        {
                                            Console.WriteLine(" Diverge = " + diverge + " " + cols[0] + "_" + cols[1]);
                                        }

                                        //raw_results.Add(cols[0] + '\t' + cols[1] + qC_1.ToString() + '\t' + qC_hat.ToString()
                                        //    + '\t' + qT_1.ToString() + '\t' + qT_hat.ToString() + '\t' + diverge.ToString() + '\n');

                                        Random rnd = new Random();
                                        if (rnd.Next(1, 3) == 1)
                                        {
                                            zraw.Add(diverge);
                                        }
                                        else
                                        {
                                            zraw.Add(-diverge);
                                        }

                                        double t_C_pop = Math.Asin(Math.Sqrt(qC_hat));
                                        double t_T_pop = Math.Asin(Math.Sqrt(qT_hat));
                                        //SNP accepted.
                                        accepted_snps.Add(cols[0] + "_" + cols[1] + '\t' + var_C + '\t' + var_T + '\t' + t_C_pop.ToString() + '\t' + t_T_pop.ToString());
                                    }
                                }
                            }
                        }
                    }
                }
                //wrap up
                Console.WriteLine("SNP processing complete, starting analysis...");
                Console.WriteLine("Accepted " + accepted_snps.Count + " SNPs");

                Console.WriteLine("Sampling/genotyping varience " + (Var_snp_specific / Convert.ToDouble(accepted_snps.Count())));

                //make the z, rankings
                //zraw UNSORTED is used for b so you have to do the ranked_z in another array
                //confirm through console the length is the same
                zraw.Sort();

                var n25 = zraw[accepted_snps.Count() / 4];
                var n50 = zraw[accepted_snps.Count() / 2];
                var n75 = zraw[3 * accepted_snps.Count() / 4];
                Console.WriteLine("Z percentiles (without direction) " + n25 + " " + n50 + " " + n75);
                Console.WriteLine("Total variance in Z (based on IQR) " + Math.Pow((n75 - n25 / 1.349), 2));

                double Var_neutral = (Math.Pow(((n75 - n25) / 1.349), 2)) - Var_snp_specific / accepted_snps.Count();
                Console.WriteLine("Bulk sampling and library variance " + Var_neutral);

                //make an unsorted B list - the output of this is the last step if just doing B
                List<double> b_std = new List<double>();
                for (int k = 0; k < zraw.Count; k++)
                {
                    string unparsed_line = accepted_snps[k].ToString();
                    string[] parsedarray = unparsed_line.Split('\t').ToArray();
                    double vdiv = Var_neutral + Convert.ToDouble(parsedarray[1]) + Convert.ToDouble(parsedarray[2]); // 0 is snnffold_X, 1 is var_C and 2 is var_T
                    double b = Math.Pow(Convert.ToDouble(zraw[k]), 2) / vdiv;
                    b_std.Add(b);
                    bList.Add("snnffold_" + parsedarray[0].ToString() + '\t' + b.ToString() + '\n');
                }


                return bList;
            }

        }
    }
}
