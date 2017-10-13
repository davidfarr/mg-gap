﻿using System;
using System.Data;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace mg_gap
{
    class VCF_Analyzer
    {
        //read the VCF... create a list of SNP's with their info in a method and return the "list" of SNPs back to the main program file

        public static List<SNP> SNP_list(int window, string vcfpath, char bs_go, string chisq_path)
        {
            //arrays
            List<double> zraw = new List<double>(); //this holds Z scores from divergence
            List<double> zranked = new List<double>(); //this will be just for reporting initial stats
            int skip_counter = 0; //to keep track of how many SNPs are skipped
            int total_counter = 0; //keep track of how many SNPs there are total

            //output list
            List<SNP> output_snps = new List<SNP>();
            List<SNP> trimmed_snps = new List<SNP>();//this one is for after B*

            // values
            int minimum_reads = 20;
            double Var_snp_specific = 0.0;
            int num_snps = 0; //keep log of all the snps not just the accepted

            //evaluate each line of the file
            using (var fs = File.OpenRead(vcfpath))
            using (var sr = new StreamReader(fs))
            {
                String line; //the raw text of the line
                while ((line = sr.ReadLine()) != null) //while there is a line to read
                {
                    string[] cols = line.Replace("\n", "").Split('\t');//get rid of newlines and empty spaces, then use the tab delimiter for array
                    if (cols.Length > 2) //sorts out any 2-col contigs
                    {
                        if (cols[0] == "#CHROM") //this is the header row (but not the first row)
                        {
                            foreach (string entity in cols)
                            {
                                Console.WriteLine("Found column {0}", entity);
                            }
                        }
                        else if (!cols[0].Contains("##"))//all contigs and non-SNP lines have this
                        {
                            //increment the total counter, it's a real snp
                            total_counter++;

                            //get the actual chromosome number out of the first column
                            string scaffold = cols[0].Split('_').ToString();
                            string chromosome = cols[0].ToString();
                            chromosome = chromosome.Remove(0, chromosome.IndexOf("_") + 1);//safely remove the non-standard id before the actual #

                            if (Convert.ToInt32(chromosome) < 15)//we are only using up to chr 14 for this data
                            {
                                int position = Convert.ToInt32(cols[1]);//grab the base pair
                                string ref_base = cols[3]; //the reference base allele
                                string alt_base = cols[4]; //the alternate base allele

                                //check for multiple bases
                                if (alt_base.Length > 1)
                                {
                                    skip_counter++;
                                }
                                else
                                {
                                    //set up more counters
                                    double C_count = 0;
                                    List<double> C_dat = new List<double>();
                                    double T_count = 0;
                                    List<double> T_dat = new List<double>();

                                    //for (int i = 10; i < cols.Length; i++) //skip 767 which is col 10
                                    for (int i = 9;i< cols.Length;i++)//for ali.vcf
                                    {
                                        string[] info = cols[i].Split(':'); //split the info column that has AD on col[1]
                                        //if (info.Length == 5 && i < 13)//ali_w_767.vcf
                                        if (info.Length == 5 && i < 12)//ali.vcf
                                        {
                                            string[] AD = info[1].Split(',');
                                            if (Convert.ToInt32(AD[0]) + Convert.ToInt32(AD[1]) > 0)
                                            {
                                                C_count++;
                                                C_dat.Add(Convert.ToInt32(AD[0]));
                                                C_dat.Add(Convert.ToInt32(AD[1]));
                                            }
                                        }
                                        //else if (info.Length == 5 && i > 12)//ali_w_767.vcf
                                        else if (info.Length == 5 && i > 11)//ali.vcf
                                        {
                                            string[] ad = info[1].Split(',');
                                            if (Convert.ToInt32(ad[0]) + Convert.ToInt32(ad[1]) > 0)
                                            {
                                                T_count++;
                                                T_dat.Add(Convert.ToInt32(ad[0]));
                                                T_dat.Add(Convert.ToInt32(ad[1]));
                                            }
                                        }
                                    }

                                    if (C_count > 0 && T_count > 0)
                                    {
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

                                        if (qT_1 >= minimum_reads && qC_1 >= minimum_reads)
                                        {
                                            double qC_hat = qC_0 / qC_1;
                                            double qT_hat = qT_0 / qT_1;
                                            if (qC_hat != qT_hat)
                                            {
                                                num_snps++;
                                                double var_C = 1.0 / qC_1; //transformed variance of C
                                                double var_T = 1.0 / qT_1; //for T
                                                Var_snp_specific += (var_C + var_T);

                                                //calculate divergence
                                                double diverge = (2.0 * ((Math.Asin(Math.Pow(qT_hat, 0.5))) - (Math.Asin(Math.Pow(qC_hat, 0.5)))));
                                                Random rnd = new Random();
                                                if (rnd.Next(1, 3) == 1)
                                                {
                                                    zraw.Add(diverge);
                                                    zranked.Add(diverge);
                                                }
                                                else
                                                {
                                                    zraw.Add(-diverge);
                                                    zranked.Add(-diverge);
                                                }

                                                double t_C_pop = Math.Asin(Math.Sqrt(qC_hat)); //transformed variance for the C population
                                                double t_T_pop = Math.Asin(Math.Sqrt(qT_hat)); //transformed variance for the T population

                                                //we have the data we need, create a snp and add to the output list
                                                SNP acceptedsnp = new SNP();
                                                acceptedsnp.Chromosome = Convert.ToInt16(chromosome);
                                                acceptedsnp.Basepair = position;
                                                acceptedsnp.C_variance = var_C;
                                                acceptedsnp.T_variance = var_T;
                                                acceptedsnp.Transformed_c_variance = t_C_pop;
                                                acceptedsnp.Transformed_t_variance = t_T_pop;
                                                acceptedsnp.Old_identifier = cols[0];
                                                output_snps.Add(acceptedsnp);
                                                //old/new output the same up until here at least
                                            }
                                            else { skip_counter++; }
                                        }
                                        else { skip_counter++; }
                                    }
                                }
                            }
                            else { skip_counter++; }
                        }
                    }
                }
            }

            //now process the list
            Console.WriteLine("SNP processing complete, starting analysis...");
            Console.WriteLine("Skipped " + skip_counter + " SNPs");
            Console.WriteLine("Accepted " + output_snps.Count + " SNPs");
            Console.WriteLine("Sampling/genotyping varience " + (Var_snp_specific / Convert.ToDouble(num_snps)));

            //sort Z list (only the one actually for ranking)
            zranked.Sort();
            var n25 = zranked[num_snps / 4];
            var n50 = zranked[num_snps / 2];
            var n75 = zranked[3 * num_snps / 4];

            //report
            Console.WriteLine("Z percentiles (without direction) " + n25 + " " + n50 + " " + n75);
            Console.WriteLine("Total variance in Z (based on IQR) " + Math.Pow((n75 - n25 / 1.349), 2));
            double Var_neutral = (Math.Pow(((n75 - n25) / 1.349), 2)) - Var_snp_specific / num_snps;
            Console.WriteLine("Bulk sampling and library variance " + Var_neutral);
            for (int i = 0; i < zraw.Count; i++)
            {
                double vdiv = Var_neutral + output_snps[i].C_variance + output_snps[i].T_variance;
                double b = Math.Pow(Convert.ToDouble(zraw[i]), 2) / vdiv;
                output_snps[i].B_standard = b;
                if ((output_snps[i].Old_identifier + "_" + output_snps[i].Basepair) == "sNNffold_1_15730")
                {
                    Console.WriteLine(output_snps[i].Old_identifier + "_" + output_snps[i].Basepair + '\t' + b.ToString());
                }
            }

            //end unless it's a go for doing b*
            if (bs_go == 'Y')
            {
                //set things up from the chisq file
                Console.WriteLine("Assembling chi square apparatus...");
                List<double> df = new List<double>();
                List<double> bs = new List<double>();
                ArrayList[] percentiles = new ArrayList[100];
                for (int i = 0; i < percentiles.Length; i++)
                {
                    percentiles[i] = new ArrayList();
                }
                int line_idx = 0;
                foreach (string line in File.ReadLines(chisq_path))
                {
                    string[] cols = line.Replace("\n", "").Split('\t');
                    df.Add(float.Parse(cols[0]));
                    for (int j = 0; j < 9; j++)
                    {
                        percentiles[line_idx].Add(float.Parse(cols[j + 1]));
                    }
                    double rx = (Convert.ToDouble(percentiles[line_idx][2]) + Convert.ToDouble(percentiles[line_idx][0]) - 2 *
                              Convert.ToDouble(percentiles[line_idx][1])) /
                        (Convert.ToDouble(percentiles[line_idx][2]) - Convert.ToDouble(percentiles[line_idx][0]));
                    bs.Add(rx);
                    line_idx++;
                }

                //now start on using the chi square to the B values
                //set up a list for ranking again
                List<double> b_sorted = new List<double>();
                List<double> braw = new List<double>();
                for (int i = window; i < num_snps; i++)
                {
                    if (i % window == 0 || i % window == window / 2)
                    {
                        double b = 0.0;
                        for (int j = 0; j < window; j++)
                        {
                            b += Math.Pow(output_snps[i - j].B_standard, 2);
                        }
                        braw.Add(b);
                        b_sorted.Add(b);
                    }
                }
                b_sorted.Sort();

                //percentiles
                n25 = b_sorted[braw.Count / 4];
                n50 = b_sorted[braw.Count / 2];
                n75 = b_sorted[3 * braw.Count / 4];

                double bowley_skew = (n75 + n25 - 2 * n50) / (n75 - n25);
                Console.WriteLine("B Bowley skew " + bowley_skew);
                double m2 = -1.0;
                int jstar = 0;
                if (bowley_skew > bs[0])
                {
                    Console.WriteLine("Too much skew.");
                }
                else
                {
                    for (int j = 1; j < bs.Count; j++)
                    {
                        if (bowley_skew > bs[j])
                        {
                            m2 = df[j];
                            jstar = j;
                            break;
                        }
                    }
                }
                double cIQR_1 = Convert.ToDouble(percentiles[jstar][2]);
                double cIQR_2 = Convert.ToDouble(percentiles[jstar][0]);
                double cIQR = cIQR_1 - cIQR_2;
                double sigB = (n75 - n25) * Math.Pow((2 * m2), 0.5) / cIQR;

                double p = 0.0;
                for (int j = 0; j < b_sorted.Count; j++)
                {
                    double bx = braw[j];
                    double b_star = m2 + (bx - window) * (Math.Pow((2 * m2), 0.5)) / sigB;
                    if (b_star < Convert.ToDouble(percentiles[jstar][3]))
                    {
                        p = 0.5;
                    }
                    else if (b_star > Convert.ToDouble(percentiles[jstar][8]))
                    {
                        p = 5.0 * Math.Pow(10, (1.0 - 8.5)); //p < table minimum
                    }
                    else
                    {
                        for (int k = 4; k < 9; k++)
                        {
                            if (b_star > Convert.ToDouble(percentiles[jstar][k - 1]) && b_star <= Convert.ToDouble(percentiles[jstar][k]))
                            {
                                double dx = (b_star - Convert.ToDouble(percentiles[jstar][k - 1])) / (Convert.ToDouble(percentiles[jstar][k]) - Convert.ToDouble(percentiles[jstar][k - 1]));
                                double x = k - 1 + dx;
                                p = 5.0 * Math.Pow(10, (1.0 - x));
                                break;
                            }
                        }
                    }
                    output_snps[j].Raw_p = p;
                    output_snps[j].B_star = b_star;
                }
                for (int i = 0; i < braw.Count(); i++)
                {
                    trimmed_snps.Add(output_snps[i]);
                }
            }

            switch (bs_go)
            {
                case 'Y': return trimmed_snps;
                case 'N': return output_snps;
                default: return output_snps;
            }
        }
    }
}
