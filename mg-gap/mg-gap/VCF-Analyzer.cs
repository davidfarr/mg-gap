using System;
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

        public static List<SNP> snp_list(int window, string vcfpath, char bs_go)
        {
            //create the empty list
            List<SNP> snp_list_raw = new List<SNP>(); //stores all SNPs and all their info in the analysis

            string in1 = (@"N:/app dev/scoville research/program files/github repo/mg-gap/mg-gap/mg-gap/mg-gap/support files/chisq.txt");
            in1 = @"C:/Users/David/Documents/GitHub/mg-gap/mg-gap/mg-gap/support files/chisq.txt"; //for mac environment only
            //these may be user defined
            int min_reads = 20;
            List<string> LineID = new List<string>(); //may not do anything useful
            List<string> Locations = new List<string>();
            List<double> ranked_z = new List<double>();

            //temporary
            List<string> bs_list = new List<string>();

            //set up counters
            int num_snps = 0; //this should be a total of all the *accepted* SNP's in the file
            double Var_snp_specific = (double)0.0;
            List<double> zraw = new List<double>();
            int skipcounter = 0;
            int diverge_0 = 0;
            int num_snp_lines = 0; //this should be a total of all the SNP's in the file

            //evaluate contents of each line of input file
            using (var fileStream = File.OpenRead(vcfpath))
            using (var streamReader = new StreamReader(fileStream))
            {
                String line;
                while ((line = streamReader.ReadLine()) != null)
                {
                    num_snp_lines++; //keep track of all lines that aren't contigs/headers
                    string[] cols = line.Replace("\n", "").Split('\t');

                    if (cols.Length < 2)
                    {
                        continue;
                    }
                    else if (cols[0] == "#CHROM")
                    {
                        foreach (string entity in cols)
                        {
                            int index = Array.IndexOf(cols, entity);
                            Console.WriteLine("Found " + entity);
                        }
                    }
                    else if (!cols[0].Contains("##")) //gets rid of contigs
                    {
                        /* guide to the cols - keep in mind the only things we are using are cols0,1,3,4 and then the AD for each data set that's not blank (S is)
                         * 0 = sNNfold_1
                         * 1 = base pair
                         * 2 = literally just a period .
                         * 3 = ref base
                         * 4 = alt base
                         * 5 = qual number
                         * 6 = also literally just a period .
                         * 7 = Example: AC=4;AF=0.500;AN=8;BaseQRankSum=-0.045;DP=23;Dels=0.00;FS=2.632;HaplotypeScore=0.2493;MLEAC=4;MLEAF=0.500;MQ=33.07;MQ0=2;MQRankSum=3.000;QD=5.11;ReadPosRankSum=-0.313
                         * 8 = GT:AD:DP:GQ:PL
                         * 9 = 767 as below
                         * 10 = 0/0:5,0:5:12:0,12,135 (CA set)
                         * 11 = 0/1:8,1:9:16:16,0,185 (CB set)
                         * 12 = ./. (S set)
                         * 13 = 1/1:2,2:4:6:80,6,0 (TA set)
                         * 14 = 0/1:4,1:5:31:31,0,49 (TB set)
                         * keep in mind that not every snp has every set for data - some might have s, some might have ca, etc
                         * */

                        string scaff = cols[0].Split('_').ToString();
                        string chrom = cols[0].ToString();
                        chrom = chrom.Remove(0, chrom.IndexOf("_") + 1);

                        if (Convert.ToInt32(chrom) < 15) //customized to only go up to and including snaffold 14
                        {
                            int position = Convert.ToInt32(cols[1]);
                            string ref_base = cols[3];
                            string alt_base = cols[4];
                            //check if multiple bases
                            if (alt_base.Length > 1)
                            {
                                skipcounter++;
                            }
                            else
                            {
                                //set up counters - these will be critical to analysis
                                double C_count = 0;
                                //ArrayList C_dat = new ArrayList();
                                List<double> C_dat = new List<double> { };
                                double T_count = 0;
                                //ArrayList T_dat = new ArrayList();
                                List<double> T_dat = new List<double> { };

                                //for (int j = 9; j < 9 + 5; j++) //line 76 in python
                                for (int j = 10; j < 9 + 6; j++) //line 76 in python with 767

                                {
                                    if (cols.Length < (9 + 5))
                                    {
                                        Console.WriteLine("whoa " + LineID + " " + cols.Length);
                                    }
                                    else
                                    {
                                        string[] info = cols[j].Split(':');
                                        //output info about info
                                        /* Guide to info array total 5 columns
                                         * 0 = 0/0 (GT)
                                         * 1 = 5,0 (this is AD)
                                         * 2 = 5 (DP)
                                         * 3 = 12 (GQ)
                                         * 4 = 0,12,135 (PL) */
                                        //for (int i = 0; i < info.Length; i++)
                                        //{
                                        //    Console.WriteLine("Row info column " + i + ": " + info[i]);
                                        //}

                                        if (info.Length == 5 && j < 13) //12 for non-767 set

                                        {
                                            string[] AD = info[1].Split(',');
                                            if (Convert.ToInt32(AD[0]) + Convert.ToInt32(AD[1]) > 0)
                                            {
                                                C_count++;
                                                C_dat.Add(Convert.ToInt32(AD[0]));
                                                C_dat.Add(Convert.ToInt32(AD[1]));
                                            }
                                        }
                                        if (info.Length == 5 && j >= 13) //12 for non-767
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
                                }


                                /* Here's what the C_dat and T_dat will look like (from our ongoing example)
                                C_dat array as follows: [ 5,0,8,1 ]
                                T_dat array as follows: [ 2,2,4,1 ]
                                */

                                if (C_count > 0 && T_count > 0) //this is the problem loop edited w/JK this whole thing should be iterated a number of times per snp
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


                                    if (qT_1 >= min_reads && qC_1 >= min_reads)
                                    {
                                        //qC_hat and qT_hat should be different almost all of the time
                                        double qC_hat = qC_0 / qC_1;
                                        double qT_hat = qT_0 / qT_1;

                                        if (qC_hat != qT_hat) //this catches any time qC_hat and qT_hat are the same - which shouldn't really happen and so throw them out
                                        {
                                            //moved down one logic level to prevent num_snp's from being off in final calculations
                                            Locations.Add(cols[0] + "_" + cols[1]);
                                            num_snps++;

                                            double var_C = 1.0 / qC_1;
                                            double var_T = 1.0 / qT_1;
                                            Var_snp_specific += (var_C + var_T);

                                            //note double provides more precision but some of the Math class methods only take double - so there may be some small loss of precision past n^-16 
                                            double diverge = (2.0 * ((Math.Asin(Math.Pow(qT_hat, 0.5))) - (Math.Asin(Math.Pow(qC_hat, 0.5)))));
                                            if (diverge == 0)
                                            {
                                                Console.WriteLine(" Diverge = " + diverge + " " + cols[0] + "_" + cols[1]);
                                                diverge_0++;
                                            }

                                            Random rnd = new Random();
                                            if (rnd.Next(1, 3) == 1)
                                            {
                                                zraw.Add(diverge);
                                                ranked_z.Add(diverge);
                                            }
                                            else
                                            {
                                                zraw.Add(-diverge);
                                                ranked_z.Add(-diverge);
                                            }

                                            double t_C_pop = Math.Asin(Math.Sqrt(qC_hat));
                                            double t_T_pop = Math.Asin(Math.Sqrt(qT_hat));
                                            //SNP accepted.
                                            //create the snp and add to list
                                            SNP snp_inprogress = new SNP();
                                            snp_inprogress.Chromosome = Convert.ToInt16(chrom);
                                            snp_inprogress.Basepair = position;
                                            snp_inprogress.C_variance = var_C;
                                            snp_inprogress.T_variance = var_T;
                                            snp_inprogress.Transformed_c_variance = t_C_pop;
                                            snp_inprogress.Transformed_t_variance = t_T_pop;
                                            snp_inprogress.Old_identifier = cols[0];
                                            snp_list_raw.Add(snp_inprogress);

                                        }
                                        else
                                        {
                                            skipcounter++; //there should be no situation where the variance is 0, this undermines the point of a VCF file
                                        }

                                    }
                                    else
                                    {
                                        skipcounter++;
                                    }
                                    //Console.WriteLine("List qC has " + qC.Count + " items and qT has " + qT.Count + " items.");
                                }
                            }
                        }
                    }
                }
            }
            //wrap up
            Console.WriteLine("SNP processing complete, starting analysis...");
            Console.WriteLine("Skipped " + skipcounter + " SNPs");
            Console.WriteLine("Accepted " + snp_list_raw.Count + " SNPs");

            if (num_snps != snp_list_raw.Count)
            {
                Console.WriteLine("number of accepted snp's varies from length of accepted list - respectively " + num_snps + " vs. " + snp_list_raw.Count);
            }

            Console.WriteLine("Sampling/genotyping varience " + (Var_snp_specific / Convert.ToDouble(num_snps)));

            //make the z, rankings
            //zraw UNSORTED is used for b so you have to do the ranked_z in another array
            //confirm through console the length is the same
            ranked_z.Sort();

            var n25 = ranked_z[num_snps / 4];
            var n50 = ranked_z[num_snps / 2];
            var n75 = ranked_z[3 * num_snps / 4];

         


            Console.WriteLine("Z percentiles (without direction) " + n25 + " " + n50 + " " + n75);
            Console.WriteLine("Total variance in Z (based on IQR) " + Math.Pow((n75 - n25 / 1.349), 2));

            double Var_neutral = (Math.Pow(((n75 - n25) / 1.349), 2)) - Var_snp_specific / num_snps;
            Console.WriteLine("Bulk sampling and library variance " + Var_neutral);

            //make an unsorted B list - the output of this is the last step if just doing B
            for (int k = 0; k < snp_list_raw.Count; k++)
            {
                double vdiv = Var_neutral + snp_list_raw[k].C_variance + snp_list_raw[k].T_variance; // 0 is snnffold_X, 1 is var_C and 2 is var_T
                //double b = Math.Pow(Convert.ToDouble(zraw[k]), 2) / vdiv;
                double b = Math.Pow(Convert.ToDouble(snp_list_raw[k].Divergence), 2) / vdiv;

                snp_list_raw[k].B_standard = b;
            }

            //#####
            //all below goes to b*
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
                foreach (string line in File.ReadLines(in1)) //this is the chisq file
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
                    //Console.WriteLine("chisq enumeration: added DF " + cols[0] + " and RX " + rx); //adds DF from 0.1 to 10 with different RX per df
                }
                List<double> braw = new List<double>();
                List<string> bloc = new List<string>();
                List<double> b_sorted = new List<double>();
                for (int k = window; k < num_snps; k++)
                {
                    if (k % window == 0 || k % window == window / 2) //is end of window?
                    {
                        double b = 0.0;
                        for (int j = 0; j < window; j++)
                        {
                            b += Math.Pow(snp_list_raw[k - j].B_standard, 2);
                        }
                        bloc.Add(Locations[k]);
                        braw.Add(b);
                        b_sorted.Add(b);
                    }
                }

                b_sorted.Sort();

                //redo percentiles?
                n25 = b_sorted[braw.Count / 4];
                n50 = b_sorted[braw.Count / 2];
                n75 = b_sorted[3 * braw.Count / 4];

                double b_skew = (n75 + n25 - 2 * n50) / (n75 - n25);
                Console.WriteLine("B Bowley skew " + b_skew);
                double m2 = -1.0;
                int jstar = 0;
                if (b_skew > bs[0])
                {
                    Console.WriteLine("Too much skew.");
                }
                else
                {
                    for (int j = 1; j < bs.Count; j++)
                    {
                        if (b_skew > bs[j])
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
                    string midpoint = string.Empty;
                    if (j > 0)
                    {
                        midpoint = bloc[j - 1].ToString();
                    }
                    else
                    {
                        midpoint = bloc[j].ToString();
                    }

                    SNP searchsnp = (from x in snp_list_raw where (x.Old_identifier + "_" + x.Basepair) == midpoint select x).First();
                    searchsnp.B_star = b_star;
                    searchsnp.B_standard = bx;
                    searchsnp.Raw_p = p;

                    //bs_list.Add(midpoint + '\t' + bx.ToString() + '\t' + b_star.ToString() + '\t' + p.ToString());
                }
            }



            //build the b list
            Console.WriteLine("Building the B list...");
            Console.WriteLine("B* available for " + bs_list.Count() + " of " + snp_list_raw.Count + " SNPs.");
            


            Console.WriteLine("Analysis complete...");

            //if this is the version of the run where we want to get B*...
            return snp_list_raw;
        }
    }
}
