using System;
using System.Data;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;


namespace mg_gap
{
    class VcfParser
    {
        public static List<string> b_processing(int window, string vcfpath, char bs_go)
        {
            string in1 = (@"N:/app dev/scoville research/program files/github repo/mg-gap/mg-gap/mg-gap/mg-gap/support files/chisq.txt");
            List<string> bList = new List<string>(); //out2
            //these may be user defined
            int min_reads = 20;
            List<string> LineID = new List<string>(); //may not do anything useful
            List<string> Locations = new List<string>();
            //List<string> raw_results = new List<string>(); //out0
            //List<string> qCqTLines = new List<string>(); //for diagnostics
            //List<string> diagnostic_log = new List<string>(); //for diagnostics
            List<double> ranked_z = new List<double>();


            //set up counters
            int num_snps = 0; //this should be a total of all the *accepted* SNP's in the file
            double Var_snp_specific = (double)0.0;
            List<double> zraw = new List<double>();
            List<string> accepted_snps = new List<string>();
            int skipcounter = 0;
            //int b0_counter = 0;
            int diverge_0 = 0;
            int num_snp_lines = 0; //this should be a total of all the SNP's in the file

            //string diags
            //string variances = string.Empty;

            //b = 0 is artificial


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
                            if (index > 8)
                            {
                                LineID.Add(cols[index]); //doesn't really do anything
                            }
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
                         * 9 = 0/0:5,0:5:12:0,12,135 (CA set)
                         * 10 = 0/1:8,1:9:16:16,0,185 (CB set)
                         * 11 = ./. (S set)
                         * 12 = 1/1:2,2:4:6:80,6,0 (TA set)
                         * 13 = 0/1:4,1:5:31:31,0,49 (TB set)
                         * keep in mind that not every snp has every set for data - some might have s, some might have ca, etc
                         * */

                        string scaff = cols[0].Split('_').ToString();
                        string chrom = cols[0].ToString();
                        chrom = chrom.Remove(0, chrom.IndexOf("_") + 1);

                        if (Convert.ToInt16(chrom) < 15) //customized to only go up to and including snaffold 14
                        {
                            int position = Convert.ToInt32(cols[1]);
                            string ref_base = cols[3];
                            string alt_base = cols[4];
                            //Console.WriteLine("Ref " + ref_base + " Alt " + alt_base);
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

                                for (int j = 9; j < 9 + 5; j++) //line 76 in python
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

                                        if (info.Length == 5 && j < 12)
                                        {
                                            string[] AD = info[1].Split(',');
                                            if (Convert.ToInt16(AD[0]) + Convert.ToInt16(AD[1]) > 0)
                                            {
                                                C_count++;
                                                C_dat.Add(Convert.ToInt16(AD[0]));
                                                C_dat.Add(Convert.ToInt16(AD[1]));
                                            }
                                        }
                                        if (info.Length == 5 && j >= 12)
                                        {
                                            string[] ad = info[1].Split(',');
                                            if (Convert.ToInt16(ad[0]) + Convert.ToInt16(ad[1]) > 0)
                                            {
                                                T_count++;
                                                T_dat.Add(Convert.ToInt16(ad[0]));
                                                T_dat.Add(Convert.ToInt16(ad[1]));
                                            }
                                        }
                                    }
                                }


                                /* Here's what the C_dat and T_dat will look like (from our ongoing example)
                                C_dat array as follows: [ 5,0,8,1 ]
                                T_dat array as follows: [ 2,2,4,1 ]
                                */
                                //Console.WriteLine("C_dat: [ " + string.Join(",", C_dat.ToArray()) + " ] count " + C_dat.Count() + "\nT_dat: [ " + string.Join(",", T_dat.ToArray()) + " ] count " + T_dat.Count() + "\n");

                                //if (C_count >= 0 && T_count >= 0) //this is the problem loop
                                if (C_count > 0 && T_count > 0) //this is the problem loop edited w/JK this whole thing should be iterated a number of times per snp
                                {

                                    double qC_0 = (double)0.0;
                                    double qC_1 = (double)0.0;
                                    double qT_0 = (double)0.0;
                                    double qT_1 = (double)0.0;

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
                                        //make a qcqtlines thingy for tracking
                                        //qCqTLines.Add("1 \t" + qT_1.ToString() + "\t" + qC_1.ToString() + "\t" + "[" + string.Join(", ", C_dat.ToArray()) + "]\t C_count array size: " + C_count.ToString() + " \t" + "[" + string.Join(", ", T_dat.ToArray()) + "]\t T_count array size: " + T_count.ToString());

                                        //skipping output for yut file

                                        //qC_hat and qT_hat should be different almost all of the time
                                        double qC_hat = qC_0 / qC_1;
                                        double qT_hat = qT_0 / qT_1;

                                        if (qC_hat != qT_hat) //this catches any time qC_hat and qT_hat are the same - which shouldn't really happen and so throw them out
                                        {
                                            //moved down one logic level to prevent num_snp's from being off in final calculations
                                            Locations.Add(cols[0] + "_" + cols[1]);
                                            num_snps++;

                                            double var_C = (double)1.0 / qC_1;
                                            double var_T = (double)1.0 / qT_1;
                                            Var_snp_specific += (var_C + var_T);

                                            //note double provides more precision but some of the Math class methods only take double - so there may be some small loss of precision past n^-16 
                                            double diverge = (double)(2.0 * ((Math.Asin(Math.Pow((double)qT_hat, 0.5))) - (Math.Asin(Math.Pow((double)qC_hat, 0.5)))));
                                            if (diverge == 0)
                                            {
                                                Console.WriteLine(" Diverge = " + diverge + " " + cols[0] + "_" + cols[1]);
                                                diverge_0++;
                                            }

                                            //raw_results.Add(cols[0] + '\t' + cols[1] + qC_1.ToString() + '\t' + qC_hat.ToString()
                                            //    + '\t' + qT_1.ToString() + '\t' + qT_hat.ToString() + '\t' + diverge.ToString() + '\n');

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

                                            //SNP accepted.
                                            accepted_snps.Add(cols[0] + "_" + cols[1] + '\t' + var_C + '\t' + var_T + '\t' + Math.Asin(Math.Sqrt((double)qC_0)) + '\t' + Math.Asin(Math.Sqrt((double)qT_0)));
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
            Console.WriteLine("Accepted " + accepted_snps.Count + " SNPs");

            if (num_snps != accepted_snps.Count)
            {
                Console.WriteLine("number of accepted snp's varies from length of accepted list - respectively " + num_snps + " vs. " + accepted_snps.Count);
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
            List<double> b_std = new List<double>();
            for (int k = 0; k < zraw.Count; k++)
            {
                string unparsed_line = accepted_snps[k].ToString();
                string[] parsedarray = unparsed_line.Split('\t').ToArray();
                double vdiv = Var_neutral + Convert.ToDouble(parsedarray[1]) + Convert.ToDouble(parsedarray[2]); // 0 is snnffold_X, 1 is var_C and 2 is var_T
                double b = Math.Pow(Convert.ToDouble(zraw[k]), 2) / vdiv;
                b_std.Add(b);
                //bList.Add(parsedarray[0].ToString() + '\t' + b.ToString() + '\n'); //traditional way
                bList.Add(parsedarray[0].ToString() + '\t' + b.ToString() + '\t' + vdiv + '\t' + parsedarray[3].ToString() + '\t' + parsedarray[4].ToString()); //verbose output
            }

            //set things up from the chisq file
            Console.WriteLine("Enumerating and working with chisq file.");
            List<double> df = new List<double>();
            List<double> bs = new List<double>();
            ArrayList[] percentiles = new ArrayList[100];
            for (int i = 0; i < percentiles.Length; i++)
            {
                percentiles[i] = new ArrayList();
            }
            var line_idx = 0;
            foreach (string line in File.ReadLines(in1))
            {
                string[] cols = line.Replace(Environment.NewLine, "").Split(new[] { '\t' });
                df.Add(float.Parse(cols[0]));
                for (int j = 1; j < 9; j++)
                {
                    percentiles[line_idx].Add(float.Parse(cols[j+1]));
                }
                double rx = (Convert.ToDouble(percentiles[line_idx][2]) + Convert.ToDouble(percentiles[line_idx][0]) - 2 *
                          Convert.ToDouble(percentiles[line_idx][1])) /
                    (Convert.ToDouble(percentiles[line_idx][2]) - Convert.ToDouble(percentiles[line_idx][0]));
                bs.Add(rx);
                line_idx++;
            }
            //#####

            List<double> braw = new List<double>();
            List<string> bloc = new List<string>();
            List<double> b_sorted = new List<double>();
            for (int k = 0; k < num_snps; k++)
            {
                if (k % window == 0 || k % window == window / 2) //is end of window?
                {
                    double b = (double)0.0;
                    for (int j = 0; j < window; j++)
                    {
                        b += ((double)Math.Pow((double)b_std[k - j], 2));
                    }
                    bloc.Add(Locations[k]);
                    braw.Add(b);
                    b_sorted.Add(b);
                }
            }

            b_sorted.Sort();
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

            List<string> bs_list = new List<string>();
            double p = 0.0;
            for (int j = 0; j < b_sorted.Count; j++)
            {
                double bstar = m2 + (braw[j] - window) * Math.Pow((2 * m2), 0.5) / sigB;
                
                if (bstar < Convert.ToDouble(percentiles[jstar][3]))//p > 0.05
                {
                    p = 0.5;
                }
                else if (bstar > Convert.ToDouble(percentiles[jstar][8])) //b less than table min
                {
                    p = 5.0 * Math.Pow(10, (1.0 - 8.5));
                }
                else
                {
                    for (int k = 4; k < 9; k++)
                    {
                        if (bstar > Convert.ToDouble(percentiles[jstar][k - 1]) && bstar <= Convert.ToDouble(percentiles[jstar][k]))
                        {
                            double dx = (bstar - Convert.ToDouble(percentiles[jstar][k - 1])) / (Convert.ToDouble(percentiles[jstar][k]) - Convert.ToDouble(percentiles[jstar][k - 1]));
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
                    bs_list.Add(midpoint + '\t' + braw[j] + '\t' + bstar.ToString() + '\t' + p.ToString() + '\t');
                }
            }

            //build the b list
            Console.WriteLine("Building the B list...");


            //Report # of b = 0 and diverge of 0
            Console.WriteLine("Number of 0 divergence: " + diverge_0); //this should always be 0
            Console.WriteLine("Analysis complete...");

            //if this is the version of the run where we want to get B*...
            if (bs_go == 'Y')
            {
                return bs_list;
            }
            else
            {
                return bList;
            }   
        }
    }
}
