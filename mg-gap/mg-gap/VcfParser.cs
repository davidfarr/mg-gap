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
            string in1 = @"N:/app dev/scoville research/program files/github repo/mg-gap/mg-gap/mg-gap/mg-gap/support files/chisq.txt");
            List<string> bList = new List<string>(); //out2
            //these may be user defined
            int min_reads = 20;
            List<string> LineID = new List<string>(); //may not do anything useful
            List<string> Locations = new List<string>();
            //List<string> raw_results = new List<string>(); //out0
            //List<string> qCqTLines = new List<string>(); //for diagnostics
            //List<string> diagnostic_log = new List<string>(); //for diagnostics
            List<decimal> ranked_z = new List<decimal>();


            //set up counters
            int num_snps = 0; //this should be a total of all the *accepted* SNP's in the file
            decimal Var_snp_specific = (decimal)0.0;
            List<decimal> zraw = new List<decimal>();
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
                                decimal C_count = 0;
                                //ArrayList C_dat = new ArrayList();
                                List<decimal> C_dat = new List<decimal> { };
                                decimal T_count = 0;
                                //ArrayList T_dat = new ArrayList();
                                List<decimal> T_dat = new List<decimal> { };

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

                                    decimal qC_0 = (decimal)0.0;
                                    decimal qC_1 = (decimal)0.0;
                                    decimal qT_0 = (decimal)0.0;
                                    decimal qT_1 = (decimal)0.0;

                                    for (int i = 0; i < C_dat.Count / 2; i++)
                                    {
                                        decimal m = C_dat[2 * i] + C_dat[2 * i + 1];
                                        qC_0 += C_dat[2 * i]; //ref base
                                        qC_1 += m; //ref + alt base
                                    }

                                    for (int i = 0; i < T_dat.Count / 2; i++)
                                    {
                                        decimal m = T_dat[2 * i] + T_dat[2 * i + 1];
                                        qT_0 += T_dat[2 * i];
                                        qT_1 += m;
                                    }


                                    if (qT_1 >= min_reads && qC_1 >= min_reads)
                                    {
                                        //make a qcqtlines thingy for tracking
                                        //qCqTLines.Add("1 \t" + qT_1.ToString() + "\t" + qC_1.ToString() + "\t" + "[" + string.Join(", ", C_dat.ToArray()) + "]\t C_count array size: " + C_count.ToString() + " \t" + "[" + string.Join(", ", T_dat.ToArray()) + "]\t T_count array size: " + T_count.ToString());

                                        //skipping output for yut file

                                        //qC_hat and qT_hat should be different almost all of the time
                                        decimal qC_hat = qC_0 / qC_1;
                                        decimal qT_hat = qT_0 / qT_1;

                                        if (qC_hat != qT_hat) //this catches any time qC_hat and qT_hat are the same - which shouldn't really happen and so throw them out
                                        {
                                            //moved down one logic level to prevent num_snp's from being off in final calculations
                                            Locations.Add(cols[0] + "_" + cols[1]);
                                            num_snps++;

                                            decimal var_C = (decimal)1.0 / qC_1;
                                            decimal var_T = (decimal)1.0 / qT_1;
                                            Var_snp_specific += (var_C + var_T);

                                            //note decimal provides more precision but some of the Math class methods only take double - so there may be some small loss of precision past n^-16 
                                            decimal diverge = (decimal)(2.0 * ((Math.Asin(Math.Pow((double)qT_hat, 0.5))) - (Math.Asin(Math.Pow((double)qC_hat, 0.5)))));
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

            Console.WriteLine("Sampling/genotyping varience " + (Var_snp_specific / Convert.ToDecimal(num_snps)));

            //make the z, rankings
            //zraw UNSORTED is used for b so you have to do the ranked_z in another array
            //confirm through console the length is the same
            ranked_z.Sort();

            var n25 = ranked_z[num_snps / 4];
            var n50 = ranked_z[num_snps / 2];
            var n75 = ranked_z[3 * num_snps / 4];
            Console.WriteLine("Z percentiles (without direction) " + n25 + " " + n50 + " " + n75);
            Console.WriteLine("Total variance in Z (based on IQR) " + Math.Pow((double)(Convert.ToDecimal(n75) - Convert.ToDecimal(n25) / (decimal)1.349), 2));

            double Var_neutral = (Math.Pow(((double)(Convert.ToDecimal(n75) - Convert.ToDecimal(n25)) / 1.349), 2)) - (double)Var_snp_specific / num_snps;
            Console.WriteLine("Bulk sampling and library variance " + Var_neutral);

            /*variances = "Z percentiles (without direction) " + n25 + " " + n50 + " " + n75 + "\n" +
                "Total variance in Z (based on IQR) " + Math.Pow((Convert.ToDouble(n75) - Convert.ToDouble(n25) / 1.349), 2) + "\n" +
                "Bulk sampling and library variance " + Var_neutral;*/


            //build the b list
            List<decimal> ranked_b = new List<decimal>();
            List<decimal> bRaw = new List<decimal>();
            Console.WriteLine("Building the B list...");
            for (int k = 0; k < zraw.Count; k++)
            {
                string unparsed_line = accepted_snps[k].ToString();
                string[] parsedarray = unparsed_line.Split('\t').ToArray();
                decimal vdiv = (decimal)Var_neutral + Convert.ToDecimal(parsedarray[1]) + Convert.ToDecimal(parsedarray[2]); // 0 is snnffold_X, 1 is var_C and 2 is var_T
                decimal b = (decimal)Math.Pow(Convert.ToDouble(zraw[k]), 2) / vdiv;
                ranked_b.Add(b);
                bRaw.Add(b);
                //bList.Add(parsedarray[0].ToString() + '\t' + b.ToString() + '\n'); //traditional way
                bList.Add(parsedarray[0].ToString() + '\t' + b.ToString() + '\t' + vdiv + '\t' + parsedarray[3].ToString() + '\t' + parsedarray[4].ToString()); //verbose output

                if (b == 0)
                {
                    Console.WriteLine("B=0");
                }
            }

            //if this is the version of the run where we want to get B*...
            if (bs_go == 'Y')
            {
                //make a sorted B list
                ranked_b.Sort();
                //Get percentiles array
                //set things up from the chisq file
                Console.WriteLine("Enumerating and working with chisq file.");
                List<double> df = new List<double>();
                List<decimal> bs = new List<decimal>();
                ArrayList[] percentiles = new ArrayList[100];
                for (int i = 0; i < percentiles.Length; i++)
                {
                    percentiles[i] = new ArrayList();
                }
                var line_idx = 0;
                foreach (var line in File.ReadLines(in1))
                {
                    var cols = line.Replace(Environment.NewLine, "").Split(new[] { '\t' });
                    df.Add(float.Parse(cols[0]));
                    for (int j = 1; j < 9; j++)
                    {
                        percentiles[line_idx].Add(float.Parse(cols[j]));
                    }
                    var rx = (Convert.ToDouble(percentiles[line_idx][2]) + Convert.ToDouble(percentiles[line_idx][0]) - 2 *
                              Convert.ToDouble(percentiles[line_idx][1])) /
                        (Convert.ToDouble(percentiles[line_idx][2]) - Convert.ToDouble(percentiles[line_idx][0]));
                    bs.Add((decimal)rx);
                    line_idx++;
                }
                decimal b_skew = (n75 + n25 - 2 * n50) / (n75 - n25);
                Console.WriteLine("B Bowley skew " + b_skew);
                double m = -1.0;
                int jstar = 0;
                if (b_skew > bs[0])
                {
                    Console.WriteLine("Too much skew.");
                }
                else
                {
                    for (int j = 0; j < bs.Count; j++)
                    {
                        if (b_skew > bs[j])
                        {
                            m = df[j];
                            jstar = j;
                            break;
                        }
                    }
                }
                decimal cIQR = (decimal)percentiles[jstar][2] - (decimal)percentiles[jstar][0];
                decimal sigB = (n75 - n25) * (decimal)Math.Pow((2 * m), 0.5) / cIQR;

                List<string> bs_list = new List<string>();
                double p = 0.0;
                for (int j = 0;j < ranked_b.Count;j++)
                {
                    decimal bstar = (decimal)m + (bRaw[j] - window) * (decimal)Math.Pow((2 * m), 0.5) / sigB;
                    if (bstar < (decimal)percentiles[jstar][3])//p > 0.05
                    {
                        p = 0.5;
                    }
                    else if (bstar > (decimal)percentiles[jstar][8]) //b less than table min
                    {
                        p = 5.0 * Math.Pow(10, (1.0 - 8.5));
                    }
                    else
                    {
                        for (int k = 4;k < 9; k++)
                        {
                            if (bstar > (decimal)percentiles[jstar][k-1] && bstar <= (decimal)percentiles[jstar][k])
                            {
                                decimal dx = (bstar - (decimal)percentiles[jstar][k - 1]) / ((decimal)percentiles[jstar][k] - (decimal)percentiles[jstar][k - 1]);
                                decimal x = k - 1 + dx;
                                p = 5.0 * Math.Pow(10, (1.0 - (double)x));
                                break;
                            }
                        }
                    }
                }
                //add to bs list to output
                //bstar is the value to put in
                foreach (string line in bList)
                {
                    string[] parsedline = line.Split('\t');
                    if ()
                }

                return bs_list;

            }
            else if (bs_go == 'N')
            {
                //Report # of b = 0 and diverge of 0
                Console.WriteLine("Number of 0 divergence: " + diverge_0); //this should always be 0
                Console.WriteLine("Analysis complete...");
                return bList;
            }

           
            //return bList;
        }
    }
}
