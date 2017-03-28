using System;
using System.Data;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Text.RegularExpressions;
using System.Diagnostics;
using System.Linq;
using RDotNet;
using RDotNet.NativeLibrary;

namespace mg_gap
{
    class VcfParser
    {
        public static ArrayList b_processing(int window, string vcfpath)
        {
            ArrayList bList = new ArrayList(); //out2
            //these may be user defined
            int min_reads = 20;
            ArrayList LineID = new ArrayList();
            ArrayList Locations = new ArrayList();
            ArrayList raw_results = new ArrayList(); //out0
            ArrayList divergeissues = new ArrayList();

            //set up counters
            int num_snps = 0;
            double Var_snp_specific = 0.0;
            ArrayList zraw = new ArrayList();
            ArrayList accepted_snps = new ArrayList();
            int skipcounter = 0;
            int b0_counter = 0;
            int diverge0_counter = 0;


            //evaluate contents of each line of input file
            using (var fileStream = File.OpenRead(vcfpath))
            using (var streamReader = new StreamReader(fileStream))
            {
                String line;
                while ((line = streamReader.ReadLine()) != null)
                {
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
                                LineID.Add(cols[index]);
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
                                float C_count = 0;
                                //ArrayList C_dat = new ArrayList();
                                List<double> C_dat = new List<double> { };
                                float T_count = 0;
                                //ArrayList T_dat = new ArrayList();
                                List<double> T_dat = new List<double> { };

                                for (int j = 9; j < 9 + 5; j++)
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
                                //status check
                                //Console.WriteLine("For " + cols[0] + " bp# " + cols[1] + " the C_dat is [" + string.Join(",", C_dat.ToArray()) + @"] and the T_dat is [" + string.Join(",", T_dat.ToArray()) + @"]");

                                /* Here's what the C_dat and T_dat will look like (from our ongoing example)
                                C_dat array as follows: [ 5,0,8,1 ]
                                T_dat array as follows: [ 2,2,4,1 ]
                                */

                                if (C_count >= 0 && T_count >= 0) //this is the problem loop
                                {
                                    double qC_0 = 0.0;
                                    double qC_1 = 0.0;
                                    double qT_0 = 0.0;
                                    double qT_1 = 0.0;

                                    for (int i = 0; i < C_dat.Count/2; i++)
                                    {
                                        double m = C_dat[2 * i] + C_dat[2 * i + 1];
                                        qC_0 += C_dat[2 * i];
                                        qC_1 += m;
                                    }

                                    for (int i = 0; i < T_dat.Count / 2; i++)
                                    {
                                        double m = T_dat[2 * i] + T_dat[2 * i + 1];
                                        qT_0 += T_dat[2 * i];
                                        qT_1 += m;
                                    }


                                    if (qT_1 >= min_reads && qC_1 >= min_reads)
                                    {
                                        //skipping output for yut file
                                        Locations.Add(cols[0] + "_" + cols[1]);
                                        //num_snps++;
                                        float qC_hat = (float)qC_0 / (float)qC_1;
                                        float qT_hat = (float)qT_0 / (float)qT_1;
                                        double var_C = 1.0 / qC_1;
                                        double var_T = 1.0 / qT_1;
                                        Var_snp_specific += (var_C + var_T);

                                        float diverge = (float)(2.0 * (Math.Asin(Math.Pow(qT_hat, 0.5)) - Math.Asin(Math.Pow(qC_hat, 0.5))));
                                        //double diverge = 2.0 * (Math.Asin(Math.Pow(qT_hat, 0.5)) - (Math.Asin(Math.Pow(qC_hat, 0.5))));
                                        //if (diverge == 0)
                                        //{
                                        //    divergeissues.Add(cols[0] + "_" + cols[1] + "\t" + "qC_hat = qC[0]/qC[1] = " + qC_0 + " / " + qC_1 +
                                        //        " = " + qC_hat + "\t" + "qt_hat = " + qT_0 + " / " + qT_1 +
                                        //        " = " + qT_hat + "\t" + " var_C = 1.0/qC[1] = " + var_C + "\t" + " var_T = " + var_T + "\t" + "Diverge = 0 \n");
                                        //}

                                        //raw_results.Add(cols[0] + '\t' + cols[1] + qC_1.ToString() + '\t' + qC_hat.ToString()
                                        //    + '\t' + qT_1.ToString() + '\t' + qT_hat.ToString() + '\t' + diverge.ToString() + '\n');

                                        Random rnd = new Random();
                                        if (rnd.Next(1, 3) == 1 && diverge > 0) //also put the diverge filter here
                                        {
                                            zraw.Add(diverge);
                                        }
                                        else if (diverge > 0)
                                        {
                                            zraw.Add(-diverge);
                                        }

                                        //this might not be what we want - throw out diverge = 0
                                        if (diverge > 0)
                                        {
                                            //SNP accepted.
                                            accepted_snps.Add(cols[0] + "_" + cols[1] + '\t' + var_C + '\t' + var_T);
                                            num_snps++;
                                        }
                                        else
                                        {
                                            skipcounter++;
                                            diverge0_counter++;
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


                //wrap up
                Console.WriteLine("Skipped " + skipcounter + " SNPs");
                Console.WriteLine("Accepted " + accepted_snps.Count + " SNPs");
                Console.WriteLine("Sampling/genotyping varience " + (Var_snp_specific / Convert.ToDouble(num_snps)));

                //make the z, rankings
                zraw.Sort();
                var n25 = zraw[num_snps / 4];
                var n50 = zraw[num_snps / 2];
                var n75 = zraw[3 * num_snps / 4];
                Console.WriteLine("Z percentiles (without direction) " + n25 + " " + n50 + " " + n75);
                Console.WriteLine("Total variance in Z (based on IQR) " + Math.Pow((Convert.ToDouble(n75) - Convert.ToDouble(n25) / 1.349), 2));

                double Var_neutral = (Math.Pow(((Convert.ToDouble(n75) - Convert.ToDouble(n25)) / 1.349), 2)) - Var_snp_specific / num_snps;
                Console.WriteLine("Bulk sampling and library variance " + Var_neutral);

                //build the b list
                for (int k = 0; k < zraw.Count; k++)
                {
                    string unparsed_line = accepted_snps[k].ToString();
                    string[] parsedarray = unparsed_line.Split('\t').ToArray();
                    var vdiv = Var_neutral + Convert.ToDouble(parsedarray[1]) + Convert.ToDouble(parsedarray[2]);
                    double b = Math.Pow(Convert.ToDouble(zraw[k]), 2) / vdiv;
                    bList.Add(parsedarray[0].ToString() + '\t' + b.ToString() + '\n');
                    if(b == 0)
                    {
                        b0_counter++;
                    }
                }
            }

            //Report # of b = 0
            Console.WriteLine("Number of B = 0: " + b0_counter);
            //Report # of SNP's that had a diverge = 0
            Console.WriteLine("Number of skipped SNPs with divergence of 0: " + diverge0_counter);

            //return divergeissues;
            return bList;
        }
    }
}
