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
        public static List<string> b_processing(int window, string vcfpath)
        {
            List<string> bList = new List<string>(); //out2
            //these may be user defined
            int min_reads = 20;
            List<string> LineID = new List<string>(); //may not do anything useful
            List<string> Locations = new List<string>();
            List<string> raw_results = new List<string>(); //out0
            List<string> qCqTLines = new List<string>(); //for diagnostics
            List<string> diagnostic_log = new List<string>(); //for diagnostics
            List<decimal> ranked_z = new List<decimal>();
            List<string> qcqthat_same = new List<string>();


            //set up counters
            int num_snps = 0;
            decimal Var_snp_specific = (decimal)0.0;
            List<decimal> zraw = new List<decimal>();
            List<string> accepted_snps = new List<string>();
            int skipcounter = 0;
            int b0_counter = 0;
            int diverge_0 = 0;
            int num_snp_lines = 0;

            //string diags
            string variances = string.Empty;


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
                                        Locations.Add(cols[0] + "_" + cols[1]);
                                        num_snps++;
                                        //qC_hat and qT_hat should be different almost all of the time
                                        decimal qC_hat = qC_0 / qC_1;
                                        decimal qT_hat = qT_0 / qT_1;

                                        //if (qC_hat == qT_hat) //this catches any time qC_hat and qT_hat are the same - which shouldn't really happen
                                        //{
                                        //    qcqthat_same.Add(cols[0] + "_" + cols[1] + " \nqC_0 / qC_1 = " + qC_0 + " / " + qC_1 + " = " + qC_hat + " = qC_hat" + "\n" +
                                        //        "qT_0 / qT_1 = " + qT_0 + " / " + qT_1 + " = " + qT_hat + " = qT_hat \n");
                                        //    Console.WriteLine(cols[0] + "_" + cols[1] + " \nqC_0 / qC_1 = " + qC_0 + " / " + qC_1 + " = " + qC_hat + " = qC_hat" + "\n" +
                                        //        "qT_0 / qT_1 = " + qT_0 + " / " + qT_1 + " = " + qT_hat + " = qT_hat \n");
                                        //}

                                        decimal var_C = (decimal)1.0 / qC_1;
                                        decimal var_T = (decimal)1.0 / qT_1;
                                        Var_snp_specific += (var_C + var_T);

                                        //converted to float for more precision
                                        decimal diverge = (decimal)(2.0 * ((Math.Asin(Math.Pow((double)qT_hat, 0.5))) - (Math.Asin(Math.Pow((double)qC_hat, 0.5)))));
                                        //double diverge = 2.0 * (Math.Asin(Math.Pow(qT_hat, 0.5)) - (Math.Asin(Math.Pow(qC_hat, 0.5))));
                                        if (diverge == 0)
                                        {
                                            //Console.WriteLine(" Diverge = " + diverge + " " + cols[0] + "_" + cols[1]);
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
                                        accepted_snps.Add(cols[0] + "_" + cols[1] + '\t' + var_C + '\t' + var_T);
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

            variances = "Z percentiles (without direction) " + n25 + " " + n50 + " " + n75 + "\n" +
                "Total variance in Z (based on IQR) " + Math.Pow((Convert.ToDouble(n75) - Convert.ToDouble(n25) / 1.349), 2) + "\n" +
                "Bulk sampling and library variance " + Var_neutral;


            //build the b list
            Console.WriteLine("Building the B list...");
            for (int k = 0; k < zraw.Count; k++)
            {
                string unparsed_line = accepted_snps[k].ToString();
                string[] parsedarray = unparsed_line.Split('\t').ToArray();
                decimal vdiv = (decimal)Var_neutral + Convert.ToDecimal(parsedarray[1]) + Convert.ToDecimal(parsedarray[2]); // 0 is snnffold_X, 1 is var_C and 2 is var_T
                decimal b = (decimal)Math.Pow(Convert.ToDouble(zraw[k]), 2) / vdiv;
                bList.Add(parsedarray[0].ToString() + '\t' + b.ToString() + '\n');
                if (b == 0)
                {
                    b0_counter++;
                }
            }

            //make diag file
            Console.WriteLine("Building the diagnostic file...");
            diagnostic_log.Add("-----------------------------------------");
            diagnostic_log.Add("Task: B Processing at SNP window " + window);
            diagnostic_log.Add("Num of actual SNP's: " + num_snp_lines);
            diagnostic_log.Add("Processed SNP lines: " + num_snps);
            diagnostic_log.Add("Skipped SNP's: " + skipcounter);
            diagnostic_log.Add("Accepted SNP's: " + accepted_snps.Count);
            diagnostic_log.Add("Num of divergence = 0: " + diverge_0);
            diagnostic_log.Add("Num of B = 0: " + b0_counter);
            diagnostic_log.Add(variances);
            diagnostic_log.Add("Num of items in qCqTLines file: " + qCqTLines.Count);
            diagnostic_log.Add("-----------------------------------------");
            using (StreamWriter diagwriter = File.CreateText("diagnostics.txt"))
            {
                foreach (var line in diagnostic_log)
                {
                    diagwriter.WriteLine(line);
                }
            }


            //Report # of b = 0 and diverge of 0
            Console.WriteLine("Number of 0 divergence: " + diverge_0);
            Console.WriteLine("Number of 0 B: " + b0_counter);

            //Console.WriteLine("Writing qCqT_Lines.txt file...");
            //for diagnostics make a qcqtlines file
            using (StreamWriter qcqtwriter = File.CreateText("qCqT_same_hats.txt"))
            {
                foreach (var line in qcqthat_same)
                {
                    qcqtwriter.WriteLine(line);
                }
            }

            Console.WriteLine("Analysis complete...");
            //return divergeissues;
            return bList;
        }
    }
}
