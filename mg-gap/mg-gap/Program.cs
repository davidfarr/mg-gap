﻿using System;
using System.Collections;
using System.Data;
using System.Diagnostics;
using System.IO;
using System.Text.RegularExpressions;


namespace mg_gap
{
    class Program
    {
        static void Main(string[] args)
        {
            //set up the filepath - in this version it's hard-coded
            string vcf_path = "N:/app dev/scoville research/program files/dev migration for windows/vcf files/ali.vcf";
            //string vcf_path = "N:/app dev/scoville research/program files/dev migration for windows/vcf files/chrom4.vcf";
            ArrayList passedLines = new ArrayList();
            ArrayList qCqTLines = new ArrayList();

            //more set-up
            //in1 is the chisq file
            string in1 = "N:/app dev/scoville research/program files/dev migration for windows/python/chisq.txt";
            //choose the first to compare, in this case it's S.rmdup.bam
            int initialWindow = 1;
            var extraSpecial_B = (6.0 * initialWindow);
            var Min_reads = 20;
            ArrayList LineID = new ArrayList();
            ArrayList Locations = new ArrayList();
            var Var_neutral = 0.0844026674093 - 0.0515598702709;
            //focus is not to make a file out of the results of the vcf parse but do make an arraylist that has the info
            ArrayList resultsList = new ArrayList(); //out0 in old script
            ArrayList vzList = new ArrayList(); //out1 in old script
            ArrayList zList = new ArrayList(); //out2 in old script
            ArrayList bList = new ArrayList(); //out3 "
            ArrayList zSpecialList = new ArrayList(); //out4 "

            //status update
            Console.WriteLine("Filepath set.");

            //this will be the main datatable
            DataTable vcfTable = new DataTable();
            ArrayList bams = new ArrayList();

            //"counters"
            Console.WriteLine("Preparing counters...");
            int depth_cc = 0;
            int raw_read_count = 0;
            int num_snps = 0;
            double Var_snp_specific = 0.0;
            ArrayList zraw = new ArrayList();
            ArrayList z_std = new ArrayList();
            ArrayList q_C = new ArrayList();
            ArrayList q_T = new ArrayList();

            var C_count = 0;
            ArrayList C_dat = new ArrayList();
            var T_count = 0;
            ArrayList T_dat = new ArrayList();

            //start stopwatch on operation
            Stopwatch stopwatch = new Stopwatch();
            Console.WriteLine("Table build started at {0}", DateTime.Now);
            stopwatch.Start();


            using (var fileStream = File.OpenRead(vcf_path))
            using (var streamReader = new StreamReader(fileStream))
            {
                String line;
                while ((line = streamReader.ReadLine()) != null)
                {
                    if (!line.ToString().Contains("##") && line != null)
                    {
                        if (line.ToString().Contains("#CHROM"))
                        {
                            //this is a header row
                            string[] headerrows = line.Split('\t');
                            foreach (string header in headerrows)
                            {
                                vcfTable.Columns.Add(header.ToString());
                                Console.WriteLine("New header column: " + header.ToString());
                                if (header.ToString().Contains(".bam"))
                                {
                                    bams.Add(header.ToString());
                                }
                            }
                        }
                        else
                        {
                            try
                            {
                                if (line.ToString().Contains("sNNffold") && line != null)
                                {
                                    //Console.WriteLine("sNNffold line is present");
                                    line.Replace("\n", "");
                                    string[] rowdata = line.Split('\t');
                                    //we don't want things than have multiple ALT bases or that don't have REF/ALT pairs
                                    string scaffstring = rowdata[0];
                                    scaffstring = Regex.Replace(scaffstring, "[^$0-9.]", "");
                                    //Console.WriteLine(scaffstring);
                                    int scaffnum = int.Parse(scaffstring);
                                    if (scaffnum < 5 && scaffnum > 3) //modified to only allow scaff 4+
                                    {
                                        //Console.WriteLine("scaffnum < 15 (" + scaffnum + ")"); 
                                        //in the right scaffold range
                                        if (!string.IsNullOrEmpty(rowdata[3].ToString()) || !string.IsNullOrEmpty(rowdata[4].ToString()))
                                        {
                                            //Console.WriteLine("meets having ref and alt base not empty");
                                            //REF and ALT have something in them, good
                                            if (rowdata[4] != null && rowdata[4].ToString().Length <= 1) //checking to make sure only one alt base
                                            {
                                                //cols = rowdata
                                                //var range1 = Enumerable.Range(0,rowdata.Length);
                                                //foreach (int i in range1)
                                                //{
                                                //Console.WriteLine(i + " " + rowdata[i].ToString());
                                                //if (i>8)
                                                //{
                                                //LineID.Add(rowdata[i]);
                                                //}
                                                //else
                                                //{
                                                //reset counters
                                                C_count = 0;
                                                C_dat.Clear();
                                                T_count = 0;
                                                T_dat.Clear();

                                                //looking at CA.rmdup.bam column forward
                                                //var jrange = Enumerable.Range(9,13); //this makes 14 but there's 13 cols
                                                for (int j = 9; j < rowdata.Length; j++)
                                                //foreach (int j in jrange)
                                                //for (int j = 9; j <= 9+5; j++)
                                                {
                                                    //Console.WriteLine("Starting on index " + j.ToString());

                                                    if (rowdata.Length < 9 + 5)
                                                    {
                                                        Console.WriteLine("whoa ");
                                                    }
                                                    else
                                                    {
                                                        var info = rowdata[j].Split(':');
                                                        //Console.WriteLine("Info is  " + info.Length);
                                                        //This will be the control populations CA/CB/S (blank)
                                                        if (info.Length == 5 && j < 12)//data present for C
                                                        {
                                                            var AD = info[1].Split(',');
                                                            if (Convert.ToInt16(AD[0]) + Convert.ToInt16(AD[1]) > 0)
                                                            {
                                                                C_count++;
                                                                C_dat.Add(Convert.ToInt16(AD[0]));
                                                                //Console.WriteLine("Added " + AD[0].ToString() + "to C_dat");
                                                                C_dat.Add(Convert.ToInt16(AD[1]));
                                                                //Console.WriteLine("Added " + AD[1].ToString() + "to C_dat");
                                                                //Console.WriteLine("added to C_dat");
                                                            }

                                                        }
                                                        //This is what the selected populations (T) populations (high/low) are being compared with
                                                        else if (info.Length == 5 && j >= 12) //data present
                                                        {
                                                            var AD = info[1].Split(',');
                                                            AD = info[1].Split(',');
                                                            if (Convert.ToInt16(AD[0]) + Convert.ToInt16(AD[1]) > 0)
                                                            {
                                                                T_count++;
                                                                T_dat.Add(Convert.ToInt16(AD[0]));
                                                                //Console.WriteLine("Added " + AD[0].ToString() + " to T_dat");
                                                                T_dat.Add(Convert.ToInt16(AD[1]));
                                                                //Console.WriteLine("Added " + AD[1].ToString() + " to T_dat");
                                                                //Console.WriteLine("added line to T_dat");
                                                            }
                                                        }
                                                        //Console.WriteLine("Working on " + j + " in jrange."); 
                                                    }
                                                }
                                                //print dip_count, tet_count
                                                //Console.WriteLine("C_count = " + C_count.ToString() + ". T_count = " + T_count.ToString());
                                                if (C_count >= 0 && T_count >= 0)
                                                {
                                                    double[] qC = new double[] { 0.0, 0.0 };
                                                    double[] qT = new double[] { 0.0, 0.0 };
                                                    for (int j = 0; j < C_dat.Count / 2; j++)
                                                    {
                                                        var m = Convert.ToDouble(C_dat[2 * j]) + Convert.ToDouble(C_dat[2 * j + 1]);
                                                        qC[0] += Convert.ToDouble(C_dat[2 * j]); //ref allele freq
                                                        qC[1] += Convert.ToDouble(m);
                                                    }
                                                    for (int j = 0; j < T_dat.Count / 2; j++)
                                                    {
                                                        var m = Convert.ToDouble(T_dat[2 * j]) + Convert.ToDouble(T_dat[2 * j + 1]);
                                                        qT[0] += Convert.ToDouble(T_dat[2 * j]); //ref allele freq
                                                        qT[1] += Convert.ToDouble(m);
                                                    }
                                                    //Console.WriteLine("qT row count: " +qT.Length.ToString());
                                                    //Console.WriteLine("qC row count: " + qC.Length.ToString());
                                                    if (qT[1] >= Min_reads && qC[1] >= Min_reads)
                                                    {
                                                        qCqTLines.Add("1 \t" + qT[1].ToString() + "\t" + qC[1].ToString() + "\t" + "[" + string.Join(", ", C_dat.ToArray()) + "]\t" + C_count.ToString() + "\t" + "[" + string.Join(", ", T_dat.ToArray()) + "]" + T_count.ToString());
                                                        //Console.WriteLine("About to increment num_snps");
                                                        Locations.Add(rowdata[0] + "_" + rowdata[1]);
                                                        num_snps++;
                                                        //Console.WriteLine(rowdata[0] + "\t" + rowdata[1]);
                                                        var qC_hat = qC[0] / qC[1];
                                                        var qT_hat = qT[0] / qT[1];

                                                        var var_C = 1.0 / qC[1];
                                                        var var_T = 1.0 / qT[1];

                                                        Var_snp_specific += (var_C + var_T); //is this really appending?
                                                        var vdiv = Var_neutral + var_C + var_T;
                                                        var diverge = 2.0 * (Math.Asin(Math.Pow(qT_hat, 0.5)) - Math.Asin(Math.Pow(qC_hat, 0.5)));
                                                        resultsList.Add(rowdata[0] + '\t' + rowdata[1] + '\t' + qC[1].ToString() +
                                                        '\t' + qC_hat.ToString() + '\t' + qT[1].ToString() +
                                                        '\t' + qT_hat.ToString() + '\t' + diverge.ToString() +
                                                        '\t' + vdiv.ToString() + '\t' +
                                                        (diverge / (Math.Pow(vdiv, 0.5))).ToString() + '\n');
                                                        Random rnd = new Random();
                                                        if (rnd.Next(1, 3) == 1)
                                                        {
                                                            zraw.Add(diverge);
                                                        }
                                                        else
                                                        {
                                                            zraw.Add(-diverge);
                                                        }
                                                        z_std.Add(diverge / (Math.Pow(vdiv, 0.5)));
                                                        q_T.Add(qT_hat);
                                                        q_C.Add(qC_hat);
                                                        zList.Add(diverge / (Math.Pow(vdiv, 0.5)) + '\n');
                                                    }
                                                    else
                                                    {
                                                        qCqTLines.Add("0 \t" + qT[1].ToString() + "\t" + qC[1].ToString() + "\t" + "[" + string.Join(", ", C_dat.ToArray()) + "]\t" + C_count.ToString() + "\t" + "[" + string.Join(", ", T_dat.ToArray()) + "]" + T_count.ToString());
                                                        passedLines.Add("qT < Min_reads and/or qC < Min_reads");
                                                    }
                                                }
                                                //}
                                                //}



                                                vcfTable.Rows.Add(line.Split('\t'));
                                            }
                                            else
                                            {
                                                passedLines.Add("Col 4 is empty or null or has more than 1 alt_base) + " + line.ToString());
                                            }
                                        }
                                        else
                                        {
                                            passedLines.Add("(Col3 or 4 is empty or null) " + line.ToString());
                                        }
                                    }
                                    else
                                    {
                                        passedLines.Add("(line sNNfold # too high/low) " + line.ToString());
                                    }
                                }
                                else
                                {
                                    passedLines.Add("(line null/doesn't contain sNNfold) " + line.ToString());
                                }
                            }
                            catch (Exception ex)
                            {
                                Console.WriteLine("Exception: " + ex.Message);
                                Console.WriteLine("With line " + line.ToString());
                                //Example:
                                //Exception: Object reference not set to an instance of an object
                                //this occurs at this point and breaks the memory
                                continue; //leaving off because it will break the program
                                          //break;
                            }
                        }
                    }
                }
            }
            //make text file of lines that were passed
            //create timestamp MMDDYYHHMM
            using (StreamWriter passedFile = File.CreateText("Passed_Lines_" + DateTime.Now.ToString("yyMMddHHmm") + ".txt"))
            {
                passedFile.WriteLine("Passed a total of " + passedLines.Count);
                foreach (string line in passedLines)
                {
                    passedFile.WriteLine(line);
                }
            }
            using (StreamWriter qcqtFile = File.CreateText("qCqT_" + DateTime.Now.ToString("yyMMddHHmm") + ".txt"))
            {
                qcqtFile.WriteLine(qCqTLines.Count + " lines in file.");
                foreach (string line in qCqTLines)
                {
                    qcqtFile.WriteLine(line.ToString());
                }
            }

            //stop stopwatch
            stopwatch.Stop();
            Console.WriteLine("Table build time: {0}", stopwatch.Elapsed);
            Console.WriteLine("Table includes " + num_snps.ToString() + " SNPs.");
            Console.WriteLine("Sampling/genotyping variance " + Var_snp_specific / num_snps + ".");
            Console.WriteLine("Identified the following as sets of unique genetic sets:");

            //display scaffolds of interest
            ArrayList distinctScaffolds = new ArrayList();
            foreach (DataRow row in vcfTable.Rows)
            {
                if (!distinctScaffolds.Contains(row["#CHROM"].ToString()))
                {
                    distinctScaffolds.Add(row["#CHROM"].ToString());
                }
            }
            distinctScaffolds.Sort();
            Console.WriteLine("You have the following range of distinct scaffolds present in your data:");
            foreach (var item in distinctScaffolds)
            {
                Console.WriteLine(item.ToString());
            }

            //set things up from the chisq file
            Console.WriteLine("Enumerating and working with chisq file.");
            var df = new ArrayList();
            var bs = new ArrayList();
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
                bs.Add(rx);
                line_idx++;
            }

            //continue
            zraw.Sort();
            Console.WriteLine(num_snps.ToString());
            var n25 = zraw[num_snps / 4];
            var n50 = zraw[num_snps / 2];
            var n75 = zraw[3 * num_snps / 4];
            Console.WriteLine("Z Percentiles " + n25 + ", " + n50 + ", " + n75);
            Console.WriteLine("Total varience in Z (based on IQR) " + (Math.Pow((Convert.ToDouble(n75) / 1.349 - Convert.ToDouble(n25) / 1.349), 2)));
            //conflict with scope for cols addint to vzllist (out1)

            //wait for key press to end
            Console.ReadKey();
        }
    }
}
