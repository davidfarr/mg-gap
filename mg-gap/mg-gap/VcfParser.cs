using System;
using System.Data;
using System.Collections;
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

            //set up counters
            int depth_cc = 0;
            int raw_read_count = 0;
            int num_snps = 0;
            double Var_snp_specific = 0.0;
            ArrayList zraw = new ArrayList();
            ArrayList accepted_snps = new ArrayList();
            int skipcounter = 0;

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
                    else
                    {
                        string scaff = cols[0].Split('_').ToString();
                        string chrom = cols[0].ToString();
                        chrom = chrom.Remove(0, chrom.IndexOf("_") + 1);

                        if (Convert.ToInt16(chrom) < 15) //customized to only go up to snaffold 15
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
                                float C_count = 0;
                                ArrayList C_dat = new ArrayList();
                                float T_count = 0;
                                ArrayList T_dat = new ArrayList();

                                for (int j = 9; j < 9 + 5; j++)
                                {
                                    if (cols.Length < (9 + 5))
                                    {
                                        Console.WriteLine("whoa " + LineID + " " + cols.Length);
                                    }
                                    else
                                    {
                                        string[] info = cols[j].Split(':');
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
                                //this is probably where things are getting messed up - this loop is not filling qC and qT properly
                                double[] qC = { 0.0, 0.0 };
                                double[] qT = { 0.0, 0.0 };
                                if (C_count >= 0 && T_count >= 0)
                                {
                                    //double[] qC = { 0.0, 0.0 };
                                    //double[] qT = { 0.0, 0.0 };
                                    Array.Clear(qC, 0, qC.Length);
                                    Array.Clear(qT, 0, qT.Length);


                                    for (int j = 0; j < C_dat.Count / 2; j++)
                                    {
                                        double m = Convert.ToDouble(C_dat[2 * j]) + Convert.ToDouble(C_dat[2 * j + 1]);
                                        qC[0] = Convert.ToDouble(C_dat[2 * j]);
                                        qC[1] = m;
                                    }
                                    for (int j = 0; j < T_dat.Count / 2; j++)
                                    {
                                        double m = Convert.ToDouble(T_dat[2 * j]) + Convert.ToDouble(T_dat[2 * j + 1]);
                                        qT[0] = Convert.ToDouble(T_dat[2 * j]);
                                        qT[1] = m;
                                    }
                                    if (qT[1] >= min_reads && qC[1] >= min_reads)
                                    {
                                        //skipping output for yut file
                                        Locations.Add(cols[0] + "_" + cols[1]);
                                        num_snps++;
                                        double qC_hat = qC[0] / qC[1];
                                        double qT_hat = qT[0] / qT[1];
                                        double var_C = 1.0 / qC[1];
                                        double var_T = 1.0 / qT[1];
                                        Var_snp_specific += (var_C + var_T);

                                        double diverge = 2.0 * (Math.Asin(Math.Pow(qT_hat, 0.5)) - (Math.Asin(Math.Pow(qC_hat, 0.5))));
                                        raw_results.Add(cols[0] + '\t' + cols[1] + qC[1].ToString() + '\t' + qC_hat.ToString()
                                            + '\t' + qT[1].ToString() + '\t' + qT_hat.ToString() + '\t' + diverge.ToString() + '\n');

                                        Random rnd = new Random();
                                        if (rnd.Next(1, 3) == 1)
                                        {
                                            zraw.Add(diverge);
                                        }
                                        else
                                        {
                                            zraw.Add(-diverge);
                                        }

                                        //SNP accepted.
                                        accepted_snps.Add(cols[0] + "_" + cols[1] + '\t' + var_C + '\t' + var_T);
                                    }
                                    else
                                    {
                                        skipcounter++;
                                    }
                                }
                            }
                        }
                    }
                }
                        

                //wrap up
                Console.WriteLine("Skipped " + skipcounter + " SNPS");
                Console.WriteLine("Accepted " + accepted_snps.Count + "SNPS");
                Console.WriteLine("Sampling/genotyping varience " + (Var_snp_specific / Convert.ToDouble(num_snps)));

                //make the z, rankings
                zraw.Sort();
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
                }
            }

                return bList;
        }
    }
}
