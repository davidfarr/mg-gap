using System;
using System.Data;
using System.Collections;
using System.IO;
using System.Text.RegularExpressions;
using System.Diagnostics;
using System.Linq;
using System.Collections.Generic;
using RDotNet;
using RDotNet.NativeLibrary;


namespace v1_gap
{
    class MainClass
    {
        public static void Main(string[] args)
        {
            //set up the filepath - in this version it's hard-coded
            string vcf_path = "N:/app dev/scoville research/program files/dev migration for windows/vcf files/Ali_w_767.vcf";
            //vcf_path = "N:/app dev/scoville research/program files/dev migration for windows/vcf files/Ali.vcf";
            vcf_path = @"//ENTROPY/All/School/Biology Research/Ali_w_767.vcf"; //private environment
            string qtlpath = "N:/app dev/scoville research/program files/github repo/mg-gap/mg-gap/mg-gap/support files/QTL10_778_781_RNASEQ_2016.csv";
            qtlpath = "C:/Users/David/Documents/GitHub/mg-gap/mg-gap/mg-gap/support files/QTL10_778_781_RNASEQ_2016.csv";
            string chisq_path = @"N:/app dev/scoville research/program files/github repo/mg-gap/mg-gap/mg-gap/support files/chisq.txt";
            chisq_path = @"C:/Users/David/Documents/GitHub/mg-gap/mg-gap/mg-gap/support files/chisq.txt";

            //run the vcf parser for SNP window of 1
            Stopwatch methodTime = new Stopwatch();
            Console.WriteLine("Starting B processing at " + DateTime.Now + "...");
            methodTime.Start();

            // TODO add this in, but maybe pass variable instead of writing to file, the variable is a list of SNPs
            // but maybe gen win needs a file to read in the script
            //new way for R work prep
            //using (StreamWriter write_results = File.CreateText("B1_new.txt"))
            //{
            //    write_results.WriteLine("CHR\tBP\tB");

            //    foreach (mg_gap.SNP snp in mg_gap.VCF_Analyzer.SNP_list(1, vcf_path, 'N',chisq_path))
            //    {
            //        write_results.WriteLine("" + snp.Chromosome + '\t' + snp.Basepair + '\t' + snp.B_standard);
            //    }
            //}


            methodTime.Stop();
            Console.WriteLine("B processing time: " + methodTime.Elapsed.ToString());
            // TODO after b values in a window of s=1 is determined above, this feeds results into genwin
            //try
            //{
            //    //do R stuff
            //    Console.WriteLine("Beginning R execution at " + DateTime.Now);
            //    Stopwatch r_stopwatch = new Stopwatch();
            //    r_stopwatch.Start();
            //    REngine engine = REngine.GetInstance();
            //    Console.WriteLine("Successfully created R engine instance. Evaluating script...");
            //    string rscriptpath = @"N:/app dev/scoville research/program files/github repo/mg-gap/mg-gap/mg-gap/support files/GenWin_script_12_29_2016.R"; //lab env
            //    //string rscriptpath = @"C:/Users/David/Documents/GitHub/mg-gap/mg-gap/mg-gap/support files/GenWin_script_12_29_2016.R"; //personal env
            //    engine.Evaluate(@"source('" + rscriptpath + "')");
            //    r_stopwatch.Stop();
            //    Console.WriteLine("R exited successfully at " + DateTime.Now + "\nRun time " + r_stopwatch.Elapsed.ToString());
            //}
            //catch (Exception ex)
            //{
            //    Console.WriteLine(ex.Message);
            //}

            //after this then find the median out of the output file and then run the b processing at that window and then b* processing - then move to java program
            //check if we recognize that the file has now been put there
            //string splinepath = "N:/app dev/scoville research/program files/github repo/mg-gap/mg-gap/mg-gap/bin/Debug/splinewindows.txt";
            //splinepath = "C:/Users/David/Documents/GitHub/mg-gap/mg-gap/mg-gap/bin/Debug/splinewindows.txt"; //private mac env


            double median = 0.0;
            // TODO add this code. it calculates median
            //if (File.Exists(splinepath) == true)
            //{
            //    Console.WriteLine("GenWin file found, obtaining median window value...");
            //    try
            //    {
            //        List<int> windows = new List<int>();
            //        using (var fileStream = File.OpenRead(splinepath))
            //        using (var streamReader = new StreamReader(fileStream))
            //        {
            //            String line;
            //            while ((line = streamReader.ReadLine()) != null)
            //            {
            //                /*Guide to splinewindows format
            //                 * 0 = CHRcol (snnaffold #)
            //                 * 1 = Window start bp position
            //                 * 2 = window end bp position
            //                 * 3 = SNP count (number of snps in the window)
            //                 * 4 = Mean Y
            //                 * 5 = W statistic
            //                 * */
            //                //first row is headers, skip
            //                if (!line.Contains("CHRcol"))
            //                {
            //                    string[] cols = line.Replace("\n", "").Split('\t');
            //                    windows.Add(Convert.ToInt16(cols[3]));
            //                }
            //            }
            //        }
            //        //windows obtained successfully, now get the median
            //        if (windows == null || windows.Count == 0)
            //        {
            //            Console.WriteLine("No window sizes to calculate!");
            //        }
            //        else
            //        {
            //            windows.Sort();
            //            int middle = windows.Count / 2;
            //            median = (windows.Count % 2 != 0) ? (double)windows[middle] : ((double)windows[middle] + (double)windows[middle - 1]) / 2;
            //            Console.WriteLine("The median window size is " + median.ToString() + ".");
            //        }
            //    }
            //    catch (Exception ex)
            //    {
            //        Console.WriteLine(ex.Message);
            //    }
            //}
            //else
            //{
            //    Console.WriteLine("No finished file found from spline analysis!");
            //}

            //now feed the median value back through b processing, then b*
            //need to hold on to the data for FDR though

            //for speeding up testing
            median = 6;
            if (median > 0)
            {
                List<mg_gap.SNP> snpList = mg_gap.VCF_Analyzer.SNP_list(Convert.ToInt32(median), vcf_path, chisq_path);
                Console.WriteLine("Re-analyzing for B* based on median window size " + median + " @ " + DateTime.Now);

                //Start the FDR process
                //This should be a user config variable in the future if they would prefer any different settings
                double fdr_input = 0.05;
                Console.WriteLine("Running FDR analysis at " + fdr_input + "...");
                List<mg_gap.SNP> fdrlist = mg_gap.FDR.Process(snpList, fdr_input);

                //start getting the annotations
                Console.WriteLine("\nFDR Analysis complete, gathering annotations...");
                List<mg_gap.SNP> annotatedlist = mg_gap.Annotator.AnnotatedList(fdrlist, qtlpath);

                using (StreamWriter finalfile = File.CreateText("Bs_" + median + "_FDR5" + "_SCvsT_annotated" + DateTime.Now.Day + DateTime.Now.Hour + DateTime.Now.Year + ".txt"))
                {
                    finalfile.WriteLine("CHR\tBP\tB\tBs\tPraw\tRnaSeqP\tRnaSeqPadj\tFDR_ThresholdValue\tDescription\tGene");
                    foreach (mg_gap.SNP snp in annotatedlist)
                    {
                        finalfile.WriteLine(snp.Chromosome.ToString() + '\t' + snp.Basepair.ToString() + '\t' + snp.B_standard.ToString() + '\t' +
                            snp.B_star.ToString() + '\t' + snp.Raw_p.ToString() + '\t' + snp.RnaSeqPval.ToString() + '\t' + snp.Adjusted_P.ToString() + '\t' + snp.Threshold_Value.ToString() + '\t' + snp.Description.ToString() + '\t' + snp.Gene.ToString());
                    }
                }
            }

            //force to wait before closing
            Console.WriteLine("Program complete, press any key to exit.");
            Console.ReadKey();
        }
    }
}

