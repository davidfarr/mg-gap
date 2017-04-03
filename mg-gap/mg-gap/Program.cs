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
            //string vcf_path = "/Users/david/Desktop/Ali.vcf";
            string vcf_path = "N:/app dev/scoville research/program files/dev migration for windows/vcf files/ali.vcf";

            //run the vcf parser for SNP window of 1
            Stopwatch methodTime = new Stopwatch();
            Console.WriteLine("Starting B processing at " + DateTime.Now + "...");
            methodTime.Start();
            List<string> bResult = mg_gap.VcfParser.b_processing(1, vcf_path, 'Y'); //window of 1, path, do not get b*
            methodTime.Stop();
            Console.WriteLine("B processing time: " + methodTime.Elapsed.ToString());
            using (StreamWriter bfilenew = File.CreateText("B1_new.txt"))
            {
                foreach (var line in bResult)
                {
                    bfilenew.WriteLine(line);
                }
            }

            try
            {
                //do R stuff
                Console.WriteLine("Beginning R execution at " + DateTime.Now);
                Stopwatch r_stopwatch = new Stopwatch();
                r_stopwatch.Start();
                REngine engine = REngine.GetInstance();
                Console.WriteLine("Successfully created R engine instance. Evaluating script...");
                string rscriptpath = @"N:/app dev/scoville research/program files/github repo/mg-gap/mg-gap/mg-gap/mg-gap/support files/GenWin_script_12_29_2016.R";
                engine.Evaluate(@"source('" + rscriptpath + "')");
                r_stopwatch.Stop();
                Console.WriteLine("R exited successfully at " + DateTime.Now + "\nRun time " + r_stopwatch.Elapsed.ToString());
            }
            catch (Exception ex)
            {
                Console.WriteLine(ex.Message);
            }

            //after this then find the median out of the output file and then run the b processing at that window and then b* processing - then move to java program
            //check if we recognize that the file has now been put there
            string splinepath = @"N:/app dev/scoville research/program files/github repo/mg-gap/mg-gap/mg-gap/mg-gap/bin/Debug/splinewindows.txt";
            double median = 0.0;
            if (File.Exists(@"N:/app dev/scoville research/program files/github repo/mg-gap/mg-gap/mg-gap/mg-gap/bin/Debug/splinewindows.txt") == true)
            {
                Console.WriteLine("GenWin file found, obtaining median window value...");
                try
                {
                    List<int> windows = new List<int>();
                    using (var fileStream = File.OpenRead(splinepath))
                    using (var streamReader = new StreamReader(fileStream))
                    {
                        String line;
                        while ((line = streamReader.ReadLine()) != null)
                        {
                            /*Guide to splinewindows format
                             * 0 = CHRcol (snnaffold #)
                             * 1 = Window start bp position
                             * 2 = window end bp position
                             * 3 = SNP count (number of snps in the window)
                             * 4 = Mean Y
                             * 5 = W statistic
                             * */ 
                            string[] cols = line.Replace("\n", "").Split('\t');
                            windows.Add(Convert.ToInt16(cols[3]));
                        }
                    }
                    //windows obtained successfully, now get the median
                    if (windows == null || windows.Count == 0)
                    {
                        Console.WriteLine("No window sizes to calculate!");
                    }
                    else
                    {
                        windows.Sort();
                        int middle = windows.Count / 2;
                        median = (windows.Count % 2 != 0) ? (double)windows[middle] : ((double)windows[middle] + (double)windows[middle - 1]) / 2;
                        Console.WriteLine("The median window size is " + median.ToString() + ".");
                    }
                }
                catch (Exception ex)
                {
                    Console.WriteLine(ex.Message);
                }
            }
            else
            {
                Console.WriteLine("No finished file found from spline analysis!");
            }

            //now feed the median value back through b processing, then b*
            if (median > 0)
            {
                List<string> b_star_Result = mg_gap.VcfParser.b_processing((int)median, vcf_path, 'N'); //window of 1, path, do not get b*
                using (StreamWriter bsfile = File.CreateText("B" + median + ".txt"))
                {
                    foreach (var line in b_star_Result)
                    {
                        bsfile.WriteLine(line);
                    }
                }
            }

            //force to wait before closing
            Console.WriteLine("Program complete, press any key to exit.");
            Console.ReadKey();
        }
    }
}

