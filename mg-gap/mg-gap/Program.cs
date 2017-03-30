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
            Console.WriteLine("Starting B processing...");
            methodTime.Start();
            int windowsize = 1;
            List<string> bResult = mg_gap.VcfParser.b_processing(windowsize, vcf_path);
            methodTime.Stop();
            Console.WriteLine("B processing time: " + methodTime.Elapsed.ToString());
            using (StreamWriter bfilenew = File.CreateText("B" + windowsize + "_new.txt"))
            {
                foreach (var line in bResult)
                {
                    bfilenew.WriteLine(line);
                }
            }

                //try
                //{
                //    //do R stuff
                //    Console.WriteLine("Beginning R execution");
                //    REngine engine = REngine.GetInstance();
                //    Console.WriteLine("Successfully created R engine instance. Evaluating script...");
                //    string rscriptpath = @"N:/app dev/scoville research/program files/github repo/mg-gap/mg-gap/mg-gap/mg-gap/support files/GenWin_script_12_29_2016.R";
                //    engine.Evaluate(@"source('" + rscriptpath + "')");
                //    Console.WriteLine("R exited successfully.");
                //    //check if we recognize that the file has now been put there
                //    if (File.Exists(@"'N:\\app dev\\scoville research\\program files\\dev migration for windows\\genwin\\splinewindows.txt"))
                //    {
                //        Console.WriteLine("GenWin file completed, created result file.");
                //    }
                //    else
                //    {
                //        Console.WriteLine("No finished file found.");
                //    }
                //}
                //catch (Exception ex)
                //{
                //    Console.WriteLine(ex.Message);
                //}

            //force to wait before closing
            Console.WriteLine("Program complete, press any key to exit.");
            Console.ReadKey();
        }
    }
}

