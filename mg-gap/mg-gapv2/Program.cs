using System;
using System.IO;
using System.Collections.Generic;
using RDotNet;
using RDotNet.NativeLibrary;

namespace mg_gapv2
{
    class Program
    {
        static void Main(string[] args)
        {
            //set file locations
            string vcf_path = "/Users/david/Desktop/Ali.vcf"; //Mac only


            //clean up the messy VCF file to make it easier to work with
            Console.WriteLine("Prepping VCF from " + vcf_path);
            //VCF_Prep.RunPrep(vcf_path);
            Console.WriteLine("Complete - checking file...");
            string ipath = "/Users/david/Documents/visual studio projects/mg-gap/mg-gap/mg-gapv2/VCF-Intermediate.csv";
            Console.WriteLine(File.Exists(ipath) ? "Intermediate VCF file exists." : "File does not exist!");

            //ready to get B values for each snp
            List<string> b_results = Analysis.GetB(ipath);
			using (StreamWriter bfilenew = File.CreateText("B1_new.txt"))
			{
				foreach (var line in b_results)
				{
					bfilenew.WriteLine(line);
				}
			}
            b_results.Clear(); //free up that memory

            string bpath = "/Users/david/Documents/visual studio projects/mg-gap/mg-gap/mg-gapv2/B1_new.txt";
            Console.WriteLine(File.Exists(bpath) ? "B result file exists." : "B result file does not exist!");

			//Throw B into R genwin
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
			string splinepath = "N:/app dev/scoville research/program files/github repo/mg-gap/mg-gap/mg-gap/mg-gap/bin/Debug/splinewindows.txt";
			double median = 0.0;
			if (File.Exists(splinepath) == true)
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
							//first row is headers, skip
							if (!line.Contains("CHRcol"))
							{
								string[] cols = line.Replace("\n", "").Split('\t');
								windows.Add(Convert.ToInt16(cols[3]));
							}
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

        }
    }
}
