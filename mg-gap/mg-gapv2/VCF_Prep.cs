using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
namespace mg_gapv2
{
	public class VCF_Prep
{
        public static void RunPrep(string path)
        {
			string vcf_filepath = path;
										   //counter for bam files
			int bamcounter = 0;

			//Set up the list for the output file
			List<string> outputfile = new List<string>();

			//load the vcf into streamreader
			using (var filestream = File.OpenRead(vcf_filepath))
			using (var streamreader = new StreamReader(filestream))
			{
				String line; //this is the kind of generic container for every line of the file
				while ((line = streamreader.ReadLine()) != null) //as long as it's not the end of the file
				{
					//make an array out of the line
					string[] cols = line.Replace("\n", "").Split('\t'); //get rid of any accidental newline characters and split into array based on tab delimiter
					if (cols.Length < 2) //this might not be useful
					{
						continue;
					}
					else if (cols[0] == "#CHROM") //if true then this is the header row, so we will try to preserve some basic info
					{
						string header_row = "CHROM\tSNP\t";

						foreach (string bamfile in cols)
						{
							if (bamfile.Contains(".bam"))
							{
								//inform the user we caught it, then put it in the array file as header row
								//general format:
								//identifier     transformed bamname1 etc etc
								Console.WriteLine("Found BAM file: " + bamfile);
								header_row += bamfile + " Allele Depth\t";
								bamcounter++;
							}
						}
						Console.WriteLine("BAM counter: {0}", bamcounter);
						//add the header row to the output file list
						outputfile.Add(header_row);
					}

					//many VCF files have contigs, so we ignore anything with two ##
					else if (!cols[0].Contains("##"))
					{
						//guide to the cols - keep in mind the only things we are using are cols0,1,3,4 and then the AD for each data set that's not blank (S is)
						//0 = sNNfold_1
						//1 = base pair
						//2 = literally just a period .
						//3 = ref base
						//4 = alt base
						//5 = qual number
						//6 = also literally just a period .
						//7 = Example: AC=4;AF=0.500;AN=8;BaseQRankSum=-0.045;DP=23;Dels=0.00;FS=2.632;HaplotypeScore=0.2493;MLEAC=4;MLEAF=0.500;MQ=33.07;MQ0=2;MQRankSum=3.000;QD=5.11;ReadPosRankSum=-0.313
						//8 = GT:AD:DP:GQ:PL
						//9 = 0/0:5,0:5:12:0,12,135 (CA set)
						//10 = 0/1:8,1:9:16:16,0,185 (CB set)
						//11 = ./. (S set)
						//12 = 1/1:2,2:4:6:80,6,0 (TA set)
						//13 = 0/1:4,1:5:31:31,0,49 (TB set)
						//keep in mind that not every snp has every set for data - some might have s, some might have ca, etc

						//break down cols[0], (snnaffold_x) to just x since everybody does their own "snnaffold" equivalent. This leaves CHROM
						string scaff = cols[0].Split('_').ToString();
						string chrom = cols[0].ToString();
						chrom = chrom.Remove(0, chrom.IndexOf("_") + 1); //from now on chrom is just the number
						if (Convert.ToInt16(chrom) < 15) //this is a custom line - only look at SNP's in chromosomes 1-15
						{
							int basepair = Convert.ToInt32(cols[1]);
							string ref_base = cols[3];
							string alt_base = cols[4];
							//ignore anything with multiple bases 
							if (alt_base.Length < 2)
							{
								//string for each actual data row to add to until it goes to output
								string prelim_line = chrom + "\t" + basepair + "\t";

								for (int i = cols.Length - bamcounter; i < cols.Length; i++)//for every bam file assuming all actual bam data is at the end of the row
								{
									string[] info = cols[i].Split(':');
									//splits the data-rich mess into an array dividing between stats
									//Guide to info array total 5 columns
									//* 0 = 0/0 (GT)
									//* 1 = 5,0 (this is AD)
									//* 2 = 5 (DP)
									//* 3 = 12 (GQ)
									//* 4 = 0,12,135 (PL)
									//We only want AD from each set
									//potential problem - you get the ref and alt depth frequency here but next code uses CA/CB then TA/TB for SC vs T
									if (info.Count() > 1)
									{
										prelim_line += info[1] + "\t";
									}
									else
									{
										prelim_line += "\t";
									}

								}
								//this should wrap up the line; add the prelim to the output
								outputfile.Add(prelim_line);
							}
						}

					}
				}
			}


			//write the file
			using (StreamWriter outputwriter = File.CreateText("VCF-Intermediate.csv"))
				foreach (var line in outputfile)
				{
					outputwriter.WriteLine(line);
				}
        }

}
}
