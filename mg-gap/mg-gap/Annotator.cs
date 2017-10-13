using System;
using System.Data;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace mg_gap
{
    class Annotator
    {
        public static List<SNP> AnnotatedList (List<SNP> fdrList, string csvpath)
        {
            for (int i = 0; i< fdrList.Count();i++)
            {
                using (var fileStream = File.OpenRead(csvpath))
                using (var streamReader = new StreamReader(fileStream))
                {
                    String line;
                    while ((line = streamReader.ReadLine()) != null)
                    {
                        string[] cols = line.Replace("\n", "").Split(','); //true CSV file
                        if (cols[7].ToString() != "padj")
                        {
                            string chrom = cols[19].ToString();
                            chrom = chrom.Remove(0, chrom.IndexOf("_") + 1);
                                if (fdrList[i].Chromosome.ToString() == chrom)
                                {
                                    if (Convert.ToInt32(cols[20]) < fdrList[i].Basepair && fdrList[i].Basepair > Convert.ToDouble(cols[21]))
                                    {
                                        fdrList[i].Description = cols[22];
                                        fdrList[i].RnaSeqPval = cols[6];
                                        fdrList[i].Adjusted_P = cols[7];
                                        fdrList[i].Gene = cols[1];
                                        Console.Write("\r{0} of {1} Annotated. ( {2} )\t", i, fdrList.Count(), (i/fdrList.Count())*100);
                                    }
                                }
                        }
                    }
                }
            }
                return fdrList;
        }
    }
}
