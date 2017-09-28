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
            using (var fileStream = File.OpenRead(csvpath))
            using (var streamReader = new StreamReader(fileStream))
            {
                String line;
                while ((line = streamReader.ReadLine())!= null)
                {
                    string[] cols = line.Replace("\n", "").Split(','); //true CSV file
                    string chrom = cols[19].ToString();
                    chrom = chrom.Remove(0, chrom.IndexOf("_") + 1);
                    foreach (SNP snp in fdrList)
                    {
                        if (snp.Chromosome.ToString() == chrom)
                        {
                            if (Convert.ToInt32(cols[20]) <  snp.Basepair && snp.Basepair > Convert.ToDouble(cols[21]))
                            {
                                snp.Description = cols[31];
                                snp.Adjusted_P = Convert.ToDouble(cols[7]);
                            }
                        }
                    }
                }
            }

                return fdrList;
        }
    }
}
