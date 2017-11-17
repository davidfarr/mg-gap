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
            //set up linq version for speed?
            List<QTLEnumerator> qtl_selection = new List<QTLEnumerator>();
            using (var filesStream = File.OpenRead(csvpath))
            using (var streamReader = new StreamReader(filesStream))
            {
                String line;
                while ((line = streamReader.ReadLine()) != null)
                {
                    string[] cols = line.Replace("\n", "").Split(',');
                    if (cols[7].ToString() != "padj")
                    {
                        string chrom = cols[19].ToString();
                        chrom = chrom.Remove(0, chrom.IndexOf("_") + 1);

                        QTLEnumerator row = new QTLEnumerator
                        {
                            Chromosome = Convert.ToInt16(chrom),
                            Description = cols[22],
                            P_Adj = cols[7],
                            Gene = cols[1],
                            QTL_P = cols[6],
                            StartRange = Convert.ToInt32(cols[20]),
                            EndRange = Convert.ToInt32(cols[21])
                        };

                        qtl_selection.Add(row);
                    }
                }
            }

            Console.WriteLine("\n");
            for(int i = 0;i < fdrList.Count;i++)
            {
                Console.Write("\rFinding SNP {0} of {1}...\t", i + 1, fdrList.Count());
                try
                {
                    QTLEnumerator match = qtl_selection.Find(x => (fdrList[i].Chromosome == x.Chromosome) && (fdrList[i].Basepair >= x.StartRange) && (fdrList[i].Basepair <= x.EndRange));
                    if (match != null)
                    {
                        fdrList[i].Description = match.Description;
                        fdrList[i].Adjusted_P = match.P_Adj;
                        fdrList[i].RnaSeqPval = match.QTL_P;
                        fdrList[i].Gene = match.Gene;
                    }
                }
                catch
                {
                    continue;
                }
            }
            Console.WriteLine("\n");


                //for (int i = 0; i < fdrList.Count(); i++)
                //{
                //    using (var fileStream = File.OpenRead(csvpath))
                //    using (var streamReader = new StreamReader(fileStream))
                //    {
                //        String line;
                //        while ((line = streamReader.ReadLine()) != null)
                //        {
                //            string[] cols = line.Replace("\n", "").Split(','); //true CSV file
                //            if (cols[7].ToString() != "padj")
                //            {
                //                string chrom = cols[19].ToString();
                //                chrom = chrom.Remove(0, chrom.IndexOf("_") + 1);
                //                if (fdrList[i].Chromosome.ToString() == chrom)
                //                {
                //                    if (Convert.ToInt32(cols[20]) < fdrList[i].Basepair && fdrList[i].Basepair > Convert.ToDouble(cols[21]))
                //                    {
                //                        fdrList[i].Description = cols[22];
                //                        fdrList[i].RnaSeqPval = cols[6];
                //                        fdrList[i].Adjusted_P = cols[7];
                //                        fdrList[i].Gene = cols[1];
                //                    }
                //                }
                //            }
                //        }
                //    }
                //}
                return fdrList;
        }
    }
}
