using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace mg_gap
{
    class FDR
    {
        //The FDR process is:
        //1. Sort B* list by p-value
        //2. Add new property FDR where FDR = 0.1 * index of the ranked SNP / (count of SNPs in the file / 2)
        //3*. FDR can be changed... 0.1 above is FDR of 10 and 0.05 is FDR 5

        public static List<SNP> Process (List<SNP> bs_list, double fdr_selected)
        {
            List<SNP> sortedlist = bs_list.OrderBy(x => x.Raw_p).ToList();

            int rankassignment = 1;
            foreach(SNP snp in sortedlist)
            {
                snp.FDR_Rank = rankassignment;
                rankassignment++;
            }
            int numwindow = rankassignment / 2;
            foreach(SNP snp in sortedlist)
            {
                snp.Threshold_Value = (fdr_selected * snp.FDR_Rank) / numwindow;
            }

            int bscapacity = bs_list.Count();

            for (int i = 0; i < bs_list.Count(); i++)
            {
                bs_list[i].Threshold_Value = fdr_selected * (i + 1) / (bscapacity / 2);
            }

            Console.WriteLine("\n{0} SNPs removed below FDR threshold leaving {1}", sortedlist.RemoveAll(x => x.Raw_p > x.Threshold_Value), sortedlist.Count());

            //now show the sig b*
            //Console.WriteLine("\nSignificant B* = {0}\nSignificant B = {1}", bs_list.Min(x => x.B_star));
            //Console.WriteLine("\nSignificant B* = {0}", sortedlist.Min(x => x.B_star));
            return sortedlist;
        }
    }
}
