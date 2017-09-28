﻿using System;
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

        public static void Process (List<SNP> bs_list, double fdr_selected)
        {
            bs_list.OrderBy(x => x.Raw_p);
            foreach (SNP snp in bs_list)
            {
                snp.FDR = fdr_selected * bs_list.IndexOf(snp) / (bs_list.Count / 2);
                //get the significant b. the sig b is the last snp where the p val < FDR
                if (snp.Raw_p > snp.FDR)
                {
                    bs_list.Remove(snp); //did not make FDR cutoff
                }
            }
            //now show the sig b*
            int lastsnp = bs_list.Count();
            Console.WriteLine("Significant B* = " + bs_list[lastsnp].B_star);
        }
    }
}
