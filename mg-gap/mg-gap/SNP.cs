﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace mg_gap
{
    class SNP
    {
        //set up the properties that each SNP should have
        private int basepair = 0;
        private int chromosome = 0;
        private double b_standard = 0.00;
        private double b_star = 0.00;
        private double raw_p = 0.00;
        private string adjusted_p = string.Empty;
        private string description = string.Empty;
        private string old_identifier = string.Empty; //snnfold_X
        private double c_variance = 0.00;
        private double t_variance = 0.00;
        private double transformed_c_variance = 0.0;
        private double transformed_t_variance = 0.0;
        private double fdr = 0.0;
        private string gene = string.Empty;

        public int Basepair { get { return basepair; } set { basepair = value; } }
        public int Chromosome { get { return chromosome; } set { chromosome = value; } }
        public double B_standard { get { return b_standard; } set { b_standard = value; } }
        public double B_star { get { return b_star; } set { b_star = value; } }
        public double Raw_p { get { return raw_p; } set { raw_p = value; } }
        public string Adjusted_P { get { return adjusted_p; } set { adjusted_p = value; } }
        public double C_variance { get { return c_variance; } set { c_variance = value; } }
        public double T_variance { get { return t_variance; } set { t_variance = value; } }
        public double Transformed_c_variance { get { return transformed_c_variance; } set { transformed_c_variance = value; } }
        public double Transformed_t_variance { get { return transformed_t_variance; } set { transformed_t_variance = value; } }
        public string Old_identifier { get { return old_identifier; } set { old_identifier = value; } }
        public double FDR { get { return fdr; } set { fdr = value; } }
        public string Description { get { return description; } set { description = value; } }
        public string Gene { get { return gene; } set { gene = value; } }

    }
}
