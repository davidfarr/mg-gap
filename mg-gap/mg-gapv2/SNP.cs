using System;
namespace mg_gapv2
{
    public class SNP
    {
        //set up the properties of each SNP
        private int bp = 0;
        private int chr = 0;
        private string chr_asID = string.Empty;
        private string AD = string.Empty;

        public int basepair
        {
            get
            {
                return bp;
            }
            set
            {
                bp = value;
            }
        }

        public int shortCHR
        {
            get
            {
                return chr;
            }
            set
            {
                chr = value;
            }
        }

        public string longCHR
        {
            get
            {
                return chr_asID;
            }
            set
            {
                chr_asID = value;
            }
        }
        public string alleleDepth
        {
            get
            {
                return AD;
            }
            set
            {
                AD = value;
            }
        }
    }
}
