#---------------------
# Name: VCF parse v3 :: works with AD - Allele Depth
# Purpose: Takes and calls genotypes
# v3: Pools [CA, CB, S] vs. [TA, TB]
# This is a small modification of the script VCF parse v2 written by John Kelly at KU
#---------------------

import sys
import math
import random

if len(sys.argv) != 3: #argument needs to be S, input vcf
    S = 6 ## of snps per window (i.e. window size). The GenWin median was 6 so for expediency this has been manually entered.
    inpath = '//ENTROPY/All/School/Biology Research/'

    #Set up the input/output files
    try:
        print("setting up paths")
        src  =open(inpath+"Ali_w_767.vcf", "rU")
        print("opened vcf")
        in1  =open(inpath+"chisq.txt", "rU")
        print("opened chisq")
        out0 =open(inpath+"Results" + str(6) + ".txt", "w") #a more verbose out3
        out1 =open(inpath+"Vz" + str(6) + ".txt", "w")#not typically used
        out2 =open(inpath+"z" + str(6) + ".txt","w")#not typically used
        out3 =open(inpath+"B" + str(6) + ".txt","w") #this is the important one for getting B, B*, P, in one result file.
    except:
        print("error with paths")

    #Output in a standard format?
    outkk = open(inpath + "Ali_w_767.vcf" + ".yut","w") #left this in, but no current use in application

    #User may define option for:
    Min_reads = 20 #this was left the same - seems standard but may need to be adjusted if data is nanopore sequenced... will revisit.
    
    #Set up arrays for handling data *during* analysis
    LineID = []
    Locations = []
    bs = []
    df = []

    #enumerate chisq file
    percentiles=[[] for i in range(100)]
    for line_idx, line in enumerate(in1):
        cols = line.replace('\n', '').split('\t')       
        df.append(float(cols[0]))
        # bs.append(float(cols[1]))
        for j in range(9):
            percentiles[line_idx].append(float(cols[j+1]))
        rx=(percentiles[line_idx][2]+percentiles[line_idx][0]-2*percentiles[line_idx][1])/(percentiles[line_idx][2]-percentiles[line_idx][0])
        bs.append(rx)

    #Counters
    depth_cc = 0
    raw_read_count = 0
    Var_snp_specific = 0.0
    zraw = []
    accepted_snps = []
    outcomes = [0,0,0]
    num_snps = 0
    z_std = [] #this is a tricky one - use to store diverge until we get a vdiv and then go back through to calculate standardized divergence

    #evaluate contents of each lien of the VCF input file
    for line_idx, line in enumerate(src):
        cols = line.replace('\n', '').split('\t')
        if len(cols) < 2: #should be contig row... these list contigs and have no data for use
            continue #skip the header row
        elif cols[0] == "#CHROM": #this is indicative of the first row of the actual data set. Print out every bam file name
            for i in range(len(cols)):
                print(i,cols[i])
                if i > 8:
                    LineID.append(cols[i])
        else:
            chromosome = cols[0] #snnafold_x
            chromosome = chromosome[9:] #x (we got rid of everything so it's just the number)
            if int(chromosome) < 15: #Many more contigs than chromosomes - we're only looking between chr 1 and 14 for meaningful data.
                        scaff = cols[0].split('_')
                        position = int(cols[1])
                        ref_base = cols[3]
                        alt_base = cols[4]
                        if len(alt_base) > 1: #too many bases at site
                            outcomes[0] += 1 #pass this SNP
                        else:
                            C_count = 0
                            C_dat = []
                            T_count = 0
                            T_dat = []


                            for j in range(10,len(cols)): #ali_w_767.vcf... keep in mind that the data *structure* difference between ali.vcf and w/767 is the 767 bam file is the first bam column index, so must add +1 to all the indexes that referenced bam file cols
                                info = cols[j].split(':')
                                if len(info) == 5 and (j < 13): #ali_w_767.vcf
                                    AD = info[1].split(",")
                                    if int(AD[0]) + int(AD[1]) > 0:
                                        C_count += 1
                                        C_dat.append(int(AD[0]))
                                        C_dat.append(int(AD[1]))
                                if len(info) == 5 and (j > 12): #ali_w_767.vcf
                                    AD = info[1].split(",")
                                    if int(AD[0]) + int(AD[1]) > 0:
                                        T_count += 1
                                        T_dat.append(int(AD[0]))
                                        T_dat.append(int(AD[1]))
                            if C_count >= 0 and T_count >= 0:
                                qC = [0.0,0.0]
                                qT = [0.0,0.0]
                                
                                for j in range(len(C_dat)//2):
                                    m = C_dat[2*j]+C_dat[2*j+1]
                                    qC[0] += float(C_dat[2*j])
                                    qC[1] += float(m)
                                for j in range(len(T_dat)//2):
                                    m = T_dat[2*j] + T_dat[2*j+1]
                                    qT[0] += float(T_dat[2*j])
                                    qT[1] += float(m)
                                if (qT[1] >= Min_reads and qC[1] >= Min_reads):
                                    #get this first so we can tell if it actually works to take the SNP
                                    qC_hat = qC[0]/qC[1]
                                    qT_hat = qT[0]/qT[1]
                                    if (qC_hat != 0 and qC_hat != 1 and qT_hat != 0 and qT_hat != 1):
                                        # columns: scaffold, position, ref allele count 1, total count 1,ref allele count 2, total count 2
                                        outkk.write(cols[0]+"\t"+cols[1]+"\t"+str(qC[0])+"\t"+str(qC[1])+"\t"+str(qT[0])+"\t"+str(qT[1])+"\n")
                                        outcomes[2] += 1
                                        Locations.append(cols[0] + "_" + cols[1])
                                        num_snps += 1
                                        var_C = 1.0/qC[1]
                                        var_T = 1.0/qT[1]

                                        Var_snp_specific += (var_C + var_T)
                                        #vdiv = Var_neutral + var_C + var_T
                                        diverge = 2.0 * (math.asin(qT_hat**0.5)-math.asin(qC_hat**0.5))
                                        out0.write(cols[0] + '\t' + cols[1] + '\t' + str(qC[1]) + '\t' + str(qC_hat) + '\t' + str(qT[1]) + '\t' + str(qT_hat) + '\t' + str(diverge) + '\n')
                                        if random.randrange(1,3) == 1:
                                            zraw.append(diverge)
                                        else:
                                            zraw.append(-diverge)
                                        z_std.append(diverge)
                                        #Major change:
                                        #vdiv used to be calculated here, but it was based off the var_neutral that was highly specific and precalculated in
                                        #the original python script - changing the window size would have resulted in an incorrect var_neutral.
                                        #The solution was to keep a seperate list of only (+) diverge and calculate it programmatically after we get the z percentiles
                                        accepted_snps.append([cols[0] + "_" + cols[1], var_C, var_T])
                                    else:
                                        outcomes[1] += 1

                                if line_idx % 100000 == 0:
                                    print(cols[0],line_idx)



    print("Outcomes ",outcomes)
    print("Included:Number of snps ",num_snps)
    print("Sampling/genotyping variance ", Var_snp_specific/float(num_snps))
    ranked_z=sorted(zraw)
    n25=ranked_z[int(num_snps/4)]
    n50=ranked_z[int(num_snps/2)]
    n75=ranked_z[int(3*num_snps/4)]
    print("Z percentiles (without direction) ",n25,n50,n75)
    print("total variance in Z (based on IQR) ",((n75-n25)/1.349)**2)

    Var_neutral=((n75-n25)/1.349)**2 - Var_snp_specific/float(num_snps)  # this is the bulk variance
    print("Bulk sampling and library variance ",Var_neutral)

# Run through all accepted SNPs to calculate standardized divergence ^ 2 :: B for 1 snp
#Major Change:
#since we only just got var neutral now we need to go back through the z_std list with only the (+) divergence and
# finish the calculation to standardize and keep going
    for k in range(0,len(z_std)):
        vdiv = Var_neutral + accepted_snps[k][1] + accepted_snps[k][2]
        z_std[k] = z_std[k]/(vdiv**0.5)

    ranked_z = sorted(z_std)
    n25 = ranked_z[int(num_snps/4)]
    n50 = ranked_z[int(num_snps/2)]
    n75 = ranked_z[int(3*num_snps/4)]
    print("Zs percentiles ",n25,n50,n75)

    Braw=[]
    Bloc=[]
    for k in range(S,num_snps):
        if k % S == 0 or k % S == S/2:
            vdiv = Var_neutral + accepted_snps[k][1] + accepted_snps[k][2]
            b=0.0
            for j in range(S):
                b+=( z_std[k-j]**2 )
            Bloc.append(Locations[k])
            Braw.append(b)

    ranked_B=sorted(Braw)
    n25=ranked_B[int(len(Braw)/4)]
    n50=ranked_B[int(len(Braw)/2)]
    n75=ranked_B[int(3*len(Braw)/4)]
    print("B percentiles ",n25,n50,n75)
    b_skew=(n75+n25-2*n50)/(n75-n25)
    print("B Bowley skew ",b_skew)

    m = -1
    if b_skew > bs[0]:
        print("Too much skew")
    else:
        for j in range(1,len(bs)):
            if b_skew > bs[j]:
                m = df[j]
                jstar=j
                break

    print ("Degrees of freedom ",m)
    cIQR=percentiles[jstar][2]-percentiles[jstar][0]
    sigB=(n75-n25)*(2*m)**0.5/cIQR

    for j in range(len(ranked_B)):
        bx=Braw[j]
        bs=m+(bx-S)*((2*m)**0.5)/sigB
        
        if bs < percentiles[jstar][3]: # p greater than 0.05
            p=0.5
        elif bs > percentiles[jstar][8]:
            p = 5.0* 10**(1.0-8.5) # p less than table min
        else:
            for k in range(4,9):
                if bs > percentiles[jstar][k-1] and bs <= percentiles[jstar][k]:
                    dx = (bs-percentiles[jstar][k-1])/(percentiles[jstar][k]-percentiles[jstar][k-1])
                    x = k-1+dx
                    p = 5.0* 10**(1.0-x)
                    break

        if j > 0:
            midpoint=str(Bloc[j-1])
        else:
            midpoint=str(Bloc[j])
        out3.write(midpoint+'\t'+str(bx)+'\t'+str(bs)+'\t'+str(p)+'\n')