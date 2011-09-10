#!/usr/bin/env python
#
# $File: test_import_vcf.py $
# $LastChangedDate: 2011-06-16 20:10:41 -0500 (Thu, 16 Jun 2011) $
# $Rev: 4234 $
#
# This file is part of variant_tools, a software application to annotate,
# summarize, and filter variants for next-gen sequencing ananlysis.
# Please visit http://variant_tools.sourceforge.net # for details.
#
# Copyright (C) 2004 - 2010 Bo Peng (bpeng@mdanderson.org)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#

import os
import glob
import unittest
import subprocess
from testUtils import ProcessTestCase, runCmd, numOfVariant

ANNOVAR_data = '''\
1 161003087 161003087 C T comments: rs1000050, a SNP in Illumina SNP arrays
1 84647761 84647761 C T comments: rs6576700 or SNP_A-1780419, a SNP in Affymetrix SNP arrays
1 13133880 13133881 TC - comments: rs59770105, a 2-bp deletion
1 11326183 11326183 - AT comments: rs35561142, a 2-bp insertion
1 105293754 105293754 A ATAAA comments: rs10552169, a block substitution
1 67478546 67478546 G A comments: rs11209026 (R381Q), a SNP in IL23R associated with Crohn's disease
2 233848107 233848107 T C comments: rs2241880 (T300A), a SNP in the ATG16L1 associated with Crohn's disease
16 49303427 49303427 C T comments: rs2066844 (R702W), a non-synonymous SNP in NOD2
16 49314041 49314041 G C comments: rs2066845 (G908R), a non-synonymous SNP in NOD2
16 49321279 49321279 - C comments: rs2066847 (c.3016_3017insC), a frameshift SNP in NOD2
13 19661686 19661686 G - comments: rs1801002 (del35G), a frameshift mutation in GJB2, associated with hearing loss
13 19695176 20003944 0 - comments: a 342kb deletion encompassing GJB6, associated with hearing loss
'''

Illumina_SNP_data	=	'''\
#	**	CASAVA	depth-filtered	snp	calls	**
#$	CMDLINE	/CASAVA-1.8.0a19/filterSmallVariants.pl	--chrom=chr1
#$	SEQ_MAX_DEPTH	chr1	142.345267150165
#
#$	COLUMNS	seq_name	pos	bcalls_used	bcalls_filt	ref	Q(snp)	max_gt	Q(max_gt)	max_gt|poly_site	Q(max_gt|poly_site)	A_used	C_used	G_used	T_used
chr1	10231	5	9	C	28	AC	28	AC	59	3	2	0	0
chr1	10255	14	29	A	1	AA	9	AT	25	12	0	0	2
chr1	10264	15	19	C	18	AC	18	AC	51	4	11	0	0
chr1	10291	2	16	C	1	CC	10	CT	21	0	1	0	1
chr1	10330	3	14	C	2	CC	5	AC	28	2	1	0	0
chr1	13273	9	0	G	58	CG	54	CG	57	0	6	3	0
chr1	14464	18	0	A	60	AT	60	AT	93	12	0	0	6
chr1	14673	19	0	G	63	CG	63	CG	96	0	8	11	0
chr1	14699	23	0	C	72	CG	72	CG	105	0	14	9	0
chr1	14907	13	0	A	118	AG	65	AG	65	4	0	9	0
chr1	14930	14	2	A	119	AG	68	AG	68	5	0	9	0
chr1	14933	14	2	G	78	AG	78	AG	110	6	0	8	0
chr1	14976	4	0	G	18	AG	18	AG	47	2	0	2	0
chr1	15211	2	0	T	37	GG	5	GG	5	0	0	2	0
chr1	15817	1	0	G	11	GT	3	GT	3	0	0	0	1
chr1	15820	1	0	G	11	GT	3	GT	3	0	0	0	1
chr1	16487	12	0	T	62	CT	62	CT	94	0	6	0	6
chr1	17538	64	0	C	88	AC	88	AC	121	18	46	0	0
chr1	17746	53	1	A	22	AG	22	AG	55	39	0	14	0
chr1	17765	47	1	G	26	AG	26	AG	59	13	0	34	0
chr1	20131	1	0	G	8	CG	2	CG	3	0	1	0	0
chr1	20144	1	0	G	9	AG	2	AG	3	1	0	0	0
chr1	20206	2	0	C	4	CT	4	CT	30	0	1	0	1
chr1	20245	3	0	G	4	AG	4	AG	34	1	0	2	0
chr1	20304	2	0	G	2	GG	5	CG	27	0	1	1	0
'''

Illumina_INDEL_data	=	'''\
#	**	CASAVA	depth-filtered	indel	calls	**
#$	CMDLINE	/filterSmallVariants.pl	--chrom=chr1
#$	SEQ_MAX_DEPTH	chr1	143.988337798483
#
#$	COLUMNS	seq_name	pos	type	ref_upstream	ref/indel	ref_downstream	Q(indel)	max_gtype	Q(max_gtype)	depth	alt_reads	indel_reads	other_reads	repeat_unit	ref_repeat_count	indel_repeat_count
chr1	10147	1D	CTAACCCTAA	C/-	CCCTAACCCT	70	het	70	6	2	3	1	C	4	3
chr1	10231	1D	CTAACCCTAA	C/-	CCCTAACCCT	1203	het	284	53	7	30	17	C	4	3
chr1	10353	1I	CCCTAACCCT	-/A	ACCCTAACCC	434	het	118	17	3	8	9	A	1	2
chr1	10390	1D	CTAACCCTAA	C/-	CCCTAACCCC	765	het	399	39	9	19	12	C	4	3
chr1	10397	1D	TAACCCCTAA	C/-	CCCTAACCCT	730	het	496	38	11	20	9	C	4	3
chr1	10440	1D	CTAACCCTAA	C/-	CCCTAACCCT	774	het	302	31	7	21	3	C	4	3
chr1	28327	1D	AAGCCTGTAG	T/-	TGCTCATCTG	3	het	3	2	1	1	0	T	2	1
chr1	54711	1I	AAACCTTGTA	-/T	TTTTTCTTTC	37	het	37	21	8	2	12	T	5	6
chr1	62240	2D	AGACACACAT	AC/--	ACACACACAC	100	het	100	22	16	4	2	AC	8	7
chr1	83830	8D	AGAAAGAAAG	AGAAAGAA/--------	AGAAAGAAAG	273	het	161	13	3	6	4	AGAA	11	9
chr1	108546	BP_RIGHT	N/A	------/CTATCA	AAAAAAAAAA	28	het	28	13	9	2	2	N/A	0	0
chr1	123089	2D	TGTGGACATG	TA/--	TATATATATA	142	het	142	13	9	4	0	TA	6	5
chr1	128590	1D	CTTCAAGTTC	A/-	CCCCCTTTTT	220	het	220	13	4	5	8	A	1	0
chr1	129011	3D	GGGATGTAGA	ATG/---	ATAAGGCTCT	258	het	258	12	5	6	1	ATG	1	0
chr1	136743	1I	GGTGAGGCAA	-/C	GGGCTCACAC	76	het	76	80	66	6	12	C	0	1
chr1	136889	1D	TGTGAGGCAA	G/-	GGGCTCGGGC	205	het	205	41	29	8	8	G	4	3
chr1	237577	1I	AAAGGGGGTT	-/C	ATTATCTAGG	60	het	60	51	45	4	2	C	0	1
chr1	247917	3D	ACCCAACCTC	AGG/---	AGTTCAGGGC	69	hom	5	2	0	2	0	AGG	1	0
chr1	255910	2I	TGTGTGTGTA	--/TG	TGTGTGTGTG	257	het	28	7	1	5	1	TG	10	11
chr1	531809	2D	CACACTTATG	CA/--	CACATTCACA	327	het	327	25	17	8	1	CA	3	2
chr1	532239	2D	TGTTCACATT	CA/--	CACTCATACA	325	het	325	64	53	10	2	CA	2	1
chr1	532259	3D	CACAGCCCAA	AAT/---	AATATACACA	303	het	303	61	43	9	10	AAT	2	1
chr1	537252	2D	AGCCACATGT	GG/--	GACAGGGCAG	6	hom	2	1	0	1	0	G	3	1
chr1	537494	5I	CAGCGTCCAT	-----/GCCCA	GCCGGCCTCC	23	het	3	2	0	1	1	GCCCA	0	1
chr1	537641	50D	ATCCCCCTCT	CCATCCCCCTCTCCATCTCCCTCTCCTTTCTCCTCTCTAGCCCCCTCTCC/--------------------------------------------------	TTTCTCCTCT	66	het	66	22	18	3	10	CCATCCCCCTCTCCATCTCCCTCTCCTTTCTCCTCTCTAGCCCCCTCTCC	1	0
'''

class TestImportTXT(ProcessTestCase):
    def setUp(self):
        'Create a project'
        runCmd('vtools init test -f')

    def removeProj(self):
        runCmd('vtools remove project')

    def testImportTXT(self):
        'Test command import_txt'
        self.assertFail('vtools import_txt')
        self.assertFail('vtools import_txt non_existing.txt')
        # help information
        self.assertSucc('vtools import_txt -h')
        # no format information, fail
        self.assertFail('vtools import_txt input.tsv')
        # no build information, fail
        self.assertFail('vtools import_txt --format ../input_fmt/ANNOVAR input.tsv')

    def testANNOVAR(self):
        'Testing the annovar input format'
        input = open('ANNOVAR.txt', 'w')
        input.write(ANNOVAR_data)
        input.close()
        self.assertSucc('vtools import_txt --build hg18 --format ../input_fmt/ANNOVAR ANNOVAR.txt')
        # one of the variant cannot be imported.
        self.assertEqual(numOfVariant(), 11)

    def testIllumina_SNP(self):
        'Testing the illumina SNP input format'
        input = open('Illumina_SNP.txt', 'w')
        input.write(Illumina_SNP_data)
        input.close()
        self.assertSucc('vtools import_txt --build hg18 --format ../input_fmt/Illumina_SNP Illumina_SNP.txt')
        self.assertEqual(numOfVariant(), 12)

    def testIllumina_INDEL(self):
        'Testing the illumina INDEL input format'
        input = open('Illumina_INDEL.txt', 'w')
        input.write(Illumina_INDEL_data)
        input.close()
        self.assertSucc('vtools import_txt --build hg18 --format ../input_fmt/Illumina_INDEL Illumina_INDEL.txt')
        self.assertEqual(numOfVariant(), 12)

if __name__ == '__main__':
    unittest.main()
