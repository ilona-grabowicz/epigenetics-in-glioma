from matplotlib_venn import venn3
from matplotlib import pyplot as plt
from numpy import random
def compare(fdr, what_to_compare):
	q = open('DeSeq2_'+what_to_compare+'.csv', 'r').readlines()
	t = open('DeSeq2_H3K27ac_'+what_to_compare+'_promoters.csv', 'r').readlines()
	u = open('DeSeq2_H3K4me3_'+what_to_compare+'_promoters.csv', 'r').readlines()
		
	DE_gradeX_gradeY = set([line.strip().split('\t')[0][0:15] for line in q[1:] if (line.strip().split('\t')[6]!='NA' and float(line.strip().split('\t')[6])<fdr) and float(line.strip().split('\t')[1])>=10])
	DA_gradeX_gradeY = set([line.strip().split('\t')[0][0:15] for line in t[1:] if (line.strip().split('\t')[6]!='NA' and float(line.strip().split('\t')[6])<fdr) and float(line.strip().split('\t')[1])>=10] )
	DM_gradeX_gradeY = set([line.strip().split('\t')[0][0:15] for line in u[1:] if (line.strip().split('\t')[6]!='NA' and float(line.strip().split('\t')[6])<fdr) and float(line.strip().split('\t')[1])>=10] )
	
	### Figures
	plt.figure()
	out = venn3([DE_gradeX_gradeY, DA_gradeX_gradeY, DM_gradeX_gradeY], ('Differential\nexpression', 'Differential\nH3K27ac', 'Differential\nH3K4me3'))
	label = out.get_label_by_id('111')
	for text in out.set_labels:
		text.set_fontsize(15)
	for x in range(len(out.subset_labels)):
		if out.subset_labels[x] is not None:
			out.subset_labels[x].set_fontsize(15)
	print(what_to_compare)
	if what_to_compare=='grade23_grade4':
		plt.title('DA vs GB/PG', fontsize=18)
	if what_to_compare=='grade1_grade4':
		plt.title('PA vs GB/PG', fontsize=18)
	if what_to_compare=='grade1_grade23':
		plt.title('PA vs DA', fontsize=18)
	# H3K4me3 color:
	out.get_patch_by_id('001').set_color('#DB8D40')
	out.get_patch_by_id('001').set_alpha(0.5)
	# H3K27ac color:
	out.get_patch_by_id('010').set_color('#3B8EDB')
	out.get_patch_by_id('010').set_alpha(0.5)
	# Gene expression color:
	out.get_patch_by_id('100').set_color('#2B9900')
	out.get_patch_by_id('100').set_alpha(0.5)
	#Intercolors:
	out.get_patch_by_id('110').set_color('#3B7C79')
	out.get_patch_by_id('110').set_alpha(0.5)
	if what_to_compare=='grade1_grade4':
		out.get_patch_by_id('011').set_color('#AA938D')
		out.get_patch_by_id('011').set_alpha(0.5)
	out.get_patch_by_id('101').set_color('#897840')
	out.get_patch_by_id('101').set_alpha(0.5)
	
	# Edges color :
	if what_to_compare=='grade1_grade4':
		out.get_patch_by_id('111').set_edgecolor('none')
		out.get_patch_by_id('011').set_edgecolor('none')
	out.get_patch_by_id('001').set_edgecolor('none')
	out.get_patch_by_id('010').set_edgecolor('none')
	out.get_patch_by_id('100').set_edgecolor('none')
	out.get_patch_by_id('110').set_edgecolor('none')
	

	
	print(DE_gradeX_gradeY)
	plt.tight_layout()
	plt.savefig('venn_fdr_'+what_to_compare+'_'+str(fdr)+'.svg') ### svg  i test
	plt.show()
	
	
	'''### Permutation DE and DO: H3K27ac
	overlaps, percentages = [], []
	M_nazwy = [line.strip().split('\t')[0][0:15] for line in q[1:] if float(line.strip().split('\t')[1])>=10] # selecting genes with coverage >10
	M = len(M_nazwy) 
	n = len(DE_gradeX_gradeY)
	N = len(DA_gradeX_gradeY)
	x = len(DE_gradeX_gradeY&DA_gradeX_gradeY)
	c = open('random_DEvsDA_'+what_to_compare+'_overlaps.csv', 'w')
	c.write('counts\tpercentage\n')
	for i in range(0, 100):
		a_gradeX_gradeY_random = random.choice(M_nazwy, n, replace = False)
		b_gradeX_gradeY_random = random.choice(M_nazwy, N, replace = False)
		overlap = len( set(a_gradeX_gradeY_random) & set(b_gradeX_gradeY_random) )
		percent = len( set(a_gradeX_gradeY_random) & set(b_gradeX_gradeY_random) ) / float(M) *100
		c.write(str(overlap)+'\t'+str(percent)+'\n')
		overlaps.append(overlap)
		percentages.append(percent)
	
	c.write('\n'+str(sum(i>=x for i in overlaps))+'\nx='+str(x))
	
	### Permutation DE and DO: H3K4me3
	overlaps, percentages = [], []
	M_nazwy = [line.strip().split('\t')[0][0:15] for line in q[1:] if float(line.strip().split('\t')[1])>=10] # selecting genes with coverage >10
	M = len(M_nazwy) 
	n = len(DE_gradeX_gradeY)
	N = len(DM_gradeX_gradeY)
	x = len(DE_gradeX_gradeY&DM_gradeX_gradeY)
	d = open('random_DEvsDM_'+what_to_compare+'_overlaps.csv', 'w')
	d.write('counts\tpercentage\n')
	for i in range(0, 100):
		a_gradeX_gradeY_random = random.choice(M_nazwy, n, replace = False)
		b_gradeX_gradeY_random = random.choice(M_nazwy, N, replace = False)
		overlap = len( set(a_gradeX_gradeY_random) & set(b_gradeX_gradeY_random) )
		percent = len( set(a_gradeX_gradeY_random) & set(b_gradeX_gradeY_random) ) / float(M) *100
		d.write(str(overlap)+'\t'+str(percent)+'\n')
		overlaps.append(overlap)
		percentages.append(percent)
	
	d.write('\n'+str(sum(i>=x for i in overlaps))+'\nx='+str(x))'''
compare(0.01, 'grade23_grade4')
compare(0.01, 'grade1_grade4')
compare(0.01, 'grade1_grade23')


