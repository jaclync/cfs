import re
import operator
import pickle
import numpy

list_sampleinfo_fields = ['sampleid', 'barcodeindex', 'run', 'source', 'project', 'patientid', 'date', 'sampletype', 'comment']

def GetSampleInfo(filename):
	f = open(filename, 'r')
	list_samples = []
	
	# only one line
	for line in f:
		list_tuples = re.findall(r'\(([^)]+)\)', line)
		break
	
	#print list_tuples
	
	for tuple in list_tuples:
		list_fields_int = re.findall(r'(\d+),', tuple)
		tuple = re.sub(r'\d+,', '', tuple)
		list_fields_str = [str for str in re.findall(r'\'([^\']+)\'', tuple)]
		list_samples.append(list_fields_int + list_fields_str)
			
	return list_samples

def CountVJPairs(filename):
	f = open(filename, 'r')
	dict_count = {}
	dict_count_perc = {}
	
	# count VJ pairs
	for line in f:
		if re.findall(r'^>', line):
			id_v = re.findall(r'v=([^ ]+)', line)[0]
			id_j = re.findall(r'j=([^ ]+)', line)[0]
			id_pair = id_v + '#' + id_j
			if (id_pair not in dict_count):
				dict_count[id_pair] = 0
			dict_count[id_pair] += 1
			
	# calculate percentage
	sum_count = sum(dict_count.values())
	for (pair, count) in dict_count.items():
		dict_count_perc[pair] = count * 1.0 / sum_count
			
	return dict_count, dict_count_perc			

def GetCountDictionaries():

	filename_sampleinfo = '../Data/SampleInfo/sampleInfo.txt'
	filename_out_count_normal = '../Output/count_normal.txt'
	filename_out_count_patient = '../Output/count_patient.txt'
	
	list_samples = GetSampleInfo(filename_sampleinfo)
	
	dict_count_normal = {}
	dict_count_patient = {}
	
	dict_all = {}
	
	# output count dictionary to pickle
	fid_dict_count_normal = open('../Pickle/dict_count_normal.pck', 'w')
	fid_dict_count_patient = open('../Pickle/dict_count_patient.pck', 'w')
	fid_dict_all = open('../Pickle/dict_all.pck', 'w')
	
	# get count dictionary of normal and patient samples
	for sample in list_samples:
		print sample	
		filename_fasta = '../Data/JGM_deliverable/fasta/' + sample[list_sampleinfo_fields.index('source')] + '_clusters.fasta'
		dict_count, dict_count_perc = CountVJPairs(filename_fasta)

		# update global dictionary
		for (pair, count) in dict_count.items():
			label = 'NOR'
			# normal sample
			if (sample[list_sampleinfo_fields.index('project')] == 'NOR'):
				dict_count_all = dict_count_normal
				if (pair not in dict_count_all):
					dict_count_all[pair] = 0
				dict_count_all[pair] += count
				dict_count_normal = dict_count_all
				label = 'NOR'
				
			# patient sample
			elif (sample[list_sampleinfo_fields.index('project')] == 'CFS'):
				dict_count_all = dict_count_patient
				if (pair not in dict_count_all):
					dict_count_all[pair] = 0
				dict_count_all[pair] += count
				dict_count_patient = dict_count_all
				label = 'CFS'
				
			# add to dict_all
			if (pair not in dict_all.keys()):
				dict_all[pair] = {'NOR':{'count':[], 'perc':[]}, 'CFS':{'count':[], 'perc':[]}}
			dict_all[pair][label]['count'].append(count)
			dict_all[pair][label]['perc'].append(dict_count_perc[pair])
	
	# for normal/patient, sort by count
	sorted_count_normal = sorted(dict_count_normal.iteritems(), key=operator.itemgetter(1))
	sorted_count_normal.reverse()
	sorted_count_patient = sorted(dict_count_patient.iteritems(), key=operator.itemgetter(1))
	sorted_count_patient.reverse()
	
	# write sorted counts to files
	# normal
	f_out = open(filename_out_count_normal, 'w')
	for (key, val) in sorted_count_normal:
		f_out.write('%s\t%s\t%d\n' % (key.split('#')[0], key.split('#')[1], val))
	f_out.close()
	# patient
	f_out = open(filename_out_count_patient, 'w')
	for (key, val) in sorted_count_patient:
		f_out.write('%s\t%s\t%d\n' % (key.split('#')[0], key.split('#')[1], val))
	f_out.close()
	
	# dump dictionary to file
	# normal
	pickle.dump(dict_count_normal, fid_dict_count_normal)
	fid_dict_count_normal.close()
	del fid_dict_count_normal
	# patient
	pickle.dump(dict_count_patient, fid_dict_count_patient)
	fid_dict_count_patient.close()
	del fid_dict_count_patient
	# all
	pickle.dump(dict_all, fid_dict_all)
	fid_dict_all.close()
	del fid_dict_all

def ReadPickleDictionary(filename):
	fid_dict = open(filename, 'r')
	return pickle.load(fid_dict)
	
def PrintAvgPercentage(dict_all):
	fid_out = open('../Output/avgPercentage.txt', 'w')
	fid_out.write('V\tJ\tNormal\tCFS\n')
	for (pair, dict) in dict_all.items():
		list_labels = ['NOR', 'CFS']
		dict_pair_perc = {}
		for label in list_labels:
			list_perc = dict[label]['perc']
			if (len(list_perc) > 0):
				dict_pair_perc[label] = numpy.mean(list_perc)
			else:
				dict_pair_perc[label] = 0.0
		fid_out.write('%s\t%s\t%1.6f\t%1.6f\n' % (pair.split('#')[0], pair.split('#')[1], dict_pair_perc['NOR'], dict_pair_perc['CFS']))
			
	fid_out.close()
	
if __name__ == '__main__':	
	# only run this line once
	GetCountDictionaries()
	
	dict_all = ReadPickleDictionary('../Pickle/dict_all.pck')
	PrintAvgPercentage(dict_all)
	
	
	
