import argparse
import copy
import re
import sys
from datetime import datetime
from tqdm import tqdm
from intervaltree import Interval, IntervalTree


def get_args(argv = None):
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--filename", help="take filename.gtf as input")
	parser.add_argument("-s", "--stringtie_output", help="GTF output from stringtie")
	parser.add_argument("-o", "--output", help="name of the gtf output file with extended genes", default = False)
	return parser.parse_args(argv)


def find_pos(filename_ref):
	chromosomes = []
	chromosomes_set = set()
	with open(filename_ref, "r") as filin:
		for line in filin:
			# Pass the header.
			if line.startswith("#!"):
				continue
			else:
			# Formate line as list demilited by tab.			
				line = line.strip()
				line = re.split("\t|;", line)[:-1]
			# Save the name of the first chromosome.	
			if line[2] == "exon":				
				if line[0] not in chromosomes_set:				
					chromosomes.append(line[0])
					chromosomes_set.add(line[0])
				# Find which index of the list correspond to gene_id and transcript_id.
				for n, attributes in enumerate(line):
					if attributes.startswith(" gene_id ") or attributes.startswith("gene_id "):
						gene_col = n
					if attributes.startswith(" transcript_id ") or attributes.startswith("transcript_id "):
						transcript_col = n
				final_transcript = line[transcript_col]	
	first_chromosome = chromosomes[0]
	last_chromosome = chromosomes[-1]
	return [gene_col, transcript_col, first_chromosome, final_transcript]
	
	
def find_exon(filename):
	dict_exon_signal = {}	
	with open(filename, "r") as filin:
		for line in filin:
			# Pass the header.
			if line.startswith("#"):
				continue
			else:
			# Formate line as list demilited by tab.
				line = line.strip()
				line = re.split("\t|;", line)[:-1]
				# Create a new key in the dictionnary for each new chromosome, with values as Intervaltree object
				if line[0] not in dict_exon_signal.keys():
					dict_exon_signal[str(line[0])] = IntervalTree()
				# Add the exon to the IntervalTree				
				if line[2] == "exon":					
					dict_exon_signal[str(line[0])][int(line[3]):int(line[4])] = [line[6], int(round(float(line[-1][6:-1])))]
	return dict_exon_signal
				

def gtf_to_dict(filename_ref, first_chromosome, final_transcript, gene_col, transcript_col):
	dict_transcript = {first_chromosome:IntervalTree()}
	dict_all = {first_chromosome:[]} 
	header = []
	exons_wrong_interval = [] 
	transcript = []
	current_chromosome = first_chromosome
	with open(filename_ref, "r") as filin:
		for line in filin:
			# Save the header.
			if line.startswith("#!"):
				header.append(line)
			else:
			# Formate line as list demilited by tab.			
				line = line.strip()
				line = re.split("\t|;", line)[:-1]
				# If new chromosome is encounter, save the last transcript of the previous one.
				if line[0] not in dict_transcript.keys():
					dict_transcript[str(current_chromosome)][min_transcript:max_transcript] = transcript
					current_chromosome = line[0]
					transcript = []
					# And convert the new chromosome in an intervaltree object.
					dict_transcript[str(current_chromosome)] = IntervalTree()
					dict_all[str(line[0])] = []					
				# Check for errors in the original GTF file (exon's start finishing after its end).					
				if int(line[3]) < int(line[4]):
					dict_all[line[0]].append(line)
					# When it's the first exon of the chromosome.
					if line[2] == "exon" and dict_transcript[str(current_chromosome)].is_empty() and len(transcript) == 0:
						last_transcript = line[transcript_col]
						min_transcript = int(line[3])
						max_transcript = int(line[4])
					# If the current exon is part of the same transcript than the previous exon.
					if line[2] == "exon" and line[transcript_col] == last_transcript:							
						# Add the new exon to the transcript and iterate start and end of the transcript.
						transcript.append(line)
						if int(line[3]) < min_transcript:
							min_transcript = int(line[3])
						if int(line[4]) > max_transcript:
							max_transcript = int(line[4])
					# When the new exon is part of a new transcript, create a new one and add the previous one as value.
					if line[2] == "exon" and line[transcript_col] != last_transcript:
						dict_transcript[str(current_chromosome)][int(min_transcript):int(max_transcript)] = transcript
						transcript = []
						last_transcript = line[transcript_col]
						transcript.append(line)
						min_transcript = int(line[3])
						max_transcript = int(line[4])
					# When it's the last transcript of the file.					
					if line[2] == "exon" and line[transcript_col] == final_transcript:
						# First exon of the final transcript.
						if line[transcript_col] == last_transcript:
							iv = Interval(int(line[3]), int(line[4]), [line])
							dict_transcript[str(current_chromosome)].add(iv)
							last_transcript = 0
							min_transcript = int(line[3])
							max_transcript = int(line[4])
							transcript = copy.deepcopy(line)
						# Other exons of the final transcript.
						else:
							transcript.append(line)
							if int(line[3]) < min_transcript:
								tree.removei(iv)
								iv = Interval(int(line[3]), int(max_transcript), line)
								dict_transcript[str(current_chromosome)].add(iv)
								min_transcript = int(line[3])
							if int(line[4]) > max_transcript:
								tree.removei(iv)
								iv = Interval(int(min_transcript), int(line[4]), line)
								dict_transcript[str(current_chromosome)].add(iv)
								max_transcript = int(line[4])		
				else:
					exons_wrong_interval.append(line)
		dict_transcript = copy.deepcopy(dict_transcript)
	return [dict_all, dict_transcript, header, exons_wrong_interval]


def precise_extension(dict_transcript, dict_exon_signal, gene_col):
	precisely_extended_dict = {}
	overlapped_transcripts = []
	coverage = 5000 # Average length of an exon = 200pb.
	# Boolean if the introns of a gene car be the exon of an other one.
	intron_exon = False
	for chromosome in dict_transcript:
		# Create a new dictionnary with the same model than dict_transcript.
		precisely_extended_dict[str(chromosome)] = IntervalTree()
		for transcript in sorted(dict_transcript[chromosome]):
			overlap_start = 0
			# Introduce the boolean extension with false as default for each transcript.
			extension = False
			# Case where the transcript is from the positive strand.
			if transcript.data[0][6] == "+":
				# Check if there is others transcripts in the area to extend.
				if len(dict_transcript[chromosome][transcript.end + 1:transcript.end + 5001]) != 0:
					exons_it = IntervalTree()
					introns_it = IntervalTree()
					max_extension = 0
					for transcript_in_iv in sorted(dict_transcript[chromosome][transcript.end + 1:transcript.end + 5001]):
						# If others transcripts are from the same strand but not the same gene, store in an IV the exons and the overlapping start.
						if transcript_in_iv.data[0][gene_col] != transcript.data[0][gene_col] and transcript_in_iv.data[0][6] == "+":
							if overlap_start == 0:
								if transcript_in_iv.begin > transcript.end:
									overlap_start = transcript_in_iv.begin
								# If transcripts are already overlapping before extension, error in the original GTF.
								else:
									overlap_start = transcript.end + 1
									overlapped_transcripts.append(transcript)
							for exon_in_transcript in transcript_in_iv.data:
								if int(exon_in_transcript[3]) > transcript.end:
									exons_it[int(exon_in_transcript[3])+ 1:int(exon_in_transcript[4])] = "exon"
								else:
									continue							
					# Comeback to the case where there is an overlapping issue.	
					if len(exons_it) > 1:
						# If there is a signal in the area where overlapping start in the stringtie output.
						if chromosome in dict_exon_signal:
							if len(dict_exon_signal[chromosome][overlap_start:transcript.end + 5001]) != 0 and intron_exon == True:
								exons_it.merge_overlaps()
								# Convert the exon intervaltree in a intron one.
								for exon_number, exons in enumerate(sorted(exons_it)):
									if exon_number == 0:
										previous_end = exons.end
									else:
										introns_it[previous_end + 1:exons.begin] = "intron"
										previous_end = exons.end
								# Check if signal overlap introns and assign max extension in consequence. 
								for signal in sorted(dict_exon_signal[chromosome][overlap_start:introns_it.end()], reverse = True):
									if signal.data[0] == "+":
										for intron in sorted(introns_it, reverse = True):
											if signal.end > intron.begin and signal.begin < intron.end and signal.end <= transcript.end + 5001:
												if signal.end < intron.end:
													max_extension = signal.end
												else:
													max_extension = intron.end
												extension = True
												break
										if max_extension != 0:
											new_transcript_end = max_extension
											break
									else:
										continue
								# Case where no signal overlap introns.
								if max_extension == 0:
									if len(dict_exon_signal[chromosome][transcript.end + 1: overlap_start]) != 0:
										for signal in sorted(dict_exon_signal[chromosome][transcript.end + 1:overlap_start], reverse = True):
											if signal.data[0] == "+" and signal.end <= transcript.end + 5001:
												new_transcript_end = signal.end
												extension = True
												break
									else:
										extension = False		
							# Case where no signal overlap transcripts.						
							else:
								if len(dict_exon_signal[chromosome][transcript.end + 1: overlap_start]) != 0:
									for signal in sorted(dict_exon_signal[chromosome][transcript.end + 1:overlap_start], reverse = True):
										if signal.data[0] == "+" and signal.end <= transcript.end + 5001 and signal.end < overlap_start:
											if signal.data[1] * (signal.end - signal.begin) > coverage:
												new_transcript_end = signal.end
												extension = True
												break
						else:
							extension = False
					elif len(exons_it) == 1:
						if chromosome in dict_exon_signal:
							if len(dict_exon_signal[chromosome][transcript.end + 1:exons_it.begin()]) != 0:
								for signal in sorted(dict_exon_signal[chromosome][transcript.end:exons_it.begin()], reverse = True):
									if signal.data[0] == "+" and signal.end <= exons_it.begin() - 1:
										if signal.data[1] * (signal.end - signal.begin) > coverage:
											new_transcript_end = signal.end
											extension = True
											break
							else:
								extension = False
						else:
							extension = False	
						
					else:
						# If there is a signal present from the stringtie output overlapping from the end of the transcript to an inputted value, save the signal's end.
						if chromosome in dict_exon_signal:
							if len(dict_exon_signal[chromosome][transcript.end + 1:transcript.end + 5001]) != 0:
								for signal in sorted(dict_exon_signal[chromosome][transcript.end:transcript.end + 5001], reverse = True):
									if signal.data[0] == "+" and signal.end <= transcript.end + 5001:
										if signal.data[1] * (signal.end - signal.begin) > coverage:
											new_transcript_end = signal.end
											extension = True
											break
							else:
								extension = False
						else:
							extension = False														
				# When extension is true, end of the transcript is changed with the signal's end and added to the new dict.
				if extension is True:
					modified_transcript = copy.deepcopy(transcript)
					modified_transcript.data[-1][4] = str(new_transcript_end)
					modified_transcript.data[-1][1] = "BestScriptEver"
					modified_transcript.data[-1].append("extension +" + str(new_transcript_end - transcript.end))
					precisely_extended_dict[chromosome][int(transcript.begin):int(new_transcript_end)] = modified_transcript.data
				# Otherwise, unmodified transcript is added.
				else:
					precisely_extended_dict[chromosome][int(transcript.begin):int(transcript.end)] = transcript.data
			
			# Case where the transcript is from the negative strand.
			if transcript.data[0][6] == "-":
				# Check if there is others transcripts in the area to extend.
				if len(dict_transcript[chromosome][transcript.begin - 5000:transcript.begin]) != 0:
					exons_it = IntervalTree()
					introns_it = IntervalTree()
					max_extension = 0
					for transcript_in_iv in sorted(dict_transcript[chromosome][transcript.begin - 5000:transcript.begin], reverse = True):
						# If others transcripts are from the same strand but not the same gene, store in an IV the exons and the overlapping start.
						if transcript_in_iv.data[0][gene_col] != transcript.data[0][gene_col] and transcript_in_iv.data[0][6] == "-":
							if overlap_start == 0:
								if transcript_in_iv.begin < transcript.begin:
									overlap_start = transcript_in_iv.end
								# If transcripts are already overlapping before extension, error in the original GTF.
								else:
									overlap_start = transcript.begin - 1
									overlapped_transcripts.append(transcript)
							for exon_in_transcript in transcript_in_iv.data:
								if int(exon_in_transcript[4]) < transcript.begin:
									exons_it[int(exon_in_transcript[3])+ 1:int(exon_in_transcript[4])] = "exon"
								else:
									continue
									
					# Comeback to the case where there is an overlapping issue.	
					if len(exons_it) > 1:
						# If there is a signal in the area where overlapping start in the stringtie output.
						if chromosome in dict_exon_signal:
							if len(dict_exon_signal[chromosome][transcript.begin - 5000:overlap_start + 1]) != 0 and intron_exon == True:
								exons_it.merge_overlaps()						
								# Convert the exon intervaltree in a intron one.
								for exon_number, exons in enumerate(sorted(exons_it)):
									if exon_number == 0:
										previous_end = exons.end
									else:
										introns_it[previous_end + 1:exons.begin] = "intron"
										previous_end = exons.end										
								# Check if signal overlap introns and assign max extension in consequence. 
								for signal in sorted(dict_exon_signal[chromosome][introns_it.begin():overlap_start + 1]):
									if signal.data[0] == "-":								
										for intron in sorted(introns_it):
											if signal.begin < intron.end and signal.end > intron.begin and signal.begin >= transcript.begin - 5000:
												if signal.begin > intron.begin:
													max_extension = signal.begin
												else:
													max_extension = intron.begin
												extension = True
												break
										if max_extension != 0:
											new_transcript_end = max_extension
											break
									else:
										continue
								# Case where no signal overlap introns.
								if max_extension == 0:
									if len(dict_exon_signal[chromosome][overlap_start: transcript.begin]) != 0:							
										for signal in sorted(dict_exon_signal[chromosome][overlap_start:transcript.begin]):
											if signal.data[0] == "-" and signal.begin >= transcript.begin - 5001:
												new_transcript_end = signal.begin
												extension = True
												break
									else:
										extension = False		
							# Case where no signal overlap transcripts.						
							else:					
								if len(dict_exon_signal[chromosome][overlap_start: transcript.begin]) != 0:
									for signal in sorted(dict_exon_signal[chromosome][overlap_start:transcript.begin]):
										if signal.data[0] == "-" and signal.begin >= transcript.begin - 5001 and signal.begin > overlap_start:									
											if signal.data[1] * (signal.end - signal.begin) > coverage:
												new_transcript_end = signal.begin
												extension = True
												break
						else:
							extension = False
					elif len(exons_it) == 1:
						if chromosome in dict_exon_signal:				
							if len(dict_exon_signal[chromosome][exons_it.end():transcript.begin]) != 0:
								for signal in sorted(dict_exon_signal[chromosome][exons_it.end():transcript.begin]):
									if signal.data[0] == "-" and signal.begin >= exons_it.end() + 1:
										if signal.data[1] * (signal.end - signal.begin) > coverage:
											new_transcript_end = signal.begin
											extension = True
											break
							else:
								extension = False
						else:
							extension = False	
					else:
						# If there is a signal present from the stringtie output overlapping from the end of the transcript to an inputted value, save the signal's end.
						if chromosome in dict_exon_signal:
							if len(dict_exon_signal[chromosome][transcript.begin - 5000:transcript.begin]) != 0:
								for signal in sorted(dict_exon_signal[chromosome][transcript.begin - 5000:transcript.begin]):
									if signal.data[0] == "-" and signal.begin >= transcript.begin - 5000:
										if signal.data[1] * (signal.end - signal.begin) > coverage:
											new_transcript_end = signal.begin
											extension = True
											break
							else:
								extension = False
						else:
							extension = False														
				# When extension is true, end of the transcript is changed with the signal's end and added to the new dict.
				if extension is True:
					modified_transcript = copy.deepcopy(transcript)
					modified_transcript.data[0][3] = str(new_transcript_end)
					modified_transcript.data[0][1] = "BestScriptEver"
					modified_transcript.data[0].append("extension " + str(new_transcript_end - transcript.begin))
					precisely_extended_dict[chromosome][int(new_transcript_end):int(transcript.end)] = modified_transcript.data
				# Otherwise, unmodified transcript is added.
				else:
					precisely_extended_dict[chromosome][int(transcript.begin):int(transcript.end)] = transcript.data
	with open("errors_file.txt", "w") as filout:
		for ovlp_transcript in overlapped_transcripts:
			filout.write("{}\n".format(ovlp_transcript.data[0]))
	return precisely_extended_dict


def dict_to_gtf(dict_all, filenameout):
	"""Create a output file similar to the gtf input with the gene extended."""
	with open(filenameout, "w") as filout:
		for chromosome in dict_all:		
			for transcript in dict_all[chromosome]:
				for exon in transcript.data:
					exon_gtf = "\t".join(exon[0:9]) + "; " + "; ".join(exon[9:]) + ";"
					filout.write("{}\n".format(exon_gtf))
				

if __name__ == "__main__":

	
	#print(str(datetime.now()))
	argvals = None
	args = get_args(argvals)
	gene_col, transcript_col, first_chromosome, final_transcript = find_pos(args.filename)
	dict_exon_signal = find_exon(args.stringtie_output)
	dict_all, dict_transcript, header, exons_wrong_interval = gtf_to_dict(args.filename, first_chromosome, final_transcript, gene_col, transcript_col)
	precisely_extended_dict = precise_extension(dict_transcript, dict_exon_signal, gene_col)
	if args.output is False:
		dict_to_gtf(precisely_extended_dict, args.filename + ".extended")
	else:	
		dict_to_gtf(precisely_extended_dict, args.output)
	
	
	


